"""Coarse discovery pass for the DECam pilot (coarse stage of
Pipeline A, plan §3.3) — port of surveys/panstarrs/scripts/
coarse_discovery.py to the NOIRLab Astro Archive adapter.

For each endpoint x role: build a conservative discovery cone with
sglseti (CTIO observer), advanced-search every instcal exposure whose
focal plane can touch the cone (raw JSON snapshotted; image/dqmask/
wtmap joined on EXPNUM), then evaluate the adaptive locus at
exposure-clustered epochs against the static focal-plane layout placed
at each exposure's pointing centre. Emits Observation and coarse-stage
IntersectionEvaluation records (hits AND misses).

Usage:
    uv run python surveys/decam/scripts/coarse_discovery.py [endpoint ...]
Defaults to the 3 pilot corridors; records land in runs/decam/coarse_v1/.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time as _time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, load_target_registry, stable_hash

from sglsurvey.adapters.base import MjdRange
from sglsurvey.adapters.noirlab_decam import DecamInstcalAdapter
from sglsurvey.geometry import GeometryContext, discovery_cone, locus_radec
from sglsurvey.records import (IntersectionEvaluation, append_records,
                               read_records)
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from decam_corridors import PILOT_ENDPOINTS  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
RUN_ID = "coarse_v1"
RUN_DIR = REPO / "runs" / "decam" / RUN_ID
LAYOUT_PATH = (REPO / "surveys" / "decam" / "configs"
               / "decam_focal_plane_v1.json")
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "decam" / "hypotheses.md"
HYPOTHESIS_VERSION = "decam-hypotheses-v0.1-draft"

# DECam science operations: 2012-09 (commissioning/SV) through today;
# generous on both ends. Yearly refresh extends the stop epoch.
TIME_RANGE = MjdRange(start_mjd_utc=56000.0, stop_mjd_utc=62000.0)
ROLES = (Role.RX, Role.TX)
# Corridor drift at z_min is ~6.5"/day; a 0.05 d bin moves the locus
# < 0.4", covered by the cluster pad.
CLUSTER_BIN_DAYS = 0.05
CLUSTER_PAD_ARCSEC = 0.5
# Static-layout residual vs per-era TPV solutions: worst observed
# corner mismatch 10.2" (2012 era; configs/decam_focal_plane_v1.json
# validation block) — 15" covers it.
LAYOUT_PAD_ARCSEC = 15.0


def main(endpoint_ids: list[str]) -> None:
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.decam_default()
    adapter = DecamInstcalAdapter(layout_path=LAYOUT_PATH)
    store = SnapshotStore(RUN_DIR)
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()
    layout_doc = json.loads(LAYOUT_PATH.read_text())

    obs_path = RUN_DIR / "records" / "observation.jsonl"
    ixn_path = RUN_DIR / "records" / "intersection_evaluation.jsonl"
    known_obs = ({r["observation_id"] for r in read_records(obs_path)}
                 if obs_path.exists() else set())
    known_ixn = ({r["intersection_id"] for r in read_records(ixn_path)}
                 if ixn_path.exists() else set())

    run_config = {
        "run_id": RUN_ID, "stage": "coarse",
        "endpoints": endpoint_ids, "roles": [r.value for r in ROLES],
        "time_range_mjd": [TIME_RANGE.start_mjd_utc,
                           TIME_RANGE.stop_mjd_utc],
        "cluster_bin_days": CLUSTER_BIN_DAYS,
        "cluster_pad_arcsec": CLUSTER_PAD_ARCSEC,
        "layout_pad_arcsec": LAYOUT_PAD_ARCSEC,
        "focal_plane_layout": {
            "path": str(LAYOUT_PATH.relative_to(REPO)),
            "version": layout_doc["version"],
            "reference_dqmask_md5": layout_doc["reference_dqmask_md5"]},
        "registry": str(REGISTRY_PATH.relative_to(REPO)),
        "registry_source_hash": registry.source_hash,
        "hypothesis_version": HYPOTHESIS_VERSION,
        "hypothesis_hash": hyp_hash,
        "geometry": ctx.identities(),
        "archive": {"adapter": "DecamInstcalAdapter",
                    "collection": adapter.collection,
                    "listing": "adv_search/find box, EXPNUM-joined"},
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    cfg_path = RUN_DIR / f"run_config_{run_config['started_utc'][:19]}.json"
    cfg_path.write_text(json.dumps(run_config, indent=2))

    summary = {}
    for endpoint_id in endpoint_ids:
        target = registry[endpoint_id]
        t0 = _time.monotonic()
        cone = discovery_cone(ctx, target, ROLES,
                              TIME_RANGE.start_mjd_utc,
                              TIME_RANGE.stop_mjd_utc)
        print(f"[{endpoint_id}] cone ra={cone.ra_deg:.4f} "
              f"dec={cone.dec_deg:.4f} r={cone.radius_deg:.4f} deg "
              f"({_time.monotonic() - t0:.0f}s)", flush=True)

        t0 = _time.monotonic()
        observations = list(adapter.discover(cone, TIME_RANGE, store))
        bands = sorted({o.band for o in observations})
        print(f"[{endpoint_id}] discovered {len(observations)} instcal "
              f"exposures, bands {bands} "
              f"({_time.monotonic() - t0:.0f}s)", flush=True)

        new_obs = [o for o in observations
                   if o.observation_id not in known_obs]
        known_obs.update(o.observation_id for o in new_obs)
        append_records(obs_path, new_obs)

        clusters = defaultdict(list)
        for o in observations:
            clusters[round(o.t_mid_mjd_utc / CLUSTER_BIN_DAYS)].append(o)

        t0 = _time.monotonic()
        for role in ROLES:
            hits = defaultdict(int)
            total = defaultdict(int)
            new_ixn = []
            for key, cluster in sorted(clusters.items()):
                t_center = Time(np.mean([o.t_mid_mjd_utc for o in cluster]),
                                format="mjd")
                pts, warns = locus_radec(ctx, target, role, t_center)
                pad = (ctx.padding_arcsec + CLUSTER_PAD_ARCSEC
                       + LAYOUT_PAD_ARCSEC)
                for o in cluster:
                    fp = adapter.nominal_footprint(o, pad_arcsec=pad)
                    hit = fp.contains_any(pts)
                    total[o.band] += 1
                    hits[o.band] += hit
                    ixn = IntersectionEvaluation.build(
                        observation_id=o.observation_id,
                        hypothesis_version=HYPOTHESIS_VERSION,
                        hypothesis_hash=hyp_hash,
                        endpoint_id=endpoint_id, role=role.value,
                        registry_source_hash=registry.source_hash,
                        target_source_hash=stable_hash(target),
                        model_id=ctx.model.model_id,
                        model_version=ctx.model.model_version,
                        ephemeris_id=ctx.identities()["ephemeris_id"],
                        tolerance_arcsec=ctx.tolerance_arcsec,
                        confidence_level=ctx.confidence_level,
                        padding_arcsec=ctx.padding_arcsec,
                        stage="coarse", hit=bool(hit), warnings=warns,
                        extra={"cluster_pad_arcsec": CLUSTER_PAD_ARCSEC,
                               "layout_pad_arcsec": LAYOUT_PAD_ARCSEC,
                               "cluster_mjd": float(t_center.mjd),
                               "expnum": o.native_key["expnum"]},
                    )
                    if ixn.intersection_id not in known_ixn:
                        known_ixn.add(ixn.intersection_id)
                        new_ixn.append(ixn)
            append_records(ixn_path, new_ixn)
            summary[(endpoint_id, role.value)] = {
                b: (hits[b], total[b]) for b in sorted(total)}
            print(f"[{endpoint_id}] {role.value}: "
                  + "  ".join(f"{b}: {hits[b]}/{total[b]}"
                              for b in sorted(total))
                  + f"  ({_time.monotonic() - t0:.0f}s)", flush=True)

    print("\n=== coarse summary (hits/evaluated per filter) ===")
    for (eid, role), bands in summary.items():
        print(f"{eid:16s} {role}: "
              + "  ".join(f"{b}: {h}/{t}" for b, (h, t) in bands.items()))


if __name__ == "__main__":
    main(sys.argv[1:] or PILOT_ENDPOINTS)
