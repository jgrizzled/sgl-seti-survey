"""Catalog screening for the Pan-STARRS1 pilot (plan §3.4 layer 1) —
port of surveys/ztf/scripts/catalog_screen.py to the DR2 ``detection``
table.

Per corridor: one cone query of the DR2 per-epoch detection table (all
filters, paged, snapshotted) plus one of the ``mean`` object table as
recurrence context. Every usable/partial precise-pass warp defines an
exposure (filter, MJD); detections with |obsTime - MJD_start| < 120 s in
that filter are matched against the endpoint x role locus at
mid-exposure. Detections within the screen radius become ScreenMatch
records (distance to track, implied relay distance, AB magnitude from
psfFlux).

No catalog is treated as complete (plan §3.4).

Usage:
    uv run python surveys/panstarrs/scripts/catalog_screen.py [corridor ...]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time as _time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.mast_ps1 import (DETECTION_COLS, FILTER_ID,
                                         MEAN_COLS, catalog_cone)
from sglsurvey.geometry import GeometryContext, discovery_cone
from sglsurvey.records import (Observation, ScreenMatch, append_records,
                               read_records)
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ps1_corridors import MEMBERS as CORRIDORS  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "panstarrs" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "panstarrs" / "precise_v1"
RUN_DIR = REPO / "runs" / "panstarrs" / "screen_v1"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "panstarrs" / "hypotheses.md"
HYPOTHESIS_VERSION = "ps1-hypotheses-v1.0"

SCREEN_RADIUS_ARCSEC = 10.0   # frozen padding (hypotheses v1.0 section 5)
LOCUS_TOLERANCE_ARCSEC = 1.0
#: detection.obsTime sits ~50-60 s after the warp MJD-OBS (exposure end
#: rather than start, measured 2026-08-20); TTI pairs are >= 10 min apart.
EXPOSURE_MATCH_DAYS = 120.0 / 86400.0
CONTEXT_RADIUS_DEG = 0.2
TIME_RANGE_MJD = (54900.0, 57300.0)
ROLES = (Role.RX, Role.TX)


def seg_min_dist_vec(det_ra, det_dec, pts):
    """Vectorised min distance (arcsec) from many points to a polyline,
    returning (dist, nearest-segment-index) arrays."""
    cosd = np.cos(np.deg2rad(det_dec))[:, None]
    x = (pts[None, :, 0] - det_ra[:, None]) * cosd * 3600.0
    y = (pts[None, :, 1] - det_dec[:, None]) * 3600.0
    ax, ay, bx, by = x[:, :-1], y[:, :-1], x[:, 1:], y[:, 1:]
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = np.clip(np.where(L2 > 0, -(ax * dx + ay * dy)
                         / np.where(L2 > 0, L2, 1.0), 0.0), 0.0, 1.0)
    d = np.hypot(ax + t * dx, ay + t * dy)
    i = np.argmin(d, axis=1)
    return d[np.arange(len(d)), i], i


def _f(row, key, default=np.nan):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("corridors", nargs="*", default=[])
    args = ap.parse_args()
    corridors = args.corridors or list(CORRIDORS)

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.ps1_default()
    store = SnapshotStore(RUN_DIR)
    session = requests.Session()
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()

    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs_by_id[r["observation_id"]] = Observation(**r)
    # usable exposures per (endpoint, role): {(band, mjd_start): [obs_ids]}
    usable = defaultdict(lambda: defaultdict(list))
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            o = obs_by_id[r["observation_id"]]
            usable[(r["endpoint_id"], r["role"])][
                (o.band, round(o.t_start_mjd_utc, 5))].append(o.observation_id)

    match_path = RUN_DIR / "records" / "screen_match.jsonl"
    known = ({r["screen_match_id"] for r in read_records(match_path)}
             if match_path.exists() else set())
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    run_config = {
        "run_id": "screen_v1", "stage": "catalog-screen",
        "source_runs": ["coarse_v1", "precise_v1"],
        "corridors": corridors,
        "screen_radius_arcsec": SCREEN_RADIUS_ARCSEC,
        "locus_tolerance_arcsec": LOCUS_TOLERANCE_ARCSEC,
        "exposure_match_seconds": EXPOSURE_MATCH_DAYS * 86400.0,
        "source_catalog": "PS1 DR2 detection (MAST catalogs API)",
        "context_table": "PS1 DR2 mean",
        "context_radius_deg": CONTEXT_RADIUS_DEG,
        "registry_source_hash": registry.source_hash,
        "hypothesis_version": HYPOTHESIS_VERSION,
        "hypothesis_hash": hyp_hash,
        "geometry": ctx.identities(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    (RUN_DIR / f"run_config_{run_config['started_utc'][:19]}.json"
     ).write_text(json.dumps(run_config, indent=2))

    summary = defaultdict(lambda: [0, 0, 0])  # matches, exposures, no-det
    catalog_stats = {}
    for corridor in corridors:
        members = list(CORRIDORS[corridor])
        # cone covering every member's discovery envelope
        cones = [discovery_cone(ctx, registry[e], ROLES, *TIME_RANGE_MJD)
                 for e in members]
        cra = float(np.mean([c.ra_deg for c in cones]))
        cdec = float(np.mean([c.dec_deg for c in cones]))
        sep = max(np.hypot((c.ra_deg - cra) * np.cos(np.deg2rad(cdec)),
                           c.dec_deg - cdec) + c.radius_deg for c in cones)
        radius = float(sep + 0.02)
        t0 = _time.monotonic()
        dets, det_snaps = catalog_cone("detection", cra, cdec, radius,
                                       DETECTION_COLS, store, session)
        means, mean_snaps = catalog_cone("mean", cra, cdec,
                                         max(radius, CONTEXT_RADIUS_DEG),
                                         MEAN_COLS, store, session)
        print(f"[{corridor}] cone ({cra:.4f},{cdec:.4f}) r={radius:.3f}: "
              f"{len(dets)} detections, {len(means)} mean objects "
              f"({_time.monotonic() - t0:.0f}s)", flush=True)
        catalog_stats[corridor] = {
            "center": [cra, cdec], "radius_deg": radius,
            "n_detections": len(dets), "n_mean_objects": len(means),
            "detection_snapshots": det_snaps, "mean_snapshots": mean_snaps}
        det_t = np.array([_f(d, "obsTime") for d in dets])
        det_f = np.array([int(_f(d, "filterID", 0)) for d in dets])
        det_ra = np.array([_f(d, "ra") for d in dets])
        det_dec = np.array([_f(d, "dec") for d in dets])
        order = np.argsort(det_t)

        new = []
        for e in members:
            for role in ROLES:
                exposures = usable.get((e, role.value), {})
                for (band, mjd0), oids in sorted(exposures.items(),
                                                 key=lambda kv: kv[0][1]):
                    fid = {v: k for k, v in FILTER_ID.items()}[band]
                    lo = np.searchsorted(det_t, mjd0 - EXPOSURE_MATCH_DAYS,
                                         sorter=order)
                    hi = np.searchsorted(det_t, mjd0 + EXPOSURE_MATCH_DAYS,
                                         sorter=order)
                    idx = order[lo:hi]
                    idx = idx[det_f[idx] == fid]
                    summary[(e, role.value)][1] += 1
                    if len(idx) == 0:
                        summary[(e, role.value)][2] += 1
                        continue
                    obs = obs_by_id[oids[0]]
                    t_mid = Time(obs.t_mid_mjd_utc, format="mjd")
                    al = adaptive_locus(
                        target=registry[e], role=role, observation_time=t_mid,
                        observer=ctx.observer, relay_range=ctx.relay_range,
                        tolerance_arcsec=LOCUS_TOLERANCE_ARCSEC,
                        ephemeris=ctx.ephemeris, model=ctx.model)
                    pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                                    for p in al.points])
                    zs = np.array([p.z_au for p in al.points])
                    m = SCREEN_RADIUS_ARCSEC / 3600.0 * 1.5
                    cosd = np.cos(np.deg2rad(pts[:, 1].mean()))
                    sel = ((det_dec[idx] > pts[:, 1].min() - m)
                           & (det_dec[idx] < pts[:, 1].max() + m)
                           & (det_ra[idx] > pts[:, 0].min() - m / cosd)
                           & (det_ra[idx] < pts[:, 0].max() + m / cosd))
                    idx = idx[sel]
                    if len(idx) == 0:
                        continue
                    d, seg = seg_min_dist_vec(det_ra[idx], det_dec[idx], pts)
                    for k, dist, s in zip(idx, d, seg):
                        if dist > SCREEN_RADIUS_ARCSEC:
                            continue
                        row = dets[k]
                        flux = _f(row, "psfFlux")
                        ferr = _f(row, "psfFluxErr")
                        mag = (-2.5 * np.log10(flux / 3631.0)
                               if flux > 0 else None)
                        sm = ScreenMatch.build(
                            endpoint_id=e, role=role.value,
                            hypothesis_version=HYPOTHESIS_VERSION,
                            registry_source_hash=registry.source_hash,
                            source_table="ps1-dr2-detection",
                            source_cntr=str(row["detectID"]),
                            source_designation=str(row["objID"]),
                            ra_deg=float(det_ra[k]), dec_deg=float(det_dec[k]),
                            sigra_mas=None, sigdec_mas=None,
                            mjd=float(det_t[k]), scan_id=None, frame_num=None,
                            photometry={
                                "filter": band,
                                "psf_flux_jy": flux, "psf_flux_err_jy": ferr,
                                "mag_ab": mag,
                                "snr": (flux / ferr if ferr > 0 else None),
                                "zp": _f(row, "zp"),
                                "psf_qf_perfect": _f(row, "psfQfPerfect")},
                            flags={"infoFlag": int(_f(row, "infoFlag", 0)),
                                   "infoFlag2": int(_f(row, "infoFlag2", 0)),
                                   "infoFlag3": int(_f(row, "infoFlag3", 0)),
                                   "imageID": row.get("imageID")},
                            dist_arcsec=round(float(dist), 3),
                            implied_z_au=round(float(zs[s]), 2),
                            z_segment_au=(round(float(zs[s]), 2),
                                          round(float(zs[s + 1]), 2)),
                            screen_radius_arcsec=SCREEN_RADIUS_ARCSEC,
                            locus_mjd=obs.t_mid_mjd_utc, snapshot_id=None,
                            extra={"observation_ids": oids,
                                   "exposure_mjd_start": mjd0},
                        )
                        if sm.screen_match_id not in known:
                            known.add(sm.screen_match_id)
                            new.append(sm)
                        summary[(e, role.value)][0] += 1
        append_records(match_path, new)
        print(f"[{corridor}] {len(new)} new matches", flush=True)

    (RUN_DIR / "catalog_stats.json").write_text(
        json.dumps(catalog_stats, indent=2))
    print("\n=== screen summary ===")
    for (e, r), (n, nexp, nodet) in sorted(summary.items()):
        print(f"{e:16s} {r}: {n} matches over {nexp} exposures "
              f"({nodet} exposures with no catalogued detections)")


if __name__ == "__main__":
    main()
