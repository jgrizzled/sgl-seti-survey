"""Precise intersection pass for the Pan-STARRS1 pilot (port of
surveys/ztf/scripts/precise_pass.py).

Consumes coarse-stage hits from runs/panstarrs/coarse_v1, fetches the
FULL skycell mask per hit warp (3.3 MB fpack; exact TAN WCS + IPP mask
bits), evaluates sglseti's covered_z_intervals against the usable-pixel
footprint, and records a per-warp cutout spec (centre = union of the Rx
and Tx loci at mid-exposure, size from the locus extent) for the later
image/weight cutout fetch. Header EXPTIME/MJD-OBS/seeing/ZP are read
from the mask header and carried in the cutout index. Emits
precise-stage IntersectionEvaluation records.

Usage:
    uv run python surveys/panstarrs/scripts/precise_pass.py [endpoint ...]
                  [--limit N] [--workers N]
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
from astropy.time import Time

from sglseti import (Role, adaptive_locus, covered_z_intervals,
                     load_target_registry, stable_hash)

from sglsurvey.adapters.base import CutoutSpec
from sglsurvey.adapters.mast_ps1 import (PIX_ARCSEC, Ps1ExactFootprint,
                                         Ps1WarpAdapter)
from sglsurvey.geometry import GeometryContext, enclosing_cone
from sglsurvey.records import (IntersectionEvaluation, Observation,
                               append_records, read_records)

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "panstarrs" / "coarse_v1"
RUN_DIR = REPO / "runs" / "panstarrs" / "precise_v1"
PRODUCT_DIR = REPO / "runs" / "panstarrs" / "products" / "msk"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "panstarrs" / "hypotheses.md"
HYPOTHESIS_VERSION = "ps1-hypotheses-v1.0"

PRECISE_TOLERANCE_ARCSEC = 1.0
SEED_STEP_ARCSEC = 10.0  # usable stretches narrower than this may be missed
#: Cutout size per warp from the locus extent (<= ~375" at z=550 AU):
#: 2 x (radius + margin) at 0.25"/pix, rounded up to a multiple of 100.
CUTOUT_MARGIN_ARCSEC = 15.0
CUTOUT_MIN_PIX, CUTOUT_MAX_PIX = 400, 3600


def cutout_size_pix(radius_deg: float) -> int:
    px = 2.0 * (radius_deg * 3600.0 + CUTOUT_MARGIN_ARCSEC) / PIX_ARCSEC
    px = int(np.ceil(px / 100.0) * 100)
    return int(np.clip(px, CUTOUT_MIN_PIX, CUTOUT_MAX_PIX))
ROLE_BY_NAME = {r.value: r for r in Role}


def load_coarse(endpoints: set[str] | None):
    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs_by_id[r["observation_id"]] = Observation(**r)
    hits_by_obs = defaultdict(list)
    for r in read_records(COARSE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] != "coarse" or not r["hit"]:
            continue
        if endpoints and r["endpoint_id"] not in endpoints:
            continue
        hits_by_obs[r["observation_id"]].append(
            (r["endpoint_id"], r["role"]))
    return obs_by_id, hits_by_obs


def cutout_index_path() -> Path:
    return RUN_DIR / "records" / "cutout_index.jsonl"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("endpoints", nargs="*")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    endpoints = set(args.endpoints) or None

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.ps1_default()
    adapter = Ps1WarpAdapter()
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()
    target_hashes = {tid: stable_hash(registry[tid])
                     for tid in registry.ids}

    obs_by_id, hits_by_obs = load_coarse(endpoints)
    obs_ids = sorted(hits_by_obs,
                     key=lambda o: obs_by_id[o].t_mid_mjd_utc)
    if args.limit:
        obs_ids = obs_ids[:args.limit]
    n_evals = sum(len(hits_by_obs[o]) for o in obs_ids)
    print(f"{len(obs_ids)} hit warps, {n_evals} precise "
          f"evaluations", flush=True)

    ixn_path = RUN_DIR / "records" / "intersection_evaluation.jsonl"
    known_ixn = ({r["intersection_id"] for r in read_records(ixn_path)}
                 if ixn_path.exists() else set())
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    run_config = {
        "run_id": "precise_v1", "stage": "precise",
        "source_run": "coarse_v1",
        "endpoints": sorted(endpoints) if endpoints else "all",
        "limit": args.limit,
        "tolerance_arcsec": PRECISE_TOLERANCE_ARCSEC,
        "seed_step_arcsec": SEED_STEP_ARCSEC,
        "cutout_margin_arcsec": CUTOUT_MARGIN_ARCSEC,
        "cutout_pix_range": [CUTOUT_MIN_PIX, CUTOUT_MAX_PIX],
        "mask_product": "full skycell .mask.fits (fpack)",
        "fatal_mask_template": Ps1ExactFootprint.FATAL_MASK,
        "registry_source_hash": registry.source_hash,
        "hypothesis_version": HYPOTHESIS_VERSION,
        "hypothesis_hash": hyp_hash,
        "geometry": ctx.identities(),
        "interval_model": ("mid-exposure evaluation with nominal per-filter "
                           "exposure times; <= 60 s frames move the locus "
                           "< 5 mas"),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    (RUN_DIR / f"run_config_{run_config['started_utc'][:19]}.json"
     ).write_text(json.dumps(run_config, indent=2))

    summary = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0.0, 0]))
    t0 = _time.monotonic()
    new_records, cut_records = [], []

    # Phase 1 (CPU): loci + cutout spec per exposure.
    def plan(oid):
        obs = obs_by_id[oid]
        t_mid = Time(obs.t_mid_mjd_utc, format="mjd")
        loci = {}
        for endpoint_id, role_name in hits_by_obs[oid]:
            al = adaptive_locus(
                target=registry[endpoint_id], role=ROLE_BY_NAME[role_name],
                observation_time=t_mid, observer=ctx.observer,
                relay_range=ctx.relay_range,
                tolerance_arcsec=PRECISE_TOLERANCE_ARCSEC,
                ephemeris=ctx.ephemeris, model=ctx.model)
            loci[(endpoint_id, role_name)] = al
        all_pts = np.concatenate([
            np.array([[p.icrs_ra_deg, p.icrs_dec_deg] for p in al.points])
            for al in loci.values()])
        cone = enclosing_cone(all_pts, 0.0)
        return loci, CutoutSpec(ra_deg=cone.ra_deg, dec_deg=cone.dec_deg,
                                size_pix=cutout_size_pix(cone.radius_deg))

    # Phase 2 (IO): threaded cutout fetch.
    def fetch(oid, cutout):
        import requests as _rq
        ad = Ps1WarpAdapter(session=_rq.Session())
        for attempt in (1, 2, 3):
            try:
                return ad.fetch(obs_by_id[oid], ["msk"], PRODUCT_DIR)
            except FileNotFoundError:
                return None
            except Exception as exc:
                if attempt == 3:
                    print(f"  fetch failed {oid}: {exc}", flush=True)
                    return None
                _time.sleep(2.0 * attempt)

    from concurrent.futures import ThreadPoolExecutor
    CHUNK = 200
    for c0 in range(0, len(obs_ids), CHUNK):
        chunk = obs_ids[c0:c0 + CHUNK]
        plans = {oid: plan(oid) for oid in chunk}
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            fetched = dict(zip(chunk, ex.map(
                lambda o: fetch(o, plans[o][1]), chunk)))
        for j, oid in enumerate(chunk):
            i = c0 + j
            obs = obs_by_id[oid]
            t_mid = Time(obs.t_mid_mjd_utc, format="mjd")
            loci, cutout = plans[oid]
            pset = fetched[oid]
            available = pset is not None
            fp = (Ps1ExactFootprint(pset.products[0].path) if available
                  else None)
            hdr = fp.header if fp is not None else {}
            cut_records.append({"observation_id": oid, "kind": "msk",
                                "available": available,
                                "center_ra_deg": cutout.ra_deg,
                                "center_dec_deg": cutout.dec_deg,
                                "size_pix": cutout.size_pix,
                                "header": {
                                    "exptime": hdr.get("EXPTIME"),
                                    "mjd_obs": hdr.get("MJD-OBS"),
                                    "seeing_pix": hdr.get("CHIP.SEEING"),
                                    "airmass": hdr.get("AIRMASS"),
                                    "fpa_zp": hdr.get("FPA.ZP"),
                                    "gain": hdr.get("CELL.GAIN")},
                                "path": (str(pset.products[0].path.relative_to(REPO))
                                         if available else None),
                                "checksum": (pset.products[0].checksum
                                             if available else None)})
            # --- evaluation ---
            bad_quality = bool(obs.quality_flags.get("badflag"))
            for (endpoint_id, role_name), al in loci.items():
                role = ROLE_BY_NAME[role_name]
                target = registry[endpoint_id]
                pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                                for p in al.points])
                if fp is None:
                    hit, frac, covered = False, None, []
                    label, note = "unknown", "product not served (404)"
                else:
                    inb, usable = fp.status(pts)
                    frac = float(usable[inb].mean()) if inb.any() else None
                    covered = covered_z_intervals(
                        target=target, role=role, observation_time=t_mid,
                        observer=ctx.observer, relay_range=ctx.relay_range,
                        contains=fp.contains,
                        tolerance_arcsec=PRECISE_TOLERANCE_ARCSEC,
                        ephemeris=ctx.ephemeris, model=ctx.model, locus=al,
                        seed_step_arcsec=SEED_STEP_ARCSEC)
                    hit = len(covered) > 0
                    if bad_quality:
                        label, note = "unusable", "listing badflag set"
                    elif not hit:
                        label, note = "unusable", "no usable-pixel crossing"
                    elif frac is not None and frac > 0.99:
                        label, note = "usable", None
                    else:
                        label, note = "partial", None
                ixn = IntersectionEvaluation.build(
                    observation_id=oid,
                    hypothesis_version=HYPOTHESIS_VERSION,
                    hypothesis_hash=hyp_hash,
                    endpoint_id=endpoint_id, role=role_name,
                    registry_source_hash=registry.source_hash,
                    target_source_hash=target_hashes[endpoint_id],
                    model_id=ctx.model.model_id,
                    model_version=ctx.model.model_version,
                    ephemeris_id=ctx.identities()["ephemeris_id"],
                    tolerance_arcsec=PRECISE_TOLERANCE_ARCSEC,
                    confidence_level=ctx.confidence_level,
                    padding_arcsec=0.0,
                    stage="precise", hit=hit,
                    covered_z_intervals_au=tuple(
                        (zi.z_min_au, zi.z_max_au) for zi in covered),
                    envelope_pad_arcsec=0.0,
                    warnings=al.warnings,
                    usable=label, usable_fraction=frac,
                    usability_notes=note,
                    extra={"seed_step_arcsec": SEED_STEP_ARCSEC,
                           "fatal_mask_template": Ps1ExactFootprint.FATAL_MASK,
                           "bad_quality": bad_quality,
                           "product_available": available,
                           "cutout_pix": cutout.size_pix,
                           "skycell": obs.extra["skycell"]},
                )
                if ixn.intersection_id not in known_ixn:
                    known_ixn.add(ixn.intersection_id)
                    new_records.append(ixn)
                s = summary[(endpoint_id, role_name)][obs.band]
                s[0] += 1
                s[1] += hit
                s[2] += label in ("usable", "partial")
                if frac is not None:
                    s[3] += frac
                s[4] += (not available)
            if (i + 1) % 50 == 0:
                rate = (i + 1) / (_time.monotonic() - t0)
                print(f"  {i + 1}/{len(obs_ids)} exposures ({rate:.2f}/s)",
                      flush=True)
            if len(new_records) >= 200:
                append_records(ixn_path, new_records)
                new_records = []
                with cutout_index_path().open("a") as fh:
                    for rec in cut_records:
                        fh.write(json.dumps(rec) + "\n")
                cut_records = []
    append_records(ixn_path, new_records)
    with cutout_index_path().open("a") as fh:
        for rec in cut_records:
            fh.write(json.dumps(rec) + "\n")

    print("\n=== precise summary (evaluated / geometric-hit / usable, "
          "mean usable frac, missing products) ===")
    for (eid, role), bands in sorted(summary.items()):
        parts = []
        for b in sorted(bands):
            n, h, u, fsum, miss = bands[b]
            parts.append(f"{b}: {n}/{h}/{u} f={fsum / max(n, 1):.2f} "
                         f"miss={miss}")
        print(f"{eid:16s} {role}: " + "  ".join(parts))


if __name__ == "__main__":
    main()
