"""Precise intersection pass for the WISE/NEOWISE shakedown (plan §4.4,
precise stage of Pipeline A).

Consumes the coarse-stage hits from runs/wise/coarse_v1, downloads the
-msk product for every hit frame-band (md5-verified; ~130 KB each, and
the mask header carries the same full WCS as the -int frame, so masks
alone provide exact astrometry AND valid-pixel data), then evaluates
sglseti's covered_z_intervals against the exact usable-pixel footprint.
Emits precise-stage IntersectionEvaluation records with the covered
relay-distance intervals (possibly disjoint) and usable-pixel fractions.

Usage:
    uv run python surveys/wise/scripts/precise_pass.py [endpoint ...]
                  [--limit N]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time as _time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import (Role, adaptive_locus, covered_z_intervals,
                     load_target_registry, stable_hash)

from sglsurvey.adapters.irsa_wise import WiseExactFootprint, WiseMergeL1bAdapter
from sglsurvey.geometry import GeometryContext
from sglsurvey.records import (IntersectionEvaluation, Observation,
                               append_records, read_records)

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "wise" / "coarse_v1"
RUN_DIR = REPO / "runs" / "wise" / "precise_v1"
PRODUCT_DIR = REPO / "runs" / "wise" / "products" / "msk"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "wise" / "hypotheses.md"
HYPOTHESIS_VERSION = "wise-hypotheses-v1.0"

PRECISE_TOLERANCE_ARCSEC = 2.0
SEED_STEP_ARCSEC = 30.0  # usable stretches narrower than this may be missed
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


def msk_path(obs: Observation) -> Path:
    name = obs.products["msk"]["url"].rsplit("/", 1)[-1]
    return PRODUCT_DIR / name


def download_masks(adapter, obs_by_id, obs_ids) -> None:
    todo = [oid for oid in obs_ids if not msk_path(obs_by_id[oid]).exists()]
    print(f"downloading {len(todo)} msk products "
          f"({len(obs_ids) - len(todo)} cached)", flush=True)
    t0 = _time.monotonic()
    for i, oid in enumerate(todo):
        obs = obs_by_id[oid]
        for attempt in (1, 2, 3):
            try:
                adapter.fetch(obs, ["msk"], PRODUCT_DIR)
                break
            except Exception as exc:
                if attempt == 3:
                    raise
                print(f"  retry {attempt} for {oid}: {exc}", flush=True)
                _time.sleep(2.0 * attempt)
        if (i + 1) % 200 == 0:
            rate = (i + 1) / (_time.monotonic() - t0)
            print(f"  {i + 1}/{len(todo)} ({rate:.1f}/s)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("endpoints", nargs="*")
    ap.add_argument("--limit", type=int, default=None,
                    help="max observations to evaluate (testing)")
    args = ap.parse_args()
    endpoints = set(args.endpoints) or None

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.wise_coarse_default()
    adapter = WiseMergeL1bAdapter()
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
    print(f"{len(obs_ids)} hit frame-bands, {n_evals} precise evaluations",
          flush=True)

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
        "fatal_mask_bits": list(WiseExactFootprint.FATAL_BITS),
        "registry_source_hash": registry.source_hash,
        "hypothesis_version": HYPOTHESIS_VERSION,
        "hypothesis_hash": hyp_hash,
        "geometry": ctx.identities(),
        "interval_model": ("midpoint evaluation; 7.7 s frames move the "
                           "locus < 1 mas over the exposure"),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    (RUN_DIR / f"run_config_{run_config['started_utc'][:19]}.json"
     ).write_text(json.dumps(run_config, indent=2))

    download_masks(adapter, obs_by_id, obs_ids)

    summary = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0.0]))
    t0 = _time.monotonic()
    new_records = []
    for i, oid in enumerate(obs_ids):
        obs = obs_by_id[oid]
        fp = WiseExactFootprint(msk_path(obs))
        qual = obs.quality_flags.get("qual_frame")
        quality_bad = qual is not None and qual == 0
        t_mid = Time(obs.t_mid_mjd_utc, format="mjd")
        for endpoint_id, role_name in hits_by_obs[oid]:
            role = ROLE_BY_NAME[role_name]
            target = registry[endpoint_id]
            al = adaptive_locus(
                target=target, role=role, observation_time=t_mid,
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=PRECISE_TOLERANCE_ARCSEC,
                ephemeris=ctx.ephemeris, model=ctx.model)
            pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                            for p in al.points])
            inb, usable = fp.status(pts)
            frac = (float(usable[inb].mean()) if inb.any() else None)
            covered = covered_z_intervals(
                target=target, role=role, observation_time=t_mid,
                observer=ctx.observer, relay_range=ctx.relay_range,
                contains=fp.contains,
                tolerance_arcsec=PRECISE_TOLERANCE_ARCSEC,
                ephemeris=ctx.ephemeris, model=ctx.model, locus=al,
                seed_step_arcsec=SEED_STEP_ARCSEC)
            hit = len(covered) > 0
            if quality_bad:
                label, note = "unusable", f"qual_frame={qual}"
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
                       "fatal_mask_bits": list(
                           WiseExactFootprint.FATAL_BITS),
                       "qual_frame": qual},
            )
            if ixn.intersection_id not in known_ixn:
                known_ixn.add(ixn.intersection_id)
                new_records.append(ixn)
            s = summary[(endpoint_id, role_name)][obs.band]
            s[0] += 1
            s[1] += hit
            s[2] += (label == "usable" or label == "partial")
            if frac is not None:
                s[3] += frac
        if (i + 1) % 100 == 0:
            rate = (i + 1) / (_time.monotonic() - t0)
            print(f"  {i + 1}/{len(obs_ids)} frames ({rate:.2f}/s)",
                  flush=True)
        if len(new_records) >= 500:
            append_records(ixn_path, new_records)
            new_records = []
    append_records(ixn_path, new_records)

    print("\n=== precise summary "
          "(evaluated / geometric-hit / usable, mean usable frac) ===")
    for (eid, role), bands in sorted(summary.items()):
        parts = []
        for b in sorted(bands):
            n, h, u, fsum = bands[b]
            parts.append(f"{b}: {n}/{h}/{u} f={fsum / max(n, 1):.2f}")
        print(f"{eid:16s} {role}: " + "  ".join(parts))


if __name__ == "__main__":
    main()
