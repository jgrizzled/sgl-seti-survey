"""Catalog screening for the ZTF pilot (plan §5 step 4, layer 1).

For every usable/partial precise-pass exposure, fetch that quadrant's
PSF-fit source catalog (``psfcat.fits``, the deepest single-epoch list)
and match every detection against the endpoint x role locus at that
exposure's mid-time. Detections within the screen radius become
ScreenMatch records (distance to track, implied relay distance,
calibrated magnitude). Also snapshots the DR24 objects table around each
corridor as recurrence context.

No catalog is treated as complete (plan §3.4).

Usage:
    uv run python surveys/ztf/scripts/catalog_screen.py [corridor ...]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import time as _time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.irsa_ztf import ZtfSciAdapter
from sglsurvey.geometry import GeometryContext
from sglsurvey.records import (Observation, ScreenMatch, append_records,
                               read_records)
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ztf_corridors import MEMBERS as CORRIDORS  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "ztf" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "ztf" / "precise_v1"
RUN_DIR = REPO / "runs" / "ztf" / "screen_v1"
PRODUCT_DIR = REPO / "runs" / "ztf" / "products" / "psfcat"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "ztf" / "hypotheses.md"
HYPOTHESIS_VERSION = "ztf-hypotheses-v1.0"

TAP_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"
SCREEN_RADIUS_ARCSEC = 10.0   # frozen padding (hypotheses v1.0 §5)
LOCUS_TOLERANCE_ARCSEC = 1.0
OBJECTS_TABLE = "ztf_objects_dr24"
OBJECTS_COLS = ["oid", "ra", "dec", "fid", "field", "ccdid", "qid",
                "nobs", "ngoodobs", "medianmag", "medmagerr",
                "minmag", "maxmag", "refmag", "astrometricrms"]
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("corridors", nargs="*", default=[])
    args = ap.parse_args()
    corridors = args.corridors or list(CORRIDORS)

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.ztf_default()
    adapter = ZtfSciAdapter()
    store = SnapshotStore(RUN_DIR)
    session = requests.Session()
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()

    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs_by_id[r["observation_id"]] = Observation(**r)
    usable = defaultdict(list)  # obs_id -> [(endpoint, role)]
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[r["observation_id"]].append((r["endpoint_id"], r["role"]))

    match_path = RUN_DIR / "records" / "screen_match.jsonl"
    known = ({r["screen_match_id"] for r in read_records(match_path)}
             if match_path.exists() else set())
    progress_path = RUN_DIR / "progress.json"
    progress = (set(json.loads(progress_path.read_text()))
                if progress_path.exists() else set())
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    run_config = {
        "run_id": "screen_v1", "stage": "catalog-screen",
        "source_runs": ["coarse_v1", "precise_v1"],
        "corridors": corridors,
        "screen_radius_arcsec": SCREEN_RADIUS_ARCSEC,
        "locus_tolerance_arcsec": LOCUS_TOLERANCE_ARCSEC,
        "source_catalog": "psfcat (per-exposure PSF-fit, ZSDS §10.6)",
        "context_table": OBJECTS_TABLE,
        "registry_source_hash": registry.source_hash,
        "hypothesis_version": HYPOTHESIS_VERSION,
        "hypothesis_hash": hyp_hash,
        "geometry": ctx.identities(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    (RUN_DIR / f"run_config_{run_config['started_utc'][:19]}.json"
     ).write_text(json.dumps(run_config, indent=2))

    summary = defaultdict(lambda: [0, 0, 0])  # matches, exposures, missing
    for corridor in corridors:
        members = set(CORRIDORS[corridor])
        obs_ids = sorted(
            (o for o, pairs in usable.items()
             if any(e in members for e, _ in pairs)),
            key=lambda o: obs_by_id[o].t_mid_mjd_utc)
        print(f"[{corridor}] {len(obs_ids)} usable exposures", flush=True)
        t0 = _time.monotonic()
        new = []
        for i, oid in enumerate(obs_ids):
            if f"exp:{oid}" in progress:
                continue
            obs = obs_by_id[oid]
            try:
                pset = adapter.fetch(obs, ["psfcat"], PRODUCT_DIR)
            except FileNotFoundError:
                for e, r in usable[oid]:
                    summary[(e, r)][2] += 1
                progress.add(f"exp:{oid}")
                continue
            with fits.open(pset.products[0].path) as hdul:
                cat = hdul[1].data
                magzp = float(hdul[0].header.get("MAGZP", np.nan))
            if len(cat) == 0:
                progress.add(f"exp:{oid}")
                continue
            det_ra = np.asarray(cat["ra"], dtype=float)
            det_dec = np.asarray(cat["dec"], dtype=float)
            t_mid = Time(obs.t_mid_mjd_utc, format="mjd")
            for e, r in usable[oid]:
                if e not in members:
                    continue
                al = adaptive_locus(
                    target=registry[e], role=Role(r), observation_time=t_mid,
                    observer=ctx.observer, relay_range=ctx.relay_range,
                    tolerance_arcsec=LOCUS_TOLERANCE_ARCSEC,
                    ephemeris=ctx.ephemeris, model=ctx.model)
                pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                                for p in al.points])
                zs = np.array([p.z_au for p in al.points])
                # cheap prefilter: bounding box + margin
                m = SCREEN_RADIUS_ARCSEC / 3600.0 * 1.5
                cosd = np.cos(np.deg2rad(pts[:, 1].mean()))
                sel = ((det_dec > pts[:, 1].min() - m)
                       & (det_dec < pts[:, 1].max() + m)
                       & (det_ra > pts[:, 0].min() - m / cosd)
                       & (det_ra < pts[:, 0].max() + m / cosd))
                if not sel.any():
                    summary[(e, r)][1] += 1
                    continue
                idx = np.flatnonzero(sel)
                d, seg = seg_min_dist_vec(det_ra[idx], det_dec[idx], pts)
                for k, dist, s in zip(idx, d, seg):
                    if dist > SCREEN_RADIUS_ARCSEC:
                        continue
                    row = cat[k]
                    sm = ScreenMatch.build(
                        endpoint_id=e, role=r,
                        hypothesis_version=HYPOTHESIS_VERSION,
                        registry_source_hash=registry.source_hash,
                        source_table="ztf-psfcat",
                        source_cntr=f"{obs.native_key['pid']}:{int(row['sourceid'])}",
                        source_designation=None,
                        ra_deg=float(row["ra"]), dec_deg=float(row["dec"]),
                        sigra_mas=None, sigdec_mas=None,
                        mjd=obs.t_mid_mjd_utc, scan_id=None, frame_num=None,
                        photometry={
                            "filter": obs.band,
                            "mag_inst": float(row["mag"]),
                            "mag_cal": float(row["mag"]) + magzp,
                            "sigmag": float(row["sigmag"]),
                            "snr": float(row["snr"]),
                            "magzp": magzp},
                        flags={"flags": int(row["flags"]),
                               "chi": float(row["chi"]),
                               "sharp": float(row["sharp"])},
                        dist_arcsec=round(float(dist), 3),
                        implied_z_au=round(float(zs[s]), 2),
                        z_segment_au=(round(float(zs[s]), 2),
                                      round(float(zs[s + 1]), 2)),
                        screen_radius_arcsec=SCREEN_RADIUS_ARCSEC,
                        locus_mjd=obs.t_mid_mjd_utc, snapshot_id=None,
                        extra={"observation_id": oid,
                               "psfcat_sha256": pset.products[0].checksum},
                    )
                    if sm.screen_match_id not in known:
                        known.add(sm.screen_match_id)
                        new.append(sm)
                    summary[(e, r)][0] += 1
                summary[(e, r)][1] += 1
            progress.add(f"exp:{oid}")
            if (i + 1) % 100 == 0:
                append_records(match_path, new)
                new = []
                progress_path.write_text(json.dumps(sorted(progress)))
                print(f"  {i + 1}/{len(obs_ids)} "
                      f"({(i + 1) / (_time.monotonic() - t0):.2f}/s)",
                      flush=True)
        append_records(match_path, new)
        progress_path.write_text(json.dumps(sorted(progress)))

        # corridor context: DR24 objects within 0.2 deg of corridor centre
        ckey = f"{corridor}:objects"
        if ckey not in progress and obs_ids:
            e0 = next(iter(members))
            al = adaptive_locus(
                target=registry[e0], role=Role.RX,
                observation_time=Time(60000.0, format="mjd"),
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=30.0, ephemeris=ctx.ephemeris,
                model=ctx.model)
            cra = float(np.mean([p.icrs_ra_deg for p in al.points]))
            cdec = float(np.mean([p.icrs_dec_deg for p in al.points]))
            query = (f"SELECT {', '.join(OBJECTS_COLS)} FROM {OBJECTS_TABLE} "
                     f"WHERE CONTAINS(POINT('ICRS',ra,dec),"
                     f"CIRCLE('ICRS',{cra:.6f},{cdec:.6f},0.2))=1")
            request_utc = datetime.now(timezone.utc).isoformat()
            try:
                resp = session.get(TAP_SYNC, params={"QUERY": query,
                                                     "FORMAT": "CSV"},
                                   timeout=900)
                resp.raise_for_status()
                rows = list(csv.DictReader(io.StringIO(resp.text)))
                store.store(service_url=TAP_SYNC, query=query,
                            request_utc=request_utc,
                            response_bytes=resp.content, row_count=len(rows),
                            http_status=resp.status_code)
                print(f"  {ckey}: {len(rows)} DR24 objects snapshotted",
                      flush=True)
                progress.add(ckey)
                progress_path.write_text(json.dumps(sorted(progress)))
            except Exception as exc:
                print(f"  {ckey}: FAILED {exc}", flush=True)

    print("\n=== screen summary ===")
    for (e, r), (n, nexp, miss) in sorted(summary.items()):
        print(f"{e:16s} {r}: {n} matches over {nexp} exposures "
              f"({miss} psfcat missing)")


if __name__ == "__main__":
    main()
