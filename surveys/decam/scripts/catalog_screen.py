"""Catalog screening for the DECam pilot (plan §3.4 layer 1) — port of
surveys/panstarrs/scripts/catalog_screen.py to NSC DR2 via the
`astro-datalab` query client (anonymous sync; recon addendum
2026-08-24).

NSC DR2 gotchas this script encodes:

- a direct q3c cone on ``nsc_dr2.meas`` (~34G rows) times out, so the
  screen goes object-cone -> track-proximity filter -> ``meas`` by
  ``objectid IN (...)`` in chunks;
- NSC DR2 is **time-partial** (ingest ends ~2017-2019 in the probed
  fields) — exposures after the catalog's last epoch legitimately have
  no catalogued detections, and catalog absence is never a null
  (plan §3.4); the per-corridor meas MJD range is recorded so the
  later stages can tell "not catalogued" from "not detected".

Per corridor: one object-cone query (positions, PM, per-band ndet,
mags; snapshotted CSV), objects within the union locus band feed a
chunked meas query (snapshotted); meas rows are matched to usable
precise-pass exposures by (filter, |mjd - t_mid| < 60 s) and to the
endpoint x role locus at mid-exposure. Matches within the screen
radius become ScreenMatch records.

Usage:
    uv run python surveys/decam/scripts/catalog_screen.py [corridor ...]
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
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.geometry import GeometryContext, discovery_cone
from sglsurvey.records import (Observation, ScreenMatch, append_records,
                               read_records)
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from decam_corridors import MEMBERS as ALL_MEMBERS  # noqa: E402
from decam_corridors import PILOT_CORRIDORS  # noqa: E402

ROLE_BY_NAME = {r.value: r for r in Role}

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "decam" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "decam" / "precise_v1"
RUN_DIR = REPO / "runs" / "decam" / "screen_v1"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "decam" / "hypotheses.md"
HYPOTHESIS_VERSION = "decam-hypotheses-v1.0"

SCREEN_RADIUS_ARCSEC = 10.0
LOCUS_TOLERANCE_ARCSEC = 1.0
#: band half-width for the object-level track-proximity prefilter:
#: screen radius + NSC mean-position scatter + high-PM allowance.
OBJECT_BAND_ARCSEC = 20.0
EXPOSURE_MATCH_DAYS = 60.0 / 86400.0
TIME_RANGE_MJD = (56000.0, 62000.0)
ROLES = (Role.RX, Role.TX)
QUERY_TIMEOUT_S = 170  # Data Lab sync ceiling is 300 s
ID_CHUNK = 300

OBJECT_COLS = ("id, ra, dec, raerr, decerr, pmra, pmdec, mjd, deltamjd, "
               "ndet, nphot, gmag, rmag, imag, zmag, ymag, umag, vrmag, "
               "class_star, fwhm, flags, variable10sig")
MEAS_COLS = ("measid, objectid, exposure, ccdnum, filter, mjd, ra, raerr, "
             "dec, decerr, mag_auto, magerr_auto, fwhm, class_star, flags")


def dl_query(sql: str, store: SnapshotStore, timeout: int = QUERY_TIMEOUT_S
             ) -> list[dict]:
    """Sync query via the astro-datalab client, CSV snapshotted
    verbatim; returns rows as dicts."""
    from dl import queryClient as qc

    request_utc = datetime.now(timezone.utc).isoformat()
    text = qc.query(sql=sql, fmt="csv", timeout=timeout)
    rows = list(csv.DictReader(io.StringIO(text)))
    store.store(service_url="datalab:queryClient/query", query=sql,
                request_utc=request_utc, response_bytes=text.encode(),
                row_count=len(rows), http_status=200)
    return rows


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
        v = float(row[key])
        return default if v > 99.98 and key.endswith("mag") else v
    except (KeyError, TypeError, ValueError):
        return default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("corridors", nargs="*", default=[])
    args = ap.parse_args()
    corridors = args.corridors or list(PILOT_CORRIDORS)

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.decam_default()
    store = SnapshotStore(RUN_DIR)
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()

    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs_by_id[r["observation_id"]] = Observation(**r)
    # usable exposures per (endpoint, role): {(band, t_mid): [obs_ids]}
    usable = defaultdict(lambda: defaultdict(list))
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            o = obs_by_id[r["observation_id"]]
            usable[(r["endpoint_id"], r["role"])][
                (o.band, round(o.t_mid_mjd_utc, 6))].append(
                    o.observation_id)

    match_path = RUN_DIR / "records" / "screen_match.jsonl"
    known = ({r["screen_match_id"] for r in read_records(match_path)}
             if match_path.exists() else set())
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    run_config = {
        "run_id": "screen_v1", "stage": "catalog-screen",
        "source_runs": ["coarse_v1", "precise_v1"],
        "corridors": corridors,
        "screen_radius_arcsec": SCREEN_RADIUS_ARCSEC,
        "object_band_arcsec": OBJECT_BAND_ARCSEC,
        "locus_tolerance_arcsec": LOCUS_TOLERANCE_ARCSEC,
        "exposure_match_seconds": EXPOSURE_MATCH_DAYS * 86400.0,
        "source_catalog": ("NSC DR2 meas via Data Lab queryClient "
                           "(object-cone -> objectid chunks); "
                           "TIME-PARTIAL vs the instcal archive"),
        "context_table": "nsc_dr2.object",
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
        members = [e for e in ALL_MEMBERS[corridor]]
        cones = [discovery_cone(ctx, registry[e], ROLES, *TIME_RANGE_MJD)
                 for e in members]
        cra = float(np.mean([c.ra_deg for c in cones]))
        cdec = float(np.mean([c.dec_deg for c in cones]))
        sep = max(np.hypot((c.ra_deg - cra) * np.cos(np.deg2rad(cdec)),
                           c.dec_deg - cdec) + c.radius_deg for c in cones)
        radius = float(sep + OBJECT_BAND_ARCSEC / 3600.0)

        t0 = _time.monotonic()
        objs = dl_query(
            f"SELECT {OBJECT_COLS} FROM nsc_dr2.object "
            f"WHERE q3c_radial_query(ra,dec,{cra:.6f},{cdec:.6f},"
            f"{radius:.5f})", store)
        obj_ra = np.array([_f(o, "ra") for o in objs])
        obj_dec = np.array([_f(o, "dec") for o in objs])
        print(f"[{corridor}] cone ({cra:.4f},{cdec:.4f}) r={radius:.3f}: "
              f"{len(objs)} NSC objects ({_time.monotonic() - t0:.0f}s)",
              flush=True)

        # union locus band: min distance of every object to every
        # usable epoch's locus (coarse polylines, deduped by month)
        keep = np.zeros(len(objs), dtype=bool)
        month_epochs = {}
        for e in members:
            for role in ROLES:
                for (band, t_mid), _oids in usable.get(
                        (e, role.value), {}).items():
                    month_epochs[(e, role.value,
                                  round(t_mid / 30.0))] = t_mid
        for (e, role_name, _m), t_mid in sorted(month_epochs.items()):
            al = adaptive_locus(
                target=registry[e], role=ROLE_BY_NAME[role_name],
                observation_time=Time(t_mid, format="mjd"),
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=5.0,
                ephemeris=ctx.ephemeris, model=ctx.model)
            pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                            for p in al.points])
            todo = np.flatnonzero(~keep)
            if len(todo) == 0:
                break
            d, _ = seg_min_dist_vec(obj_ra[todo], obj_dec[todo], pts)
            keep[todo[d <= OBJECT_BAND_ARCSEC]] = True
        track_objs = [o for k, o in zip(keep, objs) if k]
        print(f"[{corridor}] {len(track_objs)} objects within "
              f"{OBJECT_BAND_ARCSEC}\" of the union locus band "
              f"({len(month_epochs)} monthly loci)", flush=True)

        # meas for track objects, chunked
        meas = []
        ids = [o["id"] for o in track_objs]
        t0 = _time.monotonic()
        for c0 in range(0, len(ids), ID_CHUNK):
            idlist = ",".join(f"'{i}'" for i in ids[c0:c0 + ID_CHUNK])
            meas.extend(dl_query(
                f"SELECT {MEAS_COLS} FROM nsc_dr2.meas "
                f"WHERE objectid IN ({idlist})", store))
        mjds = [(m["mjd"]) for m in meas if m.get("mjd")]
        catalog_stats[corridor] = {
            "center": [cra, cdec], "radius_deg": radius,
            "n_objects_cone": len(objs),
            "n_objects_track": len(track_objs),
            "n_meas_track": len(meas),
            "meas_mjd_range": ([float(min(mjds)), float(max(mjds))]
                               if mjds else None)}
        print(f"[{corridor}] {len(meas)} meas rows for track objects "
              f"({_time.monotonic() - t0:.0f}s); catalogued epochs "
              f"{catalog_stats[corridor]['meas_mjd_range']}", flush=True)

        det_t = np.array([_f(m, "mjd") for m in meas])
        det_band = np.array([(m.get("filter") or "?") for m in meas])
        det_ra = np.array([_f(m, "ra") for m in meas])
        det_dec = np.array([_f(m, "dec") for m in meas])

        obj_by_id = {o["id"]: o for o in track_objs}
        new = []
        for e in members:
            for role in ROLES:
                exposures = usable.get((e, role.value), {})
                for (band, t_mid), oids in sorted(exposures.items(),
                                                  key=lambda kv: kv[0][1]):
                    sel = np.flatnonzero(
                        (np.abs(det_t - t_mid) < EXPOSURE_MATCH_DAYS)
                        & (det_band == band))
                    summary[(e, role.value)][1] += 1
                    if len(sel) == 0:
                        summary[(e, role.value)][2] += 1
                        continue
                    al = adaptive_locus(
                        target=registry[e], role=role,
                        observation_time=Time(t_mid, format="mjd"),
                        observer=ctx.observer, relay_range=ctx.relay_range,
                        tolerance_arcsec=LOCUS_TOLERANCE_ARCSEC,
                        ephemeris=ctx.ephemeris, model=ctx.model)
                    pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                                    for p in al.points])
                    zs = np.array([p.z_au for p in al.points])
                    d, seg = seg_min_dist_vec(det_ra[sel], det_dec[sel],
                                              pts)
                    for k, dist, s in zip(sel, d, seg):
                        if dist > SCREEN_RADIUS_ARCSEC:
                            continue
                        row = meas[k]
                        obj = obj_by_id.get(row["objectid"], {})
                        sm = ScreenMatch.build(
                            endpoint_id=e, role=role.value,
                            hypothesis_version=HYPOTHESIS_VERSION,
                            registry_source_hash=registry.source_hash,
                            source_table="nsc-dr2-meas",
                            source_cntr=str(row["measid"]),
                            source_designation=str(row["objectid"]),
                            ra_deg=float(det_ra[k]),
                            dec_deg=float(det_dec[k]),
                            sigra_mas=_f(row, "raerr") * 1000.0,
                            sigdec_mas=_f(row, "decerr") * 1000.0,
                            mjd=float(det_t[k]), scan_id=None,
                            frame_num=None,
                            photometry={
                                "filter": band,
                                "mag_auto_ab": _f(row, "mag_auto"),
                                "magerr_auto": _f(row, "magerr_auto"),
                                "fwhm_arcsec": _f(row, "fwhm"),
                                "class_star": _f(row, "class_star")},
                            flags={"meas_flags": row.get("flags"),
                                   "exposure": row.get("exposure"),
                                   "ccdnum": row.get("ccdnum"),
                                   "object_ndet": obj.get("ndet"),
                                   "object_pmra": obj.get("pmra"),
                                   "object_pmdec": obj.get("pmdec"),
                                   "object_variable10sig":
                                       obj.get("variable10sig")},
                            dist_arcsec=round(float(dist), 3),
                            implied_z_au=round(float(zs[s]), 2),
                            z_segment_au=(round(float(zs[s]), 2),
                                          round(float(zs[s + 1]), 2)),
                            screen_radius_arcsec=SCREEN_RADIUS_ARCSEC,
                            locus_mjd=t_mid, snapshot_id=None,
                            extra={"observation_ids": oids,
                                   "exposure_t_mid": t_mid},
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
