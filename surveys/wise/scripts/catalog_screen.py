"""Catalog screening for the WISE/NEOWISE shakedown (plan §4.5,
layer 1 of the detection strategy).

For every corridor x mission-phase x visit with usable precise-pass
geometry, query the phase-appropriate single-exposure source table in a
tight region around the union of the member endpoints' locus arcs
(snapshotting raw responses), then match each returned detection against
each endpoint x role locus polyline evaluated near the detection epoch.
Detections within the screen radius become ScreenMatch records carrying
distance-to-track and the implied relay distance.

Also pulls per-corridor context tables: CatWISE2020 + reject, and the
NEOWISE known-solar-system-object association list (veto input).

No catalog is treated as complete (plan §3.4): misses here are not null
results, and matches are inputs to candidate vetting, not detections.

Usage:
    uv run python surveys/wise/scripts/catalog_screen.py [corridor ...]
                  [--limit-visits N]
Corridors: barnard alphacen sirius ross154 lalande (default all).
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
import requests
from astropy.time import Time

from sglseti import Role, load_target_registry

from sglsurvey.geometry import GeometryContext, locus_radec
from sglsurvey.records import (Observation, ScreenMatch, append_records,
                               read_records)
from sglsurvey.snapshots import SnapshotStore

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "wise" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "wise" / "precise_v1"
RUN_DIR = REPO / "runs" / "wise" / "screen_v1"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "wise" / "hypotheses.md"
HYPOTHESIS_VERSION = "wise-hypotheses-v1.0"

TAP_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"
SCREEN_RADIUS_ARCSEC = 10.0   # frozen padding (hypotheses v1.0 §5)
LOCUS_TOLERANCE_ARCSEC = 2.0
LOCUS_BIN_DAYS = 0.5          # locus reuse granularity for matching
VISIT_GAP_DAYS = 5.0

from wise_corridors import MEMBERS as CORRIDORS

ROLES = (Role.RX, Role.TX)

# Mission-phase table selection by frame MJD.
PSD_COMMON = ["cntr", "source_id", "scan_id", "frame_num", "ra", "dec",
              "sigra", "sigdec", "mjd", "w1mpro", "w1sigmpro", "w1snr",
              "w2mpro", "w2sigmpro", "w2snr", "cc_flags", "ph_qual",
              "sso_flg", "qual_frame", "moon_masked", "det_bit", "nb",
              "na"]
W3 = ["w3mpro", "w3sigmpro", "w3snr"]
W4 = ["w4mpro", "w4sigmpro", "w4snr"]
PSD_PHASES = [
    (55203.0, 55414.5, "allsky_4band_p1bs_psd", PSD_COMMON + W3 + W4),
    (55414.5, 55469.5, "allsky_3band_p1bs_psd", PSD_COMMON + W3),
    (55469.5, 55650.0, "allsky_2band_p1bs_psd", PSD_COMMON),
    (56600.0, 61000.0, "neowiser_p1bs_psd", PSD_COMMON),
]

CATWISE_COLS = ["cntr", "source_name", "ra", "dec", "sigra", "sigdec",
                "ra_pm", "dec_pm", "pmra", "pmdec", "sigpmra", "sigpmdec",
                "par_pm", "par_pmsig", "meanobsmjd", "w1mpro", "w1sigmpro",
                "w1snr", "w2mpro", "w2sigmpro", "w2snr", "cc_flags",
                "ab_flags"]
MCH_COLS = ["cntr", "objid", "ra", "dec", "mjd", "scan_id", "frame_num",
            "dra", "ddec", "mconf", "w1mpro", "w2mpro"]

PHOT_KEYS = [c for c in PSD_COMMON + W3 + W4
             if c.startswith(("w1", "w2", "w3", "w4"))]
FLAG_KEYS = ["cc_flags", "ph_qual", "sso_flg", "qual_frame",
             "moon_masked", "det_bit", "nb", "na"]


def phase_of(mjd: float):
    for lo, hi, table, cols in PSD_PHASES:
        if lo <= mjd < hi:
            return table, cols
    return None, None


def gap_cluster(mjds: list[float], gap: float) -> list[list[float]]:
    out, cur = [], [mjds[0]]
    for m in mjds[1:]:
        if m - cur[-1] > gap:
            out.append(cur)
            cur = []
        cur.append(m)
    out.append(cur)
    return out


def tap_query(session, query: str, store: SnapshotStore):
    request_utc = datetime.now(timezone.utc).isoformat()
    resp = session.get(TAP_SYNC, params={"QUERY": query, "FORMAT": "CSV"},
                       timeout=900)
    resp.raise_for_status()
    import csv as _csv
    import io as _io
    rows = list(_csv.DictReader(_io.StringIO(resp.text)))
    snap = store.store(service_url=TAP_SYNC, query=query,
                       request_utc=request_utc,
                       response_bytes=resp.content, row_count=len(rows),
                       http_status=resp.status_code)
    return rows, snap


def seg_min_dist(det_ra, det_dec, pts):
    """Min angular distance (arcsec) from a point to a polyline, plus
    the index of the nearest segment. Local tangent-plane approx (valid
    at the sub-degree scales involved)."""
    cosd = np.cos(np.deg2rad(det_dec))
    x = (pts[:, 0] - det_ra) * cosd * 3600.0
    y = (pts[:, 1] - det_dec) * 3600.0
    ax, ay = x[:-1], y[:-1]
    bx, by = x[1:], y[1:]
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = np.clip(np.where(L2 > 0, -(ax * dx + ay * dy) / np.where(
        L2 > 0, L2, 1.0), 0.0), 0.0, 1.0)
    px, py = ax + t * dx, ay + t * dy
    d = np.hypot(px, py)
    i = int(np.argmin(d))
    return float(d[i]), i


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("corridors", nargs="*", default=[])
    ap.add_argument("--limit-visits", type=int, default=None)
    args = ap.parse_args()
    unknown = set(args.corridors) - set(CORRIDORS)
    if unknown:
        ap.error(f"unknown corridors: {sorted(unknown)} "
                 f"(choose from {sorted(CORRIDORS)})")
    corridors = args.corridors or list(CORRIDORS)

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.wise_coarse_default()
    store = SnapshotStore(RUN_DIR)
    session = requests.Session()
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()

    # Usable/partial precise intersections -> epochs per (endpoint, role).
    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        obs_by_id[r["observation_id"]] = r
    usable_mjds = defaultdict(set)
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            mjd = obs_by_id[r["observation_id"]]["t_mid_mjd_utc"]
            usable_mjds[(r["endpoint_id"], r["role"])].add(mjd)

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
        "locus_bin_days": LOCUS_BIN_DAYS,
        "visit_gap_days": VISIT_GAP_DAYS,
        "psd_phases": [[lo, hi, t] for lo, hi, t, _ in PSD_PHASES],
        "registry_source_hash": registry.source_hash,
        "hypothesis_version": HYPOTHESIS_VERSION,
        "hypothesis_hash": hyp_hash,
        "geometry": ctx.identities(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    (RUN_DIR / f"run_config_{run_config['started_utc'][:19]}.json"
     ).write_text(json.dumps(run_config, indent=2))

    locus_cache: dict = {}

    def polyline(endpoint_id, role_name, mjd):
        key = (endpoint_id, role_name, round(mjd / LOCUS_BIN_DAYS))
        if key not in locus_cache:
            t = Time(key[2] * LOCUS_BIN_DAYS, format="mjd")
            pts, _ = locus_radec(ctx, registry[endpoint_id],
                                 Role(role_name), t,
                                 tolerance_arcsec=LOCUS_TOLERANCE_ARCSEC)
            locus_cache[key] = (pts, float(t.mjd))
        return locus_cache[key]

    z_cache: dict = {}

    def polyline_z(endpoint_id, role_name, mjd):
        # z values ride along with the adaptive locus points.
        from sglseti import adaptive_locus
        key = (endpoint_id, role_name, round(mjd / LOCUS_BIN_DAYS))
        if key not in z_cache:
            t = Time(key[2] * LOCUS_BIN_DAYS, format="mjd")
            al = adaptive_locus(
                target=registry[endpoint_id], role=Role(role_name),
                observation_time=t, observer=ctx.observer,
                relay_range=ctx.relay_range,
                tolerance_arcsec=LOCUS_TOLERANCE_ARCSEC,
                ephemeris=ctx.ephemeris, model=ctx.model)
            pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                            for p in al.points])
            zs = np.array([p.z_au for p in al.points])
            z_cache[key] = (pts, zs, float(t.mjd))
        return z_cache[key]

    n_matches = 0
    summary = defaultdict(lambda: [0, set(), 0])  # matches, visits, sso

    for corridor in corridors:
        members = CORRIDORS[corridor]
        pairs = [(e, r.value) for e in members for r in ROLES]
        epochs = sorted(set().union(*(usable_mjds[p] for p in pairs)))
        if not epochs:
            print(f"[{corridor}] no usable epochs", flush=True)
            continue

        # Split by mission phase, then gap-cluster into visits.
        visits = []
        for lo, hi, table, cols in PSD_PHASES:
            in_phase = [m for m in epochs if lo <= m < hi]
            if in_phase:
                visits += [(table, cols, v)
                           for v in gap_cluster(in_phase, VISIT_GAP_DAYS)]
        if args.limit_visits:
            visits = visits[:args.limit_visits]
        print(f"[{corridor}] {len(epochs)} usable epochs, "
              f"{len(visits)} phase-visits", flush=True)

        for table, cols, visit in visits:
            vkey = f"{corridor}:{table}:{visit[0]:.3f}"
            if vkey in progress:
                continue
            v0, v1 = visit[0], visit[-1]
            # Query region: union of member locus arcs at visit
            # start/mid/end, plus screen + tolerance margin.
            pts_all = []
            for e, r in pairs:
                for m in (v0, (v0 + v1) / 2, v1):
                    pts_all.append(polyline(e, r, m)[0])
            pts_all = np.concatenate(pts_all)
            cra = float(np.mean(pts_all[:, 0]))
            cdec = float(np.mean(pts_all[:, 1]))
            cosd = np.cos(np.deg2rad(cdec))
            r_deg = float(np.max(np.hypot(
                (pts_all[:, 0] - cra) * cosd, pts_all[:, 1] - cdec)))
            radius = r_deg + (SCREEN_RADIUS_ARCSEC
                              + LOCUS_TOLERANCE_ARCSEC + 3.0) / 3600.0
            query = (
                f"SELECT {', '.join(cols)} FROM {table} WHERE "
                f"CONTAINS(POINT('ICRS',ra,dec),"
                f"CIRCLE('ICRS',{cra:.6f},{cdec:.6f},{radius:.6f}))=1 "
                f"AND mjd >= {v0 - 0.1:.4f} AND mjd <= {v1 + 0.1:.4f}")
            t0 = _time.monotonic()
            rows, snap = tap_query(session, query, store)
            new = []
            for row in rows:
                try:
                    dra, ddec = float(row["ra"]), float(row["dec"])
                    dmjd = float(row["mjd"])
                except (ValueError, KeyError):
                    continue
                for e, r in pairs:
                    pts, zs, lmjd = polyline_z(e, r, dmjd)
                    dist, seg = seg_min_dist(dra, ddec, pts)
                    if dist > SCREEN_RADIUS_ARCSEC:
                        continue
                    m = ScreenMatch.build(
                        endpoint_id=e, role=r,
                        hypothesis_version=HYPOTHESIS_VERSION,
                        registry_source_hash=registry.source_hash,
                        source_table=table, source_cntr=str(row["cntr"]),
                        source_designation=row.get("source_id") or None,
                        ra_deg=dra, dec_deg=ddec,
                        sigra_mas=(float(row["sigra"]) * 1000.0
                                   if row.get("sigra") else None),
                        sigdec_mas=(float(row["sigdec"]) * 1000.0
                                    if row.get("sigdec") else None),
                        mjd=dmjd, scan_id=row.get("scan_id"),
                        frame_num=(int(row["frame_num"])
                                   if row.get("frame_num") else None),
                        photometry={k: row.get(k) for k in PHOT_KEYS
                                    if row.get(k)},
                        flags={k: row.get(k) for k in FLAG_KEYS
                               if row.get(k) not in ("", None)},
                        dist_arcsec=round(dist, 3),
                        implied_z_au=round(float(zs[seg]), 2),
                        z_segment_au=(round(float(zs[seg]), 2),
                                      round(float(zs[seg + 1]), 2)),
                        screen_radius_arcsec=SCREEN_RADIUS_ARCSEC,
                        locus_mjd=lmjd, snapshot_id=snap.snapshot_id,
                    )
                    if m.screen_match_id not in known:
                        known.add(m.screen_match_id)
                        new.append(m)
                    s = summary[(e, r)]
                    s[0] += 1
                    s[1].add(vkey)
                    if (row.get("sso_flg") or "0") not in ("0", ""):
                        s[2] += 1
            append_records(match_path, new)
            n_matches += len(new)
            progress.add(vkey)
            progress_path.write_text(json.dumps(sorted(progress)))
            print(f"  {vkey}: {len(rows)} rows -> {len(new)} matches "
                  f"({_time.monotonic() - t0:.0f}s)", flush=True)

        # Corridor context pulls (once per corridor).
        for label, table, cols, radius in [
                ("catwise", "catwise_2020", CATWISE_COLS, 0.2),
                ("catwise-reject", "catwise_2020_reject", CATWISE_COLS,
                 0.2),
                ("sso-mch", "neowiser_p1ba_mch", MCH_COLS, 0.2)]:
            ckey = f"{corridor}:{label}"
            if ckey in progress:
                continue
            cra, cdec = (float(np.mean(pts_all[:, 0])),
                         float(np.mean(pts_all[:, 1])))
            query = (f"SELECT {', '.join(cols)} FROM {table} WHERE "
                     f"CONTAINS(POINT('ICRS',ra,dec),"
                     f"CIRCLE('ICRS',{cra:.6f},{cdec:.6f},{radius}))=1")
            try:
                rows, _ = tap_query(session, query, store)
                print(f"  {ckey}: {len(rows)} rows snapshotted",
                      flush=True)
                progress.add(ckey)
                progress_path.write_text(json.dumps(sorted(progress)))
            except Exception as exc:
                print(f"  {ckey}: FAILED {exc}", flush=True)

    print(f"\n=== screen summary ({n_matches} new matches) ===")
    for (e, r), (n, visits, sso) in sorted(summary.items()):
        print(f"{e:16s} {r}: {n} matches in {len(visits)} visits "
              f"({sso} sso-flagged)")


if __name__ == "__main__":
    main()
