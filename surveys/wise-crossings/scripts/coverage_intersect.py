"""WISE crossings coverage intersection v1 (hypotheses.md freeze v1.0).

Channel A, 1.0 AU rung only (every other channel/rung is geometrically
null for an elongation-90 observer — hypotheses v1.0). For each of the
54 elongation-pre-gate targets: one TAP discovery cone at the star
(snapshotted under runs/wise-crossings) over the wise_v1 era; then a
local intersection: an L1b frame covers an event if its t_mid lies in
t_ca +/- sqrt(r^2-b^2)/v_perp and its nominal footprint contains the
per-epoch propagated star position (linear PM fit through the era
events' tabulated positions). "Primary" = qual_frame > 0, qa_status A,
moon_sep >= 15 deg. The 34 targets with no elongation-viable event are
recorded as elongation-null without queries (the pointing law is the
archive's own). Also writes per-(target, band) primary epoch MJD lists
(epochs_v1.json) so the threshold freeze can evaluate pseudo-window
support offline. Coverage counting only - no pixels are touched.
"""

from __future__ import annotations

import json
import sys
import time as _time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.adapters.base import ConeRegion, MjdRange
from sglsurvey.adapters.irsa_wise import WiseMergeL1bAdapter, WiseNominalFootprint
from sglsurvey.snapshots import SnapshotStore

ERA = MjdRange(55203.0, 60524.0)
KM_PER_AU = 1.495978707e8
RUNG_AU = 1.0
BANDS = ("W1", "W2", "W3", "W4")
VISIT_BIN_DAYS = 0.5
CONE_MARGIN_DEG = 0.02
OUT = REPO / "surveys" / "wise-crossings" / "results"
RUN = REPO / "runs" / "wise-crossings"

# targets passing the elongation pre-gate (hypotheses v1.0 table)
GATE = json.loads((REPO / "surveys" / "wise-crossings" / "configs"
                   / "elongation_gate_v1.json").read_text())


def load_events():
    t = Table.read(REPO / "crossings" / "wise_v1" / "events.ecsv")
    mjd = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    t["t_ca_mjd"] = mjd
    side = np.asarray(t["axis_distance_au"])
    era = (mjd >= ERA.start_mjd_utc) & (mjd <= ERA.stop_mjd_utc)
    A = t[era & (t["link_direction"] == "inbound") & (side > 0)]
    return A[np.asarray(A["b_min_au"]) < RUNG_AU]


def star_track(sub):
    t = np.asarray(sub["t_ca_mjd"], float)
    ra = np.asarray(sub["star_icrs_ra_deg"], float)
    de = np.asarray(sub["star_icrs_dec_deg"], float)
    cosd = np.cos(np.radians(de.mean()))
    px = np.polyfit(t, ra * cosd, 1)
    py = np.polyfit(t, de, 1)
    return lambda m: (float(np.polyval(px, m)) / cosd,
                      float(np.polyval(py, m)))


def primary(o):
    # Amendment 2026-08-24 (pre-search, coverage stage): the frozen
    # "qa_status = A" value does not exist in the merge table (the
    # archive uses "Reviewed"); primary drops qa_status and keeps the
    # numeric gates. Distribution of qa_status recorded in the summary.
    q = o.quality_flags
    return ((q.get("qual_frame") or 0) > 0
            and (q.get("moon_sep") or 0) >= 15.0)


def load_snapshot_obs(adapter, store_dir):
    """Rebuild per-query Observation lists from the stored TAP
    snapshots (offline reparse; malformed rows skipped and counted)."""
    import csv as _csv
    import io as _io

    by_query = {}
    for line in (store_dir / "records"
                 / "query_snapshot.jsonl").read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            by_query[r["query"]] = r     # last write wins (dedup)
    out = {}
    skipped = 0
    for query, rec in by_query.items():
        raw = (store_dir / rec["response_path"]).read_bytes()
        rows = list(_csv.DictReader(_io.StringIO(raw.decode())))
        obs = []
        for row in rows:
            try:
                obs.append(adapter._observation(
                    row, rec["snapshot_id"], rec["request_utc"]))
            except Exception:
                skipped += 1
        out[query] = obs
    print(f"snapshot reparse: {len(out)} queries, {skipped} malformed "
          f"rows skipped", flush=True)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    A = load_events()
    gate_targets = sorted(GATE["viable_targets"])
    adapter = WiseMergeL1bAdapter()
    store = SnapshotStore(RUN)
    snapshot_obs = (load_snapshot_obs(adapter, RUN)
                    if "--from-snapshots" in sys.argv else None)

    partial_path = OUT / "coverage_partial_v1.jsonl"
    partial = {}
    if partial_path.exists():
        for line in partial_path.read_text().splitlines():
            if line.strip():
                d = json.loads(line)
                partial[d["target_id"]] = d
        print(f"resuming: {len(partial)} targets already done", flush=True)

    rows, epochs = [], {}
    qa_counter = Counter()
    for tid in sorted(set(str(x) for x in A["target_id"])):
        sub = A[np.asarray([str(x) == tid for x in A["target_id"]])]
        if tid not in gate_targets:
            for ev in sub:
                b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
                half = np.sqrt(RUNG_AU**2 - b*b) * KM_PER_AU / vp / 86400.0
                rows.append(_row(ev, tid, half, {}, {}, {}, [],
                                 status="elongation_null"))
            continue
        if tid in partial:
            d = partial[tid]
            epochs[tid] = d["epochs"]
            rows.extend(d["rows"])
            continue
        at = star_track(sub)
        ra0, de0 = at(float(np.mean(sub["t_ca_mjd"])))
        spread = 0.0
        for m in (ERA.start_mjd_utc, ERA.stop_mjd_utc):
            r_, d_ = at(m)
            spread = max(spread, np.hypot(
                (r_ - ra0) * np.cos(np.radians(de0)), d_ - de0))
        cone = ConeRegion(ra0, de0, spread + CONE_MARGIN_DEG)
        t0 = _time.monotonic()
        obs = None
        if snapshot_obs is not None:
            obs = snapshot_obs.get(adapter.discovery_query(cone, ERA))
            if obs is None:
                print(f"[{tid}] no stored snapshot for query", flush=True)
        else:
            for attempt in range(3):
                try:
                    obs = list(adapter.discover(cone, ERA, store))
                    break
                except Exception as exc:
                    print(f"[{tid}] discovery attempt {attempt + 1} "
                          f"failed: {exc}", flush=True)
                    _time.sleep(30)
        if obs is None:
            for ev in sub:
                b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
                half = np.sqrt(RUNG_AU**2 - b*b) * KM_PER_AU / vp / 86400.0
                rows.append(_row(ev, tid, half, {}, {}, {}, [],
                                 status="query_failed"))
            continue
        print(f"[{tid}] cone ({ra0:.4f},{de0:+.4f}) r={cone.radius_deg:.3f}"
              f" -> {len(obs)} frames ({_time.monotonic() - t0:.0f}s)",
              flush=True)
        # per-(band) primary epochs with the star in-footprint
        per_band = defaultdict(list)
        for o in obs:
            qa_counter[str(o.quality_flags.get("qa_status"))] += 1
            if not primary(o):
                continue
            r_, d_ = at(o.t_mid_mjd_utc)
            try:
                if WiseNominalFootprint(o).contains(r_, d_):
                    per_band[o.band].append(o.t_mid_mjd_utc)
            except Exception:
                continue
        epochs[tid] = {b: sorted(per_band[b]) for b in BANDS if per_band[b]}
        trows = []
        for ev in sub:
            b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
            half = np.sqrt(RUNG_AU**2 - b*b) * KM_PER_AU / vp / 86400.0
            lo, hi = ev["t_ca_mjd"] - half, ev["t_ca_mjd"] + half
            n_in, n_vis, mjds = {}, {}, []
            for bd, ts in epochs[tid].items():
                ts = np.asarray(ts)
                m = ts[(ts >= lo) & (ts <= hi)]
                n_in[bd] = int(len(m))
                n_vis[bd] = int(len(np.unique(np.round(m / VISIT_BIN_DAYS))))
                mjds.extend(m.tolist())
            trows.append(_row(ev, tid, half, n_in, n_vis,
                              {bd: len(ts)
                               for bd, ts in epochs[tid].items()},
                              mjds, status="queried"))
        rows.extend(trows)
        with partial_path.open("a") as fh:
            fh.write(json.dumps({"target_id": tid, "epochs": epochs[tid],
                                 "rows": trows}) + "\n")

    out = Table(rows=rows)
    out.write(OUT / "coverage_v1_events.ecsv", format="ascii.ecsv",
              overwrite=True)
    (OUT / "epochs_v1.json").write_text(json.dumps(epochs) + "\n")

    q = [r for r in rows if r["status"] == "queried"]
    summary = {
        "era_mjd": [ERA.start_mjd_utc, ERA.stop_mjd_utc],
        "rung_au": RUNG_AU, "n_events_total": len(rows),
        "n_events_elongation_null": sum(
            1 for r in rows if r["status"] == "elongation_null"),
        "n_events_queried": len(q),
        "n_covered_any_primary": sum(1 for r in q if r["n_primary"] > 0),
        "n_covered_ge2_visits": sum(
            1 for r in q if max((r[f"vis_{b}"] for b in BANDS),
                                default=0) >= 2),
        "per_band_covered": {b: sum(1 for r in q if r[f"n_{b}"] > 0)
                             for b in BANDS},
        "qa_status_distribution": dict(qa_counter),
        "quality_amendment": "primary drops qa_status (frozen value 'A' "
                             "does not exist in the merge table; archive "
                             "uses 'Reviewed'); numeric gates unchanged",
    }
    (OUT / "coverage_v1_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def _row(ev, tid, half, n_in, n_vis, n_all, mjds, status):
    r = {"event_id": str(ev["event_id"]), "target_id": tid,
         "t_ca_mjd": float(ev["t_ca_mjd"]),
         "b_min_au": float(ev["b_min_au"]),
         "v_perp_km_s": float(ev["v_perp_km_s"]),
         "radius_au": RUNG_AU, "window_days": 2 * half,
         "status": status}
    for b in BANDS:
        r[f"n_{b}"] = int(n_in.get(b, 0))
        r[f"vis_{b}"] = int(n_vis.get(b, 0))
        r[f"era_{b}"] = int(n_all.get(b, 0))
    r["n_primary"] = sum(n_in.values())
    r["mjd_first"] = min(mjds) if mjds else np.nan
    r["mjd_last"] = max(mjds) if mjds else np.nan
    return r


if __name__ == "__main__":
    main()
