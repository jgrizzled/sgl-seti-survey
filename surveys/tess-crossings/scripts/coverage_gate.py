"""TESS crossings coverage gate (pre-freeze, metadata only).

The WISE-style kill-switch, run before any survey design: intersect
the TESS-spacecraft crossing list (crossings/tess_v1 — mandatory:
Earth-center misplaces TESS-frame grazing events by up to 0.33 R_sun
and t_ca by ~3 h) with actual TESS sector coverage.

Inputs: TESScut sector lookup (one cone per channel x target,
snapshotted under runs/tess-crossings) and the HEASARC sector table
(id, pointing type incl. 'Ecliptic', start/end dates; snapshotted).
For each event x ladder rung: window [t_ca +/- sqrt(r^2-b^2)/v_perp]
vs the date spans of the sectors covering the source position
(channel A: star; channel B: antipode/relay). Day-level precision —
the survey stage would refine with FFI timestamps and the ~1 d
mid-sector downlink gap. Cadence class per sector: 30 min (s1-26),
10 min (s27-55), 200 s (s56+). No pixels are touched.
"""

from __future__ import annotations

import json
import re
import sys
import time as _time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.snapshots import SnapshotStore

TESSCUT = "https://mast.stsci.edu/tesscut/api/v0.1/sector"
HEASARC = "https://heasarc.gsfc.nasa.gov/docs/tess/sector.html"
KM_PER_AU = 1.495978707e8
LADDER = {"A": (0.1, 1.0), "B": (0.0056, 0.0116, 0.1)}
NOW_MJD = 61276.0                     # 2026-08-24
OUT = REPO / "surveys" / "tess-crossings" / "results"
RUN = REPO / "runs" / "tess-crossings"


def cadence_s(sector):
    return 1800 if sector <= 26 else 600 if sector <= 55 else 200


def sector_dates(session, store):
    r = session.get(HEASARC, timeout=60)
    r.raise_for_status()
    store.store(service_url=HEASARC, query="sector.html",
                request_utc=Time.now().isot, response_bytes=r.content,
                row_count=0, http_status=r.status_code)
    out = {}
    for row in re.findall(r"<tr>(.*?)</tr>", r.text, re.S):
        cells = [re.sub(r"<[^>]+>", " ", c).replace("&nbsp;", " ").strip()
                 for c in re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)]
        if len(cells) < 6 or not cells[0].startswith("s"):
            continue
        try:
            sec = int(cells[0][1:])
            t0 = Time.strptime(cells[4], "%Y %b %d").mjd
            t1 = Time.strptime(cells[5], "%Y %b %d").mjd + 1.0
        except Exception:
            continue
        # source-table sanity: a TESS sector is 24-32 d. The HEASARC
        # page has at least one end-date year typo (s046 "2022 Dec 30"
        # for 2021): clamp implausible durations to start + 28 d and
        # flag the row.
        fixed = False
        if not (20.0 <= t1 - t0 <= 40.0):
            t1 = t0 + 28.0
            fixed = True
            print(f"  ! sector {sec}: implausible span "
                  f"{cells[4]} -> {cells[5]}; clamped to start+28d",
                  flush=True)
        out[sec] = {"kind": cells[1], "start": float(t0),
                    "end": float(t1), "date_clamped": fixed}
    return out


def tesscut_sectors(session, store, ra, dec):
    params = {"ra": f"{ra:.5f}", "dec": f"{dec:.5f}"}
    for attempt in range(3):
        try:
            r = session.get(TESSCUT, params=params, timeout=120)
            r.raise_for_status()
            doc = r.json()
            break
        except Exception as exc:
            print(f"  tesscut attempt {attempt+1}: {exc}", flush=True)
            _time.sleep(15)
    else:
        return None
    store.store(service_url=TESSCUT,
                query=f"ra={params['ra']}&dec={params['dec']}",
                request_utc=Time.now().isot, response_bytes=r.content,
                row_count=len(doc.get("results", [])),
                http_status=r.status_code)
    return sorted({int(x["sector"]) for x in doc.get("results", [])})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    RUN.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    store = SnapshotStore(RUN)
    dates = sector_dates(session, store)
    print(f"{len(dates)} sectors with dates "
          f"(max observed-to-date: "
          f"{max(s for s, d in dates.items() if d['end'] <= NOW_MJD)})")

    t = Table.read(REPO / "crossings" / "tess_v1" / "events.ecsv")
    t["mjd"] = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    side = np.asarray(t["axis_distance_au"])
    chans = {"A": t[(t["link_direction"] == "inbound") & (side > 0)],
             "B": t[(t["link_direction"] == "outbound") & (side < 0)]}

    rows = []
    for ch, tab in chans.items():
        tab = tab[np.asarray(tab["b_min_au"]) < max(LADDER[ch])]
        for tid in sorted(set(str(x) for x in tab["target_id"])):
            sub = tab[np.asarray([str(x) == tid
                                  for x in tab["target_id"]])]
            if ch == "A":
                ras = np.asarray(sub["star_icrs_ra_deg"], float)
                des = np.asarray(sub["star_icrs_dec_deg"], float)
            else:
                ras = np.asarray(sub["relay_icrs_ra_deg"], float)
                des = np.asarray(sub["relay_icrs_dec_deg"], float)
            ra0, de0 = float(np.median(ras)), float(np.median(des))
            secs = tesscut_sectors(session, store, ra0, de0)
            if secs is None:
                print(f"[{ch}] {tid}: QUERY FAILED", flush=True)
                continue
            print(f"[{ch}] {tid}: sectors {secs}", flush=True)
            for ev in sub:
                b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
                for r in LADDER[ch]:
                    if b >= r:
                        continue
                    hd = (np.sqrt(r * r - b * b) * KM_PER_AU
                          / vp / 86400.0)
                    lo, hi = float(ev["mjd"]) - hd, float(ev["mjd"]) + hd
                    for sec in secs:
                        d = dates.get(sec)
                        if d is None:
                            continue
                        ov = min(hi, d["end"]) - max(lo, d["start"])
                        if ov <= 0:
                            continue
                        rows.append({
                            "channel": ch, "target_id": tid,
                            "event_id": str(ev["event_id"]),
                            "radius_au": r,
                            "t_ca_mjd": float(ev["mjd"]),
                            "window_days": 2 * hd, "sector": sec,
                            "sector_kind": d["kind"],
                            "overlap_days": round(float(ov), 2),
                            "full_window": bool(lo >= d["start"]
                                                and hi <= d["end"]),
                            "cadence_s": cadence_s(sec),
                            "observed": bool(d["end"] <= NOW_MJD),
                            "b_rsun": round(b / 0.00465047, 3)})
    if rows:
        Table(rows=rows).write(OUT / "coverage_gate_v1.ecsv",
                               format="ascii.ecsv", overwrite=True)

    summary = {}
    for ch in LADDER:
        for r in LADDER[ch]:
            sel = [x for x in rows if x["channel"] == ch
                   and x["radius_au"] == r and x["observed"]]
            summary[f"{ch}_{r}"] = {
                "rows": len(sel),
                "events": len(set(x["event_id"] for x in sel)),
                "targets": len(set(x["target_id"] for x in sel)),
                "full_window_events": len(set(
                    x["event_id"] for x in sel if x["full_window"])),
            }
    sched = [x for x in rows if not x["observed"]]
    summary["scheduled_future_rows"] = len(sched)
    (OUT / "coverage_gate_v1_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
