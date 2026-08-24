"""Coverage refinement v1 (hypotheses.md freeze v1.0).

Replaces the gate's day-level sector-date arithmetic with real FFI
cadences: for every observed narrow-rung gate row, fetch the TESScut
cube (channel B: 31x31 px at the event antipode — track +/-2 px and
arcminute rings need the margin; channel A: 15x15 px at the
per-sector propagated star position) and count in-window /
off-window in-sector cadences under the primary quality mask
(QUALITY == 0). Times converted to spacecraft UTC via TIMECORR
(tesscut_lib). Windows from the spacecraft-frame events table. The
A 1.0 AU rung stays at gate-level coverage (constraint-only, zero
trials; no cubes fetched). No signal statistic is formed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from tesscut_lib import Cube, fetch_cube, load_cube  # noqa: E402

KM_PER_AU = 1.495978707e8
SIZE_B, SIZE_A = 31, 15
OUT = REPO / "surveys" / "tess-crossings" / "results"


def load_events():
    t = Table.read(REPO / "crossings" / "tess_v1" / "events.ecsv")
    t["mjd"] = Time(list(t["t_ca_utc"]), format="isot",
                    scale="utc").mjd
    return t


def star_track(t, tid):
    side = np.asarray(t["axis_distance_au"])
    A = t[(t["link_direction"] == "inbound") & (side > 0)]
    sub = A[np.asarray([str(x) == tid for x in A["target_id"]])]
    tt = np.asarray(sub["mjd"], float)
    ra = np.asarray(sub["star_icrs_ra_deg"], float)
    de = np.asarray(sub["star_icrs_dec_deg"], float)
    cosd = np.cos(np.radians(de.mean()))
    px = np.polyfit(tt, ra * cosd, 1)
    py = np.polyfit(tt, de, 1)
    return lambda m: (float(np.polyval(px, m)) / cosd,
                      float(np.polyval(py, m)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gate = Table.read(OUT / "coverage_gate_v1.ecsv")
    narrow = gate[(np.asarray(gate["radius_au"]) < 1.0)
                  & np.asarray(gate["observed"], bool)]
    ev_tab = load_events()
    ev_by_id = {str(e["event_id"]): e for e in ev_tab}

    cubes: dict[tuple, Cube] = {}
    rows = []
    for r in narrow:
        ch, tid = str(r["channel"]), str(r["target_id"])
        eid, sec = str(r["event_id"]), int(r["sector"])
        ev = ev_by_id[eid]
        if ch == "B":
            ra = float(ev["relay_icrs_ra_deg"])
            de = float(ev["relay_icrs_dec_deg"])
            key = ("B", eid, sec)
            size = SIZE_B
            label = f"B-{tid}-{eid[-6:]}"
        else:
            at = star_track(ev_tab, tid)
            ra, de = at(float(r["t_ca_mjd"]))
            key = ("A", tid, sec)
            size = SIZE_A
            label = f"A-{tid}"
        if key not in cubes:
            path = fetch_cube(label, ra, de, sec, size)
            cubes[key] = load_cube(path, label)
            c = cubes[key]
            print(f"[{label}] s{sec} cam{c.camera} ccd{c.ccd}: "
                  f"{len(c.mjd_utc)} cadences, "
                  f"Q0 {float(np.mean(c.quality == 0)):.2f}",
                  flush=True)
        c = cubes[key]
        # window from the gate row (already spacecraft-frame)
        hd = float(r["window_days"]) / 2.0
        lo, hi = float(r["t_ca_mjd"]) - hd, float(r["t_ca_mjd"]) + hd
        okq = (c.quality == 0) & np.isfinite(c.mjd_utc)
        inw = okq & (c.mjd_utc >= lo) & (c.mjd_utc <= hi)
        offw = okq & ((c.mjd_utc < lo) | (c.mjd_utc > hi))
        allw = np.isfinite(c.mjd_utc) & (c.mjd_utc >= lo) \
            & (c.mjd_utc <= hi)
        cad = float(np.median(np.diff(
            c.mjd_utc[np.isfinite(c.mjd_utc)]))) * 86400.0
        rows.append({
            "channel": ch, "target_id": tid, "event_id": eid,
            "radius_au": float(r["radius_au"]),
            "b_rsun": float(r["b_rsun"]),
            "t_ca_mjd": float(r["t_ca_mjd"]),
            "window_days": float(r["window_days"]), "sector": sec,
            "cadence_s": round(cad, 1),
            "n_inwindow_primary": int(inw.sum()),
            "n_inwindow_all": int(allw.sum()),
            "n_offwindow_primary": int(offw.sum()),
            "inwindow_primary_fraction": round(
                float(inw.sum() / max(allw.sum(), 1)), 3),
            "full_window": bool(lo >= np.nanmin(c.mjd_utc)
                                and hi <= np.nanmax(c.mjd_utc)),
        })

    tab = Table(rows=rows)
    tab.write(OUT / "coverage_refined_v1.ecsv", format="ascii.ecsv",
              overwrite=True)
    summary = {}
    for ch in ("A", "B"):
        for rad in sorted(set(float(x["radius_au"]) for x in rows
                              if x["channel"] == ch)):
            sel = [x for x in rows if x["channel"] == ch
                   and x["radius_au"] == rad]
            summary[f"{ch}_{rad}"] = {
                "rows": len(sel),
                "with_inwindow_primary": sum(
                    1 for x in sel if x["n_inwindow_primary"] > 0),
                "inwindow_primary_cadences": [
                    x["n_inwindow_primary"] for x in sel],
                "offwindow_primary_min": min(
                    (x["n_offwindow_primary"] for x in sel),
                    default=0),
            }
    (OUT / "coverage_refined_v1_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
