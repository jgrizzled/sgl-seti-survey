"""WISPR coverage intersect: for every 0.1 AU sunward event-rung with a
WISPR-I arc (psp_census_v1 rows), count the real synoptic WISPR-I (and
-O) frames inside the in-beam in-field arc and inside the source's
off-beam field transit (the star-fixed baseline), from the NRL L2
inventory snapshot. Also the recon data-contact check.

Output: results/coverage_v1.json
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psp_geometry as G

REPO = Path(__file__).resolve().parents[3]
CENSUS = REPO / "surveys" / "wispr-crossings" / "results" / "psp_census_v1.json"
INV = REPO / "runs" / "wispr-crossings" / "coverage" / "l2_inventory.json"
OUT = REPO / "surveys" / "wispr-crossings" / "results" / "coverage_v1.json"
RECON_FRAMES = ["20241221T233014_1211", "20241224T000015_1211", "20241224T000205_2222", "20241227T003010_1211"]
STEP = 10 / 1440.0


def stamp_mjd(s):
    return Time(f"{s[:4]}-{s[4:6]}-{s[6:8]}T{s[9:11]}:{s[11:13]}:{s[13:15]}").mjd


def main():
    census = json.load(open(CENSUS))
    inv = json.load(open(INV))
    allf = [f for v in inv["days"].values() if v for f in v]
    frames = {}
    for tel, w in (("I", "1"), ("O", "2")):
        sel = [f for f in allf if f["wxyz"][0] == w and f["wxyz"][3] in "12"]
        m = np.array([stamp_mjd(f["stamp"]) for f in sel])
        o = np.argsort(m)
        frames[tel] = (m[o], np.array([f["file"] for f in sel])[o], np.array([f["wxyz"] for f in sel])[o])
    enc = G.encounters()
    events = json.load(open(REPO / "crossings" / "psp_v1" / "events.ecsv".replace("events.ecsv", "summary.json")))  # sanity only
    from astropy.table import Table
    ev_tab = Table.read(REPO / "crossings" / "psp_v1" / "events.ecsv")
    star = {str(r["event_id"]): (float(r["star_icrs_ra_deg"]), float(r["star_icrs_dec_deg"])) for r in ev_tab}
    rows_out, contact = [], []
    for r in census["visibility"]["rows"]:
        if r["rung"] != "0.1AU" or (r["visible_hours_wispr_i"] == 0 and r["visible_hours_wispr_o"] == 0):
            continue
        e = {"star_icrs_ra_deg": star[r["event_id"]][0], "star_icrs_dec_deg": star[r["event_id"]][1]}
        m_ca = r["mjd_ca"]
        mj = m_ca + np.arange(-8, 8 + 1e-9, STEP)
        b_e, along, rh = G.axis_geometry(e, mj)
        s_hat = G.source_direction(e, r["channel"])
        eps, lam, bet = G.wispr_coords(s_hat, mj)
        inbeam = b_e <= 0.1
        out = {"event_id": r["event_id"], "target_id": r["target_id"], "channel": r["channel"],
               "encounter": r["encounter"], "t_ca_utc": r["t_ca_utc"], "mjd_ca": m_ca}
        for tel, cam in (("I", G.WISPR_I), ("O", G.WISPR_O)):
            fov = G.in_wispr(lam, bet, cam)
            arc = fov & inbeam
            base = fov & ~inbeam
            fm, ff, fc = frames[tel]
            def count(mask):
                if not mask.any():
                    return 0, None, None, []
                # frames whose stamp falls within STEP/2 of a masked sample
                idx = np.searchsorted(mj, fm)
                ok = np.zeros(len(fm), bool)
                for k in np.nonzero(mask)[0]:
                    lo, hi = np.searchsorted(fm, mj[k] - STEP / 2), np.searchsorted(fm, mj[k] + STEP / 2)
                    ok[lo:hi] = True
                sel = np.nonzero(ok)[0]
                return int(ok.sum()), (float(fm[sel].min() - m_ca) * 24 if ok.any() else None), \
                    (float(fm[sel].max() - m_ca) * 24 if ok.any() else None), ff[sel].tolist()
            n_arc, a0, a1, arc_files = count(arc)
            n_base, b0, b1, base_files = count(base)
            out[f"wispr_{tel.lower()}"] = {
                "arc_geom_hours": round(float(arc.sum()) * STEP * 24, 2),
                "arc_frames": n_arc, "arc_frames_rel_tca_h": [None if a0 is None else round(a0, 2), None if a1 is None else round(a1, 2)],
                "field_transit_geom_hours": round(float(fov.sum()) * STEP * 24, 2),
                "baseline_frames": n_base, "baseline_rel_tca_h": [None if b0 is None else round(b0, 2), None if b1 is None else round(b1, 2)],
                "arc_b_e_rsun": [round(float(b_e[arc].min()) / G.RSUN_AU, 2), round(float(b_e[arc].max()) / G.RSUN_AU, 2)] if arc.any() else None,
                "arc_r_helio_au": [round(float(rh[arc].min()), 4), round(float(rh[arc].max()), 4)] if arc.any() else None,
                "arc_eps_deg": [round(float(eps[arc].min()), 1), round(float(eps[arc].max()), 1)] if arc.any() else None}
            for rf in RECON_FRAMES:
                if any(rf.split("_")[0] in f and rf.split("_")[1] in f for f in arc_files):
                    contact.append({"recon_frame": rf, "event_id": r["event_id"], "target_id": r["target_id"],
                                    "channel": r["channel"], "telescope": tel, "role": "arc"})
                elif any(rf.split("_")[0] in f and rf.split("_")[1] in f for f in base_files):
                    contact.append({"recon_frame": rf, "event_id": r["event_id"], "target_id": r["target_id"],
                                    "channel": r["channel"], "telescope": tel, "role": "baseline"})
        rows_out.append(out)
    # per unit summary
    units = {}
    for o in rows_out:
        u = units.setdefault(f"{o['target_id']}:{o['channel']}", {"events_with_geom_arc": 0, "covered_ge10": 0,
                                                                   "arc_frames": [], "baseline_frames": [], "encounters": []})
        w = o["wispr_i"]
        if w["arc_geom_hours"] > 0:
            u["events_with_geom_arc"] += 1
            if w["arc_frames"] >= 10:
                u["covered_ge10"] += 1
                u["arc_frames"].append(w["arc_frames"])
                u["baseline_frames"].append(w["baseline_frames"])
                u["encounters"].append(o["encounter"])
    for u in units.values():
        u["median_arc_frames"] = float(np.median(u["arc_frames"])) if u["arc_frames"] else 0
        u["min_arc_frames"] = int(min(u["arc_frames"])) if u["arc_frames"] else 0
        u["median_baseline_frames"] = float(np.median(u["baseline_frames"])) if u["baseline_frames"] else 0
        u["total_arc_frames"] = int(sum(u["arc_frames"]))
        del u["arc_frames"], u["baseline_frames"]
    res = {"inventory_snapshot": inv["snapshot_utc"], "n_synoptic_frames": {k: int(len(v[0])) for k, v in frames.items()},
           "recon_contact": contact, "units": units, "rows": rows_out}
    OUT.write_text(json.dumps(res, indent=1) + "\n")
    print("recon contact:", json.dumps(contact, indent=None))
    for k, u in sorted(units.items(), key=lambda kv: -kv[1]["covered_ge10"]):
        print(k, u)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
