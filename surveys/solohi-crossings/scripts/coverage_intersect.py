"""SoloHI coverage intersect: for every 0.1 AU sunward event-rung with an
inner-tile arc (solo_census_v1 rows), count the real inner-tile L2
frames of the source's nominal tile inside the in-beam in-tile arc and
inside the same-tile off-beam transit before the arc (the star-fixed
baseline), from the SOAR inventory snapshot; plus the recon
data-contact check (48 recon frames + 125 sampled headers).

Output: results/coverage_v1.json
"""
from __future__ import annotations
import glob, json, sys
from pathlib import Path
import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solo_geometry as G

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "solohi-crossings"
CENSUS = SURV / "results" / "solo_census_v1.json"
INV = REPO / "runs" / "solohi-crossings" / "coverage" / "soar_l2_inventory.json"
OUT = SURV / "results" / "coverage_v1.json"
RECON_DIRS = [REPO / "runs" / "solohi-crossings" / "recon", REPO / "runs" / "solohi-crossings" / "recon" / "sample"]
STEP = 10 / 1440.0
PRE_D, POST_D = 12.0, 1.0
MIN_ARC = 10


def main():
    census = json.load(open(CENSUS))
    inv = json.load(open(INV))
    ci = {c: i for i, c in enumerate(inv["columns"])}
    frames = {}
    for tile in "12":
        sel = [r for r in inv["rows"] if r[ci["descriptor"]] == f"solohi-{tile}ft"]
        m = np.array([Time(r[ci["begin_time"]]).mjd for r in sel])
        o = np.argsort(m)
        frames[tile] = (m[o], np.array([r[ci["filename"]] for r in sel])[o])
    recon_names = {Path(f).name for d in RECON_DIRS for f in glob.glob(str(d / "*.fits"))}
    ev_tab = Table.read(REPO / "crossings" / "solo_v1" / "events.ecsv")
    star = {str(r["event_id"]): (float(r["star_icrs_ra_deg"]), float(r["star_icrs_dec_deg"])) for r in ev_tab}
    rows_out, contact = [], []
    for r in census["visibility"]["rows"]:
        if r["rung"] != "0.1AU" or r["fov_hours"] == 0:
            continue
        e = {"star_icrs_ra_deg": star[r["event_id"]][0], "star_icrs_dec_deg": star[r["event_id"]][1]}
        m_ca = r["mjd_ca"]
        mj = m_ca + np.arange(-PRE_D, POST_D + 1e-9, STEP)
        b_e, along, rh = G.axis_geometry(e, mj)
        s_hat = G.source_direction(e, r["channel"])
        eps, lam, bet = G.solohi_coords(s_hat, mj)
        tiles = G.tile_of(lam, bet)
        inbeam = b_e <= 0.1
        arc_any = inbeam & np.isin(tiles, ["1", "2"])
        if not arc_any.any():
            continue
        tile = str(np.unique(tiles[arc_any])[0])
        intile = tiles == tile
        arc = intile & inbeam
        k0 = int(np.nonzero(arc)[0][0])
        base = intile & ~inbeam & (np.arange(len(mj)) < k0)
        fm, ff = frames[tile]

        def count(mask):
            if not mask.any():
                return 0, None, None, []
            ok = np.zeros(len(fm), bool)
            for k in np.nonzero(mask)[0]:
                lo, hi = np.searchsorted(fm, mj[k] - STEP / 2), np.searchsorted(fm, mj[k] + STEP / 2)
                ok[lo:hi] = True
            sel = np.nonzero(ok)[0]
            return int(ok.sum()), (float(fm[sel].min() - m_ca) * 24 if ok.any() else None), \
                (float(fm[sel].max() - m_ca) * 24 if ok.any() else None), ff[sel].tolist()
        n_arc, a0, a1, arc_files = count(arc)
        n_base, b0, b1, base_files = count(base)
        out = {"event_id": r["event_id"], "target_id": r["target_id"], "channel": r["channel"],
               "orbit": r["orbit"], "t_ca_utc": r["t_ca_utc"], "mjd_ca": m_ca, "tile": tile,
               "arc_geom_hours": round(float(arc.sum()) * STEP * 24, 2),
               "arc_frames": n_arc, "arc_frames_rel_tca_h": [None if a0 is None else round(a0, 2), None if a1 is None else round(a1, 2)],
               "baseline_geom_hours": round(float(base.sum()) * STEP * 24, 2),
               "baseline_frames": n_base, "baseline_rel_tca_h": [None if b0 is None else round(b0, 2), None if b1 is None else round(b1, 2)],
               "arc_b_e_rsun": [round(float(b_e[arc].min()) / G.RSUN_AU, 2), round(float(b_e[arc].max()) / G.RSUN_AU, 2)],
               "arc_r_helio_au": [round(float(rh[arc].min()), 4), round(float(rh[arc].max()), 4)],
               "arc_eps_deg": [round(float(eps[arc].min()), 1), round(float(eps[arc].max()), 1)],
               "baseline_eps_deg": [round(float(eps[base].min()), 1), round(float(eps[base].max()), 1)] if base.any() else None,
               "orbit_lat_deg": round(float(np.median(bet[arc])), 2)}
        for f in arc_files:
            if f in recon_names:
                contact.append({"recon_frame": f, "event_id": r["event_id"], "target_id": r["target_id"],
                                "channel": r["channel"], "role": "arc"})
        for f in base_files:
            if f in recon_names:
                contact.append({"recon_frame": f, "event_id": r["event_id"], "target_id": r["target_id"],
                                "channel": r["channel"], "role": "baseline"})
        rows_out.append(out)
    units = {}
    for o in rows_out:
        u = units.setdefault(f"{o['target_id']}:{o['channel']}", {"events_with_geom_arc": 0, "covered_ge10": 0,
                                                                   "arc_frames": [], "baseline_frames": [], "orbits": []})
        u["events_with_geom_arc"] += 1
        if o["arc_frames"] >= MIN_ARC:
            u["covered_ge10"] += 1
            u["arc_frames"].append(o["arc_frames"])
            u["baseline_frames"].append(o["baseline_frames"])
            u["orbits"].append(o["orbit"])
    for u in units.values():
        u["median_arc_frames"] = float(np.median(u["arc_frames"])) if u["arc_frames"] else 0
        u["min_arc_frames"] = int(min(u["arc_frames"])) if u["arc_frames"] else 0
        u["median_baseline_frames"] = float(np.median(u["baseline_frames"])) if u["baseline_frames"] else 0
        u["min_baseline_frames"] = int(min(u["baseline_frames"])) if u["baseline_frames"] else 0
        u["total_arc_frames"] = int(sum(u["arc_frames"]))
        del u["arc_frames"], u["baseline_frames"]
    res = {"inventory_snapshot": inv["snapshot_utc"], "n_inner_tile_frames": {k: int(len(v[0])) for k, v in frames.items()},
           "recon_contact": contact, "units": units, "rows": rows_out}
    OUT.write_text(json.dumps(res, indent=1) + "\n")
    print("recon contact:", len(contact), json.dumps(contact[:6], indent=None))
    good = {k: u for k, u in units.items() if u["covered_ge10"] >= 3}
    print(len(good), "units with >= 3 covered events;", sum(u["covered_ge10"] for u in good.values()), "events;",
          sum(u["total_arc_frames"] for u in good.values()), "arc frames")
    for k, u in sorted(units.items(), key=lambda kv: -kv[1]["covered_ge10"]):
        print(k, u)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
