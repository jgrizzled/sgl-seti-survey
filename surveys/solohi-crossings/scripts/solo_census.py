"""Solar Orbiter / SoloHI observer-geometry pass (plan §5.23): sunward-
channel census on the Solar Orbiter crossing list (`crossings/solo_v1`),
the Earth-center comparison, and the SoloHI field-of-view visibility of
every sunward event, with the SOAR L2 archive presence folded in.

Stages (mirroring `psp_census.py`):
1. perihelion/orbit table from the observer trajectory;
2. rung x channel census, full era and the SoloHI public-data era;
3. Earth-center (`universal_v1`) comparison (expected: no correspondence);
4. per-event SoloHI visibility: for each S1/S2 event with b_min <= 0.1 AU
   and each rung it satisfies, sample the in-beam window at 1-min steps,
   compute the source's elongation / ram longitude / orbit latitude
   (solo_geometry), flag epochs inside the four nominal tiles, and —
   new for SoloHI, which observes at every heliocentric distance — the
   hours with an L2 frame of the right tile within +-1 h in the SOAR
   inventory ("covered");
5. grazing-rung verdict: the graze-source elongation bound
   atan(b_rung / r_along) per event vs the 5.2 deg inner edge, and the
   mission floor from the smallest perihelion;
6. per-orbit ledger of every sunward event-rung;
7. S2 1 AU ledger (out of scope): tile transits of every target star.

Output: results/solo_census_v1.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solo_geometry as G

REPO = Path(__file__).resolve().parents[3]
SOLO = REPO / "crossings" / "solo_v1" / "events.ecsv"
EARTH = REPO / "crossings" / "universal_v1" / "events.ecsv"
INV = REPO / "runs" / "solohi-crossings" / "coverage" / "soar_l2_inventory.json"
OUT = REPO / "surveys" / "solohi-crossings" / "results" / "solo_census_v1.json"

RSUN_AU = G.RSUN_AU
RUNGS = {
    "S1_downlink_postlens": [("1.2Rsun", 1.2 * RSUN_AU), ("2.5Rsun", 2.5 * RSUN_AU), ("0.1AU", 0.1)],
    "S2_uplink_pastsun": [("0.1AU", 0.1), ("1AU", 1.0)],
}
COMBOS = {"S1_downlink_postlens": ("outbound", "target"),
          "S2_uplink_pastsun": ("inbound", "anti_target")}
STEP_MIN = 60.0
FINE_MIN = 1.0
SPAN_D = 45.0
COVER_TOL_H = 1.0     # an L2 frame of the tile within this of the epoch
BASE_D = 10.0         # days before the arc scanned for an off-beam in-field baseline
INNER_EDGE = G.INNER_EDGE_DEG


def combo_mask(t, direction, side):
    return (np.asarray(t["link_direction"]) == direction) & (np.asarray(t["side"]) == side)


def mjd_ca(ev) -> float:
    return float(Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd)


def load_inventory():
    """Per tile: sorted MJDs of L2 frame begin times."""
    inv = json.load(open(INV))
    ci = {c: i for i, c in enumerate(inv["columns"])}
    per = {k: [] for k in G.TILES}
    for r in inv["rows"]:
        for k, tile in G.TILES.items():
            if r[ci["descriptor"]] in tile["descriptors"]:
                per[k].append(Time(r[ci["begin_time"]]).mjd)
    per = {k: np.sort(np.array(v)) for k, v in per.items()}
    era = (min(v.min() for v in per.values()), max(v.max() for v in per.values()))
    return per, era, inv["snapshot_utc"]


def covered(mjd, tile_ids, frames) -> np.ndarray:
    """True where a frame of the epoch's tile lies within COVER_TOL_H."""
    out = np.zeros(len(mjd), bool)
    for k, fr in frames.items():
        m = tile_ids == k
        if not m.any() or len(fr) == 0:
            continue
        j = np.searchsorted(fr, mjd[m])
        lo = fr[np.clip(j - 1, 0, len(fr) - 1)]
        hi = fr[np.clip(j, 0, len(fr) - 1)]
        d = np.minimum(np.abs(mjd[m] - lo), np.abs(hi - mjd[m])) * 24
        out[m] = d <= COVER_TOL_H
    return out


def which_orbit(m, orbits):
    for o in orbits:
        if o["mjd_start"] <= m <= o["mjd_end"]:
            return o
    return None


def census(t, era):
    b = np.asarray(t["b_min_au"])
    tca = np.array([mjd_ca(ev) for ev in t])
    in_era = (tca >= era[0]) & (tca <= era[1])
    out = {}
    for name, rungs in RUNGS.items():
        m0 = combo_mask(t, *COMBOS[name])
        out[name] = {}
        for label_era, me in (("full_2020_2028", m0), ("solohi_public_L2_era", m0 & in_era)):
            row = {}
            for label, rmax in rungs:
                mr = me & (b <= rmax)
                tg = sorted(set(np.asarray(t["target_id"])[mr]))
                ent = {"events": int(mr.sum()), "targets": len(tg)}
                if 0 < mr.sum() <= 80:
                    ent["target_ids"] = tg
                    ent["min_b_rsun_by_target"] = {
                        x: round(float(b[mr & (np.asarray(t["target_id"]) == x)].min()) / RSUN_AU, 3)
                        for x in tg}
                row[label] = ent
            out[name][label_era] = row
    return out


def earth_comparison(solo, earth):
    jd = np.asarray(earth["t_ca_tdb_jd"])
    e_era = earth[(jd >= Time("2020-05-01").tdb.jd) & (jd <= Time("2028-01-01").tdb.jd)]
    d_t, d_b, n_unmatched, n_total = [], [], 0, 0
    for name, (direction, side) in COMBOS.items():
        s = solo[combo_mask(solo, direction, side)]
        e = e_era[combo_mask(e_era, direction, side)]
        for tid in sorted(set(np.asarray(s["target_id"]))):
            st = s[np.asarray(s["target_id"]) == tid]
            et = e[np.asarray(e["target_id"]) == tid]
            for row in st:
                n_total += 1
                if len(et) == 0:
                    n_unmatched += 1
                    continue
                dt = np.asarray(et["t_ca_tdb_jd"]) - row["t_ca_tdb_jd"]
                j = int(np.argmin(np.abs(dt)))
                if abs(dt[j]) > 120:
                    n_unmatched += 1
                    continue
                d_t.append(dt[j])
                d_b.append(row["b_min_au"] - et["b_min_au"][j])
    d_t, d_b = np.array(d_t), np.abs(np.array(d_b))
    return {
        "solo_sunward_events": n_total, "earth_era_events": int(len(e_era)),
        "matched_within_120d": int(len(d_t)), "unmatched": n_unmatched,
        "delta_tca_days_abs": {"median": float(np.median(np.abs(d_t))) if len(d_t) else None,
                               "max": float(np.max(np.abs(d_t))) if len(d_t) else None},
        "abs_delta_b_au": {"median": float(np.median(d_b)) if len(d_b) else None,
                           "p95": float(np.percentile(d_b, 95)) if len(d_b) else None,
                           "max": float(d_b.max()) if len(d_b) else None},
        "note": "Solar Orbiter crosses each axis twice per 150-230 d orbit at "
                "heliocentric radii 0.28-1.0 AU; Earth-center events (annual, "
                "1 AU) are a different population — nearest-in-time pairs are "
                "coincidences, not corrections. The solo list is mandatory",
    }


def event_visibility(ev, channel, rung_label, rung_au, orbits, frames):
    e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
         "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
    m_ca = mjd_ca(ev)
    mjds = m_ca + np.arange(-SPAN_D, SPAN_D + 1e-9, STEP_MIN / 1440.0)
    b_e, _, _ = G.axis_geometry(e, mjds)
    inbeam = b_e <= rung_au
    k = int(np.argmin(np.abs(mjds - m_ca)))
    if not inbeam[k]:
        return None
    lo = k
    while lo > 0 and inbeam[lo - 1]:
        lo -= 1
    hi = k
    while hi < len(mjds) - 1 and inbeam[hi + 1]:
        hi += 1
    truncated = bool(lo == 0 or hi == len(mjds) - 1)
    mw = np.arange(mjds[max(lo - 1, 0)], mjds[min(hi + 1, len(mjds) - 1)] + 1e-9,
                   FINE_MIN / 1440.0)
    bw, aw, rw = G.axis_geometry(e, mw)
    keep = bw <= rung_au
    mw, bw, aw, rw = mw[keep], bw[keep], aw[keep], rw[keep]
    kc = int(np.argmin(np.abs(mw - m_ca)))
    s_hat = G.source_direction(e, channel)
    eps, lam, bet = G.solohi_coords(s_hat, mw)
    tiles = G.tile_of(lam, bet)
    fov = tiles != ""
    cov = covered(mw, tiles, frames) & fov
    post = mw > m_ca
    anti_pre = float(np.mean(lam[~post] < 0)) if (~post).any() else float("nan")
    anti_post = float(np.mean(lam[post] < 0)) if post.any() else float("nan")
    orb = which_orbit(m_ca, orbits)
    row = {"event_id": str(ev["event_id"]), "target_id": str(ev["target_id"]),
           "channel": channel, "rung": rung_label,
           "t_ca_utc": str(ev["t_ca_utc"])[:16], "mjd_ca": round(m_ca, 4),
           "b_min_rsun": round(float(ev["b_min_au"]) / RSUN_AU, 3),
           "r_along_au_at_tca": round(float(abs(aw[kc])), 4),
           "r_helio_au_at_tca": round(float(rw[kc]), 4),
           "window_hours": round(float(mw[-1] - mw[0]) * 24, 2),
           "window_truncated_by_span": truncated,
           "eps_deg_at_tca": round(float(eps[kc]), 2),
           "eps_deg_max_in_window": round(float(eps.max()), 2),
           "eps_deg_rung_edge_bound": round(float(np.degrees(np.arctan2(rung_au, abs(aw[kc])))), 2),
           "ram_lon_deg_min_in_window": round(float(lam.min()), 2),
           "anti_ram_fraction_pre_tca": round(anti_pre, 2),
           "anti_ram_fraction_post_tca": round(anti_post, 2),
           "orbit_lat_deg_at_tca": round(float(bet[kc]), 2),
           "orbit": orb["orbit"] if orb else None,
           "days_from_perihelion": round(m_ca - orb["mjd_perihelion"], 2) if orb else None,
           "fov_hours": round(float(fov.sum()) * FINE_MIN / 60, 2),
           "fov_hours_by_tile": {k: round(float((tiles == k).sum()) * FINE_MIN / 60, 2)
                                 for k in G.TILES if (tiles == k).any()},
           "covered_hours": round(float(cov.sum()) * FINE_MIN / 60, 2)}
    if fov.any():
        vm = mw[fov]
        row["fov_rel_tca_h"] = [round(float(vm.min() - m_ca) * 24, 2),
                                round(float(vm.max() - m_ca) * 24, 2)]
        row["fov_eps_range_deg"] = [round(float(eps[fov].min()), 2), round(float(eps[fov].max()), 2)]
        row["fov_b_e_range_rsun"] = [round(float(bw[fov].min()) / RSUN_AU, 2),
                                     round(float(bw[fov].max()) / RSUN_AU, 2)]
        row["fov_r_helio_range_au"] = [round(float(rw[fov].min()), 4), round(float(rw[fov].max()), 4)]
    if cov.any():
        vm = mw[cov]
        row["covered_rel_tca_h"] = [round(float(vm.min() - m_ca) * 24, 2),
                                    round(float(vm.max() - m_ca) * 24, 2)]
        row["covered_eps_range_deg"] = [round(float(eps[cov].min()), 2), round(float(eps[cov].max()), 2)]
    # off-beam, in-field baseline preceding the arc (the HI-1 star-fixed
    # pattern): up to BASE_D days before the first in-field in-beam epoch
    if fov.any():
        t0 = mw[fov].min()
        mb = np.arange(t0 - BASE_D, t0, 10.0 / 1440.0)
        bb, _, _ = G.axis_geometry(e, mb)
        _, lb, btb = G.solohi_coords(s_hat, mb)
        tb = G.tile_of(lb, btb)
        off = (bb > rung_au) & (tb != "")
        inner = off & np.isin(tb, ["1", "2"])
        covb = covered(mb, tb, frames) & off
        row["baseline_hours_infield_offbeam"] = round(float(off.sum()) / 6, 1)
        row["baseline_hours_inner_tiles"] = round(float(inner.sum()) / 6, 1)
        row["baseline_hours_covered"] = round(float(covb.sum()) / 6, 1)
        row["baseline_eps_range_deg"] = [round(float(eps_b.min()), 1), round(float(eps_b.max()), 1)] if (eps_b := np.degrees(np.arccos(np.clip(np.sum(np.broadcast_to(s_hat, (len(mb), 3)) * G.ram_frame(mb)[0], axis=-1), -1, 1)))[off]).size else None
    return row


def visibility(solo, orbits, frames):
    rows, agg = [], {"_degraded_skipped": 0}
    b = np.asarray(solo["b_min_au"])
    valid = np.asarray(solo["validity"]) == "valid"
    for name, (direction, side) in COMBOS.items():
        ch = name[:2]
        m = combo_mask(solo, direction, side) & (b <= 0.1)
        agg["_degraded_skipped"] += int((m & ~valid).sum())
        sub = solo[m & valid]
        for ev in sub:
            for label, rung_au in RUNGS[name]:
                if label == "1AU" or float(ev["b_min_au"]) > rung_au:
                    continue
                row = event_visibility(ev, ch, label, rung_au, orbits, frames)
                if row is None:
                    continue
                rows.append(row)
                a = agg.setdefault(f"{ch}_{label}", {}).setdefault(row["target_id"], {
                    "events": 0, "fov_events": 0, "covered_events": 0,
                    "hours_fov": [], "hours_cov": [], "max_eps_edge_bound_deg": 0.0})
                a["events"] += 1
                a["fov_events"] += row["fov_hours"] > 0
                a["covered_events"] += row["covered_hours"] > 0
                a["max_eps_edge_bound_deg"] = max(a["max_eps_edge_bound_deg"],
                                                  row["eps_deg_rung_edge_bound"])
                if row["fov_hours"] > 0:
                    a["hours_fov"].append(row["fov_hours"])
                if row["covered_hours"] > 0:
                    a["hours_cov"].append(row["covered_hours"])
    for key in agg:
        if key.startswith("_"):
            continue
        for tid, a in agg[key].items():
            hf, hc = a.pop("hours_fov"), a.pop("hours_cov")
            a["median_fov_hours"] = float(np.median(hf)) if hf else 0.0
            a["median_covered_hours"] = float(np.median(hc)) if hc else 0.0
            a["total_covered_hours"] = float(np.sum(hc)) if hc else 0.0
    return rows, agg


def grazing_verdict(rows, orbits):
    qmin = min(o["q_au"] for o in orbits)
    out = {"mission_min_perihelion_au": round(qmin, 4),
           "mission_min_perihelion_rsun": round(qmin / RSUN_AU, 2),
           "inner_edge_deg": INNER_EDGE,
           "graze_source_elongation_floor_deg": {
               "1.2Rsun": round(float(np.degrees(np.arctan(1.2 * RSUN_AU / qmin))), 2),
               "2.5Rsun": round(float(np.degrees(np.arctan(2.5 * RSUN_AU / qmin))), 2)},
           "r_along_needed_for_inner_edge_rsun": {
               "1.2Rsun": round(1.2 / np.tan(np.radians(INNER_EDGE)), 2),
               "2.5Rsun": round(2.5 / np.tan(np.radians(INNER_EDGE)), 2)},
           "per_rung": {}}
    for label in ("1.2Rsun", "2.5Rsun"):
        rr = [r for r in rows if r["channel"] == "S1" and r["rung"] == label]
        if not rr:
            out["per_rung"][label] = {"events": 0}
            continue
        best = max(rr, key=lambda r: r["eps_deg_rung_edge_bound"])
        out["per_rung"][label] = {
            "events": len(rr),
            "events_with_edge_bound_ge_inner_edge": sum(
                r["eps_deg_rung_edge_bound"] >= INNER_EDGE for r in rr),
            "events_with_fov_hours": sum(r["fov_hours"] > 0 for r in rr),
            "best_event": {k: best[k] for k in (
                "target_id", "t_ca_utc", "orbit", "b_min_rsun", "r_along_au_at_tca",
                "eps_deg_at_tca", "eps_deg_max_in_window", "eps_deg_rung_edge_bound",
                "ram_lon_deg_min_in_window", "fov_hours")},
            "min_r_along_au": min(r["r_along_au_at_tca"] for r in rr)}
    return out


def per_orbit(rows, orbits):
    out = []
    for o in orbits:
        ev = [r for r in rows if r["orbit"] == o["orbit"]]
        out.append({**{k: o[k] for k in ("orbit", "utc_perihelion", "q_au", "q_rsun",
                                          "utc_start", "utc_end", "period_days")},
                    "sunward_event_rungs": len(ev),
                    "fov_rungs": sum(r["fov_hours"] > 0 for r in ev),
                    "fov_hours_total": round(sum(r["fov_hours"] for r in ev), 1),
                    "covered_rungs": sum(r["covered_hours"] > 0 for r in ev),
                    "covered_hours_total": round(sum(r["covered_hours"] for r in ev), 1),
                    "events": [{k: r[k] for k in (
                        "target_id", "channel", "rung", "t_ca_utc", "days_from_perihelion",
                        "b_min_rsun", "r_helio_au_at_tca", "eps_deg_at_tca",
                        "eps_deg_rung_edge_bound", "fov_hours", "covered_hours")} for r in ev]})
    return out


def s2_1au_ledger(solo, orbits, frames):
    stars = {}
    for ev in solo:
        stars.setdefault(str(ev["target_id"]), (float(ev["star_icrs_ra_deg"]),
                                                 float(ev["star_icrs_dec_deg"])))
    per = {}
    for tid, (ra, dec) in stars.items():
        s = G.radec_to_vec(ra, dec)
        n_o, hrs = 0, []
        for o in orbits:
            m = np.arange(o["mjd_start"], o["mjd_end"], 1 / 24)
            _, lam, bet = G.solohi_coords(s, m)
            tiles = G.tile_of(lam, bet)
            c = covered(m, tiles, frames) & (tiles != "")
            n_o += (tiles != "").any()
            if c.any():
                hrs.append(float(c.sum()))
        per[tid] = {"orbits_with_fov_transit": int(n_o),
                    "orbits_with_covered_transit": len(hrs),
                    "median_covered_hours": float(np.median(hrs)) if hrs else 0.0}
    return {"n_stars": len(per), "orbits": len(orbits),
            "stars_with_any_fov_transit": sum(1 for a in per.values() if a["orbits_with_fov_transit"]),
            "total_star_orbit_covered_pairs": sum(a["orbits_with_covered_transit"] for a in per.values()),
            "per_target": per}


def main():
    solo, earth = Table.read(SOLO), Table.read(EARTH)
    assert set(np.unique(solo["validity"])) <= {"valid", "degraded"}
    orbits = G.perihelia()
    frames, era, snap = load_inventory()
    rows, agg = visibility(solo, orbits, frames)
    out = {
        "study": "solo_census_v1",
        "input": str(SOLO.relative_to(REPO)),
        "crossings_id": str(solo.meta.get("crossings_id", "")),
        "n_events": len(solo),
        "solohi_fov_model": {"tiles": G.TILES, "inner_edge_deg": G.INNER_EDGE_DEG,
                             "outer_edge_deg": G.OUTER_EDGE_DEG,
                             "note": "header-measured tile edges in the orbit-plane "
                                     "(ram) frame, anti-ram side; see solo_geometry"},
        "soar_inventory": {"snapshot_utc": snap,
                           "l2_era_mjd": [float(era[0]), float(era[1])],
                           "l2_era_utc": [Time(era[0], format="mjd").iso[:10],
                                          Time(era[1], format="mjd").iso[:10]],
                           "frames_per_tile": {k: int(len(v)) for k, v in frames.items()},
                           "cover_tolerance_h": COVER_TOL_H},
        "orbits": orbits,
        "sunward_census": census(solo, era),
        "earth_comparison": earth_comparison(solo, earth),
        "grazing_verdict": grazing_verdict(rows, orbits),
        "visibility": {"aggregates": agg, "rows": rows},
        "per_orbit": per_orbit(rows, orbits),
        "s2_1AU_transit_ledger": s2_1au_ledger(solo, orbits, frames),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["sunward_census"], indent=1))
    print(json.dumps(out["earth_comparison"], indent=1))
    print(json.dumps(out["grazing_verdict"], indent=1))
    print(json.dumps(agg, indent=1))
    led = out["s2_1AU_transit_ledger"]
    print("S2 1AU ledger:", led["stars_with_any_fov_transit"], "stars,",
          led["total_star_orbit_covered_pairs"], "star-orbit covered pairs")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
