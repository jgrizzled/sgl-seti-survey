"""PSP/WISPR observer-geometry pass (plan §5.15 O5): sunward-channel
census on the PSP-observer crossing list (`crossings/psp_v1`), the
Earth-center comparison, and the WISPR field-of-view visibility of every
sunward event per encounter.

Stages:
1. encounter table from the observer trajectory (r < 0.25 AU arcs);
2. rung x channel census, full era and inside-encounter subsets;
3. Earth-center (`universal_v1`) comparison: nearest same-target /
   direction / side event within 120 d (expected: no correspondence —
   an 88-150 d orbit vs the annual Earth cycle);
4. per-event WISPR visibility: for each S1/S2 event with b_min <= 0.1 AU
   and each rung it satisfies, sample the in-beam window at 10-min
   steps, compute the source's elongation / ram longitude / orbit
   latitude (psp_geometry), flag epochs inside WISPR-I / WISPR-O while
   in-beam and inside an encounter; report visible hours, elongation
   spans, heliocentric distance, and the encounter;
5. grazing-rung verdict: the graze-source elongation bound
   atan(b_rung / r_along) per event vs the 13.5 deg inner edge, and the
   mission-wide floor from the smallest perihelion;
6. per-encounter ledger of every sunward event inside the arc;
7. S2 1 AU ledger: WISPR-I/O transits of every target star per encounter
   (any b), for the out-of-scope decision record.

Output: results/psp_census_v1.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psp_geometry as G

REPO = Path(__file__).resolve().parents[3]
PSP = REPO / "crossings" / "psp_v1" / "events.ecsv"
EARTH = REPO / "crossings" / "universal_v1" / "events.ecsv"
OUT = REPO / "surveys" / "wispr-crossings" / "results" / "psp_census_v1.json"

RSUN_AU = G.RSUN_AU
RUNGS = {
    "S1_downlink_postlens": [("1.2Rsun", 1.2 * RSUN_AU), ("2.5Rsun", 2.5 * RSUN_AU), ("0.1AU", 0.1)],
    "S2_uplink_pastsun": [("0.1AU", 0.1), ("1AU", 1.0)],
}
COMBOS = {"S1_downlink_postlens": ("outbound", "target"),
          "S2_uplink_pastsun": ("inbound", "anti_target")}
STEP_MIN = 60.0      # coarse pass to bracket the in-beam window
FINE_MIN = 1.0       # resampling step inside the window
SPAN_D = 45.0        # +-days sampled about t_ca for the in-beam window
INNER_EDGE = G.WISPR_I["lam_deg"][0]


def combo_mask(t, direction, side):
    return (np.asarray(t["link_direction"]) == direction) & (np.asarray(t["side"]) == side)


def mjd_ca(ev) -> float:
    return float(Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd)


def in_encounter(mjd, enc) -> np.ndarray:
    mjd = np.atleast_1d(mjd)
    out = np.zeros(len(mjd), bool)
    for e in enc:
        out |= (mjd >= e["mjd_start"]) & (mjd <= e["mjd_end"])
    return out


def which_encounter(m, enc, slack_d=0.0):
    for e in enc:
        if e["mjd_start"] - slack_d <= m <= e["mjd_end"] + slack_d:
            return e
    return None


def census(t, enc):
    b = np.asarray(t["b_min_au"])
    tca = np.array([mjd_ca(ev) for ev in t])
    inside = in_encounter(tca, enc)
    out = {}
    for name, rungs in RUNGS.items():
        m0 = combo_mask(t, *COMBOS[name])
        out[name] = {}
        for era, me in (("full_2018_2026", m0), ("inside_encounters_r<0.25AU", m0 & inside)):
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
            out[name][era] = row
    return out


def earth_comparison(psp, earth):
    jd = np.asarray(earth["t_ca_tdb_jd"])
    e_era = earth[(jd >= Time("2018-08-15").tdb.jd) & (jd <= Time("2026-12-01").tdb.jd)]
    d_t, d_b, n_unmatched, n_total = [], [], 0, 0
    for name, (direction, side) in COMBOS.items():
        s = psp[combo_mask(psp, direction, side)]
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
        "psp_sunward_events": n_total, "earth_era_events": int(len(e_era)),
        "matched_within_120d": int(len(d_t)), "unmatched": n_unmatched,
        "delta_tca_days_abs": {"median": float(np.median(np.abs(d_t))) if len(d_t) else None,
                               "max": float(np.max(np.abs(d_t))) if len(d_t) else None},
        "abs_delta_b_au": {"median": float(np.median(d_b)) if len(d_b) else None,
                           "p95": float(np.percentile(d_b, 95)) if len(d_b) else None,
                           "max": float(d_b.max()) if len(d_b) else None},
        "note": "PSP crosses each axis twice per 88-150 d orbit at heliocentric "
                "radii 0.05-0.7 AU; Earth-center events (annual, 1 AU) are a "
                "different population — nearest-in-time pairs are coincidences, "
                "not corrections. The PSP list is mandatory",
    }


def event_visibility(ev, channel, rung_label, rung_au, enc):
    e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
         "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
    m_ca = mjd_ca(ev)
    # coarse bracket of the contiguous in-beam interval containing t_ca
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
    # fine resampling (1 min) of the bracketed window, one coarse step
    # of slack each side, then the exact in-beam mask
    mw = np.arange(mjds[max(lo - 1, 0)], mjds[min(hi + 1, len(mjds) - 1)] + 1e-9,
                   FINE_MIN / 1440.0)
    bw, aw, rw = G.axis_geometry(e, mw)
    keep = bw <= rung_au
    mw, bw, aw, rw = mw[keep], bw[keep], aw[keep], rw[keep]
    kc = int(np.argmin(np.abs(mw - m_ca)))
    s_hat = G.source_direction(e, channel)
    eps, lam, bet = G.wispr_coords(s_hat, mw)
    enc_w = in_encounter(mw, enc)
    fov_i = G.in_wispr(lam, bet, G.WISPR_I)
    fov_o = G.in_wispr(lam, bet, G.WISPR_O)
    vis_i, vis_o = fov_i & enc_w, fov_o & enc_w
    post = mw > m_ca
    ram_post = float(np.mean(lam[post] > 0)) if post.any() else float("nan")
    ram_pre = float(np.mean(lam[~post] > 0)) if (~post).any() else float("nan")
    enc_hit = which_encounter(m_ca, enc, slack_d=(mw[-1] - mw[0]) / 2)
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
           "ram_lon_deg_max_in_window": round(float(lam.max()), 2),
           "ram_side_fraction_post_tca": round(ram_post, 2),
           "ram_side_fraction_pre_tca": round(ram_pre, 2),
           "orbit_lat_deg_at_tca": round(float(bet[kc]), 2),
           "in_encounter_at_tca": bool(enc_w[kc]),
           "encounter": enc_hit["encounter"] if enc_hit else None,
           "days_from_perihelion": round(m_ca - enc_hit["mjd_perihelion"], 2) if enc_hit else None,
           "fov_hours_wispr_i_any_r": round(float(fov_i.sum()) * FINE_MIN / 60, 2),
           "fov_hours_wispr_o_any_r": round(float(fov_o.sum()) * FINE_MIN / 60, 2),
           "visible_hours_wispr_i": round(float(vis_i.sum()) * FINE_MIN / 60, 2),
           "visible_hours_wispr_o": round(float(vis_o.sum()) * FINE_MIN / 60, 2)}
    for cam, vis in (("wispr_i", vis_i), ("wispr_o", vis_o)):
        if vis.any():
            vm = mw[vis]
            row[f"{cam}_vis_rel_tca_h"] = [round(float(vm.min() - m_ca) * 24, 2),
                                          round(float(vm.max() - m_ca) * 24, 2)]
            row[f"{cam}_eps_range_deg"] = [round(float(eps[vis].min()), 2),
                                           round(float(eps[vis].max()), 2)]
            row[f"{cam}_b_e_range_rsun"] = [round(float(bw[vis].min()) / RSUN_AU, 2),
                                            round(float(bw[vis].max()) / RSUN_AU, 2)]
            row[f"{cam}_r_helio_range_au"] = [round(float(rw[vis].min()), 4),
                                              round(float(rw[vis].max()), 4)]
    return row


def visibility(psp, enc):
    rows, agg = [], {"_degraded_skipped": 0}
    b = np.asarray(psp["b_min_au"])
    valid = np.asarray(psp["validity"]) == "valid"
    for name, (direction, side) in COMBOS.items():
        ch = name[:2]
        m = combo_mask(psp, direction, side) & (b <= 0.1)
        agg["_degraded_skipped"] += int((m & ~valid).sum())
        sub = psp[m & valid]
        for ev in sub:
            for label, rung_au in RUNGS[name]:
                if label == "1AU" or float(ev["b_min_au"]) > rung_au:
                    continue
                row = event_visibility(ev, ch, label, rung_au, enc)
                if row is None:
                    continue
                rows.append(row)
                a = agg.setdefault(f"{ch}_{label}", {}).setdefault(row["target_id"], {
                    "events": 0, "in_encounter": 0, "visible_events_wispr_i": 0,
                    "visible_events_wispr_o": 0, "hours_i": [], "hours_o": [],
                    "max_eps_edge_bound_deg": 0.0})
                a["events"] += 1
                a["in_encounter"] += row["in_encounter_at_tca"]
                a["visible_events_wispr_i"] += row["visible_hours_wispr_i"] > 0
                a["visible_events_wispr_o"] += row["visible_hours_wispr_o"] > 0
                a["max_eps_edge_bound_deg"] = max(a["max_eps_edge_bound_deg"],
                                                  row["eps_deg_rung_edge_bound"])
                if row["visible_hours_wispr_i"] > 0:
                    a["hours_i"].append(row["visible_hours_wispr_i"])
                if row["visible_hours_wispr_o"] > 0:
                    a["hours_o"].append(row["visible_hours_wispr_o"])
    for key in agg:
        if key.startswith("_"):
            continue
        for tid, a in agg[key].items():
            hi, ho = a.pop("hours_i"), a.pop("hours_o")
            a["median_visible_hours_wispr_i"] = float(np.median(hi)) if hi else 0.0
            a["median_visible_hours_wispr_o"] = float(np.median(ho)) if ho else 0.0
            a["total_visible_hours_wispr_i"] = float(np.sum(hi)) if hi else 0.0
    return rows, agg


def grazing_verdict(rows, enc):
    qmin = min(e["q_au"] for e in enc)
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
            "events_in_encounter": sum(r["in_encounter_at_tca"] for r in rr),
            "events_with_edge_bound_ge_inner_edge": sum(
                r["eps_deg_rung_edge_bound"] >= INNER_EDGE for r in rr),
            "events_with_ram_lon_ge_inner_edge": sum(
                r["ram_lon_deg_max_in_window"] >= INNER_EDGE for r in rr),
            "events_with_wispr_i_hours": sum(r["visible_hours_wispr_i"] > 0 for r in rr),
            "best_event": {k: best[k] for k in (
                "target_id", "t_ca_utc", "encounter", "b_min_rsun", "r_along_au_at_tca",
                "eps_deg_at_tca", "eps_deg_max_in_window", "eps_deg_rung_edge_bound",
                "ram_lon_deg_max_in_window", "visible_hours_wispr_i")},
            "min_r_along_au": min(r["r_along_au_at_tca"] for r in rr)}
    return out


def per_encounter(rows, enc):
    out = []
    for e in enc:
        ev = [r for r in rows if r["encounter"] == e["encounter"]]
        ent = {**{k: e[k] for k in ("encounter", "utc_perihelion", "q_au", "q_rsun",
                                    "utc_start", "utc_end", "duration_days")},
               "sunward_event_rungs_in_arc": len(ev),
               "wispr_i_visible_rungs": sum(r["visible_hours_wispr_i"] > 0 for r in ev),
               "wispr_i_hours_total": round(sum(r["visible_hours_wispr_i"] for r in ev), 1),
               "wispr_o_visible_rungs": sum(r["visible_hours_wispr_o"] > 0 for r in ev),
               "events": [{k: r[k] for k in (
                   "target_id", "channel", "rung", "t_ca_utc", "days_from_perihelion",
                   "b_min_rsun", "r_along_au_at_tca", "eps_deg_at_tca",
                   "eps_deg_max_in_window", "eps_deg_rung_edge_bound",
                   "ram_lon_deg_max_in_window", "visible_hours_wispr_i",
                   "visible_hours_wispr_o")} for r in ev]}
        out.append(ent)
    return out


def s2_1au_ledger(psp, enc):
    """Per star: encounters during which the star direction crosses
    WISPR-I / WISPR-O (any b) — the widest S2 rung, out of scope."""
    stars = {}
    for ev in psp:
        stars.setdefault(str(ev["target_id"]), (float(ev["star_icrs_ra_deg"]),
                                                float(ev["star_icrs_dec_deg"])))
    per = {}
    for tid, (ra, dec) in stars.items():
        s = G.radec_to_vec(ra, dec)
        n_i = n_o = 0
        hrs_i = []
        for e in enc:
            m = np.arange(e["mjd_start"], e["mjd_end"], 1 / 24)
            _, lam, bet = G.wispr_coords(s, m)
            vi, vo = G.in_wispr(lam, bet, G.WISPR_I), G.in_wispr(lam, bet, G.WISPR_O)
            n_i += vi.any()
            n_o += vo.any()
            if vi.any():
                hrs_i.append(float(vi.sum()))
        per[tid] = {"encounters_with_wispr_i_transit": int(n_i),
                    "encounters_with_wispr_o_transit": int(n_o),
                    "median_wispr_i_hours": float(np.median(hrs_i)) if hrs_i else 0.0}
    return {"n_stars": len(per), "encounters": len(enc),
            "stars_with_any_wispr_i_transit": sum(1 for a in per.values()
                                                  if a["encounters_with_wispr_i_transit"]),
            "total_star_encounter_wispr_i_pairs": sum(a["encounters_with_wispr_i_transit"]
                                                      for a in per.values()),
            "per_target": per}


def main():
    psp, earth = Table.read(PSP), Table.read(EARTH)
    assert set(np.unique(psp["validity"])) <= {"valid", "degraded"}
    enc = G.encounters()
    rows, agg = visibility(psp, enc)
    out = {
        "study": "psp_census_v1",
        "input": str(PSP.relative_to(REPO)),
        "crossings_id": str(psp.meta.get("crossings_id", "")),
        "n_events": len(psp),
        "wispr_fov_model": {"WISPR_I": G.WISPR_I, "WISPR_O": G.WISPR_O,
                            "encounter_r_au": G.ENCOUNTER_R_AU,
                            "note": "edges from the WISPR Data Users Guide v5 "
                                    "(13.5 deg sunward edge, 108.5 deg outer, 3 deg "
                                    "overlap, 40/58 deg fields); ram side and "
                                    "half-heights from Vourlidas+2016 — recon "
                                    "header-verification items"},
        "encounters": enc,
        "sunward_census": census(psp, enc),
        "earth_comparison": earth_comparison(psp, earth),
        "grazing_verdict": grazing_verdict(rows, enc),
        "visibility": {"aggregates": agg, "rows": rows},
        "per_encounter": per_encounter(rows, enc),
        "s2_1AU_transit_ledger": s2_1au_ledger(psp, enc),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["sunward_census"], indent=1))
    print(json.dumps(out["earth_comparison"], indent=1))
    print(json.dumps(out["grazing_verdict"], indent=1))
    print(json.dumps(agg, indent=1))
    led = out["s2_1AU_transit_ledger"]
    print("S2 1AU ledger:", led["stars_with_any_wispr_i_transit"], "stars,",
          led["total_star_encounter_wispr_i_pairs"], "star-encounter WISPR-I pairs")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
