"""Sunward-channel census on the STEREO-A observer crossing list
(`crossings/stereoa_v1`), the Earth-center validation, and the HI-1
field-of-view visibility of every 0.1 AU sunward event.

Stages:
1. rung x era census (the soho_census pattern);
2. Earth-center (`universal_v1`) validation: nearest same-target /
   direction / side event within 120 d -> |dt_ca|, |db| (the
   drifting-observer correction, Kepler pattern);
3. helioprojective transform validation against the recon test frame
   header WCS (RA/Dec 'A' system -> pixel -> primary HPLN/HPLT);
4. per-event HI-1 visibility: for each S1/S2 event with b_min <= 0.1 AU,
   sample the filled-cone window at 10-min steps, flag epochs whose
   source (fixed ICRS direction) lies inside the nominal HI-1A footprint
   while in-beam; report visible hours, elongation span, and which side
   of t_ca the visible arc falls on;
5. S2 1.0 AU ledger: HI-1 transits of every target star (any b), for
   the out-of-scope decision record.

Output: results/stereoa_census_v1.json
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hi_geometry as G

REPO = Path(__file__).resolve().parents[3]
STA = REPO / "crossings" / "stereoa_v1" / "events.ecsv"
EARTH = REPO / "crossings" / "universal_v1" / "events.ecsv"
TEST_FRAME = REPO / "runs" / "stereo-hi-crossings" / "recon" / "20100615_000901_14h1A.fts"
OUT = REPO / "surveys" / "stereo-hi-crossings" / "results" / "stereoa_census_v1.json"

RSUN_AU = G.RSUN_AU
RUNGS = {
    "S1_downlink_postlens": [("1.2Rsun", 1.2 * RSUN_AU), ("2.5Rsun", 2.5 * RSUN_AU), ("0.1AU", 0.1)],
    "S2_uplink_pastsun": [("0.1AU", 0.1), ("1AU", 1.0)],
}
COMBOS = {"S1_downlink_postlens": ("outbound", "target"),
          "S2_uplink_pastsun": ("inbound", "anti_target")}
ERAS = {"full_2007_2026": (None, None),
        "east_pointing_2007_2014-08": (None, "2014-08-19"),
        "gap_2014-08_2015-11": ("2014-08-19", "2015-11-16"),
        "west_pointing_2015-11_2023-08": ("2015-11-17", "2023-08-15"),
        "east_pointing_2023-08_on": ("2023-08-16", None)}
CADENCE_MIN = 40.0


def combo_mask(t, direction, side):
    return (np.asarray(t["link_direction"]) == direction) & (np.asarray(t["side"]) == side)


def era_mask(t, start, end):
    jd = np.asarray(t["t_ca_tdb_jd"])
    m = np.ones(len(t), bool)
    if start:
        m &= jd >= Time(start).tdb.jd
    if end:
        m &= jd <= Time(end).tdb.jd
    return m


def census(t):
    b = np.asarray(t["b_min_au"])
    out = {}
    for name, rungs in RUNGS.items():
        m0 = combo_mask(t, *COMBOS[name])
        out[name] = {}
        for era, (s, e) in ERAS.items():
            me = m0 & era_mask(t, s, e)
            row = {}
            for label, rmax in rungs:
                mr = me & (b <= rmax)
                tg = sorted(set(np.asarray(t["target_id"])[mr]))
                ent = {"events": int(mr.sum()), "targets": len(tg)}
                if 0 < mr.sum() <= 60:
                    ent["target_ids"] = tg
                    ent["min_b_rsun_by_target"] = {
                        x: round(float(b[mr & (np.asarray(t["target_id"]) == x)].min()) / RSUN_AU, 3)
                        for x in tg}
                row[label] = ent
            out[name][era] = row
    return out


def earth_validation(sta, earth):
    e_era = earth[era_mask(earth, "2007-01-01", "2026-12-01")]
    d_t, d_b, n_unmatched = [], [], 0
    graze = []
    for name, (direction, side) in COMBOS.items():
        s = sta[combo_mask(sta, direction, side)]
        e = e_era[combo_mask(e_era, direction, side)]
        for tid in sorted(set(np.asarray(s["target_id"]))):
            st = s[np.asarray(s["target_id"]) == tid]
            et = e[np.asarray(e["target_id"]) == tid]
            for row in st:
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
                if row["b_min_au"] <= 0.1 and name.startswith("S1"):
                    graze.append((tid, float(dt[j]), float(row["b_min_au"] / RSUN_AU),
                                  float(et["b_min_au"][j] / RSUN_AU)))
    d_t, d_b = np.array(d_t), np.abs(np.array(d_b))
    return {
        "matched": int(len(d_t)), "unmatched": n_unmatched,
        "delta_tca_days": {"median_signed": float(np.median(d_t)),
                           "abs_median": float(np.median(np.abs(d_t))),
                           "abs_max": float(np.max(np.abs(d_t)))},
        "abs_delta_b": {"median_au": float(np.median(d_b)),
                        "p95_au": float(np.percentile(d_b, 95)),
                        "max_au": float(d_b.max())},
        "note": "Earth-center events lead/lag by the spacecraft's heliocentric "
                "longitude offset (weeks-months); the STEREO-A list is mandatory",
        "S1_0.1AU_examples": graze[:12],
    }


def wcs_validation():
    """Compare hi_geometry.helioprojective with the test-frame header."""
    from astropy.io import fits
    from astropy.wcs import WCS
    warnings.simplefilter("ignore")
    hd = fits.getheader(TEST_FRAME)
    wa, wp = WCS(hd, key="A"), WCS(hd)
    mjd = Time(hd["DATE-AVG"]).mjd
    cat = np.load(REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz")
    px, py = wa.all_world2pix(cat["ra"], cat["dec"], 0)
    ok = np.isfinite(px) & (px > 20) & (px < 1003) & (py > 20) & (py < 1003)
    idx = np.nonzero(ok)[0][:: max(1, ok.sum() // 40)]
    hpln_w, hplt_w = wp.all_pix2world(px[idx], py[idx], 0)
    s = G.radec_to_vec(cat["ra"][idx], cat["dec"][idx])
    hpln_g, hplt_g = zip(*[G.helioprojective(v, mjd) for v in s])
    hpln_g, hplt_g = np.array(hpln_g).ravel(), np.array(hplt_g).ravel()
    return {"n_stars": int(len(idx)), "test_frame": TEST_FRAME.name,
            "hpln_resid_deg": {"median": float(np.median(hpln_g - hpln_w)),
                               "max_abs": float(np.max(np.abs(hpln_g - hpln_w)))},
            "hplt_resid_deg": {"median": float(np.median(hplt_g - hplt_w)),
                               "max_abs": float(np.max(np.abs(hplt_g - hplt_w)))},
            "fov_center_hpln_hplt": [float(hd["CRVAL1"]), float(hd["CRVAL2"])],
            "note": "residuals vs the header's own HPLN/HPLT include the header's "
                    "ephemeris/attitude vs ours; 1 px = 0.020 deg"}


def visibility(sta):
    """Per-event HI-1 visibility for the 0.1 AU rungs of both channels."""
    rows, agg = [], {}
    for name, (direction, side) in COMBOS.items():
        ch = name[:2]
        sub = sta[combo_mask(sta, direction, side) & (np.asarray(sta["b_min_au"]) <= 0.1)]
        for ev in sub:
            e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
                 "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
            mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
            half = G.window_half_days(float(ev["b_min_au"]), 0.1, float(ev["v_perp_km_s"]))
            mjds = mjd_ca + np.arange(-half, half + 1e-9, 10 / 1440.0)
            b_e = G.impact_parameter_au(e, mjds)
            s_hat = G.source_direction(e, ch)
            hpln, hplt = G.helioprojective(s_hat, mjds)
            eps = G.elongation_deg(s_hat, mjds)
            vis = G.in_fov_era(hpln, hplt, mjds) & (b_e <= 0.1)
            n_vis = int(vis.sum())
            row = {"event_id": str(ev["event_id"]), "target_id": str(ev["target_id"]),
                   "channel": ch, "t_ca_utc": str(ev["t_ca_utc"])[:16],
                   "mjd_ca": round(float(mjd_ca), 4),
                   "b_min_rsun": round(float(ev["b_min_au"]) / RSUN_AU, 3),
                   "window_days": round(2 * half, 2),
                   "b_e_min_check_rsun": round(float(b_e.min()) / RSUN_AU, 3),
                   "visible_hours": round(n_vis * 10 / 60.0, 1),
                   "visible_frames_40min": int(round(n_vis * 10 / CADENCE_MIN)),
                   "hplt_at_min_eps": round(float(hplt[np.argmin(eps)]), 2)}
            if n_vis:
                vm = mjds[vis]
                row.update({"side": "post_tca" if vm.min() > mjd_ca else "pre_tca",
                            "vis_start_rel_tca_d": round(float(vm.min() - mjd_ca), 2),
                            "vis_end_rel_tca_d": round(float(vm.max() - mjd_ca), 2),
                            "eps_range_deg": [round(float(eps[vis].min()), 2),
                                              round(float(eps[vis].max()), 2)],
                            "hpln_range_deg": [round(float(hpln[vis].min()), 2),
                                               round(float(hpln[vis].max()), 2)]})
            rows.append(row)
            a = agg.setdefault(ch, {}).setdefault(str(ev["target_id"]), {
                "events": 0, "visible_events": 0, "visible_hours": []})
            a["events"] += 1
            a["visible_events"] += n_vis > 0
            if n_vis:
                a["visible_hours"].append(row["visible_hours"])
    for ch in agg:
        for tid, a in agg[ch].items():
            h = a.pop("visible_hours")
            a["median_visible_hours"] = float(np.median(h)) if h else 0.0
    return rows, agg


def s2_1au_ledger(sta):
    """HI-1 transits of each target star (S2 1 AU rung: the whole
    transit is in-beam), from the star's helioprojective track."""
    out = {}
    sub = sta[combo_mask(sta, "inbound", "anti_target")]
    for ev in sub:
        e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
             "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
        mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
        mjds = mjd_ca + np.arange(-40, 40, 0.25)
        hpln, hplt = G.helioprojective(G.source_direction(e, "S2"), mjds)
        vis = G.in_fov_era(hpln, hplt, mjds)
        a = out.setdefault(str(ev["target_id"]), {"conjunctions": 0, "hi1_transits": 0,
                                                   "transit_days": []})
        a["conjunctions"] += 1
        if vis.any():
            a["hi1_transits"] += 1
            a["transit_days"].append(round(float(vis.sum() * 0.25), 1))
    for a in out.values():
        d = a.pop("transit_days")
        a["median_transit_days"] = float(np.median(d)) if d else 0.0
    n_t = sum(1 for a in out.values() if a["hi1_transits"])
    return {"targets_with_transits": n_t,
            "total_transits": sum(a["hi1_transits"] for a in out.values()),
            "per_target": out}


def main():
    sta, earth = Table.read(STA), Table.read(EARTH)
    assert set(np.unique(sta["validity"])) <= {"valid", "degraded"}
    rows, agg = visibility(sta)
    out = {
        "study": "stereoa_census_v1",
        "input": str(STA.relative_to(REPO)),
        "crossings_id": str(sta.meta.get("crossings_id", "")),
        "n_events": len(sta),
        "sunward_census": census(sta),
        "earth_validation": earth_validation(sta, earth),
        "helioprojective_validation": wcs_validation(),
        "hi1_footprint_hpln_hplt_deg": {"hpln": G.HI1_HPLN, "hplt": G.HI1_HPLT,
                                        "edge_margin_deg": G.HI1_EDGE_MARGIN_DEG,
                                        "pointing_eras_mjd_hpln": G.POINTING_ERAS},
        "visibility_0.1AU": {"aggregates": agg, "rows": rows},
        "s2_1AU_transit_ledger": s2_1au_ledger(sta),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["earth_validation"], indent=1)[:1500])
    print(json.dumps(out["helioprojective_validation"], indent=1))
    print(json.dumps(out["sunward_census"]["S1_downlink_postlens"]["full_2007_2026"], indent=1))
    print(json.dumps(agg, indent=1))
    led = out["s2_1AU_transit_ledger"]
    print("S2 1AU ledger:", led["targets_with_transits"], "targets,", led["total_transits"], "transits")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
