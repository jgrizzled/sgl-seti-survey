"""Rubin DP2 recon: geometry-only intersect of the DP2 visit snapshot
with (A) the universal corridor list and (B) the universal crossing
list. No pixel or catalog-row data touched — boresight cones only.

Inputs:
  runs/rubin/recon/dp2_visit_snapshot.ecsv   (recon_fetch_visits.py)
  targets/universal_v2.json                  (anti_ra/anti_dec per system)
  crossings/universal_v1/events.ecsv
Output:
  surveys/rubin/results/recon_intersect_v0.json

Conventions copied from surveys/atlas-asassn-crossings/scripts/era_scope.py:
channel A = link_direction 'inbound' (position = star), channel B =
'outbound' (position = relay/antipode); rung window half-width =
beam_radius / v_perp per event (reproduces the frozen +/-0.35 d and
+/-5.8 d values at v_perp ~30 km/s).
"""
import collections
import json
import pathlib

import numpy as np
from astropy.table import Table

REPO = pathlib.Path(__file__).resolve().parents[3]
VISITS = REPO / "runs" / "rubin" / "recon" / "dp2_visit_snapshot.ecsv"
TARGETS = REPO / "targets" / "universal_v2.json"
EVENTS = REPO / "crossings" / "universal_v1" / "events.ecsv"
OUT = REPO / "surveys" / "rubin" / "results" / "recon_intersect_v0.json"

# LSSTCam FOV is ~3.5 deg diameter; center-within-1.5 deg is the
# conservative discovery cone (DECam recon used the same construction).
CONE_DEG = 1.5
AU_KM = 149_597_870.7
RSUN_KM = 695_700.0
RUNGS = [
    ("B", "1.2Rsun", 1.2 * RSUN_KM / AU_KM),
    ("B", "2.5Rsun", 2.5 * RSUN_KM / AU_KM),
    ("B", "0.1AU", 0.1),
    ("A", "0.1AU", 0.1),
    ("A", "1.0AU", 1.0),
]


def sep_deg(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = (np.sin(d1) * np.sin(d2)
         + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def main():
    v = Table.read(VISITS)
    vra = np.array(v["ra"], float)
    vdec = np.array(v["dec"], float)
    vmjd = np.array(v["expMidptMJD"], float)
    vband = np.array(v["band"], str)
    era = (float(vmjd.min()), float(vmjd.max()))

    out = {
        "visits_input": str(VISITS.relative_to(REPO)),
        "n_visits": len(v),
        "era_mjd": era,
        "cone_deg": CONE_DEG,
    }

    # --- Pipeline A: corridor (antipode) coverage ---
    systems = json.load(open(TARGETS))["systems"]
    rows = {}
    for s in systems:
        m = sep_deg(vra, vdec, s["anti_ra"], s["anti_dec"]) < CONE_DEG
        if not m.any():
            rows[s["system"]] = {"n": 0}
            continue
        months = sorted(set((int(x) - 60676) // 30 for x in vmjd[m]))
        rows[s["system"]] = {
            "n": int(m.sum()),
            "anti_ra": s["anti_ra"], "anti_dec": s["anti_dec"],
            "bands": dict(sorted(collections.Counter(vband[m]).items())),
            "mjd_span": [round(float(vmjd[m].min()), 1),
                         round(float(vmjd[m].max()), 1)],
            "n_distinct_30d_bins": len(months),
        }
    covered = {k: r for k, r in rows.items() if r["n"] > 0}
    out["pipeline_a"] = {
        "n_systems": len(systems),
        "n_covered": len(covered),
        "covered": dict(sorted(covered.items(),
                               key=lambda kv: -kv[1]["n"])),
        "uncovered": sorted(k for k, r in rows.items() if r["n"] == 0),
    }

    # --- Pipeline B: crossing windows in the DP2 era ---
    e = Table.read(EVENTS)
    t_ca = np.array(e["t_ca_tdb_jd"], float) - 2400000.5
    is_a = np.array(e["link_direction"]) == "inbound"
    b_au = np.array(e["b_min_au"], float)
    vperp = np.array(e["v_perp_km_s"], float)
    pos_ra = np.where(is_a, e["star_icrs_ra_deg"], e["relay_icrs_ra_deg"])
    pos_dec = np.where(is_a, e["star_icrs_dec_deg"], e["relay_icrs_dec_deg"])
    tgt = np.array(e["target_id"], str)

    chan = {}
    for ch, rung, lim_au in RUNGS:
        hw = lim_au * AU_KM / vperp / 86400.0  # window half-width, days
        m = ((is_a if ch == "A" else ~is_a) & (b_au <= lim_au)
             & (t_ca + hw >= era[0]) & (t_ca - hw <= era[1]))
        idx = np.flatnonzero(m)
        per_event = []
        for i in idx:
            inwin = ((np.abs(vmjd - t_ca[i]) <= hw[i])
                     & (sep_deg(vra, vdec, float(pos_ra[i]),
                                float(pos_dec[i])) < CONE_DEG))
            if inwin.any():
                per_event.append({
                    "target": tgt[i],
                    "t_ca_mjd": round(float(t_ca[i]), 2),
                    "b_min_rsun": round(float(e["b_min_solar_radii"][i]), 3),
                    "b_min_au": round(float(b_au[i]), 4),
                    "half_width_d": round(float(hw[i]), 2),
                    "n_visits": int(inwin.sum()),
                    "bands": dict(sorted(
                        collections.Counter(vband[inwin]).items())),
                    "pos": [round(float(pos_ra[i]), 3),
                            round(float(pos_dec[i]), 3)],
                })
        chan[f"{ch}_{rung}"] = {
            "in_era_events": int(m.sum()),
            "in_era_targets": len(set(tgt[m])),
            "events_with_visits": len(per_event),
            "per_event": sorted(per_event, key=lambda r: r["b_min_au"]),
        }
    out["pipeline_b"] = chan

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))
    print("wrote", OUT.relative_to(REPO))
    print(f"era MJD {era[0]:.1f}-{era[1]:.1f}")
    print("pipeline A covered:", len(covered), "/", len(systems))
    for k, c in chan.items():
        print(f"{k}: {c['in_era_events']} in-era events, "
              f"{c['events_with_visits']} with in-window in-cone visits")


if __name__ == "__main__":
    main()
