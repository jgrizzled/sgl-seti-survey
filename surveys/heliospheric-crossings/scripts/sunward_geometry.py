"""Sunward-channel geometry study (project plan §5.8 item 9).

Question frozen at the ZTF crossings freeze (2026-08-23): of the four
(link direction x side-of-axis) combinations, two arrive from the solar
direction and were declared out of scope. Are any of those sunward
combinations *only* visible at small solar elongation? If yes,
heliospheric imagers (LASCO, STEREO-HI, WISPR) are the unique archival
substrate; if no, the question closes for the record.

Combination taxonomy over `crossings/universal_v1/events.ecsv`
(side == 'target' <=> axis_distance_au > 0, Earth on the star side):

  A  (inbound,  target)      uplink pre-Sun    -> star near opposition (observable)
  B  (outbound, anti_target) downlink pre-lens -> relay near anti-Sun  (observable)
  S1 (outbound, target)      downlink POST-LENS: beam has grazed the Sun;
                             apparent source = limb graze point, elongation
                             ~ b_graze (16.0'/Rsun at 1 AU) regardless of
                             Earth's in-beam offset
  S2 (inbound,  anti_target) uplink PAST the Sun: apparent source = the
                             star itself at elongation ~ arcsin(b/1 AU),
                             occulted for b < ~1 Rsun

This script (1) validates that taxonomy against the frozen list by
computing the actual solar elongation of the apparent source at every
t_ca, and (2) counts sunward events/targets per beam-radius rung
(frozen ZTF ladder) and per imager era.  Earth-center observer: exact
for ground + LEO, indicative for L1 (SOHO halo cross-track ~1 Rsun --
a spacecraft-observer list is a recon-stage item), *not* valid for
STEREO/PSP (drifting / inner-heliosphere observers).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord, get_sun
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
EVENTS = REPO / "crossings" / "universal_v1" / "events.ecsv"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "sunward_geometry_v1.json"

RSUN_AU = 0.00465047  # IAU solar radius in AU

# Frozen ZTF beam-radius ladder (hypotheses.md v1.0), reused unchanged.
RUNGS = {
    "S1_downlink_postlens": [("1.2Rsun", 1.2 * RSUN_AU), ("2.5Rsun", 2.5 * RSUN_AU), ("0.1AU", 0.1)],
    "S2_uplink_pastsun": [("0.1AU", 0.1), ("1AU", 1.0)],
}

# Imager eras (from general knowledge -- recon-stage items to verify; the
# universal list itself spans 1980-01-01 -> 2028-01-01).
ERAS = {
    "full_list_1980_2028": (None, None),
    "lasco_1996_on": ("1996-01-11", None),
    "stereo_hi_2007_on": ("2007-04-01", None),
    "wispr_2018_on": ("2018-11-01", None),
}


def combo_mask(t: Table, direction: str, side: str) -> np.ndarray:
    return (np.asarray(t["link_direction"]) == direction) & (np.asarray(t["side"]) == side)


def era_mask(t: Table, start: str | None, end: str | None) -> np.ndarray:
    jd = np.asarray(t["t_ca_tdb_jd"])
    m = np.ones(len(t), dtype=bool)
    if start:
        m &= jd >= Time(start).tdb.jd
    if end:
        m &= jd <= Time(end).tdb.jd
    return m


def main() -> None:
    t = Table.read(EVENTS)
    # degraded = window-boundary minima (kept per the universal_v1 rule)
    assert set(np.unique(t["validity"])) <= {"valid", "degraded"}, "unexpected invalid rows"

    # Side <-> sign consistency check.
    ax = np.asarray(t["axis_distance_au"])
    side = np.asarray(t["side"])
    assert np.all((side == "target") == (ax > 0)), "side/sign mapping broken"

    combos = {
        "A": combo_mask(t, "inbound", "target"),
        "B": combo_mask(t, "outbound", "anti_target"),
        "S1_downlink_postlens": combo_mask(t, "outbound", "target"),
        "S2_uplink_pastsun": combo_mask(t, "inbound", "anti_target"),
    }
    assert sum(int(m.sum()) for m in combos.values()) == len(t)

    # --- elongation validation over the whole list ---------------------
    # Apparent-source direction per combo: the star for the uplink combos
    # (A, S2), the relay/antipode for the downlink combos (B, S1).  For
    # S1 the physical apparent source is the limb graze point (analytic,
    # ~16'/Rsun); the relay elongation ~0 deg confirms arrival is sunward.
    times = Time(np.asarray(t["t_ca_tdb_jd"]), format="jd", scale="tdb")
    sun = get_sun(times)
    star = SkyCoord(np.asarray(t["star_icrs_ra_deg"]), np.asarray(t["star_icrs_dec_deg"]), unit="deg")
    relay = SkyCoord(np.asarray(t["relay_icrs_ra_deg"]), np.asarray(t["relay_icrs_dec_deg"]), unit="deg")
    elong_star = sun.separation(star).deg
    elong_relay = sun.separation(relay).deg
    b = np.asarray(t["b_min_au"])
    pred_small = np.degrees(np.arcsin(np.clip(b, 0, 1)))  # elongation if sunward

    validation = {}
    for name, src_elong, sunward in [
        ("A", elong_star, False),
        ("B", elong_relay, False),
        ("S1_downlink_postlens", elong_relay, True),
        ("S2_uplink_pastsun", elong_star, True),
    ]:
        m = combos[name]
        e = src_elong[m]
        expected = pred_small[m] if sunward else 180.0 - pred_small[m]
        resid = e - expected
        validation[name] = {
            "n": int(m.sum()),
            "elongation_deg": {
                "min": float(e.min()),
                "median": float(np.median(e)),
                "max": float(e.max()),
            },
            "residual_vs_analytic_deg": {
                "median_abs": float(np.median(np.abs(resid))),
                "max_abs": float(np.max(np.abs(resid))),
            },
        }

    # --- sunward rung x era census -------------------------------------
    census = {}
    for name, rungs in RUNGS.items():
        m0 = combos[name]
        census[name] = {}
        for era, (start, end) in ERAS.items():
            me = m0 & era_mask(t, start, end)
            row = {}
            for label, rmax in rungs:
                mr = me & (b <= rmax)
                tgts = sorted(set(np.asarray(t["target_id"])[mr]))
                entry = {"events": int(mr.sum()), "targets": len(tgts)}
                if mr.sum() and mr.sum() <= 60:
                    entry["target_ids"] = tgts
                row[label] = entry
            census[name][era] = row
    # S2 occultation core: closest approach behind the photosphere.
    m_occ = combos["S2_uplink_pastsun"] & (b < RSUN_AU)
    census["S2_uplink_pastsun"]["occulted_at_tca_full_list"] = {
        "criterion_b_lt_1Rsun": int(m_occ.sum()),
        "targets": len(set(np.asarray(t["target_id"])[m_occ])),
    }

    out = {
        "study": "sunward_geometry_v1",
        "events_input": str(EVENTS.relative_to(REPO)),
        "crossings_id": str(t.meta.get("crossings_id", "")),
        "combo_counts": {k: int(v.sum()) for k, v in combos.items()},
        "elongation_validation": validation,
        "sunward_census": census,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    json.dump(out["elongation_validation"], sys.stdout, indent=2)
    print()
    print(json.dumps(out["combo_counts"]))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
