"""Sunward-channel census on the SOHO-observer crossing list, plus the
Earth-center validation (TESS pattern).

Consumes `crossings/soho_v1` (Horizons -21 spacecraft observer,
1996-01-01 -> 2026-10-01) and re-runs the sunward rung x era census
from the geometry study on the corrected impact parameters. SOHO's
halo orbit shifts b by up to ~0.9 R_sun vs Earth center, so grazing
family membership can change — the validation matches each SOHO event
to its Earth-center counterpart (same target, direction, side; nearest
t_ca) and reports |Δb| / |Δt_ca| plus the per-rung membership churn.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
SOHO = REPO / "crossings" / "soho_v1" / "events.ecsv"
EARTH = REPO / "crossings" / "universal_v1" / "events.ecsv"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "soho_census_v1.json"

RSUN_AU = 0.00465047

RUNGS = {
    "S1_downlink_postlens": [("1.2Rsun", 1.2 * RSUN_AU), ("2.5Rsun", 2.5 * RSUN_AU), ("0.1AU", 0.1)],
    "S2_uplink_pastsun": [("0.1AU", 0.1), ("1AU", 1.0)],
}
COMBOS = {
    "S1_downlink_postlens": ("outbound", "target"),
    "S2_uplink_pastsun": ("inbound", "anti_target"),
}
# archive-era sub-windows of the SOHO list span
ERAS = {
    "full_soho_1996_2026": (None, None),
    "level1_era_to_2017-08": (None, "2017-08-31"),
    "wispr_era_2018_on": ("2018-11-01", None),
}


def combo_mask(t: Table, direction: str, side: str) -> np.ndarray:
    return (np.asarray(t["link_direction"]) == direction) & (np.asarray(t["side"]) == side)


def era_mask(t: Table, start, end) -> np.ndarray:
    jd = np.asarray(t["t_ca_tdb_jd"])
    m = np.ones(len(t), dtype=bool)
    if start:
        m &= jd >= Time(start).tdb.jd
    if end:
        m &= jd <= Time(end).tdb.jd
    return m


def census(t: Table) -> dict:
    b = np.asarray(t["b_min_au"])
    out = {}
    for name, rungs in RUNGS.items():
        m0 = combo_mask(t, *COMBOS[name])
        out[name] = {}
        for era, (start, end) in ERAS.items():
            me = m0 & era_mask(t, start, end)
            row = {}
            for label, rmax in rungs:
                mr = me & (b <= rmax)
                tgts = sorted(set(np.asarray(t["target_id"])[mr]))
                entry = {"events": int(mr.sum()), "targets": len(tgts)}
                if 0 < mr.sum() <= 60:
                    entry["target_ids"] = tgts
                row[label] = entry
            out[name][era] = row
    return out


def match_validation(soho: Table, earth: Table) -> dict:
    """Match SOHO events to Earth-center events within the SOHO era."""
    e_era = earth[era_mask(earth, "1996-01-01", "2026-10-01")]
    out = {}
    d_b_all, d_t_all = [], []
    churn = {}
    for name, (direction, side) in COMBOS.items():
        s = soho[combo_mask(soho, direction, side)]
        e = e_era[combo_mask(e_era, direction, side)]
        for tid in sorted(set(np.asarray(s["target_id"]))):
            st = s[np.asarray(s["target_id"]) == tid]
            et = e[np.asarray(e["target_id"]) == tid]
            if len(et) == 0:
                continue
            for row in st:
                dt = np.abs(np.asarray(et["t_ca_tdb_jd"]) - row["t_ca_tdb_jd"])
                j = int(np.argmin(dt))
                if dt[j] > 5.0:
                    continue
                d_b_all.append(row["b_min_au"] - et["b_min_au"][j])
                d_t_all.append(dt[j])
        # per-rung membership churn (grazing rungs only, S1)
        if name == "S1_downlink_postlens":
            for label, rmax in RUNGS[name][:2]:
                sm = np.asarray(s["b_min_au"]) <= rmax
                em = np.asarray(e["b_min_au"]) <= rmax
                churn[label] = {
                    "soho_events": int(sm.sum()),
                    "earth_events_same_era": int(em.sum()),
                    "soho_targets": sorted(set(np.asarray(s["target_id"])[sm])),
                    "earth_targets_same_era": sorted(set(np.asarray(e["target_id"])[em])),
                }
    d_b = np.abs(np.array(d_b_all))
    out["matched_events"] = len(d_b_all)
    out["abs_delta_b"] = {
        "median_rsun": float(np.median(d_b) / RSUN_AU),
        "p95_rsun": float(np.percentile(d_b, 95) / RSUN_AU),
        "max_rsun": float(d_b.max() / RSUN_AU),
    }
    out["abs_delta_tca_days"] = {
        "median": float(np.median(d_t_all)),
        "max": float(np.max(d_t_all)),
    }
    out["grazing_rung_churn"] = churn
    return out


def main() -> None:
    soho = Table.read(SOHO)
    earth = Table.read(EARTH)
    assert set(np.unique(soho["validity"])) <= {"valid", "degraded"}

    out = {
        "study": "soho_census_v1",
        "soho_input": str(SOHO.relative_to(REPO)),
        "soho_crossings_id": str(soho.meta.get("crossings_id", "")),
        "earth_input": str(EARTH.relative_to(REPO)),
        "n_soho_events": len(soho),
        "sunward_census_soho": census(soho),
        "earth_validation": match_validation(soho, earth),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out["earth_validation"]["abs_delta_b"], indent=1))
    print(json.dumps(out["sunward_census_soho"]["S1_downlink_postlens"]["full_soho_1996_2026"], indent=1))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
