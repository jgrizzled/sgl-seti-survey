"""Pre-freeze elongation gate (hypotheses v1.0 scope table).

WISE observes only at solar elongation ~90 deg; every crossing window
puts the beam source near elongation 0/180 unless the geometry is
wide-rung and/or high-ecliptic-latitude. For every wise_v1 in-era
event and ladder rung: does any in-window epoch (hibernation excluded)
put the source within |elongation - 90| < 4 deg? Writes
configs/elongation_gate_v1.json (the scope input of the freeze and of
coverage_intersect.py).
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

KM_PER_AU = 1.495978707e8
ERA = (55203.0, 60524.0)
HIBERNATION = (55593.0, 56639.0)
TOL_DEG = 4.0
LADDER = {"A": (0.1, 1.0), "B": (0.0056, 0.0116, 0.1)}
OUT = REPO / "surveys" / "wise-crossings" / "configs"


def elong(ra_deg, dec_deg, mjd):
    tt = Time(mjd, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", tt).xyz.to_value("AU")
    earth = get_body_barycentric("earth", tt).xyz.to_value("AU")
    s = sun - earth
    s /= np.linalg.norm(s)
    ra, de = np.radians(ra_deg), np.radians(dec_deg)
    u = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra),
                  np.sin(de)])
    return float(np.degrees(np.arccos(np.clip(np.dot(s, u), -1, 1))))


def main():
    t = Table.read(REPO / "crossings" / "wise_v1" / "events.ecsv")
    t["mjd"] = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    side = np.asarray(t["axis_distance_au"])
    era = (t["mjd"] >= ERA[0]) & (t["mjd"] <= ERA[1])
    chans = {"A": t[era & (t["link_direction"] == "inbound") & (side > 0)],
             "B": t[era & (t["link_direction"] == "outbound") & (side < 0)]}
    gate = {"input": "crossings/wise_v1", "tolerance_deg": TOL_DEG,
            "hibernation_mjd": list(HIBERNATION), "rungs": {},
            "viable_targets": {}, "viable_events": []}
    for ch, tab in chans.items():
        for r in LADDER[ch]:
            sub = tab[np.asarray(tab["b_min_au"]) < r]
            n_ok, tgts = 0, Counter()
            for ev in sub:
                b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
                half = np.sqrt(r * r - b * b) * KM_PER_AU / vp / 86400.0
                ts = np.linspace(ev["mjd"] - half, ev["mjd"] + half,
                                 max(int(half) + 2, 5))
                ts = [m for m in ts
                      if not (HIBERNATION[0] < m < HIBERNATION[1])]
                if ch == "A":
                    ra, de = (float(ev["star_icrs_ra_deg"]),
                              float(ev["star_icrs_dec_deg"]))
                    ok = any(abs(elong(ra, de, m) - 90.0) < TOL_DEG
                             for m in ts)
                else:
                    ok = any(abs(elong(float(ev["relay_icrs_ra_deg"]),
                                       float(ev["relay_icrs_dec_deg"]), m)
                                 - 90.0) < TOL_DEG for m in ts)
                if ok:
                    n_ok += 1
                    tgts[str(ev["target_id"])] += 1
                    if ch == "A" and r == 1.0:
                        gate["viable_events"].append(str(ev["event_id"]))
            gate["rungs"][f"{ch}_{r}"] = {"n_events": len(sub),
                                          "n_viable": n_ok}
            if ch == "A" and r == 1.0:
                gate["viable_targets"] = dict(sorted(tgts.items()))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "elongation_gate_v1.json").write_text(
        json.dumps(gate, indent=1) + "\n")
    print(json.dumps({k: gate[k] for k in ("rungs",)}, indent=1))
    print("viable targets:", len(gate["viable_targets"]))


if __name__ == "__main__":
    main()
