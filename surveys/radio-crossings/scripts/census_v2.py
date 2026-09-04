"""Census of the v2 radio coverage products (plan §5.8 item 7).

Reads results/askap_inwindow_v2.ecsv, lotss_inwindow_v2.ecsv and the
per-target-channel JSONs and prints the report tables: per rung, how
many in-era events have a window inside each archive's date range, how
many of those have >= 1 observation (image / footprint / validated), and
the narrow-rung row listing. Pure bookkeeping — no archive contact.
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from coverage_intersect_v2 import (LADDER, LOTSS_DEC_MIN, OUT, VALIDATED,  # noqa: E402
                                   load_events, pos_of, windows_of)

ASKAP_RANGE = (58594.17, 61280.0)       # first RACS-low SBID .. obscore max
LOTSS_RANGE = (56800.6, 60540.5)        # DR3 pointing observations


def main():
    events = load_events()
    askap = Table.read(OUT / "askap_inwindow_v2.ecsv")
    lotss = Table.read(OUT / "lotss_inwindow_v2.ecsv")
    askap_tc = json.loads((OUT / "askap_targets_v2.json").read_text())
    lotss_tc = json.loads((OUT / "lotss_targets_v2.json").read_text())

    def key(r):
        return (str(r["event_id"]), float(r["radius_au"]))

    a_img = defaultdict(set)
    a_foot = defaultdict(set)
    a_val = defaultdict(set)
    for r in askap:
        a_img[key(r)].add(str(r["sbid"]))
        if r["in_footprint"]:
            a_foot[key(r)].add(str(r["sbid"]))
            if str(r["quality"]) in VALIDATED:
                a_val[key(r)].add(str(r["sbid"]))
    l_img = defaultdict(set)
    l_beam = defaultdict(set)
    for r in lotss:
        l_img[key(r)].add(float(r["obs_mid_mjd"]))
        if r["in_beam"]:
            l_beam[key(r)].add(float(r["obs_mid_mjd"]))

    lotss_has_pointing = {(d["channel"], d["target_id"]): d["n_pointings_in_image"] > 0
                          for d in lotss_tc}

    print("## Window census per rung (events with the whole window inside the archive date range)")
    print("| arm | rung | events in range | targets | >=1 obs in image | in footprint/beam | validated | events with any obs |")
    print("|---|---|---|---|---|---|---|---|")
    out = {}
    for ch in LADDER:
        for r_au in LADDER[ch]:
            for arm, rng in (("askap", ASKAP_RANGE), ("lotss", LOTSS_RANGE)):
                n_ev, tg = 0, set()
                n_img = n_foot = n_val = 0
                for ev in events[ch]:
                    if arm == "lotss":
                        if pos_of(ch, ev)[1] < LOTSS_DEC_MIN:
                            continue
                        if not lotss_has_pointing.get((ch, str(ev["target_id"])), False):
                            continue
                    for rr, hd, lo, hi in windows_of(ch, ev):
                        if rr != r_au:
                            continue
                        if lo < rng[0] or hi > rng[1]:
                            continue
                        n_ev += 1
                        tg.add(str(ev["target_id"]))
                        k = (str(ev["event_id"]), r_au)
                        if arm == "askap":
                            n_img += k in a_img
                            n_foot += k in a_foot
                            n_val += k in a_val
                        else:
                            n_img += k in l_img
                            n_foot += k in l_beam
                            n_val += k in l_beam
                out[(arm, ch, r_au)] = (n_ev, len(tg), n_img, n_foot, n_val)
                print(f"| {arm} | {ch} {r_au} AU | {n_ev} | {len(tg)} | {n_img} | "
                      f"{n_foot} | {n_val} | {n_img} |")

    print("\n## Narrow-rung rows (radius < 1 AU), one line per (event, SBID/pointing)")
    seen = set()
    for r in askap:
        if float(r["radius_au"]) >= 1.0:
            continue
        k = (str(r["event_id"]), str(r["sbid"]))
        if k in seen:
            continue
        seen.add(k)
        print(f"ASKAP | {r['channel']} | {r['target_id']} | {r['radius_au']} | "
              f"b={r['b_min_au']*215.032:.3f} Rsun | t_ca {r['t_ca_mjd']:.2f} | "
              f"{r['collection']} {r['sbid']} {r['field']} | {r['freq_mhz']} MHz | "
              f"offset {r['offset_days']:+.2f} d of ±{r['window_days']/2:.2f} | "
              f"sep {r['sep_deg']:.2f}/{r['footprint_half_deg']:.2f} deg "
              f"in_footprint={bool(r['in_footprint'])} | {r['quality']}")
    for r in lotss:
        if float(r["radius_au"]) >= 1.0:
            continue
        print(f"LoTSS | {r['channel']} | {r['target_id']} | {r['radius_au']} | "
              f"b={r['b_min_au']*215.032:.3f} Rsun | t_ca {r['t_ca_mjd']:.2f} | "
              f"{r['pointing']} | offset {r['offset_days']:+.2f} d of "
              f"±{r['window_days']/2:.2f} | sep {r['sep_deg']:.2f} deg "
              f"in_beam={bool(r['in_beam'])}")

    print("\n## Visit lists (observations over the position, any date)")
    print("ASKAP per target-channel n_obs (dated): ",
          {f"{d['channel']}/{d['target_id']}": d["n_obs_dated"]
           for d in askap_tc if d.get("status") == "ok"})
    print("ASKAP collections over antipodes (B):",
          Counter(c for d in askap_tc if d["channel"] == "B" and d.get("status") == "ok"
                  for c in d["collections"]))
    print("LoTSS per target-channel obs in image:",
          {f"{d['channel']}/{d['target_id']}": d["n_obs_in_image"]
           for d in lotss_tc if d["n_obs_in_image"]})
    print("\n## Wide-rung (A 1.0 AU) collection census, in-footprint rows:",
          Counter(str(r["collection"]) for r in askap
                  if float(r["radius_au"]) == 1.0 and r["in_footprint"]))
    print("Grazing-rung (B 1.2/2.5 Rsun) in-era events:",
          {r: sum(1 for ev in events['B'] for rr, *_ in windows_of('B', ev) if rr == r)
           for r in LADDER['B'][:2]})


if __name__ == "__main__":
    main()
