"""Pool the per-window completeness bins across windows per
(unit, plate-limit stratum) — shallow plates lack bright-truth
statistics individually (n < 5 per bin), but pooled strata resolve
the curve. Pure arithmetic on completeness_v1.json; no new pulls.

Adds per-unit stratified m90 (with the APASS truth-depth cap noted:
the refcat runs out near B 16, so pooled m90 saturates at the 15.5
bin — statements there are 'm90 >= 15.5 (truth-limited)').

Products: results/completeness_pooled_v1.json
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

D = Path(__file__).resolve().parents[1]
STRATA = [(0.0, 13.0, "lim<13"), (13.0, 15.0, "lim13-15"),
          (15.0, 99.0, "lim>=15")]
TRUTH_CAP_BIN = 15.5


def wilson_lo(k, n, z=1.6449):
    if n == 0:
        return 0.0
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h)


def main() -> int:
    src = json.loads((D / "results" / "completeness_v1.json").read_text())
    out = {"generated_utc": datetime.now(timezone.utc).isoformat(
               timespec="seconds"),
           "method": "pooled k/n across windows per plate-limit "
                     "stratum; m90 = faintest bin center with pooled "
                     "eff >= 0.9 walking from the bright end "
                     "(>= 10 pooled truth stars per bin); truth cap "
                     f"at the {TRUTH_CAP_BIN} bin (APASS depth)",
           "units": {}}
    for unit, u in src["units"].items():
        strata = {}
        for lo, hi, name in STRATA:
            pool = defaultdict(lambda: [0, 0])
            n_exp = 0
            for w in u["windows"]:
                try:
                    lim = float(w.get("lim_mag") or 0)
                except (TypeError, ValueError):
                    continue
                if not (lo <= lim < hi):
                    continue
                n_exp += 1
                for c, b in w["bins"].items():
                    pool[float(c)][0] += b["k"]
                    pool[float(c)][1] += b["n"]
            m90 = None
            for c in sorted(pool):
                k, n = pool[c]
                if n < 10:
                    continue
                if k / n >= 0.9:
                    m90 = c
                else:
                    break
            strata[name] = dict(
                n_exposures=n_exp,
                m90=m90,
                truth_limited=bool(m90 is not None
                                   and m90 >= TRUTH_CAP_BIN),
                bins={str(c): dict(k=pool[c][0], n=pool[c][1],
                                   eff=round(pool[c][0] / pool[c][1], 3),
                                   wilson_lo=round(
                                       wilson_lo(*pool[c]), 3))
                      for c in sorted(pool) if pool[c][1] > 0})
        out["units"][unit] = strata
        show = {s: (v["m90"], v["n_exposures"]) for s, v in
                strata.items() if v["n_exposures"]}
        print(unit, show)
    (D / "results" / "completeness_pooled_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
