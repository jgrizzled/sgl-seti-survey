"""Channel-A saturation cut against TIC (hypotheses.md freeze v1.0).

Frozen rule: excluded T < 7.3 (single-cadence saturation ~6.8 + 0.5),
marginal 7.3-7.8, ok >= 7.8. T from the TESS Input Catalog via the
MAST invoke API: 3-arcmin cone at the propagated star position, each
TIC row proper-motion-propagated (epoch 2000 -> t_ca) before
matching; the match is the nearest propagated row within 30 arcsec.
Applies to the channel-A narrow-rung targets with refined coverage;
channel B (antipode fields) has no blended star and needs no cut.
Queries logged to runs/tess-crossings/tic_queries.jsonl.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from tesscut_lib import tic_cone  # noqa: E402
from coverage_refine import load_events, star_track  # noqa: E402

OUT = REPO / "surveys" / "tess-crossings" / "results"
E_EXC, E_MARG = 7.3, 7.8


def main():
    ref = Table.read(OUT / "coverage_refined_v1.ecsv")
    A = ref[ref["channel"] == "A"]
    ev_tab = load_events()
    out = []
    for tid in sorted(set(str(x) for x in A["target_id"])):
        rows_t = A[np.asarray([str(x) == tid for x in A["target_id"]])]
        t_ref = float(np.mean(rows_t["t_ca_mjd"]))
        at = star_track(ev_tab, tid)
        ra, de = at(t_ref)
        rows = tic_cone(ra, de, 0.05)
        best, bsep = None, 1e9
        yr = (t_ref - 51544.5) / 365.25    # years since J2000
        for r in rows:
            try:
                rra = float(r["ra"]) + (float(r["pmRA"] or 0) / 3.6e6
                                        * yr / np.cos(np.radians(de)))
                rde = float(r["dec"]) + float(r["pmDEC"] or 0) / 3.6e6 * yr
                sep = np.hypot((rra - ra) * np.cos(np.radians(de)),
                               rde - de) * 3600.0
            except (TypeError, ValueError):
                continue
            if sep < bsep and r.get("Tmag") is not None:
                best, bsep = r, sep
        if best is None or bsep > 30.0:
            status, tmag, ticid = "no_tic_match", None, None
        else:
            tmag = float(best["Tmag"])
            ticid = int(best["ID"])
            status = ("excluded" if tmag < E_EXC
                      else "marginal" if tmag < E_MARG else "ok")
        out.append({"target_id": tid, "tic_id": ticid,
                    "tmag": tmag, "match_sep_arcsec":
                    None if best is None else round(bsep, 2),
                    "status": status,
                    "n_refined_rows": len(rows_t)})
        print(out[-1])
    (OUT / "tic_cut_v1.json").write_text(
        json.dumps({"rule": {"excluded_below": E_EXC,
                             "marginal_below": E_MARG},
                    "targets": out}, indent=1) + "\n")


if __name__ == "__main__":
    main()
