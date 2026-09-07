"""Re-pull the eRODat DR1_Main / DR2_Main cones at the 88 anti-star corridor
positions and the 7 deep-family stars *with rows* (the recon_scan.py
catalogs stage stored only counts), parsing the VOTable and recording the
separation from the corridor centre. Output results/recon_erodat_rows_v0.json."""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import requests
from astropy.io.votable import parse_single_table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import recon_scan as R  # noqa: E402

OUT = R.OUT_DIR / "recon_erodat_rows_v0.json"


def main():
    store = R.SnapshotStore(R.RUN)
    ends = R.all_endpoints()
    targets = [("corridor", tid, d["anti"][0], d["anti"][1]) for tid, d in ends.items()]
    targets += [("star", tid, ends[tid]["star"][0], ends[tid]["star"][1]) for tid in R.DEEP]
    out = {"utc": R.now_utc(), "cones": []}
    for kind, tid, ra, de in targets:
        r_am = R.CORRIDOR_ARCMIN if kind == "corridor" else R.STAR_ARCMIN
        for cat in ("DR1_Main", "DR2_Main"):
            status, r = R.get(store, f"{R.ERODAT}/catalogue/SCS",
                              {"CAT": cat, "RA": f"{ra:.6f}", "DEC": f"{de:.6f}", "SR": f"{r_am / 60.0:.6f}", "VERB": 2},
                              f"erodat {cat} {kind} {tid}")
            rows = []; err = None
            if status == 200:
                try:
                    tab = parse_single_table(io.BytesIO(r.content)).to_table()
                    for row in tab:
                        d = {}
                        for k in tab.colnames:
                            v = row[k]
                            d[k] = v.item() if hasattr(v, "item") else (str(v) if not isinstance(v, (int, float)) else v)
                        d["_sep_arcmin"] = R.sep_deg(float(d["ra"]), float(d["dec"]), ra, de) * 60.0
                        rows.append(d)
                except Exception as exc:
                    err = f"votable: {exc}"
                    if b'name="Error"' in r.content:
                        err = r.content.decode("utf-8", "replace")[:300]
            else:
                err = f"HTTP {status}"
            out["cones"].append({"kind": kind, "target": tid, "cat": cat, "ra": ra, "dec": de,
                                 "radius_arcmin": r_am, "n_rows": len(rows), "error": err, "rows": rows})
            print(f"[erodat] {cat} {kind:8s} {tid:11s} -> {len(rows)} {err or ''}", flush=True)
    OUT.write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
