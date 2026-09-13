"""Flux-scale recon summary (hypotheses §4): distribution of the
per-frame PS1-star zero points, gate attrition, and the star ZP vs
header MAGZPT offset (the header value is a cross-check only).
Writes surveys/ptf/results/flux_scale_check.json.

Usage: uv run python surveys/ptf/scripts/flux_scale_check.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
ZP = REPO / "runs" / "ptf" / "zeropoints.jsonl"
OUT = REPO / "surveys" / "ptf" / "results" / "flux_scale_check.json"


def q(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return None
    return {"n": int(len(x)), "median": float(np.median(x)), "p10": float(np.percentile(x, 10)),
            "p90": float(np.percentile(x, 90)), "min": float(x.min()), "max": float(x.max())}


def main():
    recs = [json.loads(l) for l in open(ZP)]
    out = {"n_frames": len(recs)}
    for band in ("g", "R"):
        rs = [r for r in recs if r["band"] == band]
        ok = [r for r in rs if r.get("zp_star") is not None]
        hdr = [(r["zp_star"], r["magzpt_header"]) for r in ok
               if isinstance(r.get("magzpt_header"), (int, float))]
        fail = [r for r in rs if r.get("zp_star") is None]
        reasons = {"error": sum(1 for r in fail if r.get("error")),
                   "too_few_stars": sum(1 for r in fail if not r.get("error") and (r.get("zp_mad") is None)),
                   "scatter_gate": sum(1 for r in fail if not r.get("error") and r.get("zp_mad") is not None)}
        out[band] = {
            "n_frames": len(rs), "n_calibrated": len(ok), "attrition": reasons,
            "zp_star": q([r["zp_star"] for r in ok]),
            "zp_scatter_mag": q([r["zp_mad"] for r in ok]),
            "n_stars": q([r["n_stars"] for r in ok]),
            "n_candidates_in_frame": q([r.get("n_in_frame", np.nan) for r in rs]),
            "fwhm_arcsec": q([r.get("fwhm_arcsec", np.nan) for r in rs]),
            "star_minus_header_zp": q([a - b for a, b in hdr]),
            "n_header_magzpt": len(hdr),
        }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))
    for band in ("g", "R"):
        b = out[band]
        print(f"[{band}] {b['n_calibrated']}/{b['n_frames']} calibrated; scatter median "
              f"{(b['zp_scatter_mag'] or {}).get('median')}; stars median {(b['n_stars'] or {}).get('median')}; "
              f"star-header {(b['star_minus_header_zp'] or {}).get('median')} (n={b['n_header_magzpt']}); "
              f"attrition {b['attrition']}")


if __name__ == "__main__":
    main()
