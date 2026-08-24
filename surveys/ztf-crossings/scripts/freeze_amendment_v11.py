"""Threshold freeze amendment v1.1 (channel A), per dev findings A1/A2.

Writes configs/threshold_freeze_v1_1.json binding the amendment to the
v1.0 freeze hash and the dev report hash. Channel B is unchanged.
"""
import hashlib
import json
from pathlib import Path

D = Path(__file__).resolve().parents[1]
v10 = json.loads((D / "configs" / "threshold_freeze_v1.json").read_text())

amendment = {
    "freeze_version": "ztf-crossings-thresholds-v1.1",
    "amends": v10["freeze_content_hash"],
    "frozen_at": "2026-08-23",
    "dev_report_hash": "sha256:" + hashlib.sha256(
        (D / "results" / "dev_search_v1.md").read_bytes()).hexdigest(),
    "scope": "channel A only; channel B constructions unchanged (validated "
             "clean on dev)",
    "a_systematics_template": {
        "model": "f_model(t) = c0 + c1*dt_yr + cx*P_ra(t) + cy*P_dec(t) "
                 "+ cs*(seeing - median_seeing)",
        "parallax_basis": "P_ra, P_dec are the standard parallax factors "
                          "computed from the astropy Earth barycentric "
                          "position and the registry (ra, dec, parallax) - "
                          "phase and shape fully fixed by astrometry; only "
                          "the two linear couplings cx, cy are fitted "
                          "(no free-phase sinusoid)",
        "fit_epochs": "off-window only (outside every real window of either "
                      "rung), deterministic uniform-in-time subsample of "
                      "<= 150 epochs per (target, band)",
        "fit_method": "least squares with 3 iterations of 3-sigma clipping",
        "statistic_input": "per-epoch forced flux minus f_model; window "
                           "statistic construction unchanged",
        "quality_gate": "unit is searchable only if >= 30 off-window fit "
                        "epochs survive clipping AND |median of control "
                        "window statistics| <= 1; otherwise constraint-only",
    },
    "b_pseudo_window_exclusion": {
        "rule": "pseudo-windows must avoid real windows of the SAME rung "
                "only (was: either rung)",
        "rationale": "wide-beam signal leakage into narrow-rung controls "
                     "inflates T - conservative direction only; rung "
                     "verdicts on one target declared correlated in trials "
                     "accounting",
        "undefined_disposition": "units that still lack 8 valid offsets are "
                                 "constraint-only (S and completeness "
                                 "reported; no exceedance test, no candidate "
                                 "claims)",
    },
    "c_exceedance_rule": {
        "rule": "exceedance iff S > max(T, 0); report margin S - T",
        "replaces": "R = S/T > 1 (ratio unstable for T -> 0+ and "
                    "sign-invalid for T <= 0)",
    },
    "unchanged": ["controls count (8) and designated offsets/redraw sequence",
                  "ring construction (channel B)", "WEIGHT_CAP", "masks",
                  "split (dev/confirmatory)", "search-unit definitions",
                  "per-search crossing probability 1/9 (for defined units)"],
}
body = json.dumps(amendment, indent=1, sort_keys=True)
amendment["freeze_content_hash"] = ("sha256:"
                                    + hashlib.sha256(body.encode()).hexdigest())
out = D / "configs" / "threshold_freeze_v1_1.json"
out.write_text(json.dumps(amendment, indent=1, sort_keys=True) + "\n")
print(out, amendment["freeze_content_hash"])
