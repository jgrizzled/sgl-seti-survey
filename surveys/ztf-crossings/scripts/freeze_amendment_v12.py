"""Threshold freeze amendment v1.2 (channel A): empirical variance
rescaling, the v1 'confusion, not noise, sets the floor' precedent."""
import hashlib
import json
from pathlib import Path

D = Path(__file__).resolve().parents[1]
v11 = json.loads((D / "configs" / "threshold_freeze_v1_1.json").read_text())

amendment = {
    "freeze_version": "ztf-crossings-thresholds-v1.2",
    "amends": v11["freeze_content_hash"],
    "frozen_at": "2026-08-23",
    "scope": "channel A only; adds one element on top of v1.1, which is "
             "otherwise unchanged (as is channel B under v1.0)",
    "d_variance_rescale": {
        "rule": "per (target, band): after the v1.1 template fit, per-epoch "
                "variances are scaled by k = median(r_i^2 / v_i) / 0.4549 "
                "over the sigma-clipped off-window fit epochs (robust "
                "chi-square-per-dof; 0.4549 = median of chi^2_1), floored at "
                "k = 1 (never below the matched-filter background variance)",
        "effect": "window statistics and controls are standardized by the "
                  "empirical residual scatter, so S is O(1) under the null "
                  "and the v1.1 quality gate |median control| <= 1 is "
                  "meaningful; the exceedance comparison S > max(T, 0) is "
                  "unchanged (scale cancels between S and T)",
        "rationale": "dev v1.1 showed matched-filter variances understate "
                     "bright-star residual scatter by ~10^2-10^4 - the WISE "
                     "v1 lesson (confusion, not noise, sets depth) applied "
                     "to blended-star difference photometry",
    },
}
body = json.dumps(amendment, indent=1, sort_keys=True)
amendment["freeze_content_hash"] = ("sha256:"
                                    + hashlib.sha256(body.encode()).hexdigest())
out = D / "configs" / "threshold_freeze_v1_2.json"
out.write_text(json.dumps(amendment, indent=1, sort_keys=True) + "\n")
print(out, amendment["freeze_content_hash"])
