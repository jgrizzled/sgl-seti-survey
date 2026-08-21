"""Survey-wide look-elsewhere calibration for threshold exceedances.

For every endpoint x role x band cell search, compute the ratio
S_max / T for the real trajectory and, leave-one-out, for each of the 8
control trajectories (T_loo = max of the other 7). The control ratios
form the empirical null distribution of "exceedance strength" under the
same trials; the real ratios are compared against it globally. Writes
calib_v1/look_elsewhere.json.

Usage: uv run python surveys/ztf/scripts/look_elsewhere.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
CAL = REPO / "runs" / "ztf" / "calib_v1"
MIN_GOOD_FRAC, WEIGHT_CAP, CLIP, MIN_EPOCHS = 0.7, 20.0, 5.0, 5
BAND = {1: "zg", 2: "zr", 3: "zi"}


def cell_maxima(f, v, g):
    valid = np.isfinite(f) & np.isfinite(v) & (v > 0) & (g >= MIN_GOOD_FRAC)
    with np.errstate(all="ignore"):
        valid &= np.abs(f) / np.sqrt(np.where(v > 0, v, np.inf)) <= CLIP
    w = np.where(valid, 1.0 / np.where(v > 0, v, 1.0), 0.0)
    with np.errstate(all="ignore"):
        wmed = np.nanmedian(np.where(w > 0, w, np.nan), axis=1)
    w = np.minimum(w, (WEIGHT_CAP * np.nan_to_num(wmed, nan=np.inf, posinf=np.inf))[:, None])
    A, B = (np.where(valid, f, 0.0) * w).sum(1), w.sum(1)
    n = (w > 0).sum(1)
    with np.errstate(all="ignore"):
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
    S[n < MIN_EPOCHS] = np.nan
    return np.array([np.nanmax(S[t]) if np.isfinite(S[t]).any() else np.nan
                     for t in range(S.shape[0])])


def main() -> None:
    real_ratios, ctrl_ratios, per_key = [], [], {}
    for path in sorted((CAL / "tensors").glob("*.npz")):
        e, r = path.stem.split("__")
        d = np.load(path)
        for b in sorted(set(d["band_idx"].tolist())):
            eb = d["band_idx"] == b
            mx = cell_maxima(d["f"][:, eb], np.asarray(d["v"][:, eb]),
                             np.asarray(d["g"][:, eb], dtype=np.float32))
            if not np.isfinite(mx[1:]).any() or not np.isfinite(mx[0]):
                continue
            T = np.nanmax(mx[1:])
            rr = float(mx[0] / T)
            loo = []
            for t in range(1, len(mx)):
                others = np.delete(mx[1:], t - 1)
                if np.isfinite(others).any() and np.isfinite(mx[t]):
                    loo.append(float(mx[t] / np.nanmax(others)))
            real_ratios.append(rr)
            ctrl_ratios.extend(loo)
            per_key[f"{e}/{r}/{BAND[b]}"] = {"real_ratio": round(rr, 3),
                                             "control_loo_ratios": [round(x, 3) for x in loo]}
        d.close()
    real, ctrl = np.array(real_ratios), np.array(ctrl_ratios)
    summary = {
        "n_pair_bands": int(len(real)),
        "real_exceedances": int((real > 1).sum()),
        "real_exceedance_rate": round(float((real > 1).mean()), 4),
        "control_loo_exceedance_rate": round(float((ctrl > 1).mean()), 4),
        "expected_exceedances_from_controls": round(float((ctrl > 1).mean() * len(real)), 1),
        "real_ratio_percentiles": {p: round(float(np.percentile(real, p)), 3) for p in (50, 90, 99, 100)},
        "control_ratio_percentiles": {p: round(float(np.percentile(ctrl, p)), 3) for p in (50, 90, 99, 100)},
        "global_p_of_real_max": round(float((ctrl >= real.max()).mean()), 4),
        "per_candidate_global_p": {k: round(float((ctrl >= v["real_ratio"]).mean()), 4)
                                   for k, v in per_key.items() if v["real_ratio"] > 1},
    }
    (CAL / "look_elsewhere.json").write_text(json.dumps({"summary": summary, "per_key": per_key}, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
