"""Step D (part 2) / step G: per-cell null ensembles, exchangeability,
survey-wide FWER threshold, rank statements, mask sensitivity.

Reads the per-trajectory summaries written by build_tensors.py for the
cells of one hold-out set and, per cell (endpoint x role x band) under
the primary quality mask:

  R for the real trajectory and for every null trajectory (ring 48 with
  leave-one-out for the designated controls, 50 trajectory-randomised,
  200 phase-coherence scrambles where two phases are populated);
  KS exchangeability between constructions; the pooled null; the
  normaliser q95; the rank statement of the real R.

Then the family-wide pseudo-experiment null of max R~ = R / q95 over
the set's cells, R~_FWER at alpha = 0.05, the global p-value of every
cell, BH q-values (informational), the list of candidate cells
(R~ >= R~_FWER), and the held-out-epoch test's false-pass rate on null
trajectories. The strict / loose masks are re-analysed the same way
for the systematic spread (development set).

The confirmatory set is analysed ONCE, blind: this script computes the
threshold from the nulls before it looks at any real R, and writes
both in one pass; it refuses to run on the confirmatory set if a
previous confirmatory output exists (--force to override, which the
report must then disclose).

Usage: uv run python surveys/wise/scripts/null_ensemble.py --set dev|confirmatory
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

from sglsurvey import nulls  # noqa: E402
from sglsurvey.vetting import holdout_prediction_test  # noqa: E402

MASKS = ("primary", "strict", "loose")
F = {"S_max": 0, "iz": 1, "imu": 2, "S_phase0": 3, "S_phase1": 4, "n_epochs": 5, "top_share": 6}
H = {"S_early": 0, "f_early": 1, "S_late": 2, "f_late": 3, "n_late": 4, "iz_early": 5, "imu_early": 6}


def xt_variants(endpoint: str, role: str) -> list[Path]:
    """Cross-track tensor variants of a pair (hypotheses v2.0 §1.4)."""
    return sorted(C.TENSOR_DIR.glob(f"{endpoint}__{role}__xt*.npz"))


def combine_xt(d, variants: list[Path]):
    """Dict of the summary arrays with the statistic maximised over the
    cross-track offsets: per trajectory / mask / band the variant with
    the larger S_max supplies the whole summary row; scramble maxima
    are combined elementwise (identically seeded draws)."""
    out = {k: np.array(d[k]) for k in ("summary", "holdout", "scramble_max",
                                        "epoch_scramble_max", "n_epochs_ok")}
    for vp in variants:
        with np.load(vp) as v:
            better = np.nan_to_num(v["summary"][..., 0], nan=-np.inf) > np.nan_to_num(out["summary"][..., 0], nan=-np.inf)
            out["summary"] = np.where(better[..., None], v["summary"], out["summary"])
            out["holdout"] = np.where(better[..., None], v["holdout"], out["holdout"])
            out["scramble_max"] = np.fmax(out["scramble_max"], v["scramble_max"])
            out["epoch_scramble_max"] = np.fmax(out["epoch_scramble_max"], v["epoch_scramble_max"])
    out["n_xt"] = 1 + len(variants)
    return out


def cell_nulls(d, mi: int, bi: int, key: str, variants=()) -> nulls.CellNull | None:
    """Build the CellNull of one (mask, band) from a tensor file (and
    its cross-track variants, if any)."""
    dd = combine_xt(d, list(variants)) if variants else d
    summ = dd["summary"][mi, :, bi]          # (99, 7)
    smax = summ[:, F["S_max"]].astype(float)
    if not np.isfinite(smax[0]) or d["n_epochs_ok"][mi, bi] < C.MIN_EPOCHS:
        return None
    designated = 1 + np.asarray(d["designated"])
    R, T = nulls.exceedance_ratios(smax, designated)
    if not np.isfinite(T) or T <= 0:
        return None
    n_ring = len(d["ring_offsets"])
    ring = R[1:1 + n_ring]
    donors = R[1 + n_ring:]
    two_phase = np.isfinite(summ[0, F["S_phase0"]]) and np.isfinite(summ[0, F["S_phase1"]])
    scr = dd["scramble_max"][mi, bi].astype(float) / T if two_phase else np.array([])
    by = {"ring": ring[np.isfinite(ring)], "trajectory": donors[np.isfinite(donors)]}
    if two_phase:
        by["phase_scramble"] = scr[np.isfinite(scr)]
    ks = nulls.ks_exchangeability(by, alpha=C.KS_ALPHA)
    # Stability (hypotheses v2.0 §3.5): the local null must not depend on
    # the offset radius — inner (20") vs outer (40") ring members.
    n_ang = len(C.RING_ANGLES)
    inner, outer = ring[:n_ang], ring[2 * n_ang:3 * n_ang]
    ks_ring = nulls.ks_exchangeability({"inner": inner, "outer": outer}, alpha=C.KS_ALPHA)
    ks["pairs"]["ring_inner|ring_outer"] = ks_ring["pairs"].get("inner|outer", {"D": None, "p": None})
    q95 = float(np.quantile(by["ring"], C.NORM_QUANTILE)) if len(by["ring"]) else np.nan
    heavy = bool(len(by["ring"]) and by["ring"].max() / q95 > C.HEAVY_TAIL_RATIO)
    unstable = bool(ks_ring["unstable"]) or heavy
    # The pooled (exchangeable) null is the spatial ring; the other two
    # constructions are reported as per-cell annotations (p_trajectory,
    # p_phase), see hypotheses v2.0 §3.4.
    pooled = by["ring"]
    cn = nulls.CellNull(key=key, R_real=float(R[0]), pooled=pooled, by_construction=by,
                        unstable=unstable, ks=ks)
    cn.heavy_tail = heavy
    cn.p_trajectory = (float((by["trajectory"] >= R[0]).mean()) if len(by["trajectory"]) else None)
    cn.p_phase = (float((by["phase_scramble"] >= R[0]).mean())
                  if two_phase and len(by["phase_scramble"]) else None)
    cn.T = float(T)
    cn.S_max = float(smax[0])
    cn.R_all = R
    cn.epoch_scramble = (dd["epoch_scramble_max"][mi, bi].astype(float) / T)
    cn.n_xt = 1 + len(variants)
    cn.two_phase = two_phase
    cn.summary = summ[0]
    cn.hold = dd["holdout"][mi, :, bi]
    cn.n_epochs = int(d["n_epochs_ok"][mi, bi])
    return cn


def holdout_rates(cells: list) -> dict:
    """False-pass rate of the held-out-epoch test on null trajectories
    (all, and those above their cell's q95)."""
    n_all = n_pass = n_hi = n_hi_pass = 0
    for c in cells:
        q = c.scale(C.NORM_QUANTILE)
        for t in range(1, len(c.R_all)):
            h = c.hold[t]
            if not np.isfinite(h[H["S_late"]]):
                continue
            res = holdout_prediction_test(h[H["S_early"]], h[H["f_early"]], h[H["S_late"]],
                                          h[H["f_late"]], int(h[H["n_late"]]),
                                          C.HOLDOUT_MIN_S_LATE, C.HOLDOUT_MIN_FLUX_RATIO)
            if not res.get("applicable"):
                continue
            n_all += 1
            n_pass += res["pass"]
            if np.isfinite(c.R_all[t]) and c.R_all[t] >= q:
                n_hi += 1
                n_hi_pass += res["pass"]
    return {"n_null_trajectories": n_all, "false_pass_rate": n_pass / n_all if n_all else None,
            "false_pass_ci68": nulls.wilson_interval(n_pass, n_all) if n_all else None,
            "n_above_q95": n_hi, "false_pass_rate_above_q95": n_hi_pass / n_hi if n_hi else None}


def analyse(set_name: str, mask: str, endpoints: list[str], verbose: bool = True) -> dict:
    mi = MASKS.index(mask)
    cells, missing = [], []
    for e in endpoints:
        for role in ("rx", "tx"):
            p = C.TENSOR_DIR / f"{e}__{role}.npz"
            if not p.exists():
                missing.append(f"{e}__{role}")
                continue
            d = np.load(p)
            variants = xt_variants(e, role)
            for bi in range(4):
                key = f"{e}/{role}/{C.BAND_NAME[bi + 1]}"
                cn = cell_nulls(d, mi, bi, key, variants)
                if cn is not None:
                    cells.append(cn)
            d.close()
    # family-wide threshold from the nulls alone (before any real R is used)
    fw = nulls.fwer_threshold(cells, alpha=C.FWER_ALPHA, n_draws=C.N_PSEUDO,
                              seed=C.SPLIT_SEED, norm_quantile=C.NORM_QUANTILE)
    per_cell = {}
    pvals = {}
    for c in cells:
        rs = nulls.rank_statement(c.R_real, c.pooled)
        pvals[c.key] = rs["p_cell"]
        q = fw["scale"].get(c.key)
        per_cell[c.key] = {
            "S_max": c.S_max, "T": c.T, "R": c.R_real, "q95": q,
            "R_norm": (c.R_real / q) if q else None,
            "candidate": bool(q and c.R_real / q >= fw["R_fwer"] and not c.unstable),
            "global_p": fw["global_p"].get(c.key), "rank": rs,
            "n_epochs": c.n_epochs, "two_phase": bool(c.two_phase), "n_cross_track": c.n_xt,
            "null_unstable": c.unstable, "null_heavy_tail": c.heavy_tail, "ks": c.ks,
            "n_null": {k: int(len(v)) for k, v in c.by_construction.items()},
            "null_R_quantiles": {k: [float(np.quantile(v, qq)) for qq in (0.5, 0.95, 0.99)]
                                 for k, v in c.by_construction.items() if len(v)},
            "epoch_scramble_R_median": float(np.nanmedian(c.epoch_scramble)),
            "p_trajectory": c.p_trajectory, "p_phase": c.p_phase,
            "iz": int(c.summary[F["iz"]]), "imu": int(c.summary[F["imu"]]),
            "S_phase": [float(c.summary[F["S_phase0"]]), float(c.summary[F["S_phase1"]])],
            "top_share": float(c.summary[F["top_share"]]),
        }
    q = nulls.bh_qvalues(pvals)
    for k in per_cell:
        per_cell[k]["bh_q"] = q.get(k)
    cands = sorted([k for k, v in per_cell.items() if v["candidate"]],
                   key=lambda k: -per_cell[k]["R_norm"])
    n_unstable = sum(1 for c in cells if c.unstable)
    n_heavy = sum(1 for c in cells if c.heavy_tail)
    ks_phase_fail = sum(1 for c in cells if c.two_phase and any(
        v["p"] is not None and v["p"] < C.KS_ALPHA for kk, v in c.ks["pairs"].items() if "phase" in kk))
    ks_traj_fail = sum(1 for c in cells if (c.ks["pairs"].get("ring|trajectory") or {}).get("p", 1) < C.KS_ALPHA)
    out = {
        "set": set_name, "mask": mask, "n_endpoints": len(endpoints),
        "n_cells": len(cells), "missing_tensors": missing,
        "n_null_unstable": n_unstable, "n_null_heavy_tail": n_heavy,
        "n_two_phase_cells": int(sum(1 for c in cells if c.two_phase)),
        "n_phase_scramble_ks_fail": int(ks_phase_fail),
        "n_ring_vs_trajectory_ks_fail": int(ks_traj_fail),
        "fwer": {k: v for k, v in fw.items() if k not in ("scale", "global_p", "R_norm_real")},
        "n_candidates": len(cands), "candidates": cands,
        "real_R_quantiles": {qq: float(np.quantile([c.R_real for c in cells], qq))
                             for qq in (0.5, 0.9, 0.99, 1.0)},
        "n_real_R_gt_1": int(sum(1 for c in cells if c.R_real > 1)),
        "expected_R_gt_1_from_ring": float(np.mean([
            (c.by_construction["ring"] > 1).mean() for c in cells]) * len(cells)),
        "holdout_null_calibration": holdout_rates(cells),
        "cells": per_cell,
    }
    if verbose:
        print(f"[{set_name}/{mask}] {len(cells)} cells, {n_unstable} null-unstable ({n_heavy} heavy-tail), "
              f"{ks_traj_fail} ring-vs-trajectory KS fails, {out['n_two_phase_cells']} two-phase "
              f"({ks_phase_fail} phase-scramble KS fails); "
              f"R~_FWER = {fw['R_fwer']:.3f}; real R > 1: {out['n_real_R_gt_1']} "
              f"(ring expectation {out['expected_R_gt_1_from_ring']:.1f}); "
              f"candidates: {len(cands)} {cands[:5]}", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["dev", "confirmatory"], required=True)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    freeze = C.load_freeze()
    key = "development" if a.set == "dev" else "confirmatory"
    endpoints = freeze["split"][key]["endpoints"]
    C.NULL_DIR.mkdir(parents=True, exist_ok=True)
    out_path = C.NULL_DIR / f"{a.set}_null_ensemble.json"
    if a.set == "confirmatory" and out_path.exists() and not a.force:
        raise SystemExit("confirmatory analysis already exists; it is run once (use --force "
                         "only with disclosure in the report)")
    results = {"freeze_hash": C.freeze_hash(), "freeze_content_hash": freeze["freeze_content_hash"]}
    masks = MASKS if a.set == "dev" else ("primary",)
    for mask in masks:
        results[mask] = analyse(a.set, mask, endpoints)
    if a.set == "dev":
        # systematic spread between masks
        prim, strict, loose = (results[m]["cells"] for m in MASKS)
        common = set(prim) & set(strict) & set(loose)
        dR = {m: [results[m]["cells"][k]["R"] - prim[k]["R"] for k in common] for m in ("strict", "loose")}
        results["mask_sensitivity"] = {
            "n_common_cells": len(common),
            "R_fwer": {m: results[m]["fwer"]["R_fwer"] for m in MASKS},
            "n_candidates": {m: results[m]["n_candidates"] for m in MASKS},
            "delta_R_vs_primary": {m: {"median": float(np.median(v)), "p16": float(np.percentile(v, 16)),
                                       "p84": float(np.percentile(v, 84))} for m, v in dR.items()},
            "n_cells": {m: results[m]["n_cells"] for m in MASKS},
        }
    out_path.write_text(json.dumps(results, indent=1, default=float))
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
