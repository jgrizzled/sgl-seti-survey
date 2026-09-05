"""SPHEREx six-detector joint cell (plan §5.15 item O3 (i); hypotheses
v2.1 amendment; notes/learnings.md §10).

The v2 SPHEREx survey searched one cell per detector D1–D6. The build
stage stored the per-trajectory stack accumulators (A = sum f w,
B = sum w, n) per detector, so the declared flat-Fnu SED (the same AB
flux density in every detector) admits one more cell per endpoint x
role without touching an image: the coherent six-detector stack

    S_J(z, mu) = sum_b A_b / sqrt(sum_b B_b)      over D1..D6,

on the common grid (96 x 3 x 3, T0 61000, ZP 23.9), with the same ring
null (48 offsets), R = S_max / T (8 designated controls), R~ = R / q95,
heavy-tail / radius-dependent void flags, and one family-wise threshold
at alpha = 0.05 over the joint family of the set. Cross-track variants
are combined as in the per-detector stage (per trajectory, the variant
with the larger S_max). This mirrors the PS1 + ZTF joint stage
(surveys/joint/joint.py) within one archive: same grid, so no
interpolation. The per-detector v2 results are untouched; the joint
family is a second, dependent test on the same data and is reported as
such.

Injections: the j-th injection of a pair is the same physical source in
all six detectors (per-injection seeded draws; flat-Fnu AB) only if the
magnitude window is shared, so the joint stage re-runs the SPHEREx
injection chain once with a common window (v1 six-detector joint m90
+/- 2, the per-detector rule applied to the joint cell) into
runs/spherex/v2/injections_joint, and combines the six per-injection
window sums exactly as the accumulators.

Rejection tests: none from catalogues (the frozen SPHEREx rule — the
template is the calibrated static-sky treatment). A surviving candidate
is `retained-ambiguous` with annotations: per-detector S and fitted flux
at the joint node with a flat-Fnu chi2, phase split, p_phase,
p_trajectory, nearest 2MASS / CatWISE source.

Usage: uv run python surveys/spherex/joint6.py <stage> [--set dev|confirmatory]
Stages: freeze, inject, nulls, completeness, adjudicate, report.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import profile as spx  # noqa: E402  (surveys/spherex/profile.py)

from sglseti import canonical_json, load_target_registry, stable_id  # noqa: E402
from sglsurvey import nulls  # noqa: E402
from sglsurvey import completeness_stage as CS  # noqa: E402
from sglsurvey.manifest import combined_hash  # noqa: E402
from sglsurvey.nulls_stage import xt_variants  # noqa: E402
from sglsurvey.records import AnalysisRun, Candidate, Constraint, append_records  # noqa: E402
from sglsurvey.vetting import parallax_phase_test  # noqa: E402

P = spx.PROFILE
SURVEY_DIR = REPO / "surveys" / "spherex"
OUT = REPO / "runs" / "spherex" / "joint6"
INJ_SUBDIR = "injections_joint"
BAND = "J6"                       # the joint cell's band label in keys / records
MASKS = ("primary", "strict", "loose")
HYPOTHESIS_VERSION = "spherex-hypotheses-v2.1-joint6"
#: dev-driven amendments (recorded in the freeze before the confirmatory run)
DEV_AMENDMENTS = [
    {"id": "a", "date": "2026-09-04", "trigger": "dev nulls: proxima-cen/rx (over-subtracted Galactic-plane field; every "
     "trajectory's joint S_max < 0, T = 0.057, ring q95 = -12.3) was not void and set R~_FWER = 5.26 for the family",
     "rule": "a cell whose ring normaliser q95 is not a positive scale (q95 <= 0) is void (null_degenerate); "
             "the heavy-tail test is only evaluated for q95 > 0",
     "note": "the per-detector v2.0 engine has no such guard; its degenerate cells happened to be caught by the "
             "heavy-tail rule (positive q95) — recorded for the QR3 re-run"}]
FREEZE_NAME = "joint6_freeze.json"
MIN_DETECTORS = 2                 # detectors with n_epochs_ok >= min_epochs for a joint cell to exist
INJ_WINDOW_HALFWIDTH = 2.0        # magnitude window = v1 joint m90 +/- this (the per-detector rule)


class JointProfile:
    """Profile-like object for the shared completeness / report code."""
    name = "spherex-joint6"
    z_grid = P.z_grid
    mu_grid = P.mu_grid
    n_z_intervals = P.n_z_intervals
    temporal_models = P.temporal_models
    zp_ref = P.zp_ref
    mag_system = P.mag_system
    holdout_min_s_late = P.holdout_min_s_late
    holdout_min_flux_ratio = P.holdout_min_flux_ratio
    min_epochs = P.min_epochs
    norm_quantile = P.norm_quantile
    heavy_tail_ratio = P.heavy_tail_ratio
    ks_alpha = P.ks_alpha
    fwer_alpha = P.fwer_alpha
    n_pseudo = P.n_pseudo
    split_seed = P.split_seed
    clip_sigma = P.clip_sigma
    recovery_window = P.recovery_window
    run_dir = OUT
    survey_dir = SURVEY_DIR
    results_dir = SURVEY_DIR / "results" / "joint6"
    null_dir = OUT / "nulls"
    v1_run_dir = REPO / "runs" / "spherex"
    spectrum = {BAND: "flat_fnu_ab"}
    extra_params = {"psf": "exposure PSF plane (10x oversampled)", "v1_calib": "calib_v4",
                    "template": "v3 static-sky template subtracted at sampling (per-detector tensors)",
                    "joint": "sum_b A_b / sqrt(sum_b B_b) over D1-D6 on the common grid"}
    hypothesis_version = HYPOTHESIS_VERSION
    registry_path = P.registry_path
    bands = (BAND,)

    @staticmethod
    def fnu_from_mag(band, mag):
        return 3631.0 * 10 ** (-0.4 * mag)

    @staticmethod
    def freeze_hash():
        return "sha256:" + hashlib.sha256((SURVEY_DIR / "configs" / FREEZE_NAME).read_bytes()).hexdigest()

    @staticmethod
    def load_freeze():
        return json.loads((SURVEY_DIR / "configs" / FREEZE_NAME).read_text())


JP = JointProfile()


def _S(A, B, n):
    with np.errstate(all="ignore"):
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
    return np.where(n >= JP.min_epochs, S, np.nan)


def joint_window(endpoint, role, band):
    """Common magnitude window of the joint injections: v1 six-detector
    joint m90 (calib_v4 m90_curves 'ALL', median over z) +/- 2."""
    m90 = spx.v1_m90(endpoint, role, "ALL")
    if m90 is None:                       # not expected (every pair has ALL curves); deepest detector as fallback
        m90 = max(m for m in (spx.v1_m90(endpoint, role, b) for b in P.bands) if m is not None)
    return (m90 - INJ_WINDOW_HALFWIDTH, m90 + INJ_WINDOW_HALFWIDTH)


def joint_profile():
    return dataclasses.replace(P, inj_subdir=INJ_SUBDIR, inj_window=joint_window)


# -- joint accumulators ---------------------------------------------------------

def joint_acc(d, mi):
    """Joint (A, B, n) over the six detectors for one mask, all
    trajectories, plus the real trajectory's per-phase parts and the
    per-detector (A, B, n) of the real trajectory (annotations)."""
    A = d["accA"][mi].astype(np.float64).sum(axis=1)     # (T, nz, nm, nm)
    B = d["accB"][mi].astype(np.float64).sum(axis=1)
    n = d["accn"][mi].sum(axis=1)
    parts = {}
    for p in (0, 1):
        Ap = d["accAp"][mi, p].astype(np.float64).sum(axis=0); Bp = d["accBp"][mi, p].astype(np.float64).sum(axis=0)
        npp = d["accnp"][mi, p].sum(axis=0)
        if npp.max() > 0:
            parts[p] = (Ap, Bp, npp)
    per_det = {"A": d["accA"][mi, 0].astype(np.float64), "B": d["accB"][mi, 0].astype(np.float64), "n": d["accn"][mi, 0]}
    return A, B, n, parts, per_det


def detector_decomposition(per_det, node, f_joint):
    """Per-detector S, fitted flux and its error at a node, and the
    flat-Fnu chi2 of the six fitted fluxes against the joint flux."""
    out = {}; chi2 = 0.0; k = 0
    for bi, b in enumerate(P.bands):
        A, B, n = per_det["A"][bi][node], per_det["B"][bi][node], int(per_det["n"][bi][node])
        if B > 0 and n >= JP.min_epochs:
            f = A / B; sig = 1.0 / np.sqrt(B)
            out[b] = {"S": float(A / np.sqrt(B)), "f": float(f), "sigma": float(sig), "n": n}
            chi2 += ((f - f_joint) / sig) ** 2; k += 1
        else:
            out[b] = {"S": None, "f": None, "sigma": None, "n": n}
    return {"per_detector": out, "flat_fnu_chi2": float(chi2), "n_detectors": k,
            "dof": max(k - 1, 0)}


def analyse_cells(endpoints, mask="primary"):
    mi = MASKS.index(mask)
    cells, meta = [], {}
    designated = 1 + np.asarray(P.designated)
    n_ring = len(P.ring)
    for e in endpoints:
        for role in ("rx", "tx"):
            tp = P.tensor_dir / f"{e}__{role}.npz"
            if not tp.exists():
                continue
            with np.load(tp) as d:
                n_ok = d["n_epochs_ok"][mi]
                n_det = int((n_ok >= JP.min_epochs).sum())
                if n_det < MIN_DETECTORS:
                    continue
                A, B, n, parts, per_det = joint_acc(d, mi)
                S = _S(A, B, n)
                mx = np.array([nulls.grid_max(S[t])[0] for t in range(S.shape[0])])
                # cross-track variants: per trajectory the larger S_max (the per-detector rule)
                variants = xt_variants(P, e, role)
                for vp in variants:
                    with np.load(vp) as v:
                        Av, Bv, nv, _, _ = joint_acc(v, mi)
                        Sv = _S(Av, Bv, nv)
                        mxv = np.array([nulls.grid_max(Sv[t])[0] for t in range(Sv.shape[0])])
                        mx = np.fmax(mx, mxv)
                R, T = nulls.exceedance_ratios(mx, designated)
                if not np.isfinite(R[0]) or not np.isfinite(T) or T <= 0:
                    continue
                ring = R[1:1 + n_ring]; ring = ring[np.isfinite(ring)]
                donors = R[1 + n_ring:]; donors = donors[np.isfinite(donors)]
                key = f"{e}/{role}/{BAND}"
                q95 = float(np.quantile(ring, JP.norm_quantile)) if len(ring) else np.nan
                # dev-driven amendment (a): a null whose normaliser is not a positive scale
                # (q95 <= 0: an over-subtracted field where every trajectory's S_max is
                # negative and T is a tiny positive number) is degenerate -> void
                degenerate = not (np.isfinite(q95) and q95 > 0)
                heavy = bool(len(ring) and not degenerate and ring.max() / q95 > JP.heavy_tail_ratio)
                ks_ring = nulls.ks_exchangeability({"inner": ring[:16], "outer": ring[32:48]}, alpha=JP.ks_alpha)
                cn = nulls.CellNull(key=key, R_real=float(R[0]), pooled=ring,
                                    by_construction={"ring": ring, "trajectory": donors},
                                    unstable=bool(ks_ring["unstable"] or heavy or degenerate), ks=ks_ring)
                cn.heavy_tail = heavy; cn.degenerate = degenerate; cn.T = float(T); cn.S_max = float(mx[0])
                s_max0, node = nulls.grid_max(S[0])
                cn.node = node; cn.n_xt = 1 + len(variants)
                cn.S_phase = {p: float(parts[p][0][node] / np.sqrt(parts[p][1][node])) if parts[p][1][node] > 0 else np.nan
                              for p in parts}
                cn.p_phase = None
                if len(parts) == 2:
                    rng = np.random.default_rng(int(hashlib.sha256(f"{key}/{JP.split_seed}".encode()).hexdigest()[:8], 16))
                    scr = nulls.phase_scrambled_maxima(parts, P.n_scramble, rng) / T
                    cn.p_phase = float((scr >= R[0]).mean())
                cn.p_trajectory = float((donors >= R[0]).mean()) if len(donors) else None
                cn.n_epochs = {b: int(n_ok[bi]) for bi, b in enumerate(P.bands)}
                cn.n_epochs_total = int(n_ok.sum()); cn.n_detectors = n_det
                f_joint = float(A[0][node] / B[0][node]) if B[0][node] > 0 else np.nan
                cn.decomp = detector_decomposition(per_det, node, f_joint)
                cn.f_joint = f_joint
                cells.append(cn)
                meta[key] = {"A": A[0], "B": B[0], "S": S[0]}
    return cells, meta


def nulls_stage(set_name, force=False):
    freeze = JP.load_freeze()
    endpoints = freeze["split"]["development" if set_name == "dev" else "confirmatory"]["endpoints"]
    out_path = OUT / "nulls" / f"{set_name}_null_ensemble.json"
    if set_name == "confirmatory" and out_path.exists() and not force:
        raise SystemExit("confirmatory joint6 analysis already exists (run once)")
    (OUT / "nulls").mkdir(parents=True, exist_ok=True)
    res = {"freeze_hash": JP.freeze_hash()}
    nm = len(P.mu_grid)
    for mask in (MASKS if set_name == "dev" else ("primary",)):
        cells, _ = analyse_cells(endpoints, mask)
        fw = nulls.fwer_threshold(cells, alpha=JP.fwer_alpha, n_draws=JP.n_pseudo, seed=JP.split_seed,
                                  norm_quantile=JP.norm_quantile)
        per = {}
        for c in cells:
            rs = nulls.rank_statement(c.R_real, c.pooled)
            q = fw["scale"].get(c.key)
            per[c.key] = {"S_max": c.S_max, "T": c.T, "R": c.R_real, "q95": q,
                          "R_norm": (c.R_real / q) if q else None,
                          "candidate": bool(q and c.R_real / q >= fw["R_fwer"] and not c.unstable),
                          "global_p": fw["global_p"].get(c.key), "rank": rs, "null_unstable": c.unstable,
                          "null_heavy_tail": c.heavy_tail, "null_degenerate": c.degenerate, "n_epochs": c.n_epochs, "n_epochs_total": c.n_epochs_total,
                          "n_detectors": c.n_detectors, "n_cross_track": c.n_xt,
                          "S_phase": [c.S_phase.get(0, float("nan")), c.S_phase.get(1, float("nan"))],
                          "p_phase": c.p_phase, "p_trajectory": c.p_trajectory,
                          "iz": c.node[0], "imu": c.node[1] * nm + c.node[2], "two_phase": len(c.S_phase) == 2,
                          "top_share": None, "ks": c.ks, "f_joint": c.f_joint, "decomposition": c.decomp}
        q = nulls.bh_qvalues({k: v["rank"]["p_cell"] for k, v in per.items()})
        for k in per:
            per[k]["bh_q"] = q.get(k)
        cands = sorted([k for k, v in per.items() if v["candidate"]], key=lambda k: -per[k]["R_norm"])
        res[mask] = {"set": set_name, "mask": mask, "n_endpoints": len(endpoints), "n_cells": len(cells),
                     "n_null_unstable": sum(c.unstable for c in cells), "n_null_heavy_tail": sum(c.heavy_tail for c in cells),
                     "n_null_degenerate": sum(c.degenerate for c in cells),
                     "n_two_phase_cells": sum(len(c.S_phase) == 2 for c in cells), "n_phase_scramble_ks_fail": None,
                     "n_detectors_hist": {str(k): sum(1 for c in cells if c.n_detectors == k) for k in range(1, 7)},
                     "fwer": {k: v for k, v in fw.items() if k not in ("scale", "global_p", "R_norm_real")},
                     "n_candidates": len(cands), "candidates": cands,
                     "real_R_quantiles": {qq: float(np.quantile([c.R_real for c in cells], qq)) for qq in (0.5, 0.9, 0.99, 1.0)} if cells else {},
                     "n_real_R_gt_1": int(sum(1 for c in cells if c.R_real > 1)),
                     "expected_R_gt_1_from_ring": float(np.mean([(c.pooled > 1).mean() for c in cells]) * len(cells)) if cells else 0.0,
                     "holdout_null_calibration": {"note": "not applicable (no epoch hold-out in QR2)"},
                     "missing_tensors": [], "cells": per}
        print(f"[joint6/{set_name}/{mask}] {len(cells)} cells, {res[mask]['n_null_unstable']} void; "
              f"R~_FWER = {fw['R_fwer']:.3f}; R > 1: {res[mask]['n_real_R_gt_1']} (ring exp. {res[mask]['expected_R_gt_1_from_ring']:.1f}); "
              f"candidates {len(cands)} {cands[:5]}", flush=True)
    if set_name == "dev" and all(m in res for m in MASKS):
        common = set(res["primary"]["cells"]) & set(res["strict"]["cells"]) & set(res["loose"]["cells"])
        res["mask_sensitivity"] = {
            "n_common_cells": len(common),
            "R_fwer": {m: res[m]["fwer"]["R_fwer"] for m in MASKS},
            "n_candidates": {m: res[m]["n_candidates"] for m in MASKS},
            "delta_R_vs_primary": {m: float(np.nanmax([abs(res[m]["cells"][k]["R"] - res["primary"]["cells"][k]["R"]) for k in common]))
                                   if common else None for m in ("strict", "loose")}}
    out_path.write_text(json.dumps(res, indent=1, default=float))


# -- injections -----------------------------------------------------------------

def inject_stage(set_name, workers=7, only_missing=True, corridors=None):
    """The SPHEREx injection chain with the common joint window into
    runs/spherex/v2/injections_joint (identical per-injection draws in
    every detector; the j-th injection is one physical source)."""
    from sglsurvey import inject_stage as IS
    Pj = joint_profile()
    spx.PROFILE = Pj      # the profile's catalog() cache dir follows PROFILE.run_dir (unchanged here)
    print(f"[joint6/inject] set {set_name} -> {Pj.inj_dir} (window = v1 joint m90 +/- {INJ_WINDOW_HALFWIDTH})", flush=True)
    IS.run(Pj, set_name, workers=workers, corridors=corridors, only_missing=only_missing)


def joint_injections(e, role):
    """Combine the six detectors' per-injection window sums of the joint
    injection run into joint injection records (primary mask)."""
    Pj = joint_profile()
    jf = Pj.inj_dir / f"{e}__{role}.json"
    if not jf.exists():
        return None
    cells = json.loads(jf.read_text())["cells"]
    wins = {}
    for b in P.bands:
        wf = Pj.inj_dir / f"{e}__{role}__{b}_windows.npz"
        if b in cells and wf.exists():
            wins[b] = np.load(wf)
    if len(wins) < MIN_DETECTORS:
        return None
    bands = list(wins)
    n = min(min(len(cells[b]["injections"]) for b in bands), min(len(wins[b]["lo"]) for b in bands))
    nz, nm = len(P.z_grid), len(P.mu_grid)
    out = []
    for j in range(n):
        recs = {b: cells[b]["injections"][j] for b in bands}
        r0 = recs[bands[0]]
        if any(abs(r["z_au"] - r0["z_au"]) > 1e-6 or abs(r["mag"] - r0["mag"]) > 1e-6 for r in recs.values()):
            continue  # draws diverged (cannot happen with per-injection seeding + the common window)
        lo = int(wins[bands[0]]["lo"][j]); wn = int(wins[bands[0]]["z_window"]) * 2 + 1
        if any(int(wins[b]["lo"][j]) != lo for b in bands):
            continue
        idx = np.arange(lo, min(lo + wn, nz)); L = len(idx)
        A = sum(wins[b]["A"][j][:L].astype(np.float64) for b in bands)
        B = sum(wins[b]["B"][j][:L].astype(np.float64) for b in bands)
        nn = sum(wins[b]["n"][j][:L].astype(int) for b in bands)
        S = _S(A, B, nn)
        iz = int(r0["iz"]); mu = r0["mu"]
        wz = np.abs(idx - iz)[:, None, None]
        wmu = np.maximum(np.abs(np.arange(nm)[:, None] - np.argmin(np.abs(P.mu_grid - mu[0]))),
                         np.abs(np.arange(nm)[None, :] - np.argmin(np.abs(P.mu_grid - mu[1]))))[None]
        near = (wz <= JP.recovery_window[0]) & (wmu <= JP.recovery_window[1])
        s_peak, pk = nulls.grid_max(np.where(near, S, np.nan))
        s_max, _ = nulls.grid_max(S)
        out.append({"j": j, "model": r0["model"], "z_au": r0["z_au"], "mu": mu, "xt_arcsec": r0.get("xt_arcsec"),
                    "mag": r0["mag"], "fnu_jy": r0["fnu_jy"], "iz": iz, "resp_median": None,
                    "n_on": {b: recs[b]["n_on"] for b in bands},
                    "masks": {"primary": {"S_peak": s_peak, "S_max": s_max, "n_detectors": len(bands)}}})
    return {"endpoint": e, "role": role, "band": BAND, "injections": out, "bands": bands,
            "mag_window": cells[bands[0]]["mag_window"]}


def completeness_stage(set_name):
    freeze = JP.load_freeze()
    endpoints = freeze["split"]["development" if set_name == "dev" else "confirmatory"]["endpoints"]
    ne = json.loads((OUT / "nulls" / f"{set_name}_null_ensemble.json").read_text())["primary"]
    r_fwer = ne["fwer"]["R_fwer"]
    rng = np.random.default_rng(JP.split_seed + 7)
    registry = load_target_registry(JP.registry_path)
    edges = CS.interval_edges(JP)
    results, constraints, hashes = {}, [], []
    for e in endpoints:
        for role in ("rx", "tx"):
            ck = f"{e}/{role}/{BAND}"
            cinfo = ne["cells"].get(ck)
            if cinfo is None:
                continue
            cell = joint_injections(e, role)
            if cell is None or not cell["injections"]:
                continue
            T = cinfo["T"]; q95 = cinfo["q95"]; unstable = cinfo["null_unstable"]
            for j in cell["injections"]:
                r = j["masks"]["primary"]; r["T"] = T
                r["R_peak"] = (r["S_peak"] / T) if (r["S_peak"] is not None and np.isfinite(r["S_peak"]) and T > 0) else None
            cls = [CS.classify(JP, j, "primary", q95, r_fwer) for j in cell["injections"]]
            # single-epoch clip limit: the brightest of the six detector limits (a source
            # clipped in one detector is not a valid joint injection either)
            with np.load(P.tensor_dir / f"{e}__{role}.npz") as dten:
                clips = [m for m in (CS.single_epoch_clip_mag(P, dten, b) for b in cell["bands"]) if m is not None]
            m_clip = max(clips) if clips else None
            thr = CS.curves_for(JP, cell["injections"], cls, "threshold", rng, mag_min=m_clip)
            fin = CS.curves_for(JP, cell["injections"], cls, "final", rng, mag_min=m_clip)
            results[ck] = {"n_injections": len(cls), "n_usable": sum(c["usable"] for c in cls), "q95": q95, "T": T,
                           "bright_limit_mag": m_clip, "mag_window": cell["mag_window"], "bands": cell["bands"],
                           "null_unstable": unstable, "n_threshold": int(sum(c["threshold"] for c in cls)),
                           "n_final": int(sum(c["final"] for c in cls)), "vetoes_fired": {},
                           "resp_median": np.nan, "threshold": thr, "final_candidate": fin, "m90_v1": spx.v1_m90(e, role, "ALL"),
                           "holdout_annotation": {"n": 0, "pass": 0, "by_model": {}}}
            hashes.append(hashlib.sha256(ck.encode()).hexdigest())
            for kind, curves in (("threshold", thr), ("final_candidate", fin)):
                for i in range(len(edges) - 1):
                    zint = (round(float(edges[i]), 1), round(float(edges[i + 1]), 1))
                    per_model = {m: res["intervals"][i] for m, res in curves.items()
                                 if res.get("status") == "ok" and res["intervals"][i]["m90"] is not None}
                    worst = min(per_model, key=lambda k: per_model[k]["m90"]) if (per_model and not unstable) else None
                    if worst is None:
                        ckind, limit, ci = "not_constrainable", None, None
                        n_inj = int(sum(res["intervals"][i]["n_injections"] for res in curves.values() if res.get("intervals")))
                    else:
                        w = per_model[worst]; ckind = "recovery_curve"
                        limit = {"value": w["m90"], "unit": "ab_mag", "band": BAND, "m50": w["m50"],
                                 "fnu_jy_at_m90": JP.fnu_from_mag(BAND, w["m90"]), "worst_temporal_model": worst,
                                 "per_model_m90": {k: v["m90"] for k, v in per_model.items()}}
                        ci = w["m90_ci68"]; n_inj = int(sum(per_model[k]["n_injections"] for k in per_model))
                    constraints.append(Constraint(
                        constraint_id=stable_id("con", {"endpoint": e, "role": role, "band": BAND, "z_interval": list(zint),
                                                        "kind": kind, "hypothesis": HYPOTHESIS_VERSION,
                                                        "freeze": freeze["freeze_content_hash"], "set": set_name}),
                        analysis_run_id="pending", endpoint_id=e, role=role, hypothesis_version=HYPOTHESIS_VERSION,
                        z_interval_au=zint, band=BAND, epoch_range_mjd=(0.0, 0.0), duty_cycle_range=(0.5, 1.0),
                        residual_motion_bound_arcsec_per_yr=1.0, morphology="point source (exposure PSF plane)",
                        kind=ckind, recovery_probability=0.9 if ckind == "recovery_curve" else None, flux_limit=limit,
                        trials_accounting_ref=f"nulls/{set_name}_null_ensemble.json",
                        extra={"set": set_name, "mask": "primary", "R_fwer_norm": r_fwer, "q95": q95,
                               "exclusion_claim": bool(kind == "final_candidate" and ckind == "recovery_curve"),
                               "detectors": cell["bands"], "joint": "six-detector flat-Fnu stack"},
                        completeness_kind=kind, ci_68=tuple(ci) if ci else None, n_injections=n_inj,
                        prf_model="exposure PSF plane (10x oversampled), 2x2 sub-pixel phases", spectrum_model="flat_fnu_ab",
                        motion_bound_norm="linf"))
    obs_hash = combined_hash(hashes)
    config = {"pipeline_id": "spherex-joint6-completeness", "set": set_name, "R_fwer_norm": r_fwer,
              "freeze_hash": JP.freeze_hash(), "spherex_freeze": P.freeze_hash()}
    run_rec = AnalysisRun(analysis_run_id=stable_id("run", {"config": json.loads(canonical_json(config)), "inputs": obs_hash}),
                          pipeline_id="spherex-joint6-completeness", pipeline_version="2.1.0", config=config,
                          observation_set_hash=obs_hash, intersection_set_hash=obs_hash,
                          registry_source_hash=registry.source_hash, hypothesis_version=HYPOTHESIS_VERSION,
                          environment={}, random_seeds={"bootstrap": JP.split_seed + 7},
                          started_utc=datetime.now(timezone.utc).isoformat(), finished_utc=datetime.now(timezone.utc).isoformat(),
                          output_files={"completeness": f"completeness/{set_name}_completeness.json"})
    constraints = [Constraint(**{**c.__dict__, "analysis_run_id": run_rec.analysis_run_id}) for c in constraints]
    append_records(OUT / "records" / "analysis_run.jsonl", [run_rec])
    append_records(OUT / "records" / "constraint.jsonl", constraints)
    (OUT / "completeness").mkdir(parents=True, exist_ok=True)

    def med(kind):
        v = [c.flux_limit["value"] for c in constraints if c.completeness_kind == kind and c.flux_limit]
        return float(np.nanmedian(v)) if v else None
    summary = {"set": set_name, "mask": "primary", "analysis_run_id": run_rec.analysis_run_id, "R_fwer_norm": r_fwer,
               "n_cells": len(results), "n_constraints": len(constraints),
               "n_recovery_curve": sum(1 for c in constraints if c.kind == "recovery_curve"),
               "median_m90": {BAND: {"threshold": med("threshold"), "final_candidate": med("final_candidate")}},
               "median_m90_v1_joint": float(np.nanmedian([r["m90_v1"] for r in results.values() if r["m90_v1"] is not None])) if results else None,
               "vetoes_fired_total": {}, "median_prf_throughput": None}
    (OUT / "completeness" / f"{set_name}_completeness.json").write_text(json.dumps({"summary": summary, "cells": results}, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


# -- adjudication ---------------------------------------------------------------

def adjudicate_stage(set_name):
    ne = json.loads((OUT / "nulls" / f"{set_name}_null_ensemble.json").read_text())["primary"]
    freeze = JP.load_freeze()
    endpoints = freeze["split"]["development" if set_name == "dev" else "confirmatory"]["endpoints"]
    cells, meta = analyse_cells(endpoints, "primary") if ne["candidates"] else ([], {})
    by_key = {c.key: c for c in cells}
    cands = []
    for key in ne["candidates"]:
        cinfo = ne["cells"][key]; c = by_key[key]
        e, role, _ = key.split("/")
        node = c.node
        z = float(P.z_grid[node[0]]); mu = (float(P.mu_grid[node[1]]), float(P.mu_grid[node[2]]))
        # nearest catalogued source annotation (never a veto under the frozen SPHEREx rule)
        nearest = {}
        with np.load(P.tensor_dir / f"{e}__{role}.npz") as d:
            centre = tuple(d["centre"]); corridor = str(d["corridor"]); cosd = np.cos(np.deg2rad(centre[1]))
            d0 = d["d0"][:, node[0], :]; dt_yr = (d["mjd"] - P.t0_mjd) / 365.25
            ra_t = centre[0] + (d0[:, 0] + mu[0] * dt_yr) / 3600.0 / cosd
            dec_t = centre[1] + (d0[:, 1] + mu[1] * dt_yr) / 3600.0
        try:
            cats = spx.catalog(corridor, centre)
        except Exception:
            cats = {}
        for b, cat in (cats or {}).items():
            if cat is None or not len(cat.ra):
                continue
            sep = np.hypot((cat.ra[None, :] - ra_t[:, None]) * cosd, cat.dec[None, :] - dec_t[:, None]) * 3600.0
            i = int(np.argmin(sep.min(axis=0)))
            nearest[b] = {"label": cat.label, "min_sep_arcsec": float(sep[:, i].min()),
                          "mag_ab": (float(cat.mag[i]) if cat.mag is not None else None)}
        cands.append(Candidate(
            candidate_id=stable_id("cnd", {"cell": key, "hypothesis": HYPOTHESIS_VERSION, "freeze": JP.freeze_hash()}),
            analysis_run_id="see-null-ensemble", endpoint_id=e, role=role, observation_ids=(),
            fitted_z_au=z, fitted_residual_motion={"mu_ra": mu[0], "mu_dec": mu[1]},
            model_comparison={"S_max": c.S_max, "T": c.T, "R": c.R_real, "q95": cinfo["q95"], "R_norm": cinfo["R_norm"],
                              "f_joint_ujy": c.f_joint, "detector_decomposition": c.decomp,
                              "tests": {"catalogue": "none (frozen rule: template is the static-sky treatment)"}},
            status="retained-ambiguous", veto_reason=None, rejection_test=None,
            global_p_value=cinfo["global_p"], rank_statement=cinfo["rank"],
            annotations={"parallax_phase": parallax_phase_test(c.S_phase), "p_phase": c.p_phase, "p_trajectory": c.p_trajectory,
                         "n_epochs": c.n_epochs, "n_detectors": c.n_detectors, "nearest_catalogue_source": nearest,
                         "decomposition": c.decomp},
            extra={"band": BAND, "set": set_name, "mask": "primary"}))
        print(f"{key}: R~={cinfo['R_norm']:.2f} -> {cands[-1].status}; flat-Fnu chi2 {c.decomp['flat_fnu_chi2']:.1f}/{c.decomp['dof']}")
    out = OUT / "records" / f"candidate_{set_name}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()
    append_records(out, cands)
    summary = {"survey": "spherex-joint6", "set": set_name, "mask": "primary", "R_fwer_norm": ne["fwer"]["R_fwer"],
               "n_candidates": len(cands), "status_counts": {s: sum(1 for c in cands if c.status == s) for s in set(c.status for c in cands)},
               "retained_ambiguous": [f"{c.endpoint_id}/{c.role}/{c.extra['band']}" for c in cands if c.status == "retained-ambiguous"],
               "written": str(out), "utc": datetime.now(timezone.utc).isoformat()}
    (OUT / "records" / f"adjudication_{set_name}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


# -- freeze ---------------------------------------------------------------------

def freeze():
    """v2.1 amendment: the per-detector v2.0 freeze unchanged (its file
    is not touched); the joint cell's rule, injection window and family
    are pinned here, with the v2.0 corridor split inherited."""
    fz = P.load_freeze()
    hyp = SURVEY_DIR / "hypotheses.md"
    prev = SURVEY_DIR / "configs" / FREEZE_NAME
    supersedes = json.loads(prev.read_text())["freeze_content_hash"] if prev.exists() else None
    fr = {"survey": "spherex-joint6", "hypothesis_version": HYPOTHESIS_VERSION, "dev_amendments": DEV_AMENDMENTS,
          "supersedes_freeze": supersedes,
          "hypotheses_hash": "sha256:" + hashlib.sha256(hyp.read_bytes()).hexdigest(),
          "spherex_v2_freeze": fz["freeze_content_hash"], "band_label": BAND,
          "rule": ("S_J(z, mu) = sum_b A_b / sqrt(sum_b B_b) over D1-D6 from the stored v2 per-trajectory "
                   "accumulators (per-frame cap, |S_e| <= 5 clip and the v3 template subtraction inherited "
                   "per detector); n_J = sum_b n_b >= 5; one cell per endpoint x role with >= "
                   f"{MIN_DETECTORS} detectors at n_epochs_ok >= 5; cross-track variants combined per trajectory "
                   "by the larger S_max; 48-offset ring null, R = S_max/T (8 designated), R~ = R/q95, heavy-tail / "
                   "inner-outer KS void flags; R~_FWER at alpha = 0.05 over the joint family of the set (a second, "
                   "dependent family on the same data; the per-detector v2.0 result is unchanged)."),
          "rejection_tests": ("none from catalogues (the frozen v2.0 rule); a surviving cell is retained-ambiguous "
                              "with annotations: per-detector S / fitted flux and flat-Fnu chi2 at the joint node, "
                              "phase split, p_phase, p_trajectory, nearest 2MASS / CatWISE source"),
          "injections": {"chain": "the v2 SPHEREx injection chain (exposure PSF plane, flat-Fnu AB, per-injection seeds)",
                         "window": f"v1 six-detector joint m90 (calib_v4 m90_curves ALL, median over z) +/- {INJ_WINDOW_HALFWIDTH}",
                         "n_per_cell": P.n_inj_per_cell, "directory": f"runs/spherex/v2/{INJ_SUBDIR}",
                         "combination": "the six per-injection window sums added exactly as the accumulators",
                         "bright_limit": "max over detectors of the single-epoch clip magnitude",
                         "temporal_models": list(P.temporal_models), "recovery_window": list(P.recovery_window)},
          "masks": {"primary": "all usable exposures (= loose)", "strict": "deep-field collection excluded",
                    "note": "dev set analysed under all three masks; confirmatory under primary only"},
          "seeds": {"fwer_pseudo": JP.split_seed, "bootstrap": JP.split_seed + 7, "phase_scramble": "sha256(cell key / split seed)"},
          "ordering": "dev set first (all masks), then the confirmatory set once (primary)",
          "split": fz["split"], "frozen_at": datetime.now(timezone.utc).date().isoformat()}
    fr["freeze_content_hash"] = "sha256:" + hashlib.sha256(canonical_json(fr).encode()).hexdigest()
    (SURVEY_DIR / "configs" / FREEZE_NAME).write_text(json.dumps(fr, indent=1, sort_keys=True) + "\n")
    print("spherex joint6 freeze", fr["freeze_content_hash"][:23])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["freeze", "inject", "nulls", "completeness", "adjudicate", "report"])
    ap.add_argument("--set", choices=["dev", "confirmatory", "all"], default="dev")
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--corridors", nargs="*")
    ap.add_argument("--no-only-missing", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if a.stage == "freeze":
        freeze()
    elif a.stage == "inject":
        inject_stage(a.set, workers=a.workers, only_missing=not a.no_only_missing, corridors=a.corridors)
    elif a.stage == "nulls":
        nulls_stage(a.set, a.force)
    elif a.stage == "completeness":
        completeness_stage(a.set)
    elif a.stage == "adjudicate":
        adjudicate_stage(a.set)
    elif a.stage == "report":
        from sglsurvey import report_stage
        raise SystemExit(report_stage.run(JP))


if __name__ == "__main__":
    main()
