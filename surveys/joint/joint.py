"""Joint PS1 + ZTF v2 stage (project plan §10.1 step 6; WISE v2_plan §12:
"the joint stage multiplies the family; needs one FWER across
archives; inject the same physical source into all archives").

Inputs are the two archives' v2 products on the common µ reference
epoch T0 = 59800 and AB zero point 25: the per-trajectory accumulated
stack sums (A, B, n per mask x trajectory x band x grid) that
`sglsurvey.build` stores in every tensor, the real trajectory's
per-phase parts, and the per-injection window sums written by
`sglsurvey.inject_stage` (the j-th injection of a cell is the same
physical source in both archives: per-injection seeded draws on the
union magnitude window).

Joint cell = endpoint x role x band pair (g: PS1 g + ZTF zg; r; i).
The ZTF sums (192 nodes) are interpolated in 1/z onto the PS1 grid
(360 nodes); trajectory index t is the same sky offset in both
archives (identical ring); S_joint = (A_z + A_p) / sqrt(B_z + B_p).
Per joint cell: R over the 8 designated controls, the 48-ring null,
q95, R~, the void flags, the joint phase scramble (annotation); the
family-wise R~_FWER over the set's joint cells; joint injections
classified with that rule; the flux-consistent static test from both
catalogues (rejection only if the catalogued sources account for the
joint S); Constraint / Candidate records; report tables.

Usage: uv run python surveys/joint/joint.py --set dev|confirmatory [--force]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "surveys" / "ztf"))
sys.path.insert(0, str(REPO / "surveys" / "panstarrs"))
import importlib.util  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ZTF = _load("ztf_profile", REPO / "surveys" / "ztf" / "profile.py").PROFILE
PS1 = _load("ps1_profile", REPO / "surveys" / "panstarrs" / "profile.py").PROFILE

from sglseti import canonical_json, load_target_registry, stable_id  # noqa: E402
from sglsurvey import nulls  # noqa: E402
from sglsurvey.manifest import combined_hash  # noqa: E402
from sglsurvey.records import AnalysisRun, Candidate, Constraint, append_records  # noqa: E402
from sglsurvey import completeness_stage as CS  # noqa: E402
from sglsurvey.adjudicate_stage import response_table, _FmLike  # noqa: E402
from sglsurvey.vetting import flux_consistency, parallax_phase_test  # noqa: E402

OUT = REPO / "runs" / "joint" / "v2"
SURVEY_DIR = REPO / "surveys" / "joint"
BAND_PAIRS = {"g": ("g", "zg"), "r": ("r", "zr"), "i": ("i", "zi")}
MASKS = ("primary", "strict", "loose")
HYPOTHESIS_VERSION = "joint-ps1-ztf-v2.0"
RECOVERY_WINDOW = (2, 1)


class JointProfile:
    """Minimal profile-like object for the shared completeness code."""
    name = "joint"
    z_grid = PS1.z_grid
    mu_grid = PS1.mu_grid
    n_z_intervals = 8
    temporal_models = PS1.temporal_models
    zp_ref = 25.0
    mag_system = "ab"
    holdout_min_s_late = 3.0
    holdout_min_flux_ratio = 0.3
    min_epochs = 5
    norm_quantile = 0.95
    heavy_tail_ratio = 2.5
    ks_alpha = 0.01
    fwer_alpha = 0.05
    n_pseudo = 10000
    split_seed = 20260822
    static_search_arcsec = 30.0
    static_factor = 2.0
    run_dir = OUT
    survey_dir = SURVEY_DIR
    spectrum = {b: "flat_fnu_ab" for b in BAND_PAIRS}
    extra_params = {"psf": "per-archive (Moffat at seeing)"}
    hypothesis_version = HYPOTHESIS_VERSION
    registry_path = PS1.registry_path
    bands = tuple(BAND_PAIRS)

    @staticmethod
    def fnu_from_mag(band, mag):
        return 3631.0 * 10 ** (-0.4 * mag)

    @staticmethod
    def freeze_hash():
        return "sha256:" + hashlib.sha256((SURVEY_DIR / "configs" / "v2_freeze.json").read_bytes()).hexdigest()

    @staticmethod
    def load_freeze():
        return json.loads((SURVEY_DIR / "configs" / "v2_freeze.json").read_text())


JP = JointProfile()


def interp_q(arr_z, q_from, q_to):
    """Interpolate an array along axis 1 (z nodes) from q_from to q_to."""
    from scipy.interpolate import interp1d
    o = np.argsort(q_from)
    f = interp1d(q_from[o], np.take(arr_z, o, axis=1), axis=1, bounds_error=False, fill_value=0.0)
    return f(q_to)


def joint_cell(dz, dp, bz, bp, mi):
    """Joint stack sums for one (mask) over all trajectories: returns
    (A, B, n) on the PS1 grid, and the real trajectory's phase parts."""
    qz, qp = 1.0 / ZTF.z_grid, 1.0 / PS1.z_grid
    Az = interp_q(dz["accA"][mi, :, bz].astype(float), qz, qp)
    Bz = interp_q(dz["accB"][mi, :, bz].astype(float), qz, qp)
    nz_ = np.rint(interp_q(dz["accn"][mi, :, bz].astype(float), qz, qp)).astype(int)
    Ap_, Bp_, np__ = dp["accA"][mi, :, bp].astype(float), dp["accB"][mi, :, bp].astype(float), dp["accn"][mi, :, bp]
    T = min(Az.shape[0], Ap_.shape[0])
    A = Az[:T] + Ap_[:T]; B = Bz[:T] + Bp_[:T]; n = nz_[:T] + np__[:T]
    Apz = interp_q(dz["accAp"][mi, :, bz].astype(float), qz, qp); Bpz = interp_q(dz["accBp"][mi, :, bz].astype(float), qz, qp)
    npz_ = np.rint(interp_q(dz["accnp"][mi, :, bz].astype(float), qz, qp)).astype(int)
    parts = {}
    for p in (0, 1):
        Ap = Apz[p] + dp["accAp"][mi, p, bp]; Bp = Bpz[p] + dp["accBp"][mi, p, bp]
        npp = npz_[p] + dp["accnp"][mi, p, bp]
        if npp.max() > 0:
            parts[p] = (Ap, Bp, npp)
    return A, B, n, parts


def _S(A, B, n):
    with np.errstate(all="ignore"):
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
    return np.where(n >= JP.min_epochs, S, np.nan)


def analyse_cells(endpoints, mask="primary"):
    mi = MASKS.index(mask)
    cells, meta = [], {}
    designated = 1 + np.asarray(PS1.designated)
    for e in endpoints:
        for role in ("rx", "tx"):
            pz, pp = ZTF.tensor_dir / f"{e}__{role}.npz", PS1.tensor_dir / f"{e}__{role}.npz"
            if not (pz.exists() and pp.exists()):
                continue
            dz, dp = np.load(pz), np.load(pp)
            for jb, (b_ps1, b_ztf) in BAND_PAIRS.items():
                bz, bp = ZTF.band_idx[b_ztf] - 1, PS1.band_idx[b_ps1] - 1
                if dz["n_epochs_ok"][mi, bz] < JP.min_epochs or dp["n_epochs_ok"][mi, bp] < JP.min_epochs:
                    continue
                A, B, n, parts = joint_cell(dz, dp, bz, bp, mi)
                S = _S(A, B, n)
                mx = np.array([nulls.grid_max(S[t])[0] for t in range(S.shape[0])])
                R, T = nulls.exceedance_ratios(mx, designated)
                if not np.isfinite(R[0]) or not np.isfinite(T) or T <= 0:
                    continue
                ring = R[1:49]; ring = ring[np.isfinite(ring)]
                donors = R[49:]; donors = donors[np.isfinite(donors)]
                key = f"{e}/{role}/{jb}"
                q95 = float(np.quantile(ring, JP.norm_quantile)) if len(ring) else np.nan
                heavy = bool(len(ring) and ring.max() / q95 > JP.heavy_tail_ratio)
                ks_ring = nulls.ks_exchangeability({"inner": ring[:16], "outer": ring[32:48]}, alpha=JP.ks_alpha)
                cn = nulls.CellNull(key=key, R_real=float(R[0]), pooled=ring,
                                    by_construction={"ring": ring, "trajectory": donors},
                                    unstable=bool(ks_ring["unstable"] or heavy), ks=ks_ring)
                cn.heavy_tail = heavy; cn.T = float(T); cn.S_max = float(mx[0])
                s_max, node = nulls.grid_max(S[0])
                cn.node = node
                cn.S_phase = {p: float(parts[p][0][node] / np.sqrt(parts[p][1][node])) if parts[p][1][node] > 0 else np.nan
                              for p in parts}
                cn.p_phase = None
                if len(parts) == 2:
                    rng = np.random.default_rng(int(hashlib.sha256(f"{key}/{JP.split_seed}".encode()).hexdigest()[:8], 16))
                    scr = nulls.phase_scrambled_maxima(parts, 200, rng) / T
                    cn.p_phase = float((scr >= R[0]).mean())
                cn.p_trajectory = float((donors >= R[0]).mean()) if len(donors) else None
                cn.n_epochs = {"ztf": int(dz["n_epochs_ok"][mi, bz]), "ps1": int(dp["n_epochs_ok"][mi, bp])}
                cells.append(cn)
                meta[key] = {"A": A[0], "B": B[0], "S": S[0]}
            dz.close(); dp.close()
    return cells, meta


def nulls_stage(set_name, force=False):
    freeze = JP.load_freeze()
    endpoints = freeze["split"]["development" if set_name == "dev" else "confirmatory"]["endpoints"]
    out_path = OUT / "nulls" / f"{set_name}_null_ensemble.json"
    if set_name == "confirmatory" and out_path.exists() and not force:
        raise SystemExit("confirmatory joint analysis already exists (run once)")
    (OUT / "nulls").mkdir(parents=True, exist_ok=True)
    res = {"freeze_hash": JP.freeze_hash()}
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
                          "null_heavy_tail": c.heavy_tail, "n_epochs": c.n_epochs,
                          "S_phase": [c.S_phase.get(0, float("nan")), c.S_phase.get(1, float("nan"))],
                          "p_phase": c.p_phase, "p_trajectory": c.p_trajectory,
                          "iz": c.node[0], "imu": c.node[1] * 5 + c.node[2], "two_phase": len(c.S_phase) == 2,
                          "top_share": None, "ks": c.ks}
        q = nulls.bh_qvalues({k: v["rank"]["p_cell"] for k, v in per.items()})
        for k in per:
            per[k]["bh_q"] = q.get(k)
        cands = sorted([k for k, v in per.items() if v["candidate"]], key=lambda k: -per[k]["R_norm"])
        res[mask] = {"set": set_name, "mask": mask, "n_endpoints": len(endpoints), "n_cells": len(cells),
                     "n_null_unstable": sum(c.unstable for c in cells), "n_null_heavy_tail": sum(c.heavy_tail for c in cells),
                     "n_two_phase_cells": sum(len(c.S_phase) == 2 for c in cells), "n_phase_scramble_ks_fail": None,
                     "fwer": {k: v for k, v in fw.items() if k not in ("scale", "global_p", "R_norm_real")},
                     "n_candidates": len(cands), "candidates": cands,
                     "real_R_quantiles": {qq: float(np.quantile([c.R_real for c in cells], qq)) for qq in (0.5, 0.9, 0.99, 1.0)} if cells else {},
                     "n_real_R_gt_1": int(sum(1 for c in cells if c.R_real > 1)),
                     "expected_R_gt_1_from_ring": float(np.mean([(c.pooled > 1).mean() for c in cells]) * len(cells)) if cells else 0.0,
                     "holdout_null_calibration": {"note": "not applicable to the joint stage"},
                     "missing_tensors": [], "cells": per}
        print(f"[joint/{set_name}/{mask}] {len(cells)} cells, {res[mask]['n_null_unstable']} void; "
              f"R~_FWER = {fw['R_fwer']:.3f}; candidates {len(cands)} {cands[:5]}", flush=True)
    out_path.write_text(json.dumps(res, indent=1, default=float))


def joint_injections(e, role, jb):
    """Combine the two archives' per-injection window sums into joint
    injection records on the PS1 grid (primary mask)."""
    b_ps1, b_ztf = BAND_PAIRS[jb]
    jz, jp = ZTF.inj_dir / f"{e}__{role}.json", PS1.inj_dir / f"{e}__{role}.json"
    wz, wp = ZTF.inj_dir / f"{e}__{role}__{b_ztf}_windows.npz", PS1.inj_dir / f"{e}__{role}__{b_ps1}_windows.npz"
    if not (jz.exists() and jp.exists() and wz.exists() and wp.exists()):
        return None
    cz = json.loads(jz.read_text())["cells"].get(b_ztf); cp = json.loads(jp.read_text())["cells"].get(b_ps1)
    if not cz or not cp:
        return None
    Wz, Wp = np.load(wz), np.load(wp)
    qz, qp = 1.0 / ZTF.z_grid, 1.0 / PS1.z_grid
    n = min(len(cz["injections"]), len(cp["injections"]), len(Wz["lo"]), len(Wp["lo"]))
    out = []
    nm = len(PS1.mu_grid)
    for j in range(n):
        iz_, ip_ = cz["injections"][j], cp["injections"][j]
        if abs(iz_["z_au"] - ip_["z_au"]) > 1e-6 or abs(iz_["mag"] - ip_["mag"]) > 1e-6:
            continue  # draws diverged (should not happen with per-injection seeding)
        loz, lop = int(Wz["lo"][j]), int(Wp["lo"][j])
        wzn = int(Wz["z_window"]) * 2 + 1; wpn = int(Wp["z_window"]) * 2 + 1
        # ZTF window nodes -> PS1 window nodes (q interpolation)
        z_idx = np.arange(loz, min(loz + wzn, len(qz)))
        p_idx = np.arange(lop, min(lop + wpn, len(qp)))
        Az = interp_q(Wz["A"][j][:len(z_idx)][None].astype(float), qz[z_idx], qp[p_idx])[0]
        Bz = interp_q(Wz["B"][j][:len(z_idx)][None].astype(float), qz[z_idx], qp[p_idx])[0]
        nzj = np.rint(interp_q(Wz["n"][j][:len(z_idx)][None].astype(float), qz[z_idx], qp[p_idx])[0]).astype(int)
        A = Az + Wp["A"][j][:len(p_idx)]; B = Bz + Wp["B"][j][:len(p_idx)]; nn = nzj + Wp["n"][j][:len(p_idx)]
        S = _S(A, B, nn)
        iz = int(ip_["iz"])
        wz_ = np.abs(p_idx - iz)[:, None, None]
        mu = ip_["mu"]
        wmu = np.maximum(np.abs(np.arange(nm)[:, None] - np.argmin(np.abs(PS1.mu_grid - mu[0]))),
                         np.abs(np.arange(nm)[None, :] - np.argmin(np.abs(PS1.mu_grid - mu[1]))))[None]
        near = (wz_ <= RECOVERY_WINDOW[0]) & (wmu <= RECOVERY_WINDOW[1])
        s_peak, pk = nulls.grid_max(np.where(near, S, np.nan))
        s_max, _ = nulls.grid_max(S)
        stz = (iz_["masks"].get("primary") or {}).get("static") or {}
        stp = (ip_["masks"].get("primary") or {}).get("static") or {}
        # joint flux-consistency: the catalogued sources of both archives
        # together must account for the JOINT peak within the factor
        s_pred = float(stz.get("S_pred") or 0.0) + float(stp.get("S_pred") or 0.0)
        consistent = bool(s_peak is not None and np.isfinite(s_peak) and s_peak > 0 and s_pred >= s_peak / JP.static_factor)
        out.append({"j": j, "model": ip_["model"], "z_au": ip_["z_au"], "mu": ip_["mu"], "mag": ip_["mag"],
                    "fnu_jy": ip_["fnu_jy"], "iz": iz, "resp_median": None,
                    "n_on": {"ztf": iz_["n_on"], "ps1": ip_["n_on"]},
                    "masks": {"primary": {"S_peak": s_peak, "S_max": s_max, "static": {"consistent": consistent}}}})
    return {"endpoint": e, "role": role, "band": jb, "injections": out, "m90_v1": None,
            "T": {"primary": None}, "n_epochs": None}


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
            for jb in BAND_PAIRS:
                ck = f"{e}/{role}/{jb}"
                cinfo = ne["cells"].get(ck)
                if cinfo is None:
                    continue
                cell = joint_injections(e, role, jb)
                if cell is None or not cell["injections"]:
                    continue
                T = cinfo["T"]; q95 = cinfo["q95"]; unstable = cinfo["null_unstable"]
                for j in cell["injections"]:
                    r = j["masks"]["primary"]; r["T"] = T
                    r["R_peak"] = (r["S_peak"] / T) if (r["S_peak"] is not None and np.isfinite(r["S_peak"]) and T > 0) else None
                cls = [CS.classify(JP, j, "primary", q95, r_fwer) for j in cell["injections"]]
                b_ps1, b_ztf = BAND_PAIRS[jb]
                clips = []
                for P_, band in ((PS1, b_ps1), (ZTF, b_ztf)):
                    with np.load(P_.tensor_dir / f"{e}__{role}.npz") as dten:
                        m = CS.single_epoch_clip_mag(P_, dten, band)
                    if m is not None:
                        clips.append(m)
                m_clip = max(clips) if clips else None
                thr = CS.curves_for(JP, cell["injections"], cls, "threshold", rng, mag_min=m_clip)
                fin = CS.curves_for(JP, cell["injections"], cls, "final", rng, mag_min=m_clip)
                vetoes = {}
                for c in cls:
                    if c.get("veto"):
                        vetoes[c["veto"]] = vetoes.get(c["veto"], 0) + 1
                results[ck] = {"n_injections": len(cls), "n_usable": sum(c["usable"] for c in cls), "q95": q95, "T": T,
                               "bright_limit_mag": m_clip,
                               "null_unstable": unstable, "n_threshold": int(sum(c["threshold"] for c in cls)),
                               "n_final": int(sum(c["final"] for c in cls)), "vetoes_fired": vetoes,
                               "resp_median": np.nan, "threshold": thr, "final_candidate": fin, "m90_v1": None,
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
                            limit = {"value": w["m90"], "unit": "ab_mag", "band": jb, "m50": w["m50"],
                                     "fnu_jy_at_m90": JP.fnu_from_mag(jb, w["m90"]), "worst_temporal_model": worst,
                                     "per_model_m90": {k: v["m90"] for k, v in per_model.items()}}
                            ci = w["m90_ci68"]; n_inj = int(sum(per_model[k]["n_injections"] for k in per_model))
                        constraints.append(Constraint(
                            constraint_id=stable_id("con", {"endpoint": e, "role": role, "band": jb, "z_interval": list(zint),
                                                            "kind": kind, "hypothesis": HYPOTHESIS_VERSION,
                                                            "freeze": freeze["freeze_content_hash"], "set": set_name}),
                            analysis_run_id="pending", endpoint_id=e, role=role, hypothesis_version=HYPOTHESIS_VERSION,
                            z_interval_au=zint, band=jb, epoch_range_mjd=(0.0, 0.0), duty_cycle_range=(0.5, 1.0),
                            residual_motion_bound_arcsec_per_yr=1.0, morphology="point source (per-archive Moffat)",
                            kind=ckind, recovery_probability=0.9 if ckind == "recovery_curve" else None, flux_limit=limit,
                            trials_accounting_ref=f"nulls/{set_name}_null_ensemble.json",
                            extra={"set": set_name, "mask": "primary", "R_fwer_norm": r_fwer, "q95": q95,
                                   "exclusion_claim": bool(kind == "final_candidate" and ckind == "recovery_curve"),
                                   "archives": ["ps1", "ztf"]},
                            completeness_kind=kind, ci_68=tuple(ci) if ci else None, n_injections=n_inj,
                            prf_model="moffat beta=3 per archive", spectrum_model="flat_fnu_ab", motion_bound_norm="linf"))
    obs_hash = combined_hash(hashes)
    config = {"pipeline_id": "joint-ps1-ztf-v2-completeness", "set": set_name, "R_fwer_norm": r_fwer,
              "freeze_hash": JP.freeze_hash(), "ztf_freeze": ZTF.freeze_hash(), "ps1_freeze": PS1.freeze_hash()}
    run_rec = AnalysisRun(analysis_run_id=stable_id("run", {"config": json.loads(canonical_json(config)), "inputs": obs_hash}),
                          pipeline_id="joint-ps1-ztf-v2-completeness", pipeline_version="2.0.0", config=config,
                          observation_set_hash=obs_hash, intersection_set_hash=obs_hash,
                          registry_source_hash=registry.source_hash, hypothesis_version=HYPOTHESIS_VERSION,
                          environment={}, random_seeds={"bootstrap": JP.split_seed + 7},
                          started_utc=datetime.now(timezone.utc).isoformat(), finished_utc=datetime.now(timezone.utc).isoformat(),
                          output_files={"completeness": f"completeness/{set_name}_completeness.json"})
    constraints = [Constraint(**{**c.__dict__, "analysis_run_id": run_rec.analysis_run_id}) for c in constraints]
    append_records(OUT / "records" / "analysis_run.jsonl", [run_rec])
    append_records(OUT / "records" / "constraint.jsonl", constraints)
    (OUT / "completeness").mkdir(parents=True, exist_ok=True)
    summary = {"set": set_name, "mask": "primary", "analysis_run_id": run_rec.analysis_run_id, "R_fwer_norm": r_fwer,
               "n_cells": len(results), "n_constraints": len(constraints),
               "n_recovery_curve": sum(1 for c in constraints if c.kind == "recovery_curve"),
               "median_m90": {b: {kind: (float(np.nanmedian([c.flux_limit["value"] for c in constraints if c.band == b and c.completeness_kind == kind and c.flux_limit]))
                                         if any(c.band == b and c.completeness_kind == kind and c.flux_limit for c in constraints) else None)
                                  for kind in ("threshold", "final_candidate")} for b in BAND_PAIRS},
               "vetoes_fired_total": {}, "median_prf_throughput": None}
    for r in results.values():
        for k, v in r["vetoes_fired"].items():
            summary["vetoes_fired_total"][k] = summary["vetoes_fired_total"].get(k, 0) + v
    (OUT / "completeness" / f"{set_name}_completeness.json").write_text(json.dumps({"summary": summary, "cells": results}, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


def adjudicate_stage(set_name):
    ne = json.loads((OUT / "nulls" / f"{set_name}_null_ensemble.json").read_text())["primary"]
    freeze = JP.load_freeze()
    endpoints = freeze["split"]["development" if set_name == "dev" else "confirmatory"]["endpoints"]
    cells, meta = analyse_cells(endpoints, "primary") if ne["candidates"] else ([], {})
    by_key = {c.key: c for c in cells}
    cands = []
    for key in ne["candidates"]:
        cinfo = ne["cells"][key]; c = by_key[key]
        e, role, jb = key.split("/")
        b_ps1, b_ztf = BAND_PAIRS[jb]
        node = c.node
        z = float(PS1.z_grid[node[0]]); mu = (float(PS1.mu_grid[node[1]]), float(PS1.mu_grid[node[2]]))
        # per-archive flux-consistency at the joint node
        S_pred = 0.0; tests = {}
        for P_, band, dname in ((PS1, b_ps1, "ps1"), (ZTF, b_ztf, "ztf")):
            with np.load(P_.tensor_dir / f"{e}__{role}.npz") as d:
                b = P_.band_idx[band]; eb = d["band_idx"] == b; ok = d["masks"][:, 0][eb]
                f0, v0, g0 = d["f"][0, eb], d["v"][0, eb], d["g"][0, eb]
                q_self = 1.0 / P_.z_grid; zi = int(np.argmin(np.abs(q_self - 1.0 / z)))
                S_, A_, B_, n_, w_ = nulls.stack_S(f0, v0, g0, ok, return_parts=True, frame_cap=d["frame_cap"][eb], clip_sigma=P_.clip_sigma)
                we = w_[:, zi, node[1], node[2]]
                mjd = d["mjd"][eb]; phase = d["phase"][eb]; centre = tuple(d["centre"]); cosd = np.cos(np.deg2rad(centre[1]))
                dt_yr = (mjd - P_.t0_mjd) / 365.25; d0 = d["d0"][eb]
                ra_t = centre[0] + (d0[:, zi, 0] + mu[0] * dt_yr) / 3600.0 / cosd
                dec_t = centre[1] + (d0[:, zi, 1] + mu[1] * dt_yr) / 3600.0
                cat = P_.catalog(str(d["corridor"]), centre)
                S_arch = float(S_[zi, node[1], node[2]])
                if cat is not None and len(cat.ra):
                    fw = float(np.nanmedian(d["fwhm_arcsec"][eb])); pix = float(np.nanmedian(d["pix_arcsec"][eb]))
                    rt = response_table(P_, band, _FmLike(fw), fw, pix)
                    fc = flux_consistency(cat, ra_t, dec_t, we * ok, phase, P_.zp_ref, rt[0], rt[1], S_arch,
                                          search_arcsec=P_.static_search_arcsec, factor=P_.static_factor, min_ndet=P_.catalog_min_ndet)
                    S_pred += fc["S_pred"]; tests[dname] = {"S_arch": S_arch, "S_pred": fc["S_pred"], "sources": fc["sources"][:3]}
                else:
                    tests[dname] = {"S_arch": S_arch, "S_pred": None}
        consistent = S_pred >= c.S_max / JP.static_factor
        rejection = "static_flux_consistent" if consistent else None
        cands.append(Candidate(
            candidate_id=stable_id("cnd", {"cell": key, "hypothesis": HYPOTHESIS_VERSION, "freeze": JP.freeze_hash()}),
            analysis_run_id="see-null-ensemble", endpoint_id=e, role=role, observation_ids=(),
            fitted_z_au=z, fitted_residual_motion={"mu_ra": mu[0], "mu_dec": mu[1]},
            model_comparison={"S_max": c.S_max, "T": c.T, "R": c.R_real, "q95": cinfo["q95"], "R_norm": cinfo["R_norm"],
                              "tests": {"static_flux_consistent": {"S_pred_joint": S_pred, "consistent": consistent, "per_archive": tests}}},
            status=f"rejected:{rejection}" if rejection else "retained-ambiguous", veto_reason=rejection, rejection_test=rejection,
            global_p_value=cinfo["global_p"], rank_statement=cinfo["rank"],
            annotations={"parallax_phase": parallax_phase_test(c.S_phase), "p_phase": c.p_phase, "p_trajectory": c.p_trajectory,
                         "n_epochs": c.n_epochs},
            extra={"band": jb, "set": set_name, "mask": "primary"}))
        print(f"{key}: R~={cinfo['R_norm']:.2f} -> {cands[-1].status}")
    out = OUT / "records" / f"candidate_{set_name}.jsonl"
    if out.exists():
        out.unlink()
    append_records(out, cands)
    summary = {"survey": "joint", "set": set_name, "mask": "primary", "R_fwer_norm": ne["fwer"]["R_fwer"],
               "n_candidates": len(cands), "status_counts": {s: sum(1 for c in cands if c.status == s) for s in set(c.status for c in cands)},
               "retained_ambiguous": [f"{c.endpoint_id}/{c.role}/{c.extra['band']}" for c in cands if c.status == "retained-ambiguous"],
               "written": str(out), "utc": datetime.now(timezone.utc).isoformat()}
    (OUT / "records" / f"adjudication_{set_name}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


def freeze():
    """The joint split is the PS1/ZTF split (identical by construction)."""
    fz, fp = ZTF.load_freeze(), PS1.load_freeze()
    assert fz["split"]["development"]["corridors"] == fp["split"]["development"]["corridors"], "ZTF/PS1 splits differ"
    fr = {"survey": "joint-ps1-ztf", "hypothesis_version": HYPOTHESIS_VERSION, "ztf_freeze": fz["freeze_content_hash"],
          "ps1_freeze": fp["freeze_content_hash"], "band_pairs": BAND_PAIRS, "t0_mjd": 59800.0, "zp_ref": 25.0,
          "rule": "joint S = (A_ztf + A_ps1)/sqrt(B_ztf + B_ps1) on the PS1 grid (ZTF interpolated in 1/z); ring-only null; R~ = R/q95; FWER alpha 0.05 over the joint family; static flux-consistency from both catalogues",
          "split": fp["split"], "frozen_at": "2026-08-22"}
    fr["freeze_content_hash"] = "sha256:" + hashlib.sha256(canonical_json(fr).encode()).hexdigest()
    (SURVEY_DIR / "configs").mkdir(parents=True, exist_ok=True)
    (SURVEY_DIR / "configs" / "v2_freeze.json").write_text(json.dumps(fr, indent=1, sort_keys=True) + "\n")
    print("joint freeze", fr["freeze_content_hash"][:23])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["freeze", "nulls", "completeness", "adjudicate", "report"])
    ap.add_argument("--set", choices=["dev", "confirmatory"], default="dev")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if a.stage == "freeze":
        freeze()
    elif a.stage == "nulls":
        nulls_stage(a.set, a.force)
    elif a.stage == "completeness":
        completeness_stage(a.set)
    elif a.stage == "adjudicate":
        adjudicate_stage(a.set)
    elif a.stage == "report":
        from sglsurvey import report_stage
        JP.results_dir = SURVEY_DIR / "results"
        JP.null_dir = OUT / "nulls"
        JP.v1_run_dir = REPO / "runs" / "joint"
        raise SystemExit(report_stage.run(JP))


if __name__ == "__main__":
    main()
