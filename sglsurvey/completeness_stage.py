"""v2 engine: completeness curves and Constraint records (profile-
parameterised port of surveys/wise/scripts/completeness.py)."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone

import numpy as np

from sglseti import canonical_json, load_target_registry, stable_id

from sglsurvey.manifest import combined_hash
from sglsurvey.records import AnalysisRun, Constraint, append_records, read_records
from sglsurvey.vetting import holdout_prediction_test

MIN_INJ_PER_MODEL = 20
N_BOOT = 200
PIPELINE_VERSION = "2.0.0"



def interval_edges(P) -> np.ndarray:
    """Reciprocal-distance midpoint edges tiling [550, 10000] (ascending z)."""
    n_intervals = P.n_z_intervals
    q = np.sort(1.0 / P.z_grid)           # ascending q = descending z
    blocks = np.array_split(np.arange(len(q)), n_intervals)
    edges_q = [q[0]] + [0.5 * (q[b[-1]] + q[b[-1] + 1]) for b in blocks[:-1]] + [q[-1]]
    z = np.sort(1.0 / np.array(edges_q))
    z[0], z[-1] = 550.0, 10000.0
    return z


def classify(P, inj: dict, mask: str, q95: float, r_fwer: float) -> dict:
    r = inj["masks"].get(mask)
    out = {"threshold": False, "final": False, "veto": None, "usable": False}
    if not r or r.get("R_peak") is None or q95 is None or not np.isfinite(q95) or q95 <= 0:
        return out
    out["usable"] = True
    rn = r["R_peak"] / q95
    out["R_norm"] = rn
    out["threshold"] = bool(rn >= r_fwer)
    if not out["threshold"]:
        return out
    veto = None
    st = r.get("static")
    if st and st.get("consistent"):
        veto = "static_flux_consistent"
    h = r.get("holdout")
    out["holdout_pass"] = None
    if h and h.get("S_late") is not None and h.get("f_late") is not None:
        res = holdout_prediction_test(h["S_early"], h["f_early"] if h["f_early"] is not None else np.nan,
                                      h["S_late"], h["f_late"], int(h["n_late"]),
                                      P.holdout_min_s_late, P.holdout_min_flux_ratio)
        if res.get("applicable"):
            out["holdout_pass"] = bool(res["pass"])   # annotation (v2.1): measured, not a veto
    out["veto"] = veto
    out["final"] = veto is None
    return out


def fit_logistic(m: np.ndarray, lnz: np.ndarray, y: np.ndarray, lnz_mid: float):
    """MLE of (a, b, w) for P = 1/(1+exp((m - a - b (lnz - lnz_mid)) / w))."""
    from scipy.optimize import minimize

    x = lnz - lnz_mid

    def nll(p):
        a, b, lw = p
        w = np.exp(lw)
        t = (m - a - b * x) / w
        t = np.clip(t, -50, 50)
        logp = -np.logaddexp(0, t)          # log P
        logq = -np.logaddexp(0, -t)         # log (1-P)
        return -(y * logp + (1 - y) * logq).sum()

    if y.sum() == 0:
        return None
    if y.sum() == len(y):
        return ("all", None, None)
    a0 = float(np.median(m))
    best = None
    for a_init in (a0 - 1, a0, a0 + 1):
        res = minimize(nll, [a_init, 0.0, np.log(0.3)], method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-4, "fatol": 1e-6})
        if best is None or res.fun < best.fun:
            best = res
    a, b, lw = best.x
    return (float(a), float(b), float(np.exp(lw)))


def m_at(fit, lnz, lnz_mid, p):
    """Magnitude at which completeness equals p (brighter = smaller)."""
    if fit is None:
        return np.nan
    if fit[0] == "all":
        return np.nan
    a, b, w = fit
    return a + b * (lnz - lnz_mid) - w * np.log(p / (1 - p))


def single_epoch_clip_mag(P, d, band):
    """Magnitude at which a source is clipped from the stack as a
    single-epoch detection (|S_e| > clip_sigma) in a typical epoch of the
    cell: ZP_REF - 2.5 log10(clip_sigma x median sigma_e at the mu = 0
    nodes). None when the profile has no clip."""
    if P.clip_sigma is None:
        return None
    b = P.band_idx[band]
    eb = d["band_idx"] == b
    v = d["v"][0, eb][:, :, len(P.mu_grid) // 2, len(P.mu_grid) // 2].astype(float)
    sig = np.sqrt(v[np.isfinite(v) & (v > 0)])
    if sig.size == 0:
        return None
    return float(P.zp_ref - 2.5 * np.log10(P.clip_sigma * np.median(sig)))


def curves_for(P, injs: list[dict], cls: list[dict], kind: str, rng: np.random.Generator,
               mag_min: float | None = None):
    """m90/m50 per z-interval with bootstrap 68 % CIs, per temporal model.
    ``mag_min``: the single-epoch clip limit — brighter injections are
    clipped out of the stack by the layered-search rule (they belong to
    the catalogue layer) and are excluded from the logistic fit."""
    edges = interval_edges(P)
    lnz_mid = 0.5 * (np.log(550.0) + np.log(10000.0))
    out = {}
    for model in P.temporal_models:
        idx = [i for i, (j, c) in enumerate(zip(injs, cls)) if j["model"] == model and c["usable"]
               and (mag_min is None or j["mag"] >= mag_min - 0.5)]
        m = np.array([injs[i]["mag"] for i in idx])
        lnz = np.log(np.array([injs[i]["z_au"] for i in idx]))
        y = np.array([cls[i][kind] for i in idx], dtype=float)
        res = {"n": int(len(idx)), "n_recovered": int(y.sum()), "intervals": []}
        if len(idx) < MIN_INJ_PER_MODEL:
            res["status"] = "too_few"
            out[model] = res
            continue
        fit = fit_logistic(m, lnz, y, lnz_mid)
        boots = []
        for _ in range(N_BOOT):
            k = rng.integers(0, len(idx), size=len(idx))
            boots.append(fit_logistic(m[k], lnz[k], y[k], lnz_mid))
        for i in range(len(edges) - 1):
            z0, z1 = edges[i], edges[i + 1]
            lz = 0.5 * (np.log(z0) + np.log(z1))
            m90 = m_at(fit, lz, lnz_mid, 0.9)
            m50 = m_at(fit, lz, lnz_mid, 0.5)
            # a fit extrapolated outside the injected magnitude range is
            # not a measurement: report it as unconstrained
            lo_m, hi_m = m.min() - 3.0, m.max() + 3.0
            if mag_min is not None:
                # the stack cannot be 90 % complete brighter than the
                # single-epoch clip: such a fit is an extrapolation into
                # the catalogue layer's regime
                lo_m = max(lo_m, mag_min - 0.5)
            if not (lo_m <= m90 <= hi_m):
                m90 = np.nan
            if not (lo_m <= m50 <= hi_m):
                m50 = np.nan
            b90 = np.array([m_at(b, lz, lnz_mid, 0.9) for b in boots if b is not None])
            b90 = b90[np.isfinite(b90) & (b90 >= lo_m) & (b90 <= hi_m)]
            in_int = (np.exp(lnz) >= z0) & (np.exp(lnz) < z1)
            res["intervals"].append({
                "z_interval_au": [round(float(z0), 1), round(float(z1), 1)],
                "m90": None if not np.isfinite(m90) else round(float(m90), 3),
                "m50": None if not np.isfinite(m50) else round(float(m50), 3),
                "m90_ci68": ([round(float(np.percentile(b90, 16)), 3), round(float(np.percentile(b90, 84)), 3)]
                             if len(b90) >= 20 else None),
                "n_injections": int(in_int.sum()),
                "n_recovered": int(y[in_int].sum()),
                "frac_recovered": (round(float(y[in_int].mean()), 3) if in_int.any() else None),
            })
        res["fit"] = None if fit is None or fit[0] == "all" else {"a": fit[0], "b": fit[1], "w": fit[2]}
        res["status"] = "ok" if res["fit"] else ("all_recovered" if fit and fit[0] == "all" else "none_recovered")
        # binned empirical curve for the report
        bins = np.linspace(m.min(), m.max(), 9)
        res["binned"] = [{"mag": round(float(0.5 * (bins[k] + bins[k + 1])), 2),
                          "n": int(((m >= bins[k]) & (m < bins[k + 1])).sum()),
                          "frac": (round(float(y[(m >= bins[k]) & (m < bins[k + 1])].mean()), 3)
                                   if ((m >= bins[k]) & (m < bins[k + 1])).any() else None)}
                         for k in range(8)]
        out[model] = res
    return out


def run(P, set_name: str, mask: str = "primary") -> None:
    freeze = P.load_freeze()
    key = "development" if set_name == "dev" else "confirmatory"
    endpoints = freeze["split"][key]["endpoints"]
    ne = json.loads((P.null_dir / f"{set_name}_null_ensemble.json").read_text())[mask]
    r_fwer = ne["fwer"]["R_fwer"]
    registry = load_target_registry(P.registry_path)
    rng = np.random.default_rng(P.split_seed + 7)
    started = datetime.now(timezone.utc).isoformat()
    # v1 constraints (latest run in the v1 ledger) for supersession links
    v1 = {}
    v1_path = P.v1_run_dir / P.extra_params.get("v1_calib", "calib_v1") / "records" / "constraint.jsonl"
    if v1_path.exists():
        recs = read_records(v1_path)
        if recs:
            last_run = recs[-1]["analysis_run_id"]
            for r in recs:
                if r["analysis_run_id"] == last_run:
                    v1.setdefault((r["endpoint_id"], r["role"], r["band"]), []).append(r)
    results, input_hashes, constraints = {}, [], []
    edges = interval_edges(P)
    threshold_only = set(P.extra_params.get("threshold_only_bands", ()))
    for e in endpoints:
        for role in ("rx", "tx"):
            p = P.inj_dir / f"{e}__{role}.json"
            if not p.exists():
                continue
            inj = json.loads(p.read_text())
            input_hashes.append(inj["tensor_input_hash"])
            for vp in sorted(P.inj_dir.glob(f"{e}__{role}__xt*.json")):
                v = json.loads(vp.read_text())
                for band, cell in v["cells"].items():
                    base = inj["cells"].get(band)
                    if not base:
                        continue
                    for j0, j1 in zip(base["injections"], cell["injections"]):
                        for m, r1 in j1["masks"].items():
                            r0 = j0["masks"].get(m)
                            if r0 is None or r1.get("S_peak") is None:
                                continue
                            if r0.get("S_peak") is None or r1["S_peak"] > r0["S_peak"]:
                                j0["masks"][m] = r1
                    base["n_cross_track"] = base.get("n_cross_track", 1) + 1
            for band, cell in inj["cells"].items():
                ck = f"{e}/{role}/{band}"
                cinfo = ne["cells"].get(ck)
                q95 = cinfo["q95"] if cinfo else None
                unstable = bool(cinfo and cinfo["null_unstable"])
                cls = [classify(P, j, mask, q95, r_fwer) for j in cell["injections"]]
                with np.load(P.tensor_dir / f"{e}__{role}.npz") as dten:
                    m_clip = single_epoch_clip_mag(P, dten, band)
                thr = curves_for(P, cell["injections"], cls, "threshold", rng, mag_min=m_clip)
                fin = curves_for(P, cell["injections"], cls, "final", rng, mag_min=m_clip)
                n_us = sum(c["usable"] for c in cls)
                vetoes = {}
                hold_n = hold_pass = 0
                for c in cls:
                    if c.get("veto"):
                        vetoes[c["veto"]] = vetoes.get(c["veto"], 0) + 1
                    if c["threshold"] and c.get("holdout_pass") is not None:
                        hold_n += 1; hold_pass += c["holdout_pass"]
                hold_by_model = {}
                for model in P.temporal_models:
                    sel = [c for j, c in zip(cell["injections"], cls) if j["model"] == model and c["threshold"] and c.get("holdout_pass") is not None]
                    hold_by_model[model] = {"n": len(sel), "pass": sum(c["holdout_pass"] for c in sel)}
                results[ck] = {
                    "n_injections": len(cls), "n_usable": n_us, "q95": q95, "T": cell["T"].get(mask),
                    "bright_limit_mag": m_clip,
                    "null_unstable": unstable, "m90_v1": cell["m90_v1"],
                    "n_threshold": int(sum(c["threshold"] for c in cls)),
                    "n_final": int(sum(c["final"] for c in cls)), "vetoes_fired": vetoes,
                    "holdout_annotation": {"n": hold_n, "pass": hold_pass, "by_model": hold_by_model},
                    "resp_median": float(np.nanmedian([j["resp_median"] or np.nan for j in cell["injections"]])),
                    "threshold": thr, "final_candidate": fin,
                }
                for kind, curves in (("threshold", thr), ("final_candidate", fin)):
                    if band in threshold_only and kind == "final_candidate":
                        continue
                    for i in range(len(edges) - 1):
                        zint = (round(float(edges[i]), 1), round(float(edges[i + 1]), 1))
                        per_model = {}
                        for model, res in curves.items():
                            if res.get("status") == "ok" and res["intervals"][i]["m90"] is not None:
                                per_model[model] = res["intervals"][i]
                        worst = min(per_model, key=lambda k: per_model[k]["m90"]) if (per_model and not unstable) else None
                        if worst is None:
                            ckind, limit, ci = "not_constrainable", None, None
                            n_inj = int(sum(res["intervals"][i]["n_injections"] for res in curves.values() if res.get("intervals")))
                            reason = "null_unstable" if unstable else "insufficient_recovery_fit"
                        else:
                            w = per_model[worst]
                            ckind, reason = "recovery_curve", None
                            limit = {"value": w["m90"], "unit": f"{P.mag_system}_mag", "band": band, "m50": w["m50"],
                                     "fnu_jy_at_m90": P.fnu_from_mag(band, w["m90"]),
                                     "worst_temporal_model": worst,
                                     "per_model_m90": {k: v["m90"] for k, v in per_model.items()}}
                            ci = w["m90_ci68"]
                            n_inj = int(sum(per_model[k]["n_injections"] for k in per_model))
                        sup = None
                        for r in v1.get((e, role, band), []):
                            lo, hi = r["z_interval_au"]
                            if lo < zint[1] and hi > zint[0]:
                                sup = r["constraint_id"]; break
                        constraints.append(Constraint(
                            constraint_id=stable_id("con", {"endpoint": e, "role": role, "band": band,
                                                            "z_interval": list(zint), "kind": kind,
                                                            "hypothesis": P.hypothesis_version,
                                                            "freeze": freeze["freeze_content_hash"], "set": set_name}),
                            analysis_run_id="pending", endpoint_id=e, role=role,
                            hypothesis_version=P.hypothesis_version, z_interval_au=zint, band=band,
                            epoch_range_mjd=(0.0, 0.0), duty_cycle_range=(0.5, 1.0),
                            residual_motion_bound_arcsec_per_yr=1.0,
                            morphology=P.extra_params.get("psf", "point source"), kind=ckind,
                            recovery_probability=0.9 if ckind == "recovery_curve" else None,
                            flux_limit=limit, false_alarm_rate=None,
                            trials_accounting_ref=f"nulls/{set_name}_null_ensemble.json",
                            extra={"set": set_name, "mask": mask, "R_fwer_norm": r_fwer, "q95": q95,
                                   "exclusion_claim": bool(kind == "final_candidate" and ckind == "recovery_curve"),
                                   "not_constrainable_reason": reason,
                                   "bright_limit_mag": m_clip,
                                   "bright_limit_note": (None if m_clip is None else
                                                         "sources brighter than this are single-epoch detections "
                                                         "clipped from the stack (catalogue-screening layer)")},
                            completeness_kind=kind, ci_68=tuple(ci) if ci else None, n_injections=n_inj,
                            prf_model=P.extra_params.get("psf", "unknown"), spectrum_model=P.spectrum.get(band),
                            motion_bound_norm="linf", supersedes=sup))
    obs_hash = combined_hash(input_hashes)
    config = {"pipeline_id": f"{P.name}-v2-completeness", "set": set_name, "mask": mask,
              "freeze_hash": P.freeze_hash(), "freeze_content_hash": freeze["freeze_content_hash"],
              "R_fwer_norm": r_fwer, "parameters": freeze["parameters"]["injections"],
              "recovery_rule": "R~_peak >= R~_FWER within the recovery window; final = and no calibrated veto",
              "logistic_model": "P = 1/(1+exp((m - a - b(ln z - ln z_mid))/w)), MLE, bootstrap 200"}
    try:
        commit = subprocess.run(["git", "-C", str(P.registry_path.parents[1].parent / "sglseti"), "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:
        commit = "unknown"
    run_rec = AnalysisRun(
        analysis_run_id=stable_id("run", {"config": json.loads(canonical_json(config)), "inputs": obs_hash,
                                          "registry": registry.source_hash}),
        pipeline_id=f"{P.name}-v2-completeness", pipeline_version=PIPELINE_VERSION, config=config,
        observation_set_hash=obs_hash, intersection_set_hash=obs_hash,
        registry_source_hash=registry.source_hash, hypothesis_version=P.hypothesis_version,
        environment={"sglseti_commit": commit, "hypotheses_hash": freeze.get("hypotheses_hash")},
        random_seeds={"bootstrap": P.split_seed + 7, "injections": "sha256(inj/endpoint/role/seed)"},
        started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
        output_files={"completeness": f"completeness/{set_name}_completeness.json"})
    constraints = [Constraint(**{**c.__dict__, "analysis_run_id": run_rec.analysis_run_id}) for c in constraints]
    rec_dir = P.run_dir / "records"
    append_records(rec_dir / "analysis_run.jsonl", [run_rec])
    append_records(rec_dir / "constraint.jsonl", constraints)
    (P.run_dir / "completeness").mkdir(parents=True, exist_ok=True)
    summary = {
        "set": set_name, "mask": mask, "analysis_run_id": run_rec.analysis_run_id, "R_fwer_norm": r_fwer,
        "n_cells": len(results), "n_constraints": len(constraints),
        "n_recovery_curve": sum(1 for c in constraints if c.kind == "recovery_curve"),
        "median_m90": {b: {kind: (float(np.nanmedian([c.flux_limit["value"] for c in constraints
                                                      if c.band == b and c.completeness_kind == kind and c.flux_limit]))
                                  if any(c.band == b and c.completeness_kind == kind and c.flux_limit for c in constraints) else None)
                           for kind in ("threshold", "final_candidate")} for b in P.bands},
        "vetoes_fired_total": {},
        "median_prf_throughput": float(np.nanmedian([r["resp_median"] for r in results.values()])) if results else None,
    }
    for r in results.values():
        for k, v in r["vetoes_fired"].items():
            summary["vetoes_fired_total"][k] = summary["vetoes_fired_total"].get(k, 0) + v
    (P.run_dir / "completeness" / f"{set_name}_completeness.json").write_text(
        json.dumps({"summary": summary, "cells": results}, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))
