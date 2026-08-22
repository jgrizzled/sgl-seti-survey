"""Step E / G: completeness curves and Constraint records from the
image-level injections (hypotheses v2.0 §4).

Per cell and temporal model, each injection is classified with the
set's frozen rule:

  threshold-recovered   R~_peak = (S_peak / T) / q95_cell >= R~_FWER, where
                        S_peak is the maximum within 2 z-nodes and 1
                        mu-node of the injection;
  final-candidate       threshold-recovered AND not rejected by a
                        calibrated veto run blind on it: the flux-
                        consistent static-source test (a catalogued
                        star accounts for the injected peak) or the
                        held-out-epoch prediction test (applicable
                        when >= 5 late epochs exist).

A logistic completeness model P(m, z) = 1 / (1 + exp((m - m50(z)) / w)),
m50(z) = a + b (ln z - ln z_mid), is fitted per cell x model x kind by
maximum likelihood; m90 and m50 per z-interval (8 reciprocal-distance
intervals tiling 550-10,000 AU with no gap) come from the fit at the
interval's log-midpoint with bootstrap 68 % intervals; the worst of
the four temporal models defines the duty >= 0.5 coverage. Cells whose
null is unstable, or with fewer than 20 usable injections per model,
are `not_constrainable`. W3/W4 records carry completeness_kind =
threshold only and exclusion_claim = false (hypotheses v2.0 §4).

Writes runs/wise-v2/records/{analysis_run,constraint}.jsonl (append),
runs/wise-v2/completeness/<set>_completeness.json, and the v1
supersession links.

Usage: uv run python surveys/wise-v2/scripts/completeness.py --set dev|confirmatory
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

from sglseti import canonical_json, load_target_registry, stable_hash, stable_id  # noqa: E402
from sglsurvey import inject, nulls  # noqa: E402
from sglsurvey.manifest import combined_hash  # noqa: E402
from sglsurvey.records import AnalysisRun, Constraint, append_records, read_records  # noqa: E402
from sglsurvey.vetting import holdout_prediction_test  # noqa: E402

MIN_INJ_PER_MODEL = 20
N_BOOT = 200
PIPELINE_VERSION = "2.0.0"


def interval_edges(n_intervals: int = C.N_Z_INTERVALS) -> np.ndarray:
    """Reciprocal-distance midpoint edges tiling [550, 10000] (ascending z)."""
    q = np.sort(1.0 / C.Z_GRID)           # ascending q = descending z
    blocks = np.array_split(np.arange(len(q)), n_intervals)
    edges_q = [q[0]] + [0.5 * (q[b[-1]] + q[b[-1] + 1]) for b in blocks[:-1]] + [q[-1]]
    z = np.sort(1.0 / np.array(edges_q))
    z[0], z[-1] = 550.0, 10000.0
    return z


def classify(inj: dict, mask: str, q95: float, r_fwer: float) -> dict:
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
                                      C.HOLDOUT_MIN_S_LATE, C.HOLDOUT_MIN_FLUX_RATIO)
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


def curves_for(injs: list[dict], cls: list[dict], kind: str, rng: np.random.Generator):
    """m90/m50 per z-interval with bootstrap 68 % CIs, per temporal model."""
    edges = interval_edges()
    lnz_mid = 0.5 * (np.log(550.0) + np.log(10000.0))
    out = {}
    for model in C.TEMPORAL_MODELS:
        idx = [i for i, (j, c) in enumerate(zip(injs, cls)) if j["model"] == model and c["usable"]]
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["dev", "confirmatory"], required=True)
    ap.add_argument("--mask", default="primary")
    a = ap.parse_args()
    freeze = C.load_freeze()
    key = "development" if a.set == "dev" else "confirmatory"
    endpoints = freeze["split"][key]["endpoints"]
    ne = json.loads((C.NULL_DIR / f"{a.set}_null_ensemble.json").read_text())[a.mask]
    r_fwer = ne["fwer"]["R_fwer"]
    registry = load_target_registry(C.REGISTRY_PATH)
    rng = np.random.default_rng(C.SPLIT_SEED + 7)
    started = datetime.now(timezone.utc).isoformat()
    # v1 constraints for supersession links
    v1 = {}
    v1_path = C.V1_CAL / "records" / "constraint.jsonl"
    if v1_path.exists():
        for r in read_records(v1_path):
            if r["analysis_run_id"] == "run-1b2e86e9219e":
                v1.setdefault((r["endpoint_id"], r["role"], r["band"]), []).append(r)
    results, input_hashes, constraints = {}, [], []
    edges = interval_edges()
    for e in endpoints:
        for role in ("rx", "tx"):
            p = C.INJ_DIR / f"{e}__{role}.json"
            if not p.exists():
                continue
            inj = json.loads(p.read_text())
            input_hashes.append(inj["tensor_input_hash"])
            # cross-track variants: the statistic is the max over offsets
            for vp in sorted(C.INJ_DIR.glob(f"{e}__{role}__xt*.json")):
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
                cls = [classify(j, a.mask, q95, r_fwer) for j in cell["injections"]]
                thr = curves_for(cell["injections"], cls, "threshold", rng)
                fin = curves_for(cell["injections"], cls, "final", rng)
                n_us = sum(c["usable"] for c in cls)
                vetoes = {}
                hold_n = hold_pass = 0
                for c in cls:
                    if c.get("veto"):
                        vetoes[c["veto"]] = vetoes.get(c["veto"], 0) + 1
                    if c["threshold"] and c.get("holdout_pass") is not None:
                        hold_n += 1; hold_pass += c["holdout_pass"]
                hold_by_model = {}
                for model in C.TEMPORAL_MODELS:
                    sel = [c for j, c in zip(cell["injections"], cls) if j["model"] == model and c["threshold"] and c.get("holdout_pass") is not None]
                    hold_by_model[model] = {"n": len(sel), "pass": sum(c["holdout_pass"] for c in sel)}
                results[ck] = {
                    "n_injections": len(cls), "n_usable": n_us, "q95": q95, "T": cell["T"].get(a.mask),
                    "null_unstable": unstable, "m90_v1": cell["m90_v1"],
                    "n_threshold": int(sum(c["threshold"] for c in cls)),
                    "n_final": int(sum(c["final"] for c in cls)), "vetoes_fired": vetoes,
                    "holdout_annotation": {"n": hold_n, "pass": hold_pass, "by_model": hold_by_model},
                    "resp_median": float(np.nanmedian([j["resp_median"] or np.nan for j in cell["injections"]])),
                    "threshold": thr, "final_candidate": fin,
                }
                # ---- Constraint records ------------------------------------
                ep_range = None
                for kind, curves in (("threshold", thr), ("final_candidate", fin)):
                    if band in ("W3", "W4") and kind == "final_candidate":
                        continue
                    for i in range(len(edges) - 1):
                        zint = (round(float(edges[i]), 1), round(float(edges[i + 1]), 1))
                        per_model = {}
                        for model, res in curves.items():
                            if res.get("status") == "ok" and res["intervals"][i]["m90"] is not None:
                                per_model[model] = res["intervals"][i]
                        worst = None
                        if per_model and not unstable:
                            worst = min(per_model, key=lambda k: per_model[k]["m90"])
                        if worst is None:
                            ckind, limit, ci, n_inj = "not_constrainable", None, None, int(sum(
                                res["intervals"][i]["n_injections"] for res in curves.values() if res.get("intervals")))
                            reason = "null_unstable" if unstable else "insufficient_recovery_fit"
                        else:
                            w = per_model[worst]
                            ckind = "recovery_curve"
                            limit = {"value": w["m90"], "unit": "wise_vega_mag", "band": band,
                                     "m50": w["m50"],
                                     "fnu_jy_at_m90": inject.fnu_from_vega_mag(band, w["m90"], C.SPECTRUM[band]),
                                     "worst_temporal_model": worst,
                                     "per_model_m90": {k: v["m90"] for k, v in per_model.items()}}
                            ci = w["m90_ci68"]
                            n_inj = int(sum(per_model[k]["n_injections"] for k in per_model))
                            reason = None
                        sup = None
                        for r in v1.get((e, role, band), []):
                            lo, hi = r["z_interval_au"]
                            if lo < zint[1] and hi > zint[0]:
                                sup = r["constraint_id"]; break
                        constraints.append(Constraint(
                            constraint_id=stable_id("con", {
                                "endpoint": e, "role": role, "band": band, "z_interval": list(zint),
                                "kind": kind, "hypothesis": C.HYPOTHESIS_VERSION,
                                "freeze": freeze["freeze_content_hash"], "set": a.set}),
                            analysis_run_id="pending", endpoint_id=e, role=role,
                            hypothesis_version=C.HYPOTHESIS_VERSION, z_interval_au=zint, band=band,
                            epoch_range_mjd=(0.0, 0.0), duty_cycle_range=(0.5, 1.0),
                            residual_motion_bound_arcsec_per_yr=1.0,
                            morphology="wise-empirical-prf-point-source", kind=ckind,
                            recovery_probability=0.9 if ckind == "recovery_curve" else None,
                            flux_limit=limit, false_alarm_rate=None,
                            trials_accounting_ref=f"nulls/{a.set}_null_ensemble.json",
                            extra={"set": a.set, "mask": a.mask, "R_fwer_norm": r_fwer, "q95": q95,
                                   "exclusion_claim": bool(kind == "final_candidate" and ckind == "recovery_curve"),
                                   "not_constrainable_reason": reason,
                                   "w34_note": ("raw threshold sensitivity, 300 K blackbody, empirical PRF; "
                                                "no exclusion claim") if band in ("W3", "W4") else None},
                            completeness_kind=kind, ci_68=tuple(ci) if ci else None, n_injections=n_inj,
                            prf_model="irsa-pass2-psf-wpro-09x09", spectrum_model=C.SPECTRUM[band],
                            motion_bound_norm="linf", supersedes=sup))
    # ---- AnalysisRun with real content hashes --------------------------------
    obs_hash = combined_hash(input_hashes)
    config = {"pipeline_id": "wise-v2-completeness", "set": a.set, "mask": a.mask,
              "freeze_hash": C.freeze_hash(), "freeze_content_hash": freeze["freeze_content_hash"],
              "R_fwer_norm": r_fwer, "parameters": freeze["parameters"]["injections"],
              "recovery_rule": "R~_peak >= R~_FWER within 2 z / 1 mu nodes; final = and no calibrated veto",
              "logistic_model": "P = 1/(1+exp((m - a - b(ln z - ln z_mid))/w)), MLE, bootstrap 200"}
    try:
        commit = subprocess.run(["git", "-C", str(C.REPO.parent / "sglseti"), "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:
        commit = "unknown"
    run = AnalysisRun(
        analysis_run_id=stable_id("run", {"config": json.loads(canonical_json(config)),
                                          "inputs": obs_hash, "registry": registry.source_hash}),
        pipeline_id="wise-v2-completeness", pipeline_version=PIPELINE_VERSION, config=config,
        observation_set_hash=obs_hash, intersection_set_hash=obs_hash,
        registry_source_hash=registry.source_hash, hypothesis_version=C.HYPOTHESIS_VERSION,
        environment={"sglseti_commit": commit, "hypotheses_hash": freeze["hypotheses_hash"]},
        random_seeds={"bootstrap": C.SPLIT_SEED + 7, "injections": "sha256(inj/endpoint/role/seed)"},
        started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
        output_files={"completeness": f"completeness/{a.set}_completeness.json"})
    constraints = [Constraint(**{**c.__dict__, "analysis_run_id": run.analysis_run_id}) for c in constraints]
    rec_dir = C.RUN_DIR / "records"
    append_records(rec_dir / "analysis_run.jsonl", [run])
    append_records(rec_dir / "constraint.jsonl", constraints)
    (C.RUN_DIR / "completeness").mkdir(parents=True, exist_ok=True)
    summary = {
        "set": a.set, "mask": a.mask, "analysis_run_id": run.analysis_run_id, "R_fwer_norm": r_fwer,
        "n_cells": len(results), "n_constraints": len(constraints),
        "n_recovery_curve": sum(1 for c in constraints if c.kind == "recovery_curve"),
        "median_m90": {
            b: {kind: (float(np.nanmedian([c.flux_limit["value"] for c in constraints
                                           if c.band == b and c.completeness_kind == kind and c.flux_limit]))
                       if any(c.band == b and c.completeness_kind == kind and c.flux_limit for c in constraints) else None)
                for kind in ("threshold", "final_candidate")} for b in ("W1", "W2", "W3", "W4")},
        "vetoes_fired_total": {},
        "median_prf_throughput": float(np.nanmedian([r["resp_median"] for r in results.values()])),
    }
    for r in results.values():
        for k, v in r["vetoes_fired"].items():
            summary["vetoes_fired_total"][k] = summary["vetoes_fired_total"].get(k, 0) + v
    (C.RUN_DIR / "completeness" / f"{a.set}_completeness.json").write_text(
        json.dumps({"summary": summary, "cells": results}, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
