"""Stage 2 v2: joint PS1 + ZTF stack over the FULL (z, mu) hypothesis
family, for tensors built on a common mu reference epoch (PS1 rebuilt
with --t0 59800 = ZTF's T0; see surveys/panstarrs/scripts/run_common_t0.sh).

Same estimator, flux-scale, band pairing, phase reference, controls and
candidate rules as joint_ps1_ztf.py (v1, mu = 0 only), extended to the
5 x 5 mu grid (|mu| <= 1"/yr, 0.5"/yr steps): thresholds from the control
maxima over the whole (z, mu) cube, injections with the WISE v0.2.0
grid + mu mismatch factors at duty 0.5, candidates from the cube maximum
with the phase, split-half, other-band and catalogued-static-source
tests at the fitted (z, mu).

Refuses to run unless every paired tensor reports the same t0_mjd.

Outputs runs/joint/ps1_ztf_v2/{threshold_report.json, m90_curves.npz,
records/...}.

Usage: uv run python surveys/joint/scripts/joint_ps1_ztf_mugrid.py
           [--ps1-tensors runs/panstarrs/calib_t0_59800/tensors]
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

from sglseti import canonical_json, load_target_registry, stable_id

from sglsurvey.geometry import GeometryContext
from sglsurvey.records import (AnalysisRun, Candidate, Constraint,
                               append_records)
from sglsurvey.vetting import (STATIC_RADIUS_ARCSEC, load_ps1_mean,
                               static_reason, static_source_test)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from joint_ps1_ztf import (ALBEDO_REF, BAND_PAIRS, CLIP_SIGMA,  # noqa: E402
                           MIN_EPOCHS, MIN_GOOD_FRAC, MIN_PHASE_EPOCHS,
                           N_REPEATS, N_Z_INTERVALS, RECOVERY_P, WEIGHT_CAP,
                           ZP_REF, ZTF_THROUGHPUT_MAG, interp_1overz,
                           phase_of, phase_ref_doy, size_limit_km, stack)

REPO = Path(__file__).resolve().parents[3]
ZTF_T = REPO / "runs" / "ztf" / "calib_v1" / "tensors"
OUT = REPO / "runs" / "joint" / "ps1_ztf_v2"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYP = {"ps1": REPO / "surveys" / "panstarrs" / "hypotheses.md",
       "ztf": REPO / "surveys" / "ztf" / "hypotheses.md"}
HYPOTHESIS_VERSION = ("joint-ps1-ztf-v2.0 (ps1-hypotheses-v1.0 + ztf-hypotheses-v1.0, "
                      "|mu| <= 1\"/yr, common T0 = 59800)")
SEED = 20260822
DUTY = 0.5


def mismatch_factor(fwhm_arcsec, halfcell_arcsec, mean_abs_dt_yr):
    sigma = fwhm_arcsec / 2.3548
    u = np.linspace(0, halfcell_arcsec, 64)
    grid_f = float(np.mean(np.exp(-u * u / (4 * sigma * sigma))))
    r_mu = 0.25 * mean_abs_dt_yr * np.sqrt(2)
    mu_f = float(np.exp(-r_mu * r_mu / (4 * sigma * sigma)))
    return grid_f * mu_f


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ps1-tensors",
                    default=str(REPO / "runs" / "panstarrs" / "calib_t0_59800" / "tensors"))
    args = ap.parse_args()
    ps1_t = Path(args.ps1_tensors)
    rng = np.random.default_rng(SEED)
    registry = load_target_registry(REGISTRY_PATH)
    OUT.mkdir(parents=True, exist_ok=True)
    hyp_hash = {k: "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()
                for k, p in HYP.items()}
    pairs = sorted({p.stem for p in ps1_t.glob("*.npz")}
                   & {p.stem for p in ZTF_T.glob("*.npz")})
    print(f"{len(pairs)} endpoint-roles in both archives ({ps1_t})", flush=True)
    ztf_scale = 10 ** (0.4 * ZTF_THROUGHPUT_MAG)
    sys.path.insert(0, str(REPO / "surveys" / "panstarrs" / "scripts"))
    from ps1_corridors import CORRIDOR_OF as _CORR
    ctx_ps1 = GeometryContext.ps1_default()
    screen_ps1 = REPO / "runs" / "panstarrs" / "screen_v1"
    catalog_cache = {}

    report, m90_store, meta, t0_common = {}, {}, {}, None
    for stem in pairs:
        e, role = stem.split("__")
        dp, dz = np.load(ps1_t / f"{stem}.npz"), np.load(ZTF_T / f"{stem}.npz")
        t0p, t0z = float(dp["t0_mjd"]), float(dz["t0_mjd"])
        if abs(t0p - t0z) > 1e-6:
            raise SystemExit(f"{stem}: mu reference epochs differ (PS1 {t0p}, ZTF {t0z})")
        t0_common = t0p
        zg, mug = dp["z_grid"], dp["mu_grid"]
        if not np.allclose(mug, dz["mu_grid"]):
            raise SystemExit(f"{stem}: mu grids differ")
        nm = len(mug)
        ref = phase_ref_doy(dp["mjd"])
        ph_p, ph_z = phase_of(dp["mjd"], ref), phase_of(dz["mjd"], ref)
        halfcell = (206265.0 / zg[0] - 206265.0 / zg[-1]) / (len(zg) - 1) / 2.0
        for band, (bp, bz) in BAND_PAIRS.items():
            ep, ez = dp["band_idx"] == bp, dz["band_idx"] == bz
            if ep.sum() + ez.sum() < MIN_EPOCHS:
                continue
            mjd = np.concatenate([dp["mjd"][ep], dz["mjd"][ez]])
            phase = np.concatenate([ph_p[ep], ph_z[ez]])
            arch = np.array([0] * int(ep.sum()) + [1] * int(ez.sum()))
            seeing = np.concatenate([dp["seeing"][ep],
                                     np.where(dz["seeing"][ez] > 0, dz["seeing"][ez], 2.0)])
            E = len(mjd)
            S = np.full((9, len(zg), nm, nm), np.nan, np.float32)
            nval = np.zeros((9, len(zg), nm, nm), int)
            W, FZ, A, B, G = ({} for _ in range(5))
            for mi in range(nm):
                for mj in range(nm):
                    fp = np.asarray(dp["f"][:, :, :, mi, mj][:, ep, :], np.float32)
                    vp = np.asarray(dp["v"][:, :, :, mi, mj][:, ep, :], np.float32)
                    gp = np.asarray(dp["g"][:, :, :, mi, mj][:, ep, :], np.float32)
                    if ez.sum():
                        fz_ = interp_1overz(np.asarray(dz["f"][:, :, :, mi, mj][:, ez, :],
                                                       np.float32) * ztf_scale, dz["z_grid"], zg)
                        vz_ = interp_1overz(np.asarray(dz["v"][:, :, :, mi, mj][:, ez, :],
                                                       np.float32) * ztf_scale ** 2, dz["z_grid"], zg)
                        gz_ = interp_1overz(np.asarray(dz["g"][:, :, :, mi, mj][:, ez, :],
                                                       np.float32), dz["z_grid"], zg)
                    else:
                        fz_ = vz_ = gz_ = np.empty((9, 0, len(zg)), np.float32)
                    f = np.concatenate([fp, fz_], axis=1)
                    v = np.concatenate([vp, vz_], axis=1)
                    g = np.concatenate([gp, gz_], axis=1)
                    S_, n_, w_, fz0, A_, B_ = stack(f, v, g, phase)
                    S[:, :, mi, mj], nval[:, :, mi, mj] = S_, n_
                    W[(mi, mj)], FZ[(mi, mj)], A[(mi, mj)], B[(mi, mj)], G[(mi, mj)] = (
                        w_[0], fz0[0], A_[0], B_[0], np.nan_to_num(g[0]))
            null_max = np.array([np.nanmax(S[t]) if np.isfinite(S[t]).any() else np.nan
                                 for t in range(1, 9)])
            if not np.isfinite(S[0]).any() or not np.isfinite(null_max).any():
                continue
            T = float(np.nanmax(null_max))
            pk = np.unravel_index(int(np.nanargmax(S[0])), S[0].shape)
            zi, mi, mj = (int(x) for x in pk)
            key = f"{e}/{role}/{band}"
            w0, f0 = W[(mi, mj)][:, zi], FZ[(mi, mj)][:, zi]
            ph_S, ph_n, ar_S, ar_n = {}, {}, {}, {}
            for p_ in (0, 1):
                sel = phase == p_
                ph_n[str(p_)] = int((w0[sel] > 0).sum())
                ph_S[str(p_)] = (float((f0[sel] * w0[sel]).sum() / np.sqrt(w0[sel].sum()))
                                 if w0[sel].sum() > 0 else -99.0)
            for a_, name in ((0, "ps1"), (1, "ztf")):
                sel = arch == a_
                ar_n[name] = int((w0[sel] > 0).sum())
                ar_S[name] = (float((f0[sel] * w0[sel]).sum() / np.sqrt(w0[sel].sum()))
                              if w0[sel].sum() > 0 else -99.0)
            order = np.argsort(mjd)
            hal = {}
            for name, idx in (("early", order[:E // 2]), ("late", order[E // 2:])):
                hal[name] = [(float((f0[idx] * w0[idx]).sum() / np.sqrt(w0[idx].sum()))
                              if w0[idx].sum() > 0 else -99.0), int((w0[idx] > 0).sum())]
            report[key] = {
                "T": T, "control_maxima": np.round(null_max, 2).tolist(),
                "real_max_S": float(S[0][pk]), "real_max_z": float(zg[zi]),
                "real_max_mu": [float(mug[mi]), float(mug[mj])],
                "exceeds": bool(S[0][pk] > T),
                "n_epochs": {"ps1": int(ep.sum()), "ztf": int(ez.sum())},
                "epoch_range_mjd": [float(mjd.min()), float(mjd.max())],
                "phase_S": ph_S, "phase_n": ph_n,
                "archive_S_at_peak": ar_S, "archive_n_at_peak": ar_n,
                "split_half": hal,
                "n_valid_at_peak": int(nval[0][pk]),
                "cells_defined_frac": float(np.isfinite(S[0]).mean()),
                "phase_split_all_epochs": [int((phase == 0).sum()), int((phase == 1).sum())],
            }
            # injections: random mu per repeat, duty 0.5, grid + mu mismatch
            mean_abs_dt = float(np.mean(np.abs(mjd - t0_common) / 365.25))
            mf_e = np.array([mismatch_factor(s_, halfcell, mean_abs_dt) for s_ in seeing])
            m90 = np.full(len(zg), np.nan)
            for zi_ in range(len(zg)):
                fmins = []
                for r_ in range(N_REPEATS):
                    ci = (int(rng.integers(0, nm)), int(rng.integers(0, nm)))
                    b0, a0 = B[ci][zi_], A[ci][zi_]
                    if b0 <= 0 or nval[0][zi_, ci[0], ci[1]] < MIN_EPOCHS:
                        continue
                    wg = W[ci][:, zi_] * G[ci][:, zi_] * mf_e
                    need = max(0.0, T * np.sqrt(b0) - a0)
                    keep = rng.random(E) < DUTY
                    den = float(wg[keep].sum())
                    if den > 0:
                        fmins.append(need / den)
                if len(fmins) >= int(0.8 * N_REPEATS):
                    f90 = float(np.percentile(fmins, 100 * RECOVERY_P))
                    if f90 > 0:
                        m90[zi_] = ZP_REF - 2.5 * np.log10(f90)
            m90_store[key] = m90
            meta[key] = {"z_grid": zg, "n_valid_z": nval[0].max(axis=(1, 2)),
                         "epoch_range": report[key]["epoch_range_mjd"],
                         "S_real_cube": S[0].copy(), "mjd": mjd, "phase": phase}
            print(f"  {key:24s} N={ep.sum():3d}+{ez.sum():3d} T={T:5.2f} S={S[0][pk]:5.2f} "
                  f"z={zg[zi]:6.0f} mu=({mug[mi]:+.1f},{mug[mj]:+.1f}) phase "
                  f"{ph_S['0']:.1f}/{ph_S['1']:.1f} (n {ph_n['0']}/{ph_n['1']}) "
                  f"m90 med={np.nanmedian(m90):5.2f}", flush=True)
        dp.close(); dz.close()

    config = {
        "pipeline_id": "joint-ps1-ztf-stack-mugrid",
        "hypothesis": f"relay on the focal line at T0={t0_common}, |mu| <= 1\"/yr (5x5), "
                      "550-10000 AU, duty >= 0.5",
        "z_grid": "PS1 360-node 1/z grid; ZTF interpolated in 1/z",
        "mu_reference_epoch_mjd": t0_common,
        "flux_scale": ("PS1 star-calibrated AB at ZP 25; ZTF MAGZP scale x 10^(0.4*0.5) "
                       "(asteroid-control throughput)"),
        "band_pairs": {k: f"ps1 {k} + ztf z{k}" for k in BAND_PAIRS},
        "threshold_rule": "max of 8 shared offset-control (z, mu)-cube maxima",
        "weight_cap": WEIGHT_CAP, "min_epochs": MIN_EPOCHS, "clip_sigma": CLIP_SIGMA,
        "injection_model": "random mu cell per repeat; grid + mu mismatch factors; duty 0.5",
        "n_repeats": N_REPEATS, "recovery_probability": RECOVERY_P, "seed": SEED,
        "candidate_rule": ("cube-maximum exceedance vetoed if either phase has < 3 epochs, "
                           "a phase S < 2, not persistent across early/late halves, absent "
                           "in the other paired bands at the same (z, mu), or a catalogued "
                           "DR2 source (>= 3 detections) within 2\" of the track at the "
                           "major-phase epochs (sglsurvey.vetting)"),
    }
    try:
        commit = subprocess.run(["git", "-C", str(REPO.parent / "sglseti"), "rev-parse",
                                 "--short", "HEAD"], capture_output=True, text=True,
                                timeout=10).stdout.strip()
    except Exception:
        commit = "unknown"
    run = AnalysisRun(
        analysis_run_id=stable_id("run", {"config": json.loads(canonical_json(config)),
                                          "registry": registry.source_hash,
                                          "hypothesis": HYPOTHESIS_VERSION, "pairs": pairs}),
        pipeline_id="joint-ps1-ztf-stack-mugrid", pipeline_version="0.1.0", config=config,
        observation_set_hash="see-ps1-and-ztf-tensor-manifests",
        intersection_set_hash="see-ps1-and-ztf-tensor-manifests",
        registry_source_hash=registry.source_hash, hypothesis_version=HYPOTHESIS_VERSION,
        environment={"sglseti_commit": commit, "hypothesis_hashes": hyp_hash,
                     "ps1_tensors": str(ps1_t.relative_to(REPO))},
        random_seeds={"numpy": SEED},
        started_utc=datetime.now(timezone.utc).isoformat(),
        finished_utc=datetime.now(timezone.utc).isoformat())

    candidates = []
    for key, r in report.items():
        if not r["exceeds"]:
            continue
        e, role, band = key.split("/")
        pn, pS, hal = r["phase_n"], r["phase_S"], r["split_half"]
        zg = meta[key]["z_grid"]
        zi = int(np.argmin(np.abs(zg - r["real_max_z"])))
        mi = int(np.argmin(np.abs(np.array([-1, -.5, 0, .5, 1]) - r["real_max_mu"][0])))
        mj = int(np.argmin(np.abs(np.array([-1, -.5, 0, .5, 1]) - r["real_max_mu"][1])))
        others = {ob: float(np.nan_to_num(meta[f"{e}/{role}/{ob}"]["S_real_cube"][zi, mi, mj],
                                          nan=-99.0))
                  for ob in BAND_PAIRS if ob != band and f"{e}/{role}/{ob}" in meta}
        r["other_bands_at_cell"] = others
        static = None
        if min(pn.values()) < MIN_PHASE_EPOCHS:
            status, reason = "vetoed", f"single-phase data (phase epochs {pn})"
        elif min(pS.values()) < 2.0:
            status, reason = "vetoed", (f"phase-split significances disagree {pS}: static "
                                        "source / artifact at one phase position")
        elif hal["early"][0] < 2.0 or hal["late"][0] < 2.0:
            status, reason = "vetoed", (f"split-half not persistent (early {hal['early']}, "
                                        f"late {hal['late']})")
        elif others and max(others.values()) < 2.0:
            status, reason = "vetoed", (f"absent in the other paired bands at the same cell "
                                        f"{ {k: round(v, 2) for k, v in others.items()} }")
        else:
            corr = _CORR[e]
            if corr not in catalog_cache:
                catalog_cache[corr] = load_ps1_mean(screen_ps1, corr, band)
            static = static_source_test(ctx_ps1, registry[e], role, r["real_max_z"],
                                        r["real_max_mu"], t0_common, meta[key]["mjd"],
                                        meta[key]["phase"], catalog_cache[corr],
                                        STATIC_RADIUS_ARCSEC["ps1"])
            r["static_source_test"] = static
            if static["static"]:
                status, reason = "vetoed", static_reason(static)
            else:
                status, reason = "retained", None
        candidates.append(Candidate(
            candidate_id=stable_id("cnd", {"endpoint": e, "role": role, "band": band,
                                           "source": "joint-ps1-ztf-v2-exceedance",
                                           "hypothesis": HYPOTHESIS_VERSION}),
            analysis_run_id=run.analysis_run_id, endpoint_id=e, role=role,
            observation_ids=(), fitted_z_au=r["real_max_z"],
            fitted_residual_motion={"mu_arcsec_yr": r["real_max_mu"], "t0_mjd": t0_common},
            model_comparison={"real_max_S": r["real_max_S"], "threshold_8_controls": r["T"],
                              "phase_S": pS, "phase_n": pn, "split_half": hal,
                              "archive_S": r["archive_S_at_peak"],
                              "archive_n": r["archive_n_at_peak"],
                              "other_bands_at_cell": others, "static_source_test": static},
            status=status, veto_reason=reason,
            extra={"band": band, "origin": "joint (z, mu)-cube exceedance census"}))

    constraints = []
    for key, m in m90_store.items():
        e, role, band = key.split("/")
        zg = meta[key]["z_grid"]
        nz = len(zg)
        edges = np.linspace(0, nz, N_Z_INTERVALS + 1, dtype=int)
        for k in range(N_Z_INTERVALS):
            lo, hi = edges[k], edges[k + 1]
            seg = m[lo:hi]
            zint = (float(min(zg[lo], zg[hi - 1])), float(max(zg[lo], zg[hi - 1])))
            nvalid = int(np.isfinite(seg).sum())
            if nvalid == 0:
                kind, limit = "not_constrainable", None
            else:
                kind = "recovery_curve"
                mlim = round(float(np.nanmin(seg)), 2)
                zmid = float(np.sqrt(zint[0] * zint[1]))
                limit = {"value": mlim, "unit": "ab_mag_ps1_scale",
                         "band": f"ps1 {band} + ztf z{band}",
                         "epochs_per_cell_max": int(meta[key]["n_valid_z"][lo:hi].max()),
                         "both_parallax_phases": bool(
                             min(report[key]["phase_split_all_epochs"]) >= MIN_PHASE_EPOCHS),
                         "interpretation": ("reflected sunlight / self-luminous optical; "
                                            "never thermal; achromatic within paired band"),
                         "reflected_light_size_km_albedo_0p1": round(
                             size_limit_km(mlim, zmid, ALBEDO_REF), 1),
                         "size_conversion_z_au": round(zmid, 1),
                         "nodes_valid": nvalid, "nodes_total": hi - lo}
            constraints.append(Constraint(
                constraint_id=stable_id("con", {"endpoint": e, "role": role, "band": band,
                                                "z_interval": [round(zint[0], 1), round(zint[1], 1)],
                                                "analysis_run": run.analysis_run_id,
                                                "hypothesis": HYPOTHESIS_VERSION}),
                analysis_run_id=run.analysis_run_id, endpoint_id=e, role=role,
                hypothesis_version=HYPOTHESIS_VERSION,
                z_interval_au=(round(zint[0], 1), round(zint[1], 1)),
                band=f"ps1 {band} + ztf z{band}",
                epoch_range_mjd=tuple(meta[key]["epoch_range"]),
                duty_cycle_range=(0.5, 1.0), residual_motion_bound_arcsec_per_yr=1.0,
                morphology="gaussian-psf-point-source", kind=kind,
                recovery_probability=RECOVERY_P if kind == "recovery_curve" else None,
                flux_limit=limit, false_alarm_rate=None,
                extra={"threshold": report[key]["T"], "n_epochs": report[key]["n_epochs"],
                       "mu_reference_epoch_mjd": t0_common}))

    rec = OUT / "records"
    append_records(rec / "analysis_run.jsonl", [run])
    append_records(rec / "constraint.jsonl", constraints)
    append_records(rec / "candidate.jsonl", candidates)
    np.savez_compressed(OUT / "m90_curves.npz",
                        **{k.replace("/", "__"): v for k, v in m90_store.items()})
    slim = {k: {kk: vv for kk, vv in v.items()} for k, v in report.items()}
    (OUT / "threshold_report.json").write_text(json.dumps(slim, indent=2, default=float))
    print(f"\nanalysis_run: {run.analysis_run_id}")
    print(f"{len(constraints)} constraints, {len(candidates)} candidates "
          f"({sum(c.status == 'retained' for c in candidates)} retained)")
    for c in candidates:
        print(f"  CANDIDATE {c.endpoint_id}/{c.role} {c.extra['band']}: {c.status}"
              + (f" ({c.veto_reason})" if c.veto_reason else ""))
    both = sum(1 for r in report.values() if min(r["phase_split_all_epochs"]) >= MIN_PHASE_EPOCHS)
    print(f"cells with both phases populated: {both}/{len(report)}")
    for band in BAND_PAIRS:
        ms = [np.nanmedian(m) for k, m in m90_store.items() if k.endswith("/" + band)
              and np.isfinite(m).any()]
        if ms:
            print(f"  {band}: joint m90 median {np.median(ms):.2f} over {len(ms)} endpoint-roles")


if __name__ == "__main__":
    main()
