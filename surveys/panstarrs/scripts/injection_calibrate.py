"""Stage-7 injection calibration and constraint generation for the
Pan-STARRS1 pilot (port of surveys/ztf/scripts/injection_calibrate.py).
Same estimator as WISE calibration v0.2.0 / ZTF v1.0 (8-offset-control
thresholds, effective-epoch weight cap, per-cell epoch floor,
single-epoch clip, duty-cycle-randomised analytic injections, 90%
recovery), with PS1 specifics: per-epoch header seeing sets the mismatch
factors; depths are AB mag on the per-warp STAR-CALIBRATED zero point
(which absorbs the Gaussian-PSF throughput, so no separate control
correction is applied; the asteroid control tests the trajectory
machinery end to end); every threshold exceedance becomes a Candidate
adjudicated by the parallax-phase split test; each constraint carries a
reflected-light size conversion (asteroid H-D convention).

Usage: uv run python surveys/panstarrs/scripts/injection_calibrate.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from sglseti import canonical_json, load_target_registry, stable_id

from sglsurvey.records import (AnalysisRun, Candidate, Constraint,
                               append_records)
from sglsurvey.vetting import (STATIC_RADIUS_ARCSEC, load_ps1_mean,
                               static_reason, static_source_test)

REPO = Path(__file__).resolve().parents[3]
CAL_DIR = REPO / "runs" / "panstarrs" / "calib_v1"
SCREEN_DIR = REPO / "runs" / "panstarrs" / "screen_v1"
TENSOR_DIR = CAL_DIR / "tensors"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "panstarrs" / "hypotheses.md"
HYPOTHESIS_VERSION = "ps1-hypotheses-v1.0"

MIN_GOOD_FRAC = 0.7
WEIGHT_CAP_FACTOR = 20.0
N_REPEATS = 32
DUTIES = [1.0, 0.5]
SEED = 20260820
RECOVERY_P = 0.90
N_Z_INTERVALS = 8
MIN_PHASE_EPOCHS = 3
#: Cells (z, mu) with fewer valid epochs than this are undefined for the
#: real trajectory AND the controls: with N=1-2 a single residual sets
#: S, and no weight cap or phase test can act (chip-gap geometry makes
#: such sparsely covered cells common in ZTF).
MIN_EPOCHS = 5
#: Layered search (plan §3.4): sources above this single-epoch
#: significance belong to the catalog-screening layer (psfcat + DR24
#: recurrence); the stack cell is defined as sub-threshold sources, and
#: epochs where |S_e| exceeds the clip are dropped from the stack for
#: real, control and injection trajectories alike. Static stars on a
#: science-image-regime track would otherwise set S ~ 100 per crossing.
CLIP_SIGMA = 5.0
#: No separate throughput correction: the per-warp zero point is
#: calibrated with catalogued stars through the identical matched filter
#: (sample_tensor.py), so estimator throughput is already in the scale.
CONTROL_THROUGHPUT_MAG = 0.0
BAND_NAME = {1: "g", 2: "r", 3: "i", 4: "z", 5: "y"}
BAND_PHYSICS = ("reflected sunlight (primary) / self-luminous optical "
                "emission; never thermal coverage")
ALBEDO_REF = 0.1


def mismatch_factor(fwhm_arcsec: float, halfcell_arcsec: float,
                    mean_abs_dt_yr: float) -> float:
    sigma = fwhm_arcsec / 2.3548
    u = np.linspace(0, halfcell_arcsec, 64)
    grid_f = float(np.mean(np.exp(-u * u / (4 * sigma * sigma))))
    r_mu = 0.25 * mean_abs_dt_yr * np.sqrt(2)
    mu_f = float(np.exp(-r_mu * r_mu / (4 * sigma * sigma)))
    return grid_f * mu_f


def size_limit_km(mag: float, z_au: float, albedo: float) -> float:
    H = mag - 5.0 * np.log10(z_au * z_au)
    return float(1329.0 / np.sqrt(albedo) * 10 ** (-H / 5.0))


def main() -> None:
    rng = np.random.default_rng(SEED)
    registry = load_target_registry(REGISTRY_PATH)
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()
    config = {
        "pipeline_id": "ps1-injection-calibration",
        "threshold_rule": "max of 8 offset-control grid maxima",
        "weight_cap": f"per-cell epoch weights capped at {WEIGHT_CAP_FACTOR}x median positive weight",
        "n_repeats": N_REPEATS, "duties": DUTIES, "seed": SEED,
        "recovery_probability": RECOVERY_P,
        "min_good_frac": MIN_GOOD_FRAC,
        "z_grid": "360 nodes uniform in 1/z, 550-10000 AU",
        "injection_model": ("gaussian-psf point source at grid nodes, "
                            "per-epoch seeing; grid/mu mismatch factors"),
        "min_epochs_per_cell": MIN_EPOCHS,
        "single_epoch_clip_sigma": CLIP_SIGMA,
        "stack_cell_definition": ("sources below the single-epoch clip; "
                                  "brighter sources are the screening "
                                  "layer's cell"),
        "candidate_rule": ("every real-trajectory exceedance of T is a "
                           "Candidate; vetoed if either parallax phase "
                           f"has < {MIN_PHASE_EPOCHS} epochs or the "
                           "phase-split significances disagree (one "
                           "phase < 2 sigma while combined > T), or a "
                           "catalogued DR2 source lies within "
                           f"{STATIC_RADIUS_ARCSEC['ps1']}\" of the track at "
                           "the major-phase epochs (sglsurvey.vetting)"),
        "size_conversion": f"asteroid H-D convention, albedo {ALBEDO_REF}",
        "control_throughput_mag": CONTROL_THROUGHPUT_MAG,
        "flux_scale": ("per-warp ZP from DR2 mean-table stars through the "
                       "same matched filter; header FPA.ZP fallback offset"),
        "control_run": "control_v1 (numbered asteroid, F51)",
    }
    try:
        commit = subprocess.run(
            ["git", "-C", str(REPO.parent / "sglseti"), "rev-parse",
             "--short", "HEAD"], capture_output=True, text=True,
            timeout=10).stdout.strip()
    except Exception:
        commit = "unknown"

    started = datetime.now(timezone.utc).isoformat()
    constraints, candidates = [], []
    threshold_report, m90_store, meta_store = {}, {}, {}
    epoch_store, catalog_cache = {}, {}
    from sglsurvey.geometry import GeometryContext
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from ps1_corridors import CORRIDOR_OF as _CORR
    ctx = GeometryContext.ps1_default()
    tensor_files = sorted(TENSOR_DIR.glob("*.npz"))

    for path in tensor_files:
        endpoint, role = path.stem.split("__")
        d = np.load(path)
        z_grid, mu_grid = d["z_grid"], d["mu_grid"]
        zp_ref = float(d["zp_ref"])
        mjd, band_idx, phase = d["mjd"], d["band_idx"], d["phase"]
        epoch_store[(endpoint, role)] = (mjd.copy(), band_idx.copy(), phase.copy(),
                                         float(d["t0_mjd"]))
        seeing = d["seeing"]
        f, v, g = d["f"], np.asarray(d["v"]), np.asarray(d["g"],
                                                        dtype=np.float32)
        nt = f.shape[0]
        nz, nm = len(z_grid), len(mu_grid)
        halfcell = (206265.0 / z_grid[0] - 206265.0 / z_grid[-1]) \
            / (nz - 1) / 2.0
        mean_abs_dt = float(np.mean(np.abs(mjd - d["t0_mjd"]) / 365.25))

        for b in sorted(set(band_idx.tolist())):
            eb = band_idx == b
            fb, vb, gb = f[:, eb], v[:, eb], g[:, eb]
            phb = phase[eb]
            valid = (np.isfinite(fb) & np.isfinite(vb) & (vb > 0)
                     & (gb >= MIN_GOOD_FRAC))
            with np.errstate(invalid="ignore", divide="ignore"):
                s_e = np.abs(fb) / np.sqrt(np.where(vb > 0, vb, np.inf))
            n_clipped = int((valid & (s_e > CLIP_SIGMA))[0].any(axis=(1, 2, 3)).sum())
            valid &= s_e <= CLIP_SIGMA
            w = np.where(valid, 1.0 / np.where(vb > 0, vb, 1.0), 0.0)
            # single-epoch 5-sigma brightness (median epoch): the bright
            # edge of the stack cell, quoted on every constraint
            with np.errstate(invalid="ignore", divide="ignore"):
                sig_e = np.sqrt(np.where(vb[0] > 0, vb[0], np.nan))
            bright_edge = float(zp_ref - 2.5 * np.log10(
                CLIP_SIGMA * np.nanmedian(sig_e)) - CONTROL_THROUGHPUT_MAG)
            with np.errstate(all="ignore"):
                wmed = np.nanmedian(np.where(w > 0, w, np.nan), axis=1)
            cap = WEIGHT_CAP_FACTOR * np.nan_to_num(wmed, nan=np.inf,
                                                    posinf=np.inf)
            w = np.minimum(w, cap[:, None])
            fz = np.where(valid, fb, 0.0)
            A = (fz * w).sum(axis=1)
            B = w.sum(axis=1)
            n_valid = (w > 0).sum(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                S = A / np.sqrt(np.where(B > 0, B, np.inf))
            S[n_valid < MIN_EPOCHS] = np.nan
            null_max = np.array([np.nanmax(S[t]) if np.isfinite(S[t]).any()
                                 else np.nan for t in range(1, nt)])
            if not np.isfinite(S[0]).any() or not np.isfinite(null_max).any():
                continue
            T = float(np.nanmax(null_max))
            i_real = np.unravel_index(np.nanargmax(S[0]), S[0].shape)
            band = BAND_NAME[b]
            key = f"{endpoint}/{role}/{band}"
            w2sum = (w[0] ** 2).sum(axis=0)
            n_eff = float(B[0][i_real] ** 2 / max(w2sum[i_real], 1e-300))
            # phase-split significances at the real peak
            ph_S, ph_n = {}, {}
            for p_ in (0, 1):
                sel = phb == p_
                wp = w[0][sel][:, i_real[0], i_real[1], i_real[2]]
                fp_ = fz[0][sel][:, i_real[0], i_real[1], i_real[2]]
                Bp = wp.sum()
                ph_n[str(p_)] = int((wp > 0).sum())
                ph_S[str(p_)] = (float((fp_ * wp).sum() / np.sqrt(Bp))
                                 if Bp > 0 else -99.0)
            threshold_report[key] = {
                "T": T, "control_maxima": np.round(null_max, 2).tolist(),
                "real_max_S": float(S[0][i_real]),
                "real_max_z": float(z_grid[i_real[0]]),
                "real_max_mu": [float(mu_grid[i_real[1]]),
                                float(mu_grid[i_real[2]])],
                "exceeds": bool(S[0][i_real] > T),
                "n_epochs": int(eb.sum()),
                "n_eff_at_peak": round(n_eff, 1),
                "phase_S": ph_S, "phase_n": ph_n,
                "median_seeing": float(np.median(seeing[eb])),
                "n_valid_at_peak": int(n_valid[0][i_real]),
                "cells_defined_frac": float(np.isfinite(S[0]).mean()),
                "epochs_with_any_clipped": n_clipped,
                "stack_bright_edge_mag": bright_edge,
                "n_epochs_star_calibrated": int((d["n_cal"][eb] >= 5).sum()),
                "median_zp_star_minus_hdr": float(np.nanmedian(
                    d["zp_star"][eb] - d["zp_hdr"][eb])),
            }
            # --- injections ----------------------------------------------
            E = int(eb.sum())
            see = np.where(seeing[eb] > 0, seeing[eb], 2.0)
            mf_e = np.array([mismatch_factor(s_, halfcell, mean_abs_dt)
                             for s_ in see])
            m90 = {duty: np.full(nz, np.nan) for duty in DUTIES}
            for zi in range(nz):
                mu_pick = rng.integers(0, nm, size=N_REPEATS)
                mu_pick2 = rng.integers(0, nm, size=N_REPEATS)
                f_min = np.full((len(DUTIES), N_REPEATS), np.nan)
                for r in range(N_REPEATS):
                    ci = (zi, mu_pick[r], mu_pick2[r])
                    wg = (w[0][:, ci[0], ci[1], ci[2]]
                          * np.nan_to_num(gb[0][:, ci[0], ci[1], ci[2]])
                          * mf_e)
                    a0, b0 = A[0][ci], B[0][ci]
                    if b0 <= 0 or n_valid[0][ci] < MIN_EPOCHS:
                        continue
                    need = max(0.0, T * np.sqrt(b0) - a0)
                    for di, duty in enumerate(DUTIES):
                        keep = (rng.random(E) < duty) if duty < 1.0 \
                            else np.ones(E, dtype=bool)
                        denom = float(wg[keep].sum())
                        if denom > 0:
                            f_min[di, r] = need / denom
                for di, duty in enumerate(DUTIES):
                    fm = f_min[di][np.isfinite(f_min[di])]
                    if len(fm) >= int(0.8 * N_REPEATS):
                        f90 = float(np.percentile(fm, 100 * RECOVERY_P))
                        if f90 > 0:
                            m90[duty][zi] = (zp_ref - 2.5 * np.log10(f90)
                                             - CONTROL_THROUGHPUT_MAG)
            m90_store[key] = {str(duty): m90[duty] for duty in DUTIES}
            meta_store[key] = {"z_grid": z_grid,
                               "epoch_range": (float(mjd[eb].min()),
                                               float(mjd[eb].max())),
                               "n_valid_z": n_valid[0].max(axis=(1, 2)),
                               "bright_edge": bright_edge}
        d.close()

    run = AnalysisRun(
        analysis_run_id=stable_id("run", {
            "config": json.loads(canonical_json(config)),
            "registry": registry.source_hash,
            "hypothesis": HYPOTHESIS_VERSION,
            "tensors": sorted(p.name for p in tensor_files)}),
        pipeline_id="ps1-injection-calibration",
        pipeline_version="0.2.0", config=config,
        observation_set_hash="see-tensor-manifest",
        intersection_set_hash="see-tensor-manifest",
        registry_source_hash=registry.source_hash,
        hypothesis_version=HYPOTHESIS_VERSION,
        environment={"sglseti_commit": commit, "hypothesis_hash": hyp_hash},
        random_seeds={"numpy": SEED},
        started_utc=started,
        finished_utc=datetime.now(timezone.utc).isoformat(),
    )

    # --- Candidates: full exceedance census --------------------------
    for key, r in threshold_report.items():
        if not r["exceeds"]:
            continue
        endpoint, role, band = key.split("/")
        pn, pS = r["phase_n"], r["phase_S"]
        if min(pn.values()) < MIN_PHASE_EPOCHS:
            status, reason = "vetoed", (
                f"single-phase data (phase epochs {pn}); a relay must "
                "appear at both parallax phases")
        elif min(pS.values()) < 2.0:
            status, reason = "vetoed", (
                f"phase-split significances disagree {pS}; consistent "
                "with a static source or artifact at one phase position")
        else:
            status, reason = "retained", None
        static = None
        if status == "retained":
            corr = _CORR[endpoint]
            if corr not in catalog_cache:
                catalog_cache[corr] = load_ps1_mean(SCREEN_DIR, corr, band if band in "gri" else "r")
            mjd_, bidx_, ph_, t0_ = epoch_store[(endpoint, role)]
            eb_ = bidx_ == {v: k for k, v in BAND_NAME.items()}[band]
            static = static_source_test(
                ctx, registry[endpoint], role, r["real_max_z"], r["real_max_mu"],
                t0_, mjd_[eb_], ph_[eb_], catalog_cache[corr],
                STATIC_RADIUS_ARCSEC["ps1"])
            if static["static"]:
                status, reason = "vetoed", static_reason(static)
            r["static_source_test"] = static
        candidates.append(Candidate(
            candidate_id=stable_id("cnd", {
                "endpoint": endpoint, "role": role, "band": band,
                "source": "calib_v1-threshold-exceedance",
                "hypothesis": HYPOTHESIS_VERSION}),
            analysis_run_id=run.analysis_run_id,
            endpoint_id=endpoint, role=role, observation_ids=(),
            fitted_z_au=r["real_max_z"],
            model_comparison={"real_max_S": r["real_max_S"],
                              "threshold_8_controls": r["T"],
                              "control_maxima": r["control_maxima"],
                              "phase_S": pS, "phase_n": pn,
                              "mu_arcsec_yr": r["real_max_mu"],
                              "static_source_test": static},
            status=status, veto_reason=reason,
            extra={"band": band, "origin": "threshold exceedance census"},
        ))

    # --- Constraints -------------------------------------------------
    for key, curves in m90_store.items():
        endpoint, role, band = key.split("/")
        rep = threshold_report[key]
        z_grid = meta_store[key]["z_grid"]
        nz = len(z_grid)
        edges = np.linspace(0, nz, N_Z_INTERVALS + 1, dtype=int)
        m = curves["0.5"]
        for k in range(N_Z_INTERVALS):
            lo, hi = edges[k], edges[k + 1]
            seg = m[lo:hi]
            zint = (float(min(z_grid[lo], z_grid[hi - 1])),
                    float(max(z_grid[lo], z_grid[hi - 1])))
            n_valid = int(np.isfinite(seg).sum())
            if n_valid == 0:
                kind, limit = "not_constrainable", None
            else:
                kind = "recovery_curve"
                mlim = round(float(np.nanmin(seg)), 2)
                zmid = float(np.sqrt(zint[0] * zint[1]))
                nvz = meta_store[key]["n_valid_z"][lo:hi]
                limit = {"value": mlim, "unit": "ps1_ab_mag", "band": band,
                         "epochs_per_cell_max": int(nvz.max()),
                         "stack_bright_edge_mag": round(
                             meta_store[key]["bright_edge"], 2),
                         "control_throughput_mag_applied": CONTROL_THROUGHPUT_MAG,
                         "flux_scale": "star-calibrated per-warp ZP",
                         "interpretation": BAND_PHYSICS,
                         "reflected_light_size_km_albedo_0p1":
                             round(size_limit_km(mlim, zmid, ALBEDO_REF), 1),
                         "size_conversion_z_au": round(zmid, 1),
                         "nodes_valid": n_valid, "nodes_total": hi - lo}
            constraints.append(Constraint(
                constraint_id=stable_id("con", {
                    "endpoint": endpoint, "role": role, "band": band,
                    "z_interval": [round(zint[0], 1), round(zint[1], 1)],
                    "analysis_run": run.analysis_run_id,
                    "hypothesis": HYPOTHESIS_VERSION}),
                analysis_run_id=run.analysis_run_id,
                endpoint_id=endpoint, role=role,
                hypothesis_version=HYPOTHESIS_VERSION,
                z_interval_au=(round(zint[0], 1), round(zint[1], 1)),
                band=band, epoch_range_mjd=meta_store[key]["epoch_range"],
                duty_cycle_range=(0.5, 1.0),
                residual_motion_bound_arcsec_per_yr=1.0,
                morphology="gaussian-psf-point-source",
                kind=kind,
                recovery_probability=(RECOVERY_P if kind == "recovery_curve"
                                      else None),
                flux_limit=limit, false_alarm_rate=None,
                extra={"threshold": rep["T"],
                       "threshold_rule": config["threshold_rule"],
                       "n_epochs": rep["n_epochs"],
                       "far_note": ("empirical: 0 of 8 control trajectories "
                                    "exceed T by construction; per-grid-"
                                    "search FAR < 1/8")},
            ))

    rec_dir = CAL_DIR / "records"
    append_records(rec_dir / "analysis_run.jsonl", [run])
    append_records(rec_dir / "constraint.jsonl", constraints)
    append_records(rec_dir / "candidate.jsonl", candidates)
    np.savez_compressed(
        CAL_DIR / "m90_curves.npz",
        **{k.replace("/", "__") + "__" + duty: arr
           for k, c in m90_store.items() for duty, arr in c.items()})
    (CAL_DIR / "threshold_report.json").write_text(
        json.dumps(threshold_report, indent=2))

    print(f"analysis_run: {run.analysis_run_id}")
    print(f"{len(constraints)} constraints, {len(candidates)} candidates")
    for c in candidates:
        print(f"  CANDIDATE {c.endpoint_id}/{c.role} {c.extra['band']}: "
              f"{c.status}" + (f" ({c.veto_reason})" if c.veto_reason else ""))
    print("\n=== recovery depths (duty>=0.5, 90% recovery, AB mag) ===")
    for key in sorted(m90_store):
        m = m90_store[key]["0.5"]
        r = threshold_report[key]
        if np.isfinite(m).any():
            print(f"  {key:22s} N={r['n_epochs']:3d} cells={r['cells_defined_frac']:.2f} "
                  f"T={r['T']:5.2f} S={r['real_max_S']:5.2f} "
                  f"m90 med={np.nanmedian(m):5.2f} best={np.nanmax(m):5.2f} "
                  f"worst={np.nanmin(m):5.2f} seeing={r['median_seeing']:.2f}\"")


if __name__ == "__main__":
    main()
