"""Injection calibration and constraint generation for the SPHEREx
pilot. Same estimator as the ZTF calibration (8-offset-control
thresholds, effective-epoch weight cap, per-cell epoch floor,
single-epoch clip, duty-cycle-randomised analytic injections exploiting
matched-filter linearity, 90% recovery), with SPHEREx specifics:
per-epoch PSF FWHM sets the mismatch factors; depths quoted in AB mag
(uJy scale, ZP 23.9) per detector with the covered wavelength range;
every threshold exceedance becomes a Candidate adjudicated by the
parallax-phase split test; each constraint carries a reflected-light
size conversion (solar-coloured reflector, albedo 0.1) and blackbody-
equivalent emitting diameters at 400 / 700 / 1000 K; the star-control
throughput offset (runs/spherex/control_v1/summary.json) is applied in
the conservative direction when present.

Usage: uv run python surveys/spherex/scripts/injection_calibrate.py
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from sglseti import canonical_json, load_target_registry, stable_id

from sglsurvey.records import (AnalysisRun, Candidate, Constraint,
                               append_records)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spherex_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
CAL_DIR = REPO / "runs" / "spherex" / os.environ.get("SPHEREX_CALIB_RUN", "calib_v1")
TENSOR_DIR = CAL_DIR / "tensors"
CONTROL_SUMMARY = REPO / "runs" / "spherex" / "control_v1" / "summary.json"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "spherex" / "hypotheses.md"
HYPOTHESIS_VERSION = "spherex-hypotheses-v1.0"

MIN_GOOD_FRAC = 0.7
WEIGHT_CAP_FACTOR = 20.0
N_REPEATS = 32
DUTIES = [1.0, 0.5]
SEED = 20260820
RECOVERY_P = 0.90
N_Z_INTERVALS = 8
MIN_PHASE_EPOCHS = 3
MIN_EPOCHS = 5
CLIP_SIGMA = 5.0
ROLE_COINCIDENCE_SIGMA = 2.0
#: Amendment 2 (2026-08-21): the stack layer's cell is "static sky
#: removed"; where the template covered fewer than this fraction of the
#: real-trajectory epochs the exceedance is not adjudicable here.
MIN_TEMPLATE_FRAC = 0.8
#: Amendment 2: a layer-1 static cluster (screen_v1) within this many
#: PSF FWHM of the peak position and brighter than the stack bright
#: edge by BRIGHT_NEIGHBOUR_DMAG leaves wing residuals after template
#: subtraction (PSF varies between epochs and zones).
BRIGHT_NEIGHBOUR_FWHM = 3.0
BRIGHT_NEIGHBOUR_DMAG = 5.0
PSF_FWHM_NOMINAL_ARCSEC = 5.3
SCREEN_SUMMARY = REPO / "runs" / "spherex" / "screen_v1" / "summary.json"
ALBEDO_REF = 0.1
THERMAL_T_K = (400.0, 700.0, 1000.0)
BAND_RANGE_UM = {"D1": (0.75, 1.12), "D2": (1.10, 1.64), "D3": (1.62, 2.42),
                 "D4": (2.40, 3.82), "D5": (3.80, 4.42), "D6": (4.40, 5.00)}
BAND_PHYSICS = {
    "D1": "reflected sunlight (primary); self-luminous emission",
    "D2": "reflected sunlight (primary); self-luminous emission",
    "D3": "reflected sunlight (primary); hot components >~ 1000 K",
    "D4": "reflected sunlight / hot components (>~ 700 K); W1 regime",
    "D5": "hot components (>~ 500 K) / reflected sunlight; W1-W2 regime",
    "D6": "hot components (>~ 400 K) / reflected sunlight; W2 regime",
    "ALL": ("joint six-detector stack for a flat F_nu spectrum (0.75-5 um); "
            "a stellar-SED source is penalised, a grey reflector or "
            "non-thermal emitter gains ~sqrt(6)"),
}
JOINT_BAND = "ALL"
#: Solar absolute AB magnitude vs wavelength (um), Willmer 2018 values
#: at V, i, z, J, H, Ks, W1, W2; linearly interpolated.
SUN_AB = (np.array([0.55, 0.75, 0.90, 1.25, 1.65, 2.20, 3.40, 4.60]),
          np.array([4.81, 4.53, 4.50, 4.54, 4.71, 5.14, 5.92, 6.58]))
H_PLANCK, K_B, C_LIGHT, AU_M = 6.62607e-34, 1.380649e-23, 2.99792458e8, 1.495978707e11


def mismatch_factor(fwhm_arcsec: float, halfcell_arcsec: float,
                    mean_abs_dt_yr: float) -> float:
    sigma = fwhm_arcsec / 2.3548
    u = np.linspace(0, halfcell_arcsec, 64)
    grid_f = float(np.mean(np.exp(-u * u / (4 * sigma * sigma))))
    r_mu = 0.5 * mean_abs_dt_yr * np.sqrt(2)  # 1"/yr node spacing
    mu_f = float(np.exp(-r_mu * r_mu / (4 * sigma * sigma)))
    return grid_f * mu_f


def reflected_size_km(mag_ab: float, wave_um: float, z_au: float,
                      albedo: float) -> float:
    m_sun = float(np.interp(wave_um, *SUN_AB))
    h_v = (mag_ab - (m_sun - 4.81)) - 5.0 * np.log10(z_au * z_au)
    return float(1329.0 / np.sqrt(albedo) * 10 ** (-h_v / 5.0))


def thermal_size_km(flux_ujy: float, wave_um: float, z_au: float,
                    t_k: float) -> float:
    nu = C_LIGHT / (wave_um * 1e-6)
    b_nu = (2 * H_PLANCK * nu ** 3 / C_LIGHT ** 2
            / np.expm1(H_PLANCK * nu / (K_B * t_k)))
    f_si = flux_ujy * 1e-32
    r_m = z_au * AU_M * np.sqrt(f_si / (np.pi * b_nu))
    return float(2 * r_m / 1e3)


def main() -> None:
    rng = np.random.default_rng(SEED)
    registry = load_target_registry(REGISTRY_PATH)
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()
    control_mag = 0.0
    control_note = "no star-control summary found; no throughput correction"
    if CONTROL_SUMMARY.exists():
        cs = json.loads(CONTROL_SUMMARY.read_text())
        control_mag = max(0.0, float(cs.get("throughput_offset_mag", 0.0)))
        control_note = cs.get("note", "")
    config = {
        "pipeline_id": "spherex-injection-calibration",
        "threshold_rule": "max of the offset-control grid maxima (16 controls since v3)",
        "weight_cap": f"per-cell epoch weights capped at {WEIGHT_CAP_FACTOR}x median positive weight",
        "n_repeats": N_REPEATS, "duties": DUTIES, "seed": SEED,
        "recovery_probability": RECOVERY_P,
        "min_good_frac": MIN_GOOD_FRAC,
        "z_grid": "96 nodes uniform in 1/z, 550-10000 AU",
        "mu_grid": "3 nodes, -1/0/+1 arcsec/yr (unresolved over 14 months)",
        "injection_model": ("exposure-PSF point source at grid nodes; "
                            "grid/mu mismatch factors from per-epoch PSF FWHM"),
        "min_epochs_per_cell": MIN_EPOCHS,
        "single_epoch_clip_sigma": CLIP_SIGMA,
        "stack_cell_definition": ("sources below the single-epoch clip; "
                                  "brighter sources are the screening layer's cell"),
        "candidate_rule": ("every real-trajectory exceedance of T is a "
                           "Candidate; vetoed if either parallax phase "
                           f"has < {MIN_PHASE_EPOCHS} epochs, or the "
                           "phase-split significances disagree, or the "
                           "other role's stack at the same cell (Rx/Tx loci "
                           "coincide to < 1 arcsec for linear endpoints) is "
                           f"below {ROLE_COINCIDENCE_SIGMA} sigma"),
        "flux_scale": "matched-filter MJy/sr x OMEGA_MEDIAN -> uJy; AB = 23.9 - 2.5 log10 F",
        "size_conversion": (f"reflected: solar-coloured H-D convention, albedo {ALBEDO_REF}; "
                            f"thermal: blackbody sphere at {THERMAL_T_K} K"),
        "control_throughput_mag": control_mag,
        "control_note": control_note,
        "variance_calibration": ("v2: per-cutout MAD variance x template per-node "
                                 "reduced chi2 (sample_tensor.py, 2026-08-21)"),
        "joint_stack": ("v2: band 'ALL' = sum of the six detector stacks at each "
                        "cell (flat F_nu injection; same controls, phase and "
                        "role-coincidence vetoes)"),
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
    threshold_report, m90_store, meta_store, s_cubes = {}, {}, {}, {}
    joint_parts = {}
    tensor_files = sorted(TENSOR_DIR.glob("*.npz"))

    for path in tensor_files:
        endpoint, role, band = path.stem.split("__")
        d = np.load(path)
        z_grid, mu_grid = d["z_grid"], d["mu_grid"]
        zp_ref = float(d["zp_ref"])
        mjd, phase = d["mjd"], d["phase"]
        fwhm = np.asarray(d["psf_fwhm"], dtype=float)
        wave = np.asarray(d["wave_um"], dtype=float)
        f, v, g = d["f"], np.asarray(d["v"]), np.asarray(d["g"], dtype=np.float32)
        tsub = np.asarray(d["tsub"], dtype=np.float32) if "tsub" in d else np.zeros_like(g)
        nt = f.shape[0]
        nz, nm = len(z_grid), len(mu_grid)
        halfcell = (206265.0 / z_grid[0] - 206265.0 / z_grid[-1]) / (nz - 1) / 2.0
        mean_abs_dt = float(np.mean(np.abs(mjd - d["t0_mjd"]) / 365.25))

        valid = (np.isfinite(f) & np.isfinite(v) & (v > 0) & (g >= MIN_GOOD_FRAC))
        with np.errstate(invalid="ignore", divide="ignore"):
            s_e = np.abs(f) / np.sqrt(np.where(v > 0, v, np.inf))
        n_clipped = int((valid & (s_e > CLIP_SIGMA))[0].any(axis=(1, 2, 3)).sum())
        valid &= s_e <= CLIP_SIGMA
        w = np.where(valid, 1.0 / np.where(v > 0, v, 1.0), 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            sig_e = np.sqrt(np.where(v[0] > 0, v[0], np.nan))
        bright_edge = float(zp_ref - 2.5 * np.log10(
            CLIP_SIGMA * np.nanmedian(sig_e)) - control_mag)
        with np.errstate(all="ignore"):
            wmed = np.nanmedian(np.where(w > 0, w, np.nan), axis=1)
        cap = WEIGHT_CAP_FACTOR * np.nan_to_num(wmed, nan=np.inf, posinf=np.inf)
        w = np.minimum(w, cap[:, None])
        fz = np.where(valid, f, 0.0)
        A = (fz * w).sum(axis=1)
        B = w.sum(axis=1)
        n_valid = (w > 0).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            S = A / np.sqrt(np.where(B > 0, B, np.inf))
        S[n_valid < MIN_EPOCHS] = np.nan
        null_max = np.array([np.nanmax(S[t]) if np.isfinite(S[t]).any()
                             else np.nan for t in range(1, nt)])
        if not np.isfinite(S[0]).any() or not np.isfinite(null_max).any():
            d.close()
            continue
        T = float(np.nanmax(null_max))
        i_real = np.unravel_index(np.nanargmax(S[0]), S[0].shape)
        key = f"{endpoint}/{role}/{band}"
        s_cubes[key] = S[0]
        w2sum = (w[0] ** 2).sum(axis=0)
        n_eff = float(B[0][i_real] ** 2 / max(w2sum[i_real], 1e-300))
        ph_S, ph_n = {}, {}
        for p_ in (0, 1):
            sel = phase == p_
            wp = w[0][sel][:, i_real[0], i_real[1], i_real[2]]
            fp_ = fz[0][sel][:, i_real[0], i_real[1], i_real[2]]
            Bp = wp.sum()
            ph_n[str(p_)] = int((wp > 0).sum())
            ph_S[str(p_)] = (float((fp_ * wp).sum() / np.sqrt(Bp)) if Bp > 0 else -99.0)
        # wavelength-resolved significance at the real peak (candidate
        # discrimination aid: a relay has no stellar SED)
        wbins = np.linspace(BAND_RANGE_UM[band][0], BAND_RANGE_UM[band][1], 5)
        spec_S = []
        for k in range(4):
            sel = (wave >= wbins[k]) & (wave < wbins[k + 1] + (1e-6 if k == 3 else 0))
            wp = w[0][sel][:, i_real[0], i_real[1], i_real[2]]
            fp_ = fz[0][sel][:, i_real[0], i_real[1], i_real[2]]
            Bp = wp.sum()
            spec_S.append(round(float((fp_ * wp).sum() / np.sqrt(Bp)), 2) if Bp > 0 else None)
        threshold_report[key] = {
            "T": T, "control_maxima": np.round(null_max, 2).tolist(),
            "n_controls": int(nt - 1),
            "real_max_S": float(S[0][i_real]),
            "real_max_z": float(z_grid[i_real[0]]),
            "real_max_mu": [float(mu_grid[i_real[1]]), float(mu_grid[i_real[2]])],
            "real_max_cell": [int(i) for i in i_real],
            "exceeds": bool(S[0][i_real] > T),
            "n_epochs": int(len(mjd)),
            "n_eff_at_peak": round(n_eff, 1),
            "phase_S": ph_S, "phase_n": ph_n,
            "spectral_S_quartiles": spec_S,
            "median_psf_fwhm": float(np.nanmedian(fwhm)),
            "wave_range_um": [float(np.nanmin(wave)), float(np.nanmax(wave))],
            "median_var_scale": float(np.median(d["var_scale"])),
            "deep_frac": float(np.mean(d["deep"])),
            "template_subtracted_frac": float(np.nanmean(tsub[0])),
            "template_frac_at_peak": float(np.nanmean(tsub[0][:, i_real[0], i_real[1], i_real[2]])),
            "n_valid_at_peak": int(n_valid[0][i_real]),
            "cells_defined_frac": float(np.isfinite(S[0]).mean()),
            "epochs_with_any_clipped": n_clipped,
            "stack_bright_edge_mag": bright_edge,
        }
        # --- injections ----------------------------------------------
        E = len(mjd)
        fw = np.where(np.isfinite(fwhm) & (fwhm > 0), fwhm, 5.3)
        mf_e = np.array([mismatch_factor(s_, halfcell, mean_abs_dt) for s_ in fw])
        m90 = {duty: np.full(nz, np.nan) for duty in DUTIES}
        for zi in range(nz):
            mu_pick = rng.integers(0, nm, size=N_REPEATS)
            mu_pick2 = rng.integers(0, nm, size=N_REPEATS)
            f_min = np.full((len(DUTIES), N_REPEATS), np.nan)
            for r in range(N_REPEATS):
                ci = (zi, mu_pick[r], mu_pick2[r])
                wg = (w[0][:, ci[0], ci[1], ci[2]]
                      * np.nan_to_num(g[0][:, ci[0], ci[1], ci[2]]) * mf_e)
                a0, b0 = A[0][ci], B[0][ci]
                if b0 <= 0 or n_valid[0][ci] < MIN_EPOCHS:
                    continue
                need = max(0.0, T * np.sqrt(b0) - a0)
                for di, duty in enumerate(DUTIES):
                    keep = (rng.random(E) < duty) if duty < 1.0 else np.ones(E, dtype=bool)
                    denom = float(wg[keep].sum())
                    if denom > 0:
                        f_min[di, r] = need / denom
            for di, duty in enumerate(DUTIES):
                fm = f_min[di][np.isfinite(f_min[di])]
                if len(fm) >= int(0.8 * N_REPEATS):
                    f90 = float(np.percentile(fm, 100 * RECOVERY_P))
                    if f90 > 0:
                        m90[duty][zi] = zp_ref - 2.5 * np.log10(f90) - control_mag
        joint_parts.setdefault((endpoint, role), []).append({
            "band": band, "A": A.astype(np.float32), "B": B.astype(np.float32),
            "n_valid": n_valid, "w0": w[0].astype(np.float32),
            "fz0": fz[0].astype(np.float32),
            "wg0": (w[0] * np.nan_to_num(g[0]) * mf_e[:, None, None, None]).astype(np.float32),
            "phase": phase, "mjd": mjd, "wave": wave, "S0": S[0],
            "tsub0": tsub[0],
            "z_grid": z_grid, "mu_grid": mu_grid, "zp_ref": zp_ref})
        m90_store[key] = {str(duty): m90[duty] for duty in DUTIES}
        meta_store[key] = {"z_grid": z_grid,
                           "epoch_range": (float(mjd.min()), float(mjd.max())),
                           "n_valid_z": n_valid[0].max(axis=(1, 2)),
                           "bright_edge": bright_edge,
                           "wave_eff": float(np.nanmedian(wave))}
        d.close()

    # --- joint six-detector stack (v2, 2026-08-21) -------------------
    for (endpoint, role), parts in sorted(joint_parts.items()):
        if len(parts) < 2:
            continue
        z_grid, mu_grid, zp_ref = parts[0]["z_grid"], parts[0]["mu_grid"], parts[0]["zp_ref"]
        nz, nm = len(z_grid), len(mu_grid)
        A = sum(p_["A"].astype(np.float64) for p_ in parts)
        B = sum(p_["B"].astype(np.float64) for p_ in parts)
        n_valid = sum(p_["n_valid"] for p_ in parts)
        nt = A.shape[0]
        with np.errstate(divide="ignore", invalid="ignore"):
            S = A / np.sqrt(np.where(B > 0, B, np.inf))
        S[n_valid < MIN_EPOCHS] = np.nan
        null_max = np.array([np.nanmax(S[t]) if np.isfinite(S[t]).any() else np.nan
                             for t in range(1, nt)])
        if not np.isfinite(S[0]).any() or not np.isfinite(null_max).any():
            continue
        T = float(np.nanmax(null_max))
        i_real = np.unravel_index(np.nanargmax(S[0]), S[0].shape)
        key = f"{endpoint}/{role}/{JOINT_BAND}"
        s_cubes[key] = S[0]
        w0 = np.concatenate([p_["w0"] for p_ in parts], axis=0)
        fz0 = np.concatenate([p_["fz0"] for p_ in parts], axis=0)
        wg0 = np.concatenate([p_["wg0"] for p_ in parts], axis=0)
        phase = np.concatenate([p_["phase"] for p_ in parts])
        mjd = np.concatenate([p_["mjd"] for p_ in parts])
        wave = np.concatenate([p_["wave"] for p_ in parts])
        E = len(mjd)
        w2sum = (w0[:, i_real[0], i_real[1], i_real[2]].astype(np.float64) ** 2).sum()
        n_eff = float(B[0][i_real] ** 2 / max(w2sum, 1e-300))
        ph_S, ph_n = {}, {}
        for p_ in (0, 1):
            sel = phase == p_
            wp = w0[sel][:, i_real[0], i_real[1], i_real[2]].astype(np.float64)
            fp_ = fz0[sel][:, i_real[0], i_real[1], i_real[2]].astype(np.float64)
            Bp = wp.sum()
            ph_n[str(p_)] = int((wp > 0).sum())
            ph_S[str(p_)] = float((fp_ * wp).sum() / np.sqrt(Bp)) if Bp > 0 else -99.0
        per_band_S = {p_["band"]: (round(float(p_["S0"][i_real]), 2)
                                   if np.isfinite(p_["S0"][i_real]) else None) for p_ in parts}
        threshold_report[key] = {
            "T": T, "control_maxima": np.round(null_max, 2).tolist(),
            "n_controls": int(nt - 1),
            "real_max_S": float(S[0][i_real]),
            "real_max_z": float(z_grid[i_real[0]]),
            "real_max_mu": [float(mu_grid[i_real[1]]), float(mu_grid[i_real[2]])],
            "real_max_cell": [int(i) for i in i_real],
            "exceeds": bool(S[0][i_real] > T),
            "n_epochs": int(E), "n_eff_at_peak": round(n_eff, 1),
            "phase_S": ph_S, "phase_n": ph_n,
            "spectral_S_quartiles": per_band_S,
            "bands_combined": [p_["band"] for p_ in parts],
            "median_psf_fwhm": None,
            "wave_range_um": [float(np.nanmin(wave)), float(np.nanmax(wave))],
            "median_var_scale": None, "deep_frac": None,
            "template_subtracted_frac": None,
            "template_frac_at_peak": float(np.mean(np.concatenate(
                [p_["tsub0"][:, i_real[0], i_real[1], i_real[2]] for p_ in parts]))),
            "n_valid_at_peak": int(n_valid[0][i_real]),
            "cells_defined_frac": float(np.isfinite(S[0]).mean()),
            "epochs_with_any_clipped": None,
            "stack_bright_edge_mag": None,
        }
        m90 = {duty: np.full(nz, np.nan) for duty in DUTIES}
        for zi in range(nz):
            mu_pick = rng.integers(0, nm, size=N_REPEATS)
            mu_pick2 = rng.integers(0, nm, size=N_REPEATS)
            f_min = np.full((len(DUTIES), N_REPEATS), np.nan)
            for r_ in range(N_REPEATS):
                ci = (zi, mu_pick[r_], mu_pick2[r_])
                wg = wg0[:, ci[0], ci[1], ci[2]].astype(np.float64)
                a0, b0 = A[0][ci], B[0][ci]
                if b0 <= 0 or n_valid[0][ci] < MIN_EPOCHS:
                    continue
                need = max(0.0, T * np.sqrt(b0) - a0)
                for di, duty in enumerate(DUTIES):
                    keep = (rng.random(E) < duty) if duty < 1.0 else np.ones(E, dtype=bool)
                    denom = float(wg[keep].sum())
                    if denom > 0:
                        f_min[di, r_] = need / denom
            for di, duty in enumerate(DUTIES):
                fm = f_min[di][np.isfinite(f_min[di])]
                if len(fm) >= int(0.8 * N_REPEATS):
                    f90 = float(np.percentile(fm, 100 * RECOVERY_P))
                    if f90 > 0:
                        m90[duty][zi] = zp_ref - 2.5 * np.log10(f90) - control_mag
        m90_store[key] = {str(duty): m90[duty] for duty in DUTIES}
        meta_store[key] = {"z_grid": z_grid,
                           "epoch_range": (float(mjd.min()), float(mjd.max())),
                           "n_valid_z": n_valid[0].max(axis=(1, 2)),
                           "bright_edge": float("nan"),
                           "wave_eff": float(np.nanmedian(wave))}
    run = AnalysisRun(
        analysis_run_id=stable_id("run", {
            "config": json.loads(canonical_json(config)),
            "registry": registry.source_hash,
            "hypothesis": HYPOTHESIS_VERSION,
            "tensors": sorted(p.name for p in tensor_files)}),
        pipeline_id="spherex-injection-calibration",
        pipeline_version="0.1.0", config=config,
        observation_set_hash="see-cutout-index",
        intersection_set_hash="see-cutout-index",
        registry_source_hash=registry.source_hash,
        hypothesis_version=HYPOTHESIS_VERSION,
        environment={"sglseti_commit": commit, "hypothesis_hash": hyp_hash},
        random_seeds={"numpy": SEED},
        started_utc=started,
        finished_utc=datetime.now(timezone.utc).isoformat(),
    )

    screen = json.loads(SCREEN_SUMMARY.read_text()) if SCREEN_SUMMARY.exists() else {}

    def bright_neighbour(endpoint, role, z_au, mjd_mid, depth_mag):
        """Nearest layer-1 static cluster to the track position at z."""
        try:
            from astropy.time import Time
            from sglseti import Role, adaptive_locus
            from sglsurvey.geometry import GeometryContext
            ctx = GeometryContext.spherex_default()
            al = adaptive_locus(target=registry[endpoint], role=Role(role),
                                observation_time=Time(mjd_mid, format="mjd"),
                                observer=ctx.observer, relay_range=ctx.relay_range,
                                tolerance_arcsec=1.0, ephemeris=ctx.ephemeris, model=ctx.model)
            pt = min(al.points, key=lambda q: abs(q.z_au - z_au))
            ra0, dec0 = pt.icrs_ra_deg, pt.icrs_dec_deg
        except Exception:
            return None
        best = None
        for c in screen.get(CORRIDOR_OF.get(endpoint, ""), {}).get("static", []):
            d = 3600 * np.hypot((c["ra"] - ra0) * np.cos(np.deg2rad(dec0)), c["dec"] - dec0)
            if d <= BRIGHT_NEIGHBOUR_FWHM * PSF_FWHM_NOMINAL_ARCSEC and \
                    (depth_mag is None or c["median_mag_ab"] <= depth_mag - BRIGHT_NEIGHBOUR_DMAG):
                if best is None or d < best["sep_arcsec"]:
                    best = {"sep_arcsec": round(float(d), 1), "mag_ab": c["median_mag_ab"],
                            "n_epochs": c["n_epochs"], "catwise": c.get("catwise")}
        return best

    def season_complement(endpoint, role, cell, p_dom):
        """Amendment 3 (2026-08-21): stack ALL detectors' epochs of the
        parallax phase the exceeding stack lacks, at the same cell.
        In the deep fields the detectors observe in different seasons,
        so the per-detector split is degenerate while the cross-
        detector one is not."""
        A_ = B_ = 0.0
        n_ = 0
        for p_ in joint_parts.get((endpoint, role), []):
            sel = p_["phase"] != p_dom
            if not sel.any():
                continue
            wp = p_["w0"][sel][:, cell[0], cell[1], cell[2]].astype(np.float64)
            fp_ = p_["fz0"][sel][:, cell[0], cell[1], cell[2]].astype(np.float64)
            A_ += float((fp_ * wp).sum()); B_ += float(wp.sum()); n_ += int((wp > 0).sum())
        return (A_ / np.sqrt(B_) if B_ > 0 else None), n_

    for key, r in threshold_report.items():
        if not r["exceeds"]:
            continue
        endpoint, role, band = key.split("/")
        pn, pS = r["phase_n"], r["phase_S"]
        tfrac = r.get("template_frac_at_peak")
        cell = tuple(r["real_max_cell"])
        p_dom = max(pn, key=lambda k_: pn[k_])
        s_comp, n_comp = season_complement(endpoint, role, cell, int(p_dom))
        r["season_complement"] = {"S": s_comp, "n": n_comp, "lacking_phase": 1 - int(p_dom)}
        depth = (float(np.nanmedian(m90_store[key]["0.5"]))
                 if key in m90_store and np.isfinite(m90_store[key]["0.5"]).any() else None)
        er = meta_store[key]["epoch_range"]
        nb = bright_neighbour(endpoint, role, r["real_max_z"], 0.5 * (er[0] + er[1]), depth)
        r["bright_neighbour"] = nb
        other = f"{endpoint}/{'tx' if role == 'rx' else 'rx'}/{band}"
        s_other = None
        if other in s_cubes:
            c = tuple(r["real_max_cell"])
            s_other = float(s_cubes[other][c])
            if not np.isfinite(s_other):
                s_other = None
        r["other_role_S_at_peak"] = s_other
        if min(pn.values()) < MIN_PHASE_EPOCHS and n_comp < MIN_PHASE_EPOCHS:
            status, reason = "vetoed", (
                f"single-season data (phase epochs {pn}, cross-detector "
                f"complement {n_comp} epochs); a relay must appear at both "
                "parallax phases")
        elif min(pn.values()) < MIN_PHASE_EPOCHS and s_comp is not None and s_comp < 2.0:
            status, reason = "vetoed", (
                f"cross-detector season complement disagrees: S={s_comp:.2f} "
                f"over {n_comp} epochs of the lacking phase in all detectors; "
                "consistent with a static source or artifact at one phase position")
        elif min(pn.values()) >= MIN_PHASE_EPOCHS and min(pS.values()) < 2.0:
            status, reason = "vetoed", (
                f"phase-split significances disagree {pS}; consistent "
                "with a static source or artifact at one phase position")
        elif s_other is not None and s_other < ROLE_COINCIDENCE_SIGMA:
            status, reason = "vetoed", (
                f"other role's stack at the same cell S={s_other:.2f} < "
                f"{ROLE_COINCIDENCE_SIGMA}; Rx/Tx loci coincide, so a real "
                "source must appear in both")
        elif tfrac is not None and np.isfinite(tfrac) and tfrac < MIN_TEMPLATE_FRAC:
            status, reason = "vetoed", (
                f"static-sky template covered only {tfrac:.2f} of epochs at "
                f"the peak cell (< {MIN_TEMPLATE_FRAC}); not adjudicable by "
                "the stack layer (unremoved static confusion)")
        elif nb is not None:
            status, reason = "vetoed", (
                f"static source {nb['mag_ab']} AB at {nb['sep_arcsec']}\" "
                f"(< {BRIGHT_NEIGHBOUR_FWHM} FWHM) from the track position; "
                "PSF-wing residual after template subtraction")
        else:
            status, reason = "retained", None
        candidates.append(Candidate(
            candidate_id=stable_id("cnd", {
                "endpoint": endpoint, "role": role, "band": band,
                "source": "calib_v1-threshold-exceedance",
                "hypothesis": HYPOTHESIS_VERSION}),
            analysis_run_id=run.analysis_run_id,
            endpoint_id=endpoint, role=role, observation_ids=(),
            fitted_z_au=r["real_max_z"],
            model_comparison={"real_max_S": r["real_max_S"],
                              "threshold_controls": r["T"], "n_controls": r.get("n_controls"),
                              "control_maxima": r["control_maxima"],
                              "phase_S": pS, "phase_n": pn,
                              "spectral_S_quartiles": r["spectral_S_quartiles"],
                              "other_role_S_at_peak": s_other,
                              "template_frac_at_peak": tfrac,
                              "bright_neighbour": nb,
                              "season_complement": r["season_complement"],
                              "mu_arcsec_yr": r["real_max_mu"]},
            status=status, veto_reason=reason,
            extra={"band": band, "origin": "threshold exceedance census"},
        ))

    for key, curves in m90_store.items():
        endpoint, role, band = key.split("/")
        rep = threshold_report[key]
        z_grid = meta_store[key]["z_grid"]
        nz = len(z_grid)
        edges = np.linspace(0, nz, N_Z_INTERVALS + 1, dtype=int)
        m = curves["0.5"]
        wave_eff = meta_store[key]["wave_eff"]
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
                flim_ujy = 10 ** ((23.9 - mlim) / 2.5)
                zmid = float(np.sqrt(zint[0] * zint[1]))
                nvz = meta_store[key]["n_valid_z"][lo:hi]
                limit = {"value": mlim, "unit": "ab_mag", "flux_ujy": round(flim_ujy, 2),
                         "band": band,
                         "wave_range_um": [round(x, 3) for x in rep["wave_range_um"]],
                         "wave_eff_um": round(wave_eff, 3),
                         "epochs_per_cell_max": int(nvz.max()),
                         "stack_bright_edge_mag": (round(meta_store[key]["bright_edge"], 2)
                                                   if np.isfinite(meta_store[key]["bright_edge"]) else None),
                         "control_throughput_mag_applied": control_mag,
                         "template_subtracted_frac": (round(rep["template_subtracted_frac"], 2)
                                                      if rep["template_subtracted_frac"] is not None else None),
                         "interpretation": BAND_PHYSICS[band],
                         "reflected_light_size_km_albedo_0p1":
                             round(reflected_size_km(mlim, wave_eff, zmid, ALBEDO_REF), 1),
                         "thermal_size_km": {
                             f"{int(t)}K": round(thermal_size_km(flim_ujy, wave_eff, zmid, t), 1)
                             for t in THERMAL_T_K},
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
                morphology="exposure-psf-point-source",
                kind=kind,
                recovery_probability=(RECOVERY_P if kind == "recovery_curve" else None),
                flux_limit=limit, false_alarm_rate=None,
                extra={"threshold": rep["T"],
                       "threshold_rule": config["threshold_rule"],
                       "n_epochs": rep["n_epochs"],
                       "n_controls": rep.get("n_controls"),
                       "far_note": (f"empirical: 0 of {rep.get('n_controls')} control "
                                    "trajectories exceed T by construction; per-grid-"
                                    f"search FAR < 1/{rep.get('n_controls')}")},
            ))

    joint_parts.clear()
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
            print(f"  {key:26s} N={r['n_epochs']:4d} cells={r['cells_defined_frac']:.2f} "
                  f"T={r['T']:5.2f} S={r['real_max_S']:5.2f} "
                  f"m90 med={np.nanmedian(m):5.2f} best={np.nanmax(m):5.2f} "
                  f"worst={np.nanmin(m):5.2f} fwhm={r['median_psf_fwhm'] or 0:.1f}")


if __name__ == "__main__":
    main()
