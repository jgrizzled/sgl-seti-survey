"""Stage 2 (plan §3.4 motion-model stage, cross-archive): joint PS1 + ZTF
stacks per endpoint x role x band on the corridors both archives cover.

Why: PS1 3pi samples every corridor at one parallax phase (97:3), ZTF
at both; a static background source at the PS1 track position can only
be excluded with the other phase. The joint stack also spans 2009-2026
(17 yr) at 1" resolution.

Hypothesis handled here: the station-kept relay, mu_resid = 0 — the
ONE trajectory that is identical in both tensor families (their mu
reference epochs differ: PS1 T0 = 56000, ZTF T0 = 59800, so a mu != 0
cell is a different track in each archive and is left to the
per-archive calibrations). The three PS1 marginal mu != 0 cells are
tested separately by direct forced photometry on the ZTF cutouts along
the PS1-extrapolated trajectory (``--marginal``).

Inputs: runs/panstarrs/calib_v1/tensors (360-node 1/z grid, star-
calibrated AB flux at ZP 25), runs/ztf/calib_v1/tensors (192-node 1/z
grid, MAGZP scale at ZP 25; the ZTF asteroid control measured a 0.5 mag
Gaussian-kernel throughput loss, applied here as x10^0.2 so both
archives are on the PS1 total-flux scale). ZTF is interpolated onto the
PS1 grid in 1/z. Bands paired g:zg, r:zr, i:zi (filter-curve
differences ignored; an achromatic-within-band source assumed — stated
on every record). Both archives use the same 9 trajectories (real + 8
offset controls), so the control thresholds are joint by construction.

Phase labels are recomputed for both archives against one reference
(the circular-median PS1 day-of-year of the corridor), so the phase
split is meaningful across archives.

Outputs runs/joint/ps1_ztf_v1/{threshold_report.json, m90_curves.npz,
records/{analysis_run,constraint,candidate}.jsonl, marginal_tests.json}.

Usage: uv run python surveys/joint/scripts/joint_ps1_ztf.py [--marginal]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from sglseti import canonical_json, load_target_registry, stable_id

from sglsurvey.records import (AnalysisRun, Candidate, Constraint,
                               append_records)

REPO = Path(__file__).resolve().parents[3]
PS1_T = REPO / "runs" / "panstarrs" / "calib_v1" / "tensors"
ZTF_T = REPO / "runs" / "ztf" / "calib_v1" / "tensors"
OUT = REPO / "runs" / "joint" / "ps1_ztf_v1"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYP = {"ps1": REPO / "surveys" / "panstarrs" / "hypotheses.md",
       "ztf": REPO / "surveys" / "ztf" / "hypotheses.md"}
HYPOTHESIS_VERSION = "joint-ps1-ztf-v1.0 (ps1-hypotheses-v1.0 + ztf-hypotheses-v1.0, mu=0)"

BAND_PAIRS = {"g": (1, 1), "r": (2, 2), "i": (3, 3)}  # name: (ps1 idx, ztf idx)
ZTF_THROUGHPUT_MAG = 0.5
MIN_GOOD_FRAC, WEIGHT_CAP, CLIP_SIGMA, MIN_EPOCHS = 0.7, 20.0, 5.0, 5
MIN_PHASE_EPOCHS = 3
N_REPEATS, RECOVERY_P, N_Z_INTERVALS, SEED = 32, 0.90, 8, 20260821
ZP_REF = 25.0
ALBEDO_REF = 0.1


def phase_ref_doy(mjd: np.ndarray) -> float:
    ang = (mjd % 365.25) / 365.25 * 2 * np.pi
    return float((np.arctan2(np.sin(ang).mean(), np.cos(ang).mean())
                  % (2 * np.pi)) / (2 * np.pi) * 365.25)


def phase_of(mjd: np.ndarray, ref: float) -> np.ndarray:
    d = (mjd % 365.25) - ref
    d = (d + 182.625) % 365.25 - 182.625
    return (np.abs(d) >= 91.3).astype(np.uint8)


def interp_1overz(arr, z_from, z_to):
    """Linear interpolation along the z axis (axis 1 of (T,E,NZ)) in 1/z."""
    q_from, q_to = 1.0 / z_from, 1.0 / z_to
    o = np.argsort(q_from)
    out = np.empty(arr.shape[:2] + (len(z_to),), dtype=np.float32)
    for t in range(arr.shape[0]):
        for e in range(arr.shape[1]):
            out[t, e] = np.interp(q_to, q_from[o], arr[t, e][o])
    return out


def stack(f, v, g, phase):
    """Weighted stack over epochs (axis 1) -> S (T,NZ), n_valid, weights."""
    valid = np.isfinite(f) & np.isfinite(v) & (v > 0) & (g >= MIN_GOOD_FRAC)
    with np.errstate(invalid="ignore", divide="ignore"):
        s_e = np.abs(f) / np.sqrt(np.where(v > 0, v, np.inf))
    valid &= s_e <= CLIP_SIGMA
    w = np.where(valid, 1.0 / np.where(v > 0, v, 1.0), 0.0)
    with np.errstate(all="ignore"):
        wmed = np.nanmedian(np.where(w > 0, w, np.nan), axis=1)
    w = np.minimum(w, WEIGHT_CAP * np.nan_to_num(wmed, nan=np.inf, posinf=np.inf)[:, None])
    fz = np.where(valid, f, 0.0)
    A, B = (fz * w).sum(axis=1), w.sum(axis=1)
    n = (w > 0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
    S[n < MIN_EPOCHS] = np.nan
    return S, n, w, fz, A, B


def size_limit_km(mag, z_au, albedo):
    H = mag - 5.0 * np.log10(z_au * z_au)
    return float(1329.0 / np.sqrt(albedo) * 10 ** (-H / 5.0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--marginal", action="store_true",
                    help="also run the direct ZTF test of the PS1 marginal cells")
    args = ap.parse_args()
    rng = np.random.default_rng(SEED)
    registry = load_target_registry(REGISTRY_PATH)
    OUT.mkdir(parents=True, exist_ok=True)
    hyp_hash = {k: "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()
                for k, p in HYP.items()}
    import time as _time
    fresh = {p.stem for p in ZTF_T.glob("*.npz")
             if _time.time() - p.stat().st_mtime < 600}
    if fresh:
        print(f"skipping {len(fresh)} ZTF tensors modified < 10 min ago "
              f"(ZTF tensor stage still running): {sorted(fresh)}")
    pairs = sorted({p.stem for p in PS1_T.glob("*.npz")}
                   & {p.stem for p in ZTF_T.glob("*.npz")} - fresh)
    print(f"{len(pairs)} endpoint-roles in both archives", flush=True)
    ztf_scale = 10 ** (0.4 * ZTF_THROUGHPUT_MAG)

    report, m90_store, meta = {}, {}, {}
    for stem in pairs:
        e, role = stem.split("__")
        dp, dz = np.load(PS1_T / f"{stem}.npz"), np.load(ZTF_T / f"{stem}.npz")
        zg = dp["z_grid"]
        mu0p = int(np.argmin(np.abs(dp["mu_grid"])))
        mu0z = int(np.argmin(np.abs(dz["mu_grid"])))
        ref = phase_ref_doy(dp["mjd"])
        ph_p, ph_z = phase_of(dp["mjd"], ref), phase_of(dz["mjd"], ref)
        halfcell = (206265.0 / zg[0] - 206265.0 / zg[-1]) / (len(zg) - 1) / 2.0
        for band, (bp, bz) in BAND_PAIRS.items():
            ep, ez = dp["band_idx"] == bp, dz["band_idx"] == bz
            if ep.sum() == 0 and ez.sum() == 0:
                continue
            fp = np.asarray(dp["f"][:, :, :, mu0p, mu0p][:, ep, :], dtype=np.float32)
            vp = np.asarray(dp["v"][:, :, :, mu0p, mu0p][:, ep, :], dtype=np.float32)
            gp = np.asarray(dp["g"][:, :, :, mu0p, mu0p][:, ep, :], dtype=np.float32)
            fz_ = np.asarray(dz["f"][:, :, :, mu0z, mu0z][:, ez, :], dtype=np.float32) * ztf_scale
            vz_ = np.asarray(dz["v"][:, :, :, mu0z, mu0z][:, ez, :], dtype=np.float32) * ztf_scale ** 2
            gz_ = np.asarray(dz["g"][:, :, :, mu0z, mu0z][:, ez, :], dtype=np.float32)
            if ez.sum():
                fz_ = interp_1overz(fz_, dz["z_grid"], zg)
                vz_ = interp_1overz(vz_, dz["z_grid"], zg)
                gz_ = interp_1overz(gz_, dz["z_grid"], zg)
            else:
                fz_ = vz_ = gz_ = np.empty((fp.shape[0], 0, len(zg)), np.float32)
            f = np.concatenate([fp, fz_], axis=1)
            v = np.concatenate([vp, vz_], axis=1)
            g = np.concatenate([gp, gz_], axis=1)
            arch = np.array([0] * ep.sum() + [1] * ez.sum())
            mjd = np.concatenate([dp["mjd"][ep], dz["mjd"][ez]])
            phase = np.concatenate([ph_p[ep], ph_z[ez]])
            seeing = np.concatenate([dp["seeing"][ep],
                                     np.where(dz["seeing"][ez] > 0, dz["seeing"][ez], 2.0)])
            S, n_valid, w, fzero, A, B = stack(f, v, g, phase)
            nt = S.shape[0]
            null_max = np.array([np.nanmax(S[t]) if np.isfinite(S[t]).any() else np.nan
                                 for t in range(1, nt)])
            if not np.isfinite(S[0]).any() or not np.isfinite(null_max).any():
                continue
            T = float(np.nanmax(null_max))
            zi = int(np.nanargmax(S[0]))
            key = f"{e}/{role}/{band}"
            ph_S, ph_n, ar_S, ar_n = {}, {}, {}, {}
            for p_ in (0, 1):
                sel = phase == p_
                wp, fp_ = w[0][sel, zi], fzero[0][sel, zi]
                ph_n[str(p_)] = int((wp > 0).sum())
                ph_S[str(p_)] = float((fp_ * wp).sum() / np.sqrt(wp.sum())) if wp.sum() > 0 else -99.0
            for a_, name in ((0, "ps1"), (1, "ztf")):
                sel = arch == a_
                wp, fp_ = w[0][sel, zi], fzero[0][sel, zi]
                ar_n[name] = int((wp > 0).sum())
                ar_S[name] = float((fp_ * wp).sum() / np.sqrt(wp.sum())) if wp.sum() > 0 else -99.0
            report[key] = {
                "T": T, "control_maxima": np.round(null_max, 2).tolist(),
                "real_max_S": float(S[0][zi]), "real_max_z": float(zg[zi]),
                "exceeds": bool(S[0][zi] > T),
                "n_epochs": {"ps1": int(ep.sum()), "ztf": int(ez.sum())},
                "epoch_range_mjd": [float(mjd.min()), float(mjd.max())],
                "phase_S": ph_S, "phase_n": ph_n,
                "archive_S_at_peak": ar_S, "archive_n_at_peak": ar_n,
                "n_valid_at_peak": int(n_valid[0][zi]),
                "cells_defined_frac": float(np.isfinite(S[0]).mean()),
                "phase_split_all_epochs": [int((phase == 0).sum()), int((phase == 1).sum())],
            }
            # injections at mu = 0 (grid mismatch only), duty 0.5
            E = f.shape[1]
            sig = seeing / 2.3548
            u = np.linspace(0, halfcell, 64)
            mf_e = np.array([float(np.mean(np.exp(-u * u / (4 * s * s)))) for s in sig])
            m90 = np.full(len(zg), np.nan)
            for zi_ in range(len(zg)):
                b0 = B[0][zi_]
                if b0 <= 0 or n_valid[0][zi_] < MIN_EPOCHS:
                    continue
                wg = w[0][:, zi_] * np.nan_to_num(g[0][:, zi_]) * mf_e
                need = max(0.0, T * np.sqrt(b0) - A[0][zi_])
                fmins = []
                for r_ in range(N_REPEATS):
                    keep = rng.random(E) < 0.5
                    den = float(wg[keep].sum())
                    if den > 0:
                        fmins.append(need / den)
                if len(fmins) >= int(0.8 * N_REPEATS):
                    f90 = float(np.percentile(fmins, 100 * RECOVERY_P))
                    if f90 > 0:
                        m90[zi_] = ZP_REF - 2.5 * np.log10(f90)
            m90_store[key] = m90
            meta[key] = {"z_grid": zg, "n_valid_z": n_valid[0],
                         "epoch_range": report[key]["epoch_range_mjd"],
                         "S_real_z": S[0].copy()}
            # split-half (early / late) at the peak, for adjudication
            order = np.argsort(mjd)
            hal = {}
            for name, idx in (("early", order[:len(order) // 2]),
                              ("late", order[len(order) // 2:])):
                wp, fp_ = w[0][idx, zi], fzero[0][idx, zi]
                hal[name] = [float((fp_ * wp).sum() / np.sqrt(wp.sum())) if wp.sum() > 0 else -99.0,
                             int((wp > 0).sum())]
            report[key]["split_half"] = hal
            print(f"  {key:24s} N={ep.sum():3d}+{ez.sum():3d} T={T:5.2f} S={S[0][zi]:5.2f} "
                  f"phase {ph_S['0']:.1f}/{ph_S['1']:.1f} (n {ph_n['0']}/{ph_n['1']}) "
                  f"m90 med={np.nanmedian(m90):5.2f}", flush=True)
        dp.close(); dz.close()

    config = {
        "pipeline_id": "joint-ps1-ztf-stack",
        "hypothesis": "station-kept relay, mu_resid = 0, 550-10000 AU, duty >= 0.5",
        "z_grid": "PS1 360-node 1/z grid; ZTF interpolated in 1/z",
        "flux_scale": ("PS1 star-calibrated AB at ZP 25; ZTF MAGZP scale x 10^(0.4*0.5) "
                       "(asteroid-control throughput)"),
        "band_pairs": {k: f"ps1 {k} + ztf z{k}" for k in BAND_PAIRS},
        "achromatic_assumption": "source colour identical in paired filters",
        "threshold_rule": "max of 8 shared offset-control grid maxima",
        "weight_cap": WEIGHT_CAP, "min_epochs": MIN_EPOCHS, "clip_sigma": CLIP_SIGMA,
        "phase_reference": "circular-median PS1 day-of-year per endpoint-role",
        "n_repeats": N_REPEATS, "recovery_probability": RECOVERY_P, "seed": SEED,
        "candidate_rule": (f"exceedance vetoed if either phase has < {MIN_PHASE_EPOCHS} "
                           "epochs or a phase S < 2; both archives' S reported"),
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
                                          "hypothesis": HYPOTHESIS_VERSION,
                                          "pairs": pairs}),
        pipeline_id="joint-ps1-ztf-stack", pipeline_version="0.1.0", config=config,
        observation_set_hash="see-ps1-and-ztf-tensor-manifests",
        intersection_set_hash="see-ps1-and-ztf-tensor-manifests",
        registry_source_hash=registry.source_hash,
        hypothesis_version=HYPOTHESIS_VERSION,
        environment={"sglseti_commit": commit, "hypothesis_hashes": hyp_hash},
        random_seeds={"numpy": SEED},
        started_utc=datetime.now(timezone.utc).isoformat(),
        finished_utc=datetime.now(timezone.utc).isoformat())

    candidates = []
    for key, r in report.items():
        if not r["exceeds"]:
            continue
        e, role, band = key.split("/")
        pn, pS = r["phase_n"], r["phase_S"]
        if min(pn.values()) < MIN_PHASE_EPOCHS:
            status, reason = "vetoed", f"single-phase data (phase epochs {pn})"
        elif min(pS.values()) < 2.0:
            status, reason = "vetoed", (f"phase-split significances disagree {pS}: "
                                        "static source / artifact at one phase position")
        else:
            # stage-7: split-half persistence and the other paired bands at
            # the same z (achromatic source expected in at least one)
            hal = r["split_half"]
            zi = int(np.argmin(np.abs(meta[key]["z_grid"] - r["real_max_z"])))
            others = {}
            for ob in BAND_PAIRS:
                ok_ = f"{e}/{role}/{ob}"
                if ob != band and ok_ in meta:
                    others[ob] = float(np.nan_to_num(meta[ok_]["S_real_z"][zi], nan=-99.0))
            r["other_bands_at_z"] = others
            if hal["early"][0] < 2.0 or hal["late"][0] < 2.0:
                status, reason = "vetoed", (f"split-half not persistent (early {hal['early']}, "
                                            f"late {hal['late']})")
            elif others and max(others.values()) < 2.0:
                status, reason = "vetoed", (f"absent in the other paired bands at the same z "
                                            f"{ {k: round(v, 2) for k, v in others.items()} }")
            else:
                status, reason = "retained", None
        candidates.append(Candidate(
            candidate_id=stable_id("cnd", {"endpoint": e, "role": role, "band": band,
                                           "source": "joint-ps1-ztf-v1-exceedance",
                                           "hypothesis": HYPOTHESIS_VERSION}),
            analysis_run_id=run.analysis_run_id, endpoint_id=e, role=role,
            observation_ids=(), fitted_z_au=r["real_max_z"],
            fitted_residual_motion={"mu_arcsec_yr": [0.0, 0.0]},
            model_comparison={"real_max_S": r["real_max_S"], "threshold_8_controls": r["T"],
                              "phase_S": pS, "phase_n": pn,
                              "archive_S": r["archive_S_at_peak"],
                              "archive_n": r["archive_n_at_peak"],
                              "split_half": r["split_half"],
                              "other_bands_at_z": r.get("other_bands_at_z")},
            status=status, veto_reason=reason,
            extra={"band": band, "origin": "joint exceedance census"}))

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
                limit = {"value": mlim, "unit": "ab_mag_ps1_scale", "band": f"ps1 {band} + ztf z{band}",
                         "epochs_per_cell_max": int(meta[key]["n_valid_z"][lo:hi].max()),
                         "both_parallax_phases": bool(
                             min(report[key]["phase_split_all_epochs"]) >= MIN_PHASE_EPOCHS),
                         "interpretation": ("reflected sunlight / self-luminous optical; "
                                            "never thermal; achromatic within paired band"),
                         "reflected_light_size_km_albedo_0p1": round(size_limit_km(mlim, zmid, ALBEDO_REF), 1),
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
                duty_cycle_range=(0.5, 1.0), residual_motion_bound_arcsec_per_yr=0.0,
                morphology="gaussian-psf-point-source", kind=kind,
                recovery_probability=RECOVERY_P if kind == "recovery_curve" else None,
                flux_limit=limit, false_alarm_rate=None,
                extra={"threshold": report[key]["T"], "n_epochs": report[key]["n_epochs"]}))

    rec = OUT / "records"
    append_records(rec / "analysis_run.jsonl", [run])
    append_records(rec / "constraint.jsonl", constraints)
    append_records(rec / "candidate.jsonl", candidates)
    np.savez_compressed(OUT / "m90_curves.npz",
                        **{k.replace("/", "__"): v for k, v in m90_store.items()})
    (OUT / "threshold_report.json").write_text(json.dumps(report, indent=2))
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

    if args.marginal:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from marginal_ztf_test import run_marginal_tests
        run_marginal_tests(OUT, run.analysis_run_id)


if __name__ == "__main__":
    main()
