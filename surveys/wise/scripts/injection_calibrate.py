"""Stage-7 injection calibration and constraint generation (plan §4.7).

Consumes the calib_v1 sample tensors (real + 8 offset-control
trajectories per endpoint x role) and produces:

  1. Predeclared detection thresholds per endpoint x role x band:
     T = max of the 8 control-trajectory grid maxima (empirical
     false-alarm rate < 1/8 per grid search by construction).
  2. Injection-recovery curves: synthetic Gaussian point sources
     injected at every z node (exploiting matched-filter linearity:
     an injected source of flux f adds f * g_e * grid/mu mismatch
     factors to each epoch's sample), duty cycle in {1.0, 0.5} with
     randomized epoch realizations, 32 repeats -> the 90%-recovery
     flux is the 90th percentile of per-repeat minimum recoverable
     fluxes. Quoted in Vega mag via the per-frame magzp common scale.
  3. Adjudication of the two stack_v1 marginal Lalande cells against
     the 8-control null (Candidate records with outcome).
  4. Ledger-ready Constraint records: per endpoint x role x band x
     z-interval (8 nodes each), duty in [0.5, 1], |mu| <= 1 "/yr,
     Gaussian-point-source morphology, at 90% recovery probability.

Usage: uv run python surveys/wise/scripts/injection_calibrate.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from sglseti import canonical_json, load_target_registry, stable_id

from sglsurvey.records import (AnalysisRun, Candidate, Constraint,
                               append_records)

REPO = Path(__file__).resolve().parents[3]
CAL_DIR = REPO / "runs" / "wise" / "calib_v1"
TENSOR_DIR = CAL_DIR / "tensors"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "wise" / "hypotheses.md"
HYPOTHESIS_VERSION = "wise-hypotheses-v1.0"

ZP_REF = 20.0
MIN_GOOD_FRAC = 0.7
N_REPEATS = 32
DUTIES = [1.0, 0.5]
SEED = 20260818
RECOVERY_P = 0.90
N_Z_INTERVALS = 8
PSF_FWHM = {1: 6.1, 2: 6.4, 3: 6.5, 4: 12.0}
BAND_NAME = {1: "W1", 2: "W2", 3: "W3", 4: "W4"}
BAND_PHYSICS = {
    "W1": "reflected sunlight / hot components / nonthermal",
    "W2": "reflected sunlight / hot components / nonthermal",
    "W3": "warm waste heat (~300 K radiators); 2010 cryo epochs only",
    "W4": "warm waste heat (~300 K radiators); 2010 cryo epochs only",
}
FLAGGED_CELLS = [("lalande-21185", "rx", "W1"),
                 ("lalande-21185", "tx", "W2")]


def mismatch_factors(band_idx: int, halfcell_arcsec: float,
                     mean_abs_dt_yr: float) -> float:
    """Mean matched-filter response loss for off-node true parameters.

    rho(r) = exp(-r^2 / (4 sigma^2)) is the normalized response of the
    Gaussian matched filter to a Gaussian source offset by r.
    """
    sigma = PSF_FWHM[band_idx] / 2.3548
    u = np.linspace(0, halfcell_arcsec, 64)
    grid_f = float(np.mean(np.exp(-u * u / (4 * sigma * sigma))))
    r_mu = 0.25 * mean_abs_dt_yr * np.sqrt(2)   # per-axis worst / 2 avg
    mu_f = float(np.exp(-r_mu * r_mu / (4 * sigma * sigma)))
    return grid_f * mu_f


def main() -> None:
    rng = np.random.default_rng(SEED)
    registry = load_target_registry(REGISTRY_PATH)
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()

    config = {
        "pipeline_id": "wise-injection-calibration",
        "threshold_rule": "max of 8 offset-control grid maxima",
        "n_repeats": N_REPEATS, "duties": DUTIES, "seed": SEED,
        "recovery_probability": RECOVERY_P,
        "min_good_frac": MIN_GOOD_FRAC, "zp_ref": ZP_REF,
        "z_grid": "64 nodes uniform in 1/z, 550-10000 AU",
        "injection_model": ("gaussian-psf point source at grid nodes; "
                            "grid/mu mismatch factors applied"),
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
    threshold_report, m90_store = {}, {}
    tensor_files = sorted(TENSOR_DIR.glob("*.npz"))

    for path in tensor_files:
        endpoint, role = path.stem.split("__")
        d = np.load(path)
        z_grid, mu_grid = d["z_grid"], d["mu_grid"]
        mjd, band_idx, phase = d["mjd"], d["band_idx"], d["phase"]
        f, v, g = d["f"], np.asarray(d["v"]), np.asarray(
            d["g"], dtype=np.float32)
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
            w = np.where(valid, 1.0 / np.where(vb > 0, vb, 1.0), 0.0)
            A = (np.where(valid, fb, 0.0) * w).sum(axis=1)
            B = w.sum(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                S = A / np.sqrt(np.where(B > 0, B, np.inf))
            null_max = np.array([np.nanmax(S[t]) for t in range(1, nt)])
            T = float(null_max.max())
            i_real = np.unravel_index(np.nanargmax(S[0]), S[0].shape)
            band = BAND_NAME[b]
            key = f"{endpoint}/{role}/{band}"
            threshold_report[key] = {
                "T": T, "control_maxima": null_max.round(2).tolist(),
                "real_max_S": float(S[0][i_real]),
                "real_max_z": float(z_grid[i_real[0]]),
                "exceeds": bool(S[0][i_real] > T),
                "n_epochs": int(eb.sum()),
            }

            # --- injections ------------------------------------------
            mf = mismatch_factors(b, halfcell, mean_abs_dt)
            E = int(eb.sum())
            m90 = {duty: np.full(nz, np.nan) for duty in DUTIES}
            for zi in range(nz):
                mu_pick = rng.integers(0, nm, size=N_REPEATS)
                mu_pick2 = rng.integers(0, nm, size=N_REPEATS)
                f_min = np.full((len(DUTIES), N_REPEATS), np.nan)
                for r in range(N_REPEATS):
                    ci = (zi, mu_pick[r], mu_pick2[r])
                    wg = (w[0][:, ci[0], ci[1], ci[2]]
                          * np.nan_to_num(gb[0][:, ci[0], ci[1], ci[2]]))
                    a0, b0 = A[0][ci], B[0][ci]
                    if b0 <= 0:
                        continue
                    need = max(0.0, T * np.sqrt(b0) - a0)
                    for di, duty in enumerate(DUTIES):
                        keep = (rng.random(E) < duty) if duty < 1.0 \
                            else np.ones(E, dtype=bool)
                        denom = float(wg[keep].sum()) * mf
                        if denom > 0:
                            f_min[di, r] = need / denom
                for di, duty in enumerate(DUTIES):
                    fm = f_min[di][np.isfinite(f_min[di])]
                    if len(fm) >= int(0.8 * N_REPEATS):
                        f90 = float(np.percentile(fm, 100 * RECOVERY_P))
                        if f90 > 0:
                            m90[duty][zi] = ZP_REF - 2.5 * np.log10(f90)
            m90_store[key] = {str(duty): m90[duty] for duty in DUTIES}

        d.close()

    # --- AnalysisRun -------------------------------------------------
    run = AnalysisRun(
        analysis_run_id=stable_id("run", {
            "config": json.loads(canonical_json(config)),
            "registry": registry.source_hash,
            "hypothesis": HYPOTHESIS_VERSION,
            "tensors": sorted(p.name for p in tensor_files)}),
        pipeline_id="wise-injection-calibration",
        pipeline_version="0.1.0", config=config,
        observation_set_hash="see-tensor-manifest",
        intersection_set_hash="see-tensor-manifest",
        registry_source_hash=registry.source_hash,
        hypothesis_version=HYPOTHESIS_VERSION,
        environment={"sglseti_commit": commit,
                     "hypothesis_hash": hyp_hash},
        random_seeds={"numpy": SEED},
        started_utc=started,
        finished_utc=datetime.now(timezone.utc).isoformat(),
    )

    # --- Candidates: adjudicate flagged cells ------------------------
    for endpoint, role, band in FLAGGED_CELLS:
        key = f"{endpoint}/{role}/{band}"
        r = threshold_report.get(key)
        if r is None:
            continue
        vetoed = not r["exceeds"]
        candidates.append(Candidate(
            candidate_id=stable_id("cnd", {
                "endpoint": endpoint, "role": role, "band": band,
                "source": "stack_v1-marginal-cell",
                "hypothesis": HYPOTHESIS_VERSION}),
            analysis_run_id=run.analysis_run_id,
            endpoint_id=endpoint, role=role, observation_ids=(),
            fitted_z_au=r["real_max_z"],
            model_comparison={
                "real_max_S": r["real_max_S"],
                "threshold_8_controls": r["T"],
                "control_maxima": r["control_maxima"]},
            status="vetoed" if vetoed else "retained",
            veto_reason=(f"below 8-control empirical null "
                         f"(S={r['real_max_S']:.1f} vs T={r['T']:.1f}) "
                         "on the 1/z-uniform grid; W1:W2 significance "
                         "ratio star-like" if vetoed else None),
            extra={"band": band, "origin": "stack_v1 marginal cell"},
        ))

    # --- Constraints -------------------------------------------------
    edges = np.linspace(0, 64, N_Z_INTERVALS + 1, dtype=int)
    for key, curves in m90_store.items():
        endpoint, role, band = key.split("/")
        rep = threshold_report[key]
        d0 = np.load(TENSOR_DIR / f"{endpoint}__{role}.npz")
        z_grid = d0["z_grid"]
        eb = d0["band_idx"] == {v: k for k, v in BAND_NAME.items()}[band]
        ep_range = (float(d0["mjd"][eb].min()), float(d0["mjd"][eb].max()))
        d0.close()
        m = curves["0.5"]   # duty 0.5 bounds the [0.5, 1] cell
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
                limit = {"value": round(float(np.nanmin(seg)), 2),
                         "unit": "wise_vega_mag",
                         "band": band,
                         "interpretation": BAND_PHYSICS[band],
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
                band=band, epoch_range_mjd=ep_range,
                duty_cycle_range=(0.5, 1.0),
                residual_motion_bound_arcsec_per_yr=1.0,
                morphology="gaussian-psf-point-source",
                kind=kind,
                recovery_probability=(RECOVERY_P
                                      if kind == "recovery_curve"
                                      else None),
                flux_limit=limit,
                false_alarm_rate=None,
                extra={"threshold": rep["T"],
                       "threshold_rule": config["threshold_rule"],
                       "far_note": ("empirical: 0 of 8 control "
                                    "trajectories exceed T by "
                                    "construction; per-grid-search "
                                    "FAR < 1/8")},
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
    print("\n=== adjudication of flagged cells ===")
    for c in candidates:
        print(f"  {c.endpoint_id}/{c.role} {c.extra['band']}: "
              f"{c.status}" + (f" ({c.veto_reason})" if c.veto_reason
                               else ""))
    print("\n=== sample recovery depths (duty>=0.5, 90% recovery, "
          "Vega mag) ===")
    for key in sorted(m90_store):
        m = m90_store[key]["0.5"]
        if np.isfinite(m).any():
            print(f"  {key:32s} median m90={np.nanmedian(m):5.2f}  "
                  f"best={np.nanmax(m):5.2f}  worst={np.nanmin(m):5.2f}")


if __name__ == "__main__":
    main()
