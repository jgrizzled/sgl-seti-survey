"""TESS crossings threshold freeze v1.0.

Declares, before any signal statistic is formed, the statistics,
control construction, threshold rule, dev/confirmatory split and
trials accounting, bound to content hashes of the frozen inputs
(hypotheses.md, coverage_gate_v1.ecsv, coverage_refined_v1.ecsv,
tic_cut_v1.json, tesscut_recon.md). Channel structure per hypotheses
v1.0: B = discovery channel (spatial ring controls, window-length
independent); A = constraint-only light-curve products (temporal
8-control family unreachable for >= 9 d windows in 27 d sectors).
Two frozen statistics per B unit — the d = 1 chord amplitude and the
per-cadence pulse maximum — each with its own controls and trial.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
D = REPO / "surveys" / "tess-crossings"
INPUTS = {
    "hypotheses": D / "hypotheses.md",
    "coverage_gate": D / "results" / "coverage_gate_v1.ecsv",
    "coverage_refined": D / "results" / "coverage_refined_v1.ecsv",
    "tic_cut": D / "results" / "tic_cut_v1.json",
    "usability": D / "results" / "usability_v1.json",
    "recon": D / "notes" / "tesscut_recon.md",
}
SEED = 20260827
MIN_INWINDOW_CADENCES = 20
MIN_OFFWINDOW_CADENCES = 200
RING_OFFSETS_ARCSEC = [(126.0, 0.0), (-126.0, 0.0), (189.0, 0.0),
                       (-189.0, 0.0), (252.0, 0.0), (-252.0, 0.0),
                       (0.0, 189.0), (0.0, -189.0)]
N_PULSE_PERIODS = 6


def sha(p: Path) -> str:
    return "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ref = Table.read(INPUTS["coverage_refined"])
    tic = json.loads(INPUTS["tic_cut"].read_text())
    tic_status = {t["target_id"]: t["status"] for t in tic["targets"]}
    usab = json.loads(INPUTS["usability"].read_text())
    usable = {(r["channel"], r["target_id"], r["event_id"],
               r["sector"]): r["status"] == "usable"
              for r in usab["rows"]}

    units, a_ref, off_science = [], [], []
    for r in ref:
        row = {
            "channel": str(r["channel"]),
            "target_id": str(r["target_id"]),
            "event_id": str(r["event_id"]),
            "radius_au": float(r["radius_au"]),
            "b_rsun": float(r["b_rsun"]),
            "t_ca_mjd": float(r["t_ca_mjd"]),
            "window_days": float(r["window_days"]),
            "sector": int(r["sector"]),
            "cadence_s": float(r["cadence_s"]),
            "n_inwindow_primary": int(r["n_inwindow_primary"]),
            "n_offwindow_primary": int(r["n_offwindow_primary"]),
        }
        row["usable"] = usable.get(
            (row["channel"], row["target_id"], row["event_id"],
             row["sector"]), False)
        if row["channel"] == "B":
            if not row["usable"]:
                row["class"] = "off_science_array"
                off_science.append(row)
            elif (row["n_inwindow_primary"] >= MIN_INWINDOW_CADENCES
                    and row["n_offwindow_primary"]
                    >= MIN_OFFWINDOW_CADENCES):
                p0 = 2 * row["cadence_s"]
                p1 = row["window_days"] * 86400.0
                row["pulse_periods_s"] = [round(float(p), 1)
                                          for p in np.geomspace(
                                              p0, p1, N_PULSE_PERIODS)]
                units.append(row)
        else:
            st = tic_status.get(row["target_id"], "no_tic_match")
            row["tic_status"] = st
            if not row["usable"]:
                row["class"] = "off_science_array"
                off_science.append(row)
            elif (st in ("ok", "marginal")
                    and row["n_inwindow_primary"] > 0):
                a_ref.append(row)

    # dev/confirmatory split (unit = target). Frozen from hypotheses:
    # ross-128 (the single wide-partial unit, no grazing exposure) is
    # dev; the grazing targets are fully confirmatory. The seed is
    # recorded for the record although the split is forced.
    targets = sorted(set(u["target_id"] for u in units))
    grazing = sorted(set(u["target_id"] for u in units
                         if u["radius_au"] < 0.1))
    dev = sorted(set(targets) - set(grazing))
    conf = grazing
    n_trials = 2 * len(units)

    freeze = {
        "freeze_version": "tess-crossings-thresholds-v1.1",
        "amendment_v1_1": "science-array usability screen added "
            "pre-search (results/usability_v1.json): TESScut serves "
            "collateral (non-science) CCD pixels with aperture=1 and "
            "~zero flux; three cutouts sat on/across the science edge "
            "(cols 45-2092). ross-128 B (the v1.0 dev unit), "
            "teegarden A s44 and gj-1276 A s42 are reclassified "
            "off_science_array (nominal-covered/unusable). With no "
            "searchable non-grazing unit left, machinery validation "
            "is reassigned to the ZERO-TRIAL channel-A reference rows "
            "(teegarden s71, van-maanen s43) before any confirmatory "
            "unit is touched; all 6 B units are confirmatory. No "
            "signal statistic existed when this amendment was made.",
        "frozen_at": "2026-08-24",
        "hypothesis_version": "tess-crossings-hypotheses-v1.0",
        "input_hashes": {k: sha(p) for k, p in INPUTS.items()},
        "scope": "channel B searchable (spatial controls); channel A "
                 "constraint-only light-curve products (temporal "
                 "8-control family unreachable, hypotheses v1.0); "
                 "A 1.0 AU rung programme-wide constraint-only, "
                 "gate-level coverage only",
        "statistics": {
            "substrate": "per-cadence matched-filter forced photometry "
                         "on TESScut FFI cubes (SPOC FLUX with FLUX_BKG "
                         "subtracted; Gaussian kernel at the effective "
                         "PSF FWHM, sub-pixel phase via stamp response); "
                         "track position per cadence (Earth-center "
                         "apparent relay, <= 0.04 px budget); flux scale "
                         "= per-cutout TIC star calibration through the "
                         "identical kernel, fallback T = 20.44 - 2.5 "
                         "log10(e-/s), verified to <= 0.2 mag before "
                         "any depth is quoted",
            "chord": "S_c = weighted least-squares amplitude of the "
                     "frozen chord top-hat over the window on the "
                     "per-cadence series; baseline = off-window "
                     "in-sector cadences (robust mean, 3x3sigma clip); "
                     "empirical variance rescale k = "
                     "median(r^2/v)/0.4549 (floor 1) from off-window "
                     "residuals; WEIGHT_CAP 20x; one-sided positive",
            "pulse": "S_p = max over in-window cadences of the "
                     "per-cadence S/N (same baseline, k, mask); "
                     "sensitive to pulses >= 1 cadence; the frozen "
                     "per-unit period grids parameterize the "
                     "*injection* models, not a search grid — S_p "
                     "itself is grid-free",
            "gates": f"unit requires >= {MIN_INWINDOW_CADENCES} "
                     f"in-window and >= {MIN_OFFWINDOW_CADENCES} "
                     "off-window primary cadences (all 7 gate rows "
                     "pass)",
        },
        "controls": {
            "kind": "spatial ring trajectories (window-length "
                    "independent — the reason channel B keeps "
                    "discovery power)",
            "n": 8, "offsets_arcsec": RING_OFFSETS_ARCSEC,
            "note": "offsets are 6/9/12 px at the 21 arcsec pixel; "
                    "each control carries both statistics",
        },
        "threshold_rule": {
            "T": "per statistic: max one-sided control value over the "
                 "8 designated rings",
            "exceedance": "S > max(T, 0); margin reported",
            "per_search_crossing_probability": "1/9 per statistic",
            "n_units": len(units), "n_statistics_per_unit": 2,
            "n_trials": n_trials,
            "expected_control_crossings": round(n_trials / 9.0, 1),
            "disposition": "every exceedance individually adjudicated; "
                           "no silent drops",
        },
        "veto_ladder": [
            "1. SkyBoT known-object census at the exceedance cadence "
            "(mandatory for every pulse-statistic exceedance: 21-arcsec "
            "pixels in the near-ecliptic asteroid stream)",
            "2. chord-shape/timing test: amplitude consistency with "
            "the predicted ingress/egress times; flare morphology "
            "(fast-rise-exponential-decay) and asteroid transit "
            "(hours-scale bump) both fail it",
            "3. straylight/scattered-light: QUALITY-block adjacency "
            "and background-correlation (annotation unless the strict "
            "mask removes the exceedance)",
            "4. recurrence at a second covered window where one "
            "exists; none does in-archive for the grazing units, so "
            "surviving exceedances are retained-ambiguous with a "
            "light-curve morphology report — promotion to candidate "
            "requires independent recurrence (a future ecliptic "
            "sector or another archive), never in-window evidence "
            "alone",
        ],
        "quality": {
            "primary": "QUALITY == 0",
            "strict": "additionally drop cadences within 0.05 d of any "
                      "straylight-flagged block and cadences with "
                      "|POS_CORR| > 0.5 px; strict re-run reported for "
                      "any exceedance",
        },
        "alpha": {"fwer_alpha": 0.05, "ks_alpha": 0.01},
        "split": {
            "seed": SEED, "unit": "target",
            "dev": dev, "confirmatory": conf,
            "note": "v1.1: no searchable dev unit exists (ross-128 "
                    "off-science). Machinery validation runs on the "
                    "zero-trial A reference rows; every searchable B "
                    "unit is confirmatory, entered only after that "
                    "validation is clean.",
        },
        "off_science_array": {"n": len(off_science),
                              "rows": off_science},
        "search_units": {"n_units": len(units), "units": units},
        "a_reference_rows": {
            "n": len(a_ref), "rows": a_ref,
            "note": "constraint-only: resolved on-star crossing light "
                    "curves + injection-calibrated reference depths "
                    "(threshold-free, labeled); zero trials",
        },
        "completeness": {
            "injections": ">= 100 draws per unit per temporal model "
                          "via the SPOC per-camera/CCD PRF through the "
                          "identical kernel (stamp response); temporal "
                          "models: d = 1 chord and d = 0.1 boxcar at "
                          "the frozen per-unit period grids; magnitude "
                          "grid T 10-18 (0.5 steps); recovery against "
                          "each statistic's frozen threshold",
        },
    }
    body = json.dumps(freeze, indent=1, sort_keys=True)
    freeze["freeze_content_hash"] = (
        "sha256:" + hashlib.sha256(body.encode()).hexdigest())
    out = D / "configs" / "threshold_freeze_v1.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(freeze, indent=1, sort_keys=True) + "\n")
    print(json.dumps({
        "B_units": len(units),
        "unit_list": [f"{u['target_id']}/r{u['radius_au']}"
                      f"/s{u['sector']}/n{u['n_inwindow_primary']}"
                      for u in units],
        "n_trials": n_trials,
        "expected_control_crossings":
            freeze["threshold_rule"]["expected_control_crossings"],
        "A_reference_rows": len(a_ref),
        "dev": dev, "confirmatory": conf,
        "hash": freeze["freeze_content_hash"]}, indent=2))


if __name__ == "__main__":
    main()
