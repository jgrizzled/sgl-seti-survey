"""DASCH crossings threshold freeze v1.0.

Emits configs/threshold_freeze_v1.json binding the frozen search
constructions/constants to content hashes of hypotheses.md and the
coverage products. No photometric data is touched. Companion prose:
thresholds.md.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import daschlib as dl  # noqa: E402

D = Path(__file__).resolve().parents[1]

FREEZE = {
    "version": "v1.0",
    "hypotheses": "hypotheses.md v1.0 (FROZEN 2026-08-26)",
    "substrate": "catalogue-level (D3): lightcurve rows for A, "
                 "platephot rows for B; APASS refcat; cutouts "
                 "vetting-only",
    # channel B (locus-track hit statistic)
    "b_construction": {
        "z_grid_au": [550.0, 1000.0, 2500.0, 5500.0, 10000.0],
        "match_radius_arcsec": 10.0,
        "match_radius_sensitivity_arcsec": [5.0, 15.0],
        "candidate_classes": ["uncatalogued (blank ref_number)",
                              "refcat-matched brighter than stdmag by "
                              ">= flux_anomaly_mag"],
        "flux_anomaly_mag": 1.0,
        "fatal_aflags": ["SUSPECTED_DEFECT", "PICKERING_WEDGE",
                          "MULT_EXP_UNMATCHED", "MULT_EXP_BLEND",
                          "REJECTED_BLEND", "UNCERTAIN_DATE"],
        "fatal_aflags_mask": (1 << 25) | (1 << 24) | (1 << 8)
                             | (1 << 10) | (1 << 27) | (1 << 9),
        "strict_extra_aflags": ["CASE_B_BLEND", "CASE_C_BLEND",
                                 "CASE_BC_BLEND", "SXT_BLEND",
                                 "BAD_PLATE_QUALITY"],
        "strict_extra_mask": (1 << 20) | (1 << 21) | (1 << 22)
                             | (1 << 26) | (1 << 7),
        "ring_offsets_arcsec": [[60, 0], [-60, 0], [90, 0], [-90, 0],
                                 [120, 0], [-120, 0], [0, 90], [0, -90]],
    },
    # channel A (on-star window excess / limits-only hits)
    "a_construction": {
        "statistic": "robust z: (median in-window magcal_magdep - "
                     "off-window median) / (MAD_off * 1.4826 / sqrt(n_in)); "
                     "brightening positive",
        "baseline": "all usable off-window rows of the unit's "
                    "lightcurve outside every covered window padded to "
                    "2x half-width",
        "min_offwindow_epochs": 20,
        "limits_only_regime": "star fainter than the p90 in-window "
                              "limiting_mag_local: hit-count statistic "
                              "at the star position (B-construction "
                              "candidate classes, r_match identical)",
        "temporal_offsets_days": [-97, -71, -47, -23, 23, 47, 71, 97],
        "temporal_redraws_days": [-127, -113, 113, 127],
        "mirror_gate": True,
    },
    # decision rule (house 1/9 budget)
    "rule": {
        "exceedance": "S > max(T over the unit's 8 controls, 0); "
                      "margin reported",
        "statistics": ["S_event (max single-window)",
                        "S_stack (all covered windows)"],
        "expected_control_crossing_rate": "1/9 per trial",
    },
    "row_quality": {
        "reject_flag": "rows with nonzero archive reject_flag excluded",
        "too_bright": "rows with TOO_BRIGHT (aflags bit 29) or "
                      "SATURATED (bflags bit 2) unusable (D6)",
    },
}


def main() -> int:
    coverage = D / "results" / "coverage_v1_windows.ecsv"
    summary_f = D / "results" / "coverage_v1_summary.json"
    summary = json.loads(summary_f.read_text())

    # searched-unit list and trial tally from the frozen coverage record
    units = {u: s for u, s in summary["units"].items()
             if s["covered_windows"] > 0}
    ledger = {u: s for u, s in summary["units"].items()
              if s["covered_windows"] == 0}
    dev_units = ["van-maanen/A_0.1AU",     # forced_dev (D1)
                 "teegarden/A_0.1AU",      # limits-only machinery
                 "wolf-359/B_2.5Rs"]       # B-chain machinery
    conf_units = sorted(u for u in units if u not in dev_units)
    assert all(u in units for u in dev_units)

    FREEZE["units"] = {
        "searched": sorted(units),
        "dev": dev_units,
        "confirmatory": conf_units,
        "coverage_without_statistic": sorted(ledger),
        "trials": {"per_unit": 2,
                   "dev": 2 * len(dev_units),
                   "confirmatory": 2 * len(conf_units),
                   "total": 2 * len(units),
                   "expected_control_crossings_dev":
                       round(2 * len(dev_units) / 9.0, 2),
                   "expected_control_crossings_confirmatory":
                       round(2 * len(conf_units) / 9.0, 2)},
    }
    FREEZE["inputs"] = {
        "hypotheses_sha256": dl.sha256_file(D / "hypotheses.md"),
        "coverage_windows_sha256": dl.sha256_file(coverage),
        "coverage_summary_sha256": dl.sha256_file(summary_f),
        "events_sha256": summary["input_events_sha256"],
    }
    FREEZE["frozen_utc"] = datetime.now(timezone.utc).isoformat(
        timespec="seconds")

    out = D / "configs" / "threshold_freeze_v1.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(FREEZE, indent=1) + "\n")
    print("units searched:", len(units), "dev:", len(dev_units),
          "confirmatory:", len(conf_units), "ledger:", len(ledger))
    print("trials:", FREEZE["units"]["trials"])
    print("wrote", out, dl.sha256_file(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
