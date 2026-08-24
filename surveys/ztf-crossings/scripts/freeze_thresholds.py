"""ZTF crossings threshold freeze v1.0.

Declares, before any pixel is touched, the search statistics, control
constructions, threshold rule, dev/confirmatory split, and trials
accounting for both channels, binding them to content hashes of the
frozen inputs (hypotheses.md, coverage_v1_events.ecsv,
saturation_cut_v1.ecsv). Conventions inherited from the v2 engine:
8 designated ring controls (radii 20/30/40 arcsec), exceedance ratio
R = S/T with T = max over controls (per-search crossing probability
1/9), WEIGHT_CAP 20x effective-epoch floor, split seed 20260822,
dev fraction 0.30, fwer_alpha 0.05 / ks_alpha 0.01.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
D = REPO / "surveys" / "ztf-crossings"
INPUTS = {
    "hypotheses": D / "hypotheses.md",
    "coverage": D / "results" / "coverage_v1_events.ecsv",
    "saturation": D / "results" / "saturation_cut_v1.ecsv",
}
SEED, DEV_FRACTION = 20260822, 0.30
BANDS = {"g": "zg", "r": "zr", "i": "zi"}


def sha(p: Path) -> str:
    return "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    cov = Table.read(INPUTS["coverage"])
    sat = Table.read(INPUTS["saturation"])
    searchable = {(str(r["target_id"]), b): r[f"status_{b}"] in ("ok", "marginal")
                  for r in sat for b in BANDS}

    # -- channel A search units: (target, band, radius) with >=1 covered
    #    window and a searchable band; primary statistic needs >=2 windows.
    A = cov[cov["channel"] == "A"]
    units_A, singles_A = [], []
    for rad in (0.1, 1.0):
        sr = A[np.asarray(A["radius_au"]) == rad]
        for tid in sorted(set(str(t) for t in sr["target_id"])):
            sub = sr[np.asarray([str(x) == tid for x in sr["target_id"]])]
            for b in ("g", "r", "i"):
                if not searchable.get((tid, b), False):
                    continue
                n_win = int((np.asarray(sub[f"n_{b}_pri"]) > 0).sum())
                if n_win >= 2:
                    units_A.append({"target_id": tid, "band": b,
                                    "radius_au": rad, "n_windows": n_win})
                elif n_win == 1:
                    singles_A.append({"target_id": tid, "band": b,
                                      "radius_au": rad})

    # -- channel B search units: (event, radius, band); the z grid is a
    #    nested family (epochs are shared across z), so the unit statistic
    #    is max over z and controls take the same max -- one trial per unit.
    B = cov[cov["channel"] == "B"]
    agg = {}
    for r in B:
        for b in ("g", "r"):          # i opportunistic, excluded from units
            n = int(r[f"n_{b}_pri"])
            if n == 0:
                continue
            k = (str(r["event_id"]), str(r["target_id"]),
                 float(r["radius_au"]), b)
            agg[k] = max(agg.get(k, 0), n)
    units_B, singles_B = [], []
    for (eid, tid, rad, b), n in sorted(agg.items()):
        u = {"event_id": eid, "target_id": tid, "radius_au": rad,
             "band": b, "n_epochs_max_z": n}
        (units_B if n >= 2 else singles_B).append(u)

    # -- dev/confirmatory split (unit = target), stratified
    rng = np.random.default_rng(SEED)
    strata_A = {}
    for u in units_A + singles_A:
        tid = u["target_id"]
        gr = ("gr" if searchable.get((tid, "g")) and searchable.get((tid, "r"))
              else "g" if searchable.get((tid, "g")) else "r")
        strata_A[tid] = gr
    dev_A = []
    for cls in sorted(set(strata_A.values())):
        ts = sorted(t for t, c in strata_A.items() if c == cls)
        n_dev = max(1, round(DEV_FRACTION * len(ts))) if len(ts) > 1 else 0
        dev_A += sorted(rng.choice(ts, size=n_dev, replace=False).tolist())
    conf_A = sorted(set(strata_A) - set(dev_A))

    wide_B = sorted({u["target_id"] for u in units_B + singles_B
                     if u["radius_au"] == 0.1})
    grazing_targets = sorted({u["target_id"] for u in units_B + singles_B
                              if u["radius_au"] < 0.1})
    # grazing-rung windows are temporal subsets of the same targets'
    # wide-rung windows: exclude grazing targets from dev eligibility so
    # no dev exposure leaks into the fully-confirmatory grazing rung
    dev_eligible = sorted(set(wide_B) - set(grazing_targets))
    n_dev_B = (max(1, round(DEV_FRACTION * len(dev_eligible)))
               if len(dev_eligible) > 1 else 0)
    dev_B = sorted(rng.choice(dev_eligible, size=n_dev_B,
                              replace=False).tolist())

    n_searches = len(units_A) + len(units_B)
    freeze = {
        "freeze_version": "ztf-crossings-thresholds-v1.0",
        "frozen_at": "2026-08-23",
        "hypothesis_version": "ztf-crossings-hypotheses-v1.0",
        "input_hashes": {k: sha(p) for k, p in INPUTS.items()},
        "statistics": {
            "A": {
                "substrate": "PSF forced photometry on ZTF difference images "
                             "at the star position (star lives in the reference; "
                             "difference flux is excess by construction)",
                "per_window": "inverse-variance weighted mean in-window excess "
                              "S_w with N_eff (WEIGHT_CAP 20x median epoch weight)",
                "primary": "S_A = weighted stack of S_w over covered windows; "
                           "requires >= 2 covered windows; one-sided positive",
                "single_window_units": "constraint-only, never candidates",
                "annulus": "<=10 arcsec station annulus searched with the "
                           "channel-B construction minus the z-track",
            },
            "B": {
                "substrate": "PSF forced photometry on ZTF difference images "
                             "along the per-(event,z) apparent relay track",
                "statistic": "S = max over the nested z family of the weighted "
                             "in-window stack along the (event, z) track; "
                             "controls take the same max; one-sided positive",
                "single_epoch_units": "reportable 'single-epoch exceedance' "
                                      "class; promotion to candidate requires "
                                      "recurrence at a second covered window",
            },
        },
        "controls": {
            "A": {"kind": "temporal pseudo-windows", "n": 8,
                  "designated_offsets_days": [-97, -71, -47, -23, 23, 47, 71, 97],
                  "notes": "same window durations/counts as the real set; "
                           "offsets avoid lunation multiples (29.5k) and the "
                           "semiannual crossing spacing; any pseudo-window "
                           "overlapping a real window of either rung is "
                           "re-drawn to the next offset in a frozen sequence "
                           "(+/-113, +/-127 d)"},
            "B": {"kind": "spatial ring trajectories", "n": 8,
                  "ring_radii_arcsec": [20.0, 30.0, 40.0],
                  "designated": "identical designated_indices construction as "
                                "ztf-v2 (perpendicular offsets follow the track)"},
        },
        "threshold_rule": {
            "T": "max one-sided control statistic over the 8 designated controls",
            "exceedance": "R = S / T > 1",
            "per_search_crossing_probability": "1/9",
            "expected_control_crossings": round(n_searches / 9.0, 1),
            "disposition": "every exceedance individually adjudicated by the "
                           "frozen veto ladder (census, chromatic, rate, "
                           "recurrence); no silent drops",
        },
        "quality": {"mask": "primary (bad_quality false, seeing <= 4 arcsec)",
                    "robustness": "strict mask re-run reported for any exceedance",
                    "reference_check": "any ZTF reference built from >= 1 "
                                       "in-window epoch is flagged and the node "
                                       "masked (frozen before search)"},
        "alpha": {"fwer_alpha": 0.05, "ks_alpha": 0.01},
        "split": {
            "seed": SEED, "dev_fraction": DEV_FRACTION, "unit": "target",
            "A": {"strata": "band-searchability class (gr/g/r)",
                  "dev": dev_A, "confirmatory": conf_A},
            "B_wide": {"targets": wide_B, "dev_eligible": dev_eligible,
                       "dev": dev_B,
                       "confirmatory": sorted(set(wide_B) - set(dev_B))},
            "B_grazing": {"targets": grazing_targets, "dev": [],
                          "note": "no dev set (5 covered events); procedures "
                                  "locked from A + B-wide dev experience; "
                                  "fully confirmatory"},
        },
        "search_units": {
            "A": {"n_units": len(units_A), "n_single_window": len(singles_A),
                  "units": units_A, "single_window": singles_A},
            "B": {"n_units": len(units_B), "n_single_epoch": len(singles_B),
                  "units": units_B, "single_epoch": singles_B},
        },
        "completeness": {
            "injections": "per search unit, 100 per event-window cell on the "
                          "difference-image substrate, v2 recovery window [2,1]",
            "spectrum": {"g": "line_532nm", "r": "line_650nm"},
            "temporal_model": "chord profile of the unit's (b, v_perp, radius)",
        },
    }
    body = json.dumps(freeze, indent=1, sort_keys=True)
    freeze["freeze_content_hash"] = ("sha256:"
                                     + hashlib.sha256(body.encode()).hexdigest())
    out = D / "configs" / "threshold_freeze_v1.json"
    out.write_text(json.dumps(freeze, indent=1, sort_keys=True) + "\n")
    print(json.dumps({
        "A_units": len(units_A), "A_single_window": len(singles_A),
        "B_units": len(units_B), "B_single_epoch": len(singles_B),
        "expected_control_crossings": freeze["threshold_rule"]
            ["expected_control_crossings"],
        "dev_A": dev_A, "dev_B_wide": dev_B,
        "grazing_targets": grazing_targets,
        "hash": freeze["freeze_content_hash"]}, indent=2))


if __name__ == "__main__":
    main()
