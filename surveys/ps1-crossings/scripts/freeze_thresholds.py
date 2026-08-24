"""PS1 crossings threshold freeze v1.0.

Declares, before any pixel is touched, the search statistics, control
constructions, threshold rule, dev/confirmatory split, and trials
accounting for both channels, binding them to content hashes of the
frozen inputs (hypotheses.md, coverage_v1_events.ecsv,
saturation_cut_v1.ecsv). Port of the ZTF crossings freeze with the
final ZTF constructions (v1.0+v1.1+v1.2 lessons) adopted at freeze
time rather than rediscovered: exceedance rule S > max(T, 0), empirical
variance rescale k = median(r^2/v)/0.4549 (floor 1) in both channels,
same-rung-only pseudo-window exclusion, offset-count constraint-only
gate. PS1 substitutions: warp-direct star-calibrated substrate (no
difference images), all five grizy bands, the A 1.0-AU rung
constraint-only by frozen declaration (window ~ observing season), the
TTI-pair motion veto, and a catalogued-static annotation rule.
Conventions inherited from the v2 engine: 8 designated controls,
WEIGHT_CAP 20x effective-epoch floor, split seed 20260824, dev
fraction 0.30, fwer_alpha 0.05 / ks_alpha 0.01.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
D = REPO / "surveys" / "ps1-crossings"
INPUTS = {
    "hypotheses": D / "hypotheses.md",
    "coverage": D / "results" / "coverage_v1_events.ecsv",
    "saturation": D / "results" / "saturation_cut_v1.ecsv",
}
SEED, DEV_FRACTION = 20260824, 0.30
BANDS = ("g", "r", "i", "z", "y")
A_SEARCH_RADIUS = 0.1     # the 1.0-AU rung is constraint-only by freeze


def sha(p: Path) -> str:
    return "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    cov = Table.read(INPUTS["coverage"])
    sat = Table.read(INPUTS["saturation"])
    searchable = {(str(r["target_id"]), b):
                  r[f"status_{b}"] in ("ok", "marginal")
                  for r in sat for b in BANDS}

    # -- channel A search units: (target, band) on the 0.1-AU rung with
    #    >=1 covered window and a searchable band; primary statistic
    #    needs >=2 windows. The 1.0-AU rung yields no units (frozen
    #    constraint-only: window ~ observing season, ZTF dev lesson 1).
    A = cov[cov["channel"] == "A"]
    A01 = A[np.asarray(A["radius_au"]) == A_SEARCH_RADIUS]
    units_A, singles_A = [], []
    for tid in sorted(set(str(t) for t in A01["target_id"])):
        sub = A01[np.asarray([str(x) == tid for x in A01["target_id"]])]
        for b in BANDS:
            if not searchable.get((tid, b), False):
                continue
            n_win = int((np.asarray(sub[f"n_{b}_pri"]) > 0).sum())
            if n_win >= 2:
                units_A.append({"target_id": tid, "band": b,
                                "radius_au": A_SEARCH_RADIUS,
                                "n_windows": n_win})
            elif n_win == 1:
                singles_A.append({"target_id": tid, "band": b,
                                  "radius_au": A_SEARCH_RADIUS})
    A10 = A[np.asarray(A["radius_au"]) == 1.0]
    a10_cov = int((np.asarray(A10["n_primary"]) > 0).sum())

    # -- channel B search units: (event, radius, band); the z grid is a
    #    nested family (epochs are shared across z), so the unit
    #    statistic is max over z and controls take the same max -- one
    #    trial per unit. All five bands are primary (hypotheses v1.0).
    B = cov[cov["channel"] == "B"]
    agg, pair_any = {}, {}
    for r in B:
        for b in BANDS:
            n = int(r[f"n_{b}_pri"])
            if n == 0:
                continue
            k = (str(r["event_id"]), str(r["target_id"]),
                 float(r["radius_au"]), b)
            agg[k] = max(agg.get(k, 0), n)
            pair_any[k] = max(pair_any.get(k, 0), int(r[f"pair_{b}"]))
    units_B, singles_B = [], []
    for (eid, tid, rad, b), n in sorted(agg.items()):
        u = {"event_id": eid, "target_id": tid, "radius_au": rad,
             "band": b, "n_epochs_max_z": n,
             "n_tti_pairs": pair_any[(eid, tid, rad, b)]}
        (units_B if n >= 2 else singles_B).append(u)

    # -- dev/confirmatory split (unit = target), stratified
    rng = np.random.default_rng(SEED)
    strata_A = {}
    for u in units_A + singles_A:
        tid = u["target_id"]
        strata_A[tid] = "".join(b for b in BANDS
                                if searchable.get((tid, b), False))
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
        "freeze_version": "ps1-crossings-thresholds-v1.0",
        "frozen_at": "2026-08-24",
        "hypothesis_version": "ps1-crossings-hypotheses-v1.0",
        "input_hashes": {k: sha(p) for k, p in INPUTS.items()},
        "statistics": {
            "substrate": "star-calibrated matched-filter forced photometry "
                         "directly on DR2 warp cutouts (per-frame ZP from "
                         "DR2 mean stars through the identical filter, "
                         "filter-median fallback; no difference images "
                         "exist). Exact warp mask (fatal template 16255) "
                         "defines usable pixels.",
            "variance": "per-epoch variances from the warp wt map, rescaled "
                        "by k = median(r^2/v)/0.4549 over clipped "
                        "off-window epochs (floor 1), per (unit, band) -- "
                        "the ZTF v1.2 / v1 confusion-floor rule, frozen "
                        "here for both channels",
            "A": {
                "primary": "S_A = weighted stack over covered windows of "
                           "the per-window mean excess vs the same star's "
                           "off-window baseline (robust mean, 3x3sigma "
                           "clip); one-sided positive; WEIGHT_CAP 20x "
                           "median epoch weight",
                "gates": "primary requires >= 2 covered windows AND >= 8 "
                         "primary off-window epochs in band (baseline "
                         "support, resolved at search time from the "
                         "snapshotted listings) AND >= 8 valid "
                         "pseudo-window offsets; any failure -> "
                         "constraint-only",
                "no_template": "no parallax-factor systematics template: "
                               "the star is measured directly (no "
                               "reference image), so the ZTF v1.1 PM-"
                               "dipole construction does not apply; the "
                               "variance rescale alone absorbs the "
                               "empirical scatter",
                "constraint_only_rung": {
                    "radius_au": 1.0,
                    "reason": "window ~115 d ~ observing season; crossing "
                              "geometry phase-locks to the sidereal year "
                              "and PS1 single-phase cadence concentrates "
                              "epochs in the same season -- temporal "
                              "controls fail by construction (ZTF dev "
                              "lesson 1, adopted at freeze)",
                    "n_covered_event_rows": a10_cov,
                    "reporting": "coverage + per-window depth only; no "
                                 "discovery statistic; zero trials",
                },
                "single_window_units": "constraint-only, never candidates",
                "annulus": "<= 10 arcsec station annulus searched with the "
                           "channel-B construction minus the z-track",
            },
            "B": {
                "statistic": "S = max over the nested z family of the "
                             "weighted in-window stack along the "
                             "(event, z) track; controls take the same "
                             "max; one-sided positive",
                "static_sky": "the static sky is in the photometry (no "
                              "differencing): track nodes within 2 arcsec "
                              "of a DR2 catalogued static source (>= 3 "
                              "detections, sglsurvey.vetting snapshot "
                              "rule) are annotated at search time; ring "
                              "controls cross static sources at the same "
                              "areal rate, so the threshold absorbs the "
                              "confusion floor",
                "single_epoch_units": "reportable 'single-epoch exceedance' "
                                      "class; promotion to candidate "
                                      "requires recurrence at a second "
                                      "covered window",
            },
        },
        "controls": {
            "A": {"kind": "temporal pseudo-windows", "n": 8,
                  "designated_offsets_days": [-97, -71, -47, -23,
                                              23, 47, 71, 97],
                  "redraw_offsets_days": [-127, -113, 113, 127],
                  "validity": "an offset is valid iff its pseudo-windows "
                              "contain >= 2 windows with >= 1 primary "
                              "epoch (mirror of the real-set gate); "
                              "overlap exclusion is same-rung-only (ZTF "
                              "v1.1); units lacking 8 valid offsets are "
                              "constraint-only"},
            "B": {"kind": "spatial ring trajectories", "n": 8,
                  "ring_radii_arcsec": [20.0, 30.0, 40.0],
                  "designated": "identical designated_indices construction "
                                "as ztf-v2/ps1-v2 (perpendicular offsets "
                                "follow the track)"},
        },
        "threshold_rule": {
            "T": "max one-sided control statistic over the 8 designated "
                 "controls",
            "exceedance": "S > max(T, 0); margin S - T reported (the R "
                          "ratio was retired by ZTF v1.1: unstable near "
                          "T -> 0+, sign-invalid for T <= 0)",
            "per_search_crossing_probability": "1/9",
            "expected_control_crossings": round(n_searches / 9.0, 1),
            "disposition": "every exceedance individually adjudicated by "
                           "the frozen veto ladder; no silent drops",
        },
        "veto_ladder": [
            "1. MPC/known-object census per in-window epoch (calibrated "
            "veto; antipodes of near-ecliptic stars sit in the opposition "
            "asteroid swarm)",
            "2. rate test: candidate must move at the predicted "
            "0.3-6.5 arcsec/day retrograde rate for its z",
            "3. TTI-pair motion test: same-night warp pairs (~15-40 min) "
            "-- an ordinary mover displaces >= several arcsec between "
            "pair members, the relay track < 0.2 arcsec (calibrated veto, "
            "new vs ZTF)",
            "4. flux-consistent catalogued-static test "
            "(sglsurvey.vetting, DR2 snapshot, >= 3 detections)",
            "5. recurrence at a second covered window on the recomputed "
            "track",
            "annotations (never sole grounds for rejection): chromatic "
            "consistency across non-simultaneous bands; single-season / "
            "phase-lock structure",
        ],
        "quality": {
            "coverage_mask": "badflag == 0 (uniformly true in-era; formal)",
            "search_mask": "exact warp mask, fatal template 16255; "
                           "~75% of nominal-footprint epochs expected to "
                           "survive (Pipeline-A measurement)",
            "strict": "PSF FWHM <= 2.5 arcsec from the warp header; "
                      "strict re-run reported for any exceedance",
            "saturation_gate": "before any confirmatory channel-A search, "
                               "the dev stage must verify the frozen "
                               "per-band exclusion levels against >= 1 "
                               "bright-star warp (CELL.SATURATION + "
                               "SAT/STARCORE mask bits); amendment "
                               "required if off by > 0.5 mag",
        },
        "alpha": {"fwer_alpha": 0.05, "ks_alpha": 0.01},
        "split": {
            "seed": SEED, "dev_fraction": DEV_FRACTION, "unit": "target",
            "A": {"strata": "searchable-band class (subset of grizy)",
                  "dev": dev_A, "confirmatory": conf_A,
                  "note": "if dev is empty or near-empty, A constructions "
                          "stand as frozen (imported from the validated "
                          "ZTF v1.0+1.1+1.2 chain); the saturation gate "
                          "still runs at dev"},
            "B_wide": {"targets": wide_B, "dev_eligible": dev_eligible,
                       "dev": dev_B,
                       "confirmatory": sorted(set(wide_B) - set(dev_B))},
            "B_grazing": {"targets": grazing_targets, "dev": [],
                          "note": "no dev set; grazing-rung windows are "
                                  "temporal subsets of wide-rung windows "
                                  "of the same targets; fully "
                                  "confirmatory with procedures locked "
                                  "from B-wide dev"},
        },
        "search_units": {
            "A": {"n_units": len(units_A),
                  "n_single_window": len(singles_A),
                  "units": units_A, "single_window": singles_A},
            "B": {"n_units": len(units_B),
                  "n_single_epoch": len(singles_B),
                  "units": units_B, "single_epoch": singles_B},
        },
        "completeness": {
            "injections": "per search unit, >= 100 per event-window cell, "
                          "stamp-response (Moffat beta=3 at the warp "
                          "header FWHM) on the warp substrate through the "
                          "identical star-calibrated filter, v2 recovery "
                          "window [2,1]",
            "spectrum": {"g": "line_532nm", "r": "line_617nm",
                         "i": "line_752nm", "z": "line_866nm",
                         "y": "line_962nm"},
            "temporal_model": "chord profile of the unit's "
                              "(b, v_perp, radius)",
        },
    }
    body = json.dumps(freeze, indent=1, sort_keys=True)
    freeze["freeze_content_hash"] = (
        "sha256:" + hashlib.sha256(body.encode()).hexdigest())
    out = D / "configs" / "threshold_freeze_v1.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(freeze, indent=1, sort_keys=True) + "\n")
    print(json.dumps({
        "A_units": len(units_A), "A_single_window": len(singles_A),
        "A_1au_constraint_only_rows": a10_cov,
        "B_units": len(units_B), "B_single_epoch": len(singles_B),
        "expected_control_crossings": freeze["threshold_rule"]
            ["expected_control_crossings"],
        "dev_A": dev_A, "dev_B_wide": dev_B,
        "grazing_targets": grazing_targets,
        "hash": freeze["freeze_content_hash"]}, indent=2))


if __name__ == "__main__":
    main()
