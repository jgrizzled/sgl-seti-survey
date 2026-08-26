"""PTF crossings threshold freeze v1.0.

Declares, before any pixel is touched, the search statistics, control
constructions, threshold rule, calibration gate, unit population and
trials accounting for both channels, binding them to content hashes
of the frozen inputs (hypotheses.md, coverage_v1_events.ecsv,
saturation_cut_v1.ecsv). Port of the PS1 crossings freeze with the
frozen hypotheses-v1.0 substitutions: scie-direct star-calibrated
substrate, g/R bands, unit = (target, radius, band) with the
S_event / S_stack statistic family (freeze D4 — S_event primary;
S_stack only where >= 2 covered events), the same-night repeat veto,
and the dev/confirmatory split frozen verbatim at D8 (no random
draw). Conventions inherited from the v2 engine: 8 designated
controls, WEIGHT_CAP 20x effective-epoch floor, fwer_alpha 0.05,
ks_alpha 0.01.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
D = REPO / "surveys" / "ptf-crossings"
INPUTS = {
    "hypotheses": D / "hypotheses.md",
    "coverage": D / "results" / "coverage_v1_events.ecsv",
    "saturation": D / "results" / "saturation_cut_v1.ecsv",
}
BANDS = ("g", "R")
A_SEARCH_RADIUS = 0.1     # the 1.0-AU rung is constraint-only by freeze

#: D8 split, frozen at the hypothesis freeze — recorded, not drawn.
DEV_UNITS = {
    ("wolf-359", "B", 0.1),
    ("gj-1276", "A", 0.1),
    ("ross-128", "A", 0.1),
}


def sha(p: Path) -> str:
    return "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    cov = Table.read(INPUTS["coverage"])
    sat = Table.read(INPUTS["saturation"])
    sat_status = {(str(r["target_id"]), b): str(r[f"status_{b}"])
                  for r in sat for b in BANDS}

    # -- unit construction: (target, radius, band); per event take the
    #    max epoch count over the nested z family (epochs are shared
    #    across z, so the statistic and controls take the same max —
    #    the z grid multiplies nothing).
    def build_units(ch):
        sub = cov[cov["channel"] == ch]
        ev = {}
        for r in sub:
            if ch == "A" and float(r["radius_au"]) != A_SEARCH_RADIUS:
                continue
            for b in BANDS:
                n = int(r[f"n_{b}"])
                if n == 0:
                    continue
                k = (str(r["target_id"]), float(r["radius_au"]), b)
                e = ev.setdefault(k, {})
                eid = str(r["event_id"])
                if n > e.get(eid, (0, 0))[0]:
                    e[eid] = (n, int(r[f"pair_{b}"]))
        units = []
        for (tid, rad, b), events in sorted(ev.items()):
            n_events = len(events)
            n_multi = sum(1 for n, _ in events.values() if n >= 2)
            stats = []
            if n_multi >= 1:
                stats.append("S_event")
            if n_events >= 2:
                stats.append("S_stack")
            cls = "searched" if stats else "single_epoch"
            if ch == "A":
                st = sat_status.get((tid, b), "no_photometry")
                if st == "excluded":
                    cls = "saturation_excluded"
                sat_cls = st
            else:
                sat_cls = None
            units.append({
                "target_id": tid, "channel": ch, "radius_au": rad,
                "band": b, "n_covered_events": n_events,
                "n_multi_epoch_events": n_multi,
                "epochs_per_event": {e: n for e, (n, _) in
                                     sorted(events.items())},
                "pairs_per_event": {e: p for e, (_, p) in
                                    sorted(events.items())},
                "statistics": stats if cls == "searched" else [],
                "n_trials": len(stats) if cls == "searched" else 0,
                "class": cls,
                "saturation_status": sat_cls,
            })
        return units

    units_B = build_units("B")
    units_A = build_units("A")

    # dev membership: D8 lists positions (target, channel, 0.1 AU rung);
    # a dev B target's grazing-rung units would also be dev, but none
    # exist (wolf-359 has no covered grazing events).
    for u in units_A + units_B:
        u["split"] = ("dev" if (u["target_id"],
                                u["channel"], 0.1) in DEV_UNITS
                      else "confirmatory")

    A10 = cov[(cov["channel"] == "A")
              & (np.asarray(cov["radius_au"]) == 1.0)]
    a10_cov = int((np.asarray(A10["n_total"]) > 0).sum())

    searched = [u for u in units_A + units_B if u["class"] == "searched"]
    n_trials = sum(u["n_trials"] for u in searched)
    trials_dev = sum(u["n_trials"] for u in searched
                     if u["split"] == "dev")

    freeze = {
        "freeze_version": "ptf-crossings-thresholds-v1.0",
        "frozen_at": "2026-08-26",
        "hypothesis_version": "ptf-crossings-hypotheses-v1.0",
        "input_hashes": {k: sha(p) for k, p in INPUTS.items()},
        "statistics": {
            "substrate": "star-calibrated forced PSF photometry directly "
                         "on level-1 science cutouts (no public "
                         "difference images); per-frame ZP from PS1 DR2 "
                         "mean stars through the frozen band transform "
                         "(g: PS1 g direct; R: Jordi 2006 "
                         "R = r - 0.153(r-i) - 0.117); exact dmask "
                         "(fatal template 65533) defines usable pixels",
            "calibration_gate": {
                "min_calibrators": 5,
                "max_zp_scatter_mag": 0.2,
                "rule": "per-frame field-star ZP from >= 5 unsaturated "
                        "PS1 DR2 calibrators through the identical-band "
                        "transform, robust scatter <= 0.2 mag; a frame "
                        "failing the gate is unusable (never "
                        "fallback-calibrated) — freeze D3; header MAGZPT "
                        "recorded as a cross-check diagnostic only",
            },
            "variance": "per-epoch variances from the PSF fit + local "
                        "background, rescaled by k = "
                        "median(r^2/v)/0.4549 over clipped off-window "
                        "epochs (floor 1), per (unit, band) — the "
                        "ZTF v1.2 / v1 confusion-floor rule",
            "unit": "(target, radius, band); channel B epochs per event "
                    "are the max over the nested z family (z multiplies "
                    "nothing; controls take the same max)",
            "S_event": "primary: max over the unit's multi-epoch covered "
                       "events of the weighted in-window mean excess "
                       "(chord-weighted, one-sided positive) vs the "
                       "off-window baseline at the same position/track; "
                       "requires >= 1 covered event with >= 2 in-window "
                       "epochs",
            "S_stack": "secondary: recurrence stack — weighted stack of "
                       "the per-event in-window excesses over all "
                       "covered events of the unit; applies only where "
                       ">= 2 covered events (freeze D4); WEIGHT_CAP 20x "
                       "median epoch weight",
            "single_epoch_units": "units whose covered events all have "
                                  "1 in-window epoch: reportable "
                                  "'single-epoch exceedance' class, "
                                  "zero trials; promotion requires "
                                  "recurrence at a second covered "
                                  "window",
            "A_gates": "channel-A S_event additionally requires >= 8 "
                       "usable off-window epochs in band and >= 8 valid "
                       "pseudo-window offsets; any failure -> "
                       "constraint-only",
            "A_constraint_only_rung": {
                "radius_au": 1.0,
                "n_covered_event_rows": a10_cov,
                "reason": "window ~ observing season under the "
                          "sidereal-year phase-lock (standing theorem, "
                          "frozen at hypotheses §2); coverage + depth "
                          "only, zero trials",
            },
            "A_saturation_excluded": "excluded-class (target, band) "
                                     "units carry no statistic and no "
                                     "contrast-limit claim; reported as "
                                     "excluded in the ledger",
        },
        "controls": {
            "A": {"kind": "temporal pseudo-windows", "n": 8,
                  "designated_offsets_days": [-97, -71, -47, -23,
                                              23, 47, 71, 97],
                  "redraw_offsets_days": [-127, -113, 113, 127],
                  "validity": "an offset is valid iff its pseudo-window "
                              "set mirrors the real gate (>= 1 "
                              "pseudo-window with >= 2 usable epochs "
                              "for S_event); same-rung-only overlap "
                              "exclusion (ZTF v1.1); < 8 valid offsets "
                              "-> constraint-only (freeze D7: assessed "
                              "per unit, not assumed)"},
            "B": {"kind": "spatial ring trajectories", "n": 8,
                  "ring_radii_arcsec": [20.0, 30.0, 40.0],
                  "designated": "identical designated_indices "
                                "construction as ztf-v2/ps1-v2 "
                                "(perpendicular offsets follow the "
                                "track)"},
        },
        "threshold_rule": {
            "T": "max one-sided control statistic over the 8 designated "
                 "controls, per statistic",
            "exceedance": "S > max(T, 0); margin S - T reported",
            "per_trial_crossing_probability": "1/9",
            "n_trials": n_trials,
            "n_trials_dev": trials_dev,
            "n_trials_confirmatory": n_trials - trials_dev,
            "expected_control_crossings": round(n_trials / 9.0, 2),
            "disposition": "every exceedance individually adjudicated "
                           "by the frozen veto ladder; no silent drops",
        },
        "veto_ladder": [
            "1. SkyBoT known-object census per in-window epoch "
            "(calibrated veto; antipodes sit in the opposition asteroid "
            "stream)",
            "2. rate test: candidate must move at the predicted "
            "0.36-6.5 arcsec/day retrograde rate for its z",
            "3. same-night repeat test: where same-night exposure pairs "
            "exist (11 of 14 covered narrow-rung events), an ordinary "
            "mover displaces arcsec-scale between pair members while "
            "the relay track moves < 0.2 arcsec; where no pair exists "
            "the ladder lacks this rung, recorded per exceedance "
            "(hypotheses §5)",
            "4. flux-consistent catalogued-static test (ptf_objects "
            "snapshot, >= 3 detections, 2 arcsec annotation radius)",
            "5. recurrence at a second covered window on the recomputed "
            "track",
            "annotations (never sole grounds): cross-epoch chromatic "
            "consistency (bands not simultaneous); campaign-cadence "
            "clustering structure",
        ],
        "quality": {
            "coverage_mask": "none (freeze D5): imgtype='object', "
                             "fid<=2 only",
            "search_mask": "exact dmask, fatal template 65533 (Laher "
                           "2014 Table 15, all bits except 2^1 object)",
            "strict": "seeing <= 2.5 arcsec (metadata) + fatal template "
                      "+ calibration gate at scatter <= 0.1 mag; strict "
                      "re-run reported for any exceedance",
            "saturation_gate": "before any confirmatory search, the dev "
                               "stage verifies the frozen exclusion "
                               "levels (E_R 14.5 / E_g 15.0) against "
                               ">= 1 bright-star frame (dmask bit-8 "
                               "extent vs magnitude); amendment "
                               "required if off by > 0.5 mag",
        },
        "alpha": {"fwer_alpha": 0.05, "ks_alpha": 0.01},
        "split": {
            "rule": "frozen verbatim at hypothesis freeze D8 — no "
                    "random draw",
            "dev_positions": sorted(f"{t} {c} {r}" for t, c, r
                                    in DEV_UNITS),
            "note": "dev exercises the machinery (wolf-359 B: "
                    "single-epoch class; ross-128 A: "
                    "saturation-excluded class; gj-1276 A: the one "
                    "searchable A unit) — confirmatory keeps both "
                    "headline B families blind",
        },
        "search_units": {
            "A": units_A,
            "B": units_B,
        },
        "completeness": {
            "injections": "per searched unit, >= 100 per event-window "
                          "cell, stamp-response (Moffat beta=3 at the "
                          "header SEEING) on the scie substrate through "
                          "the identical star-calibrated chain (C1 "
                          "rule: no depth statement without the "
                          "injection-measured response), v2 recovery "
                          "window [2,1]",
            "spectrum": {"g": "line_532nm", "R": "line_658nm "
                                                 "(generic leakage)"},
            "temporal_model": "chord profile of the unit's "
                              "(b, v_perp, radius) per event",
            "positive_control": "one numbered main-belt asteroid "
                                "crossing a scoped antipode field "
                                "(SkyBoT-selected at dev) through the "
                                "identical chain, recovered to "
                                "<= 0.1 mag against predicted "
                                "magnitudes",
        },
    }
    body = json.dumps(freeze, indent=1, sort_keys=True)
    freeze["freeze_content_hash"] = (
        "sha256:" + hashlib.sha256(body.encode()).hexdigest())
    out = D / "configs" / "threshold_freeze_v1.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(freeze, indent=1, sort_keys=True) + "\n")

    print(json.dumps({
        "searched_units": [
            f"{u['target_id']} {u['channel']} {u['radius_au']} "
            f"{u['band']}: {u['n_covered_events']} ev "
            f"({u['n_multi_epoch_events']} multi) "
            f"{'+'.join(u['statistics'])} [{u['split']}]"
            for u in searched],
        "single_epoch": [
            f"{u['target_id']} {u['channel']} {u['radius_au']} "
            f"{u['band']} [{u['split']}]"
            for u in units_A + units_B if u["class"] == "single_epoch"],
        "saturation_excluded": [
            f"{u['target_id']} {u['channel']} {u['band']}"
            for u in units_A if u["class"] == "saturation_excluded"],
        "n_trials": n_trials, "n_trials_dev": trials_dev,
        "expected_control_crossings": round(n_trials / 9.0, 2),
        "hash": freeze["freeze_content_hash"]}, indent=2))


if __name__ == "__main__":
    main()
