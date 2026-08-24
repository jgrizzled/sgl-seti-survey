"""WISE crossings threshold freeze v1.0.

Single-channel freeze (A, 1.0 AU rung — everything else is
elongation-null per hypotheses v1.0). Declares, before any pixel is
touched, the statistic, control construction, threshold rule,
dev/confirmatory split and trials accounting, bound to content hashes
of the frozen inputs (hypotheses.md, elongation_gate_v1.json,
coverage_v1_events.ecsv, epochs_v1.json, saturation_cut_v1.ecsv).
The final ZTF/PS1 constructions are adopted at freeze: S > max(T, 0),
empirical variance rescale (floor 1), same-rung pseudo-window
exclusion with epoch-supported validity, WEIGHT_CAP 20x, off-window
baseline clip <= 120 epochs uniform-in-time. Pseudo-window support is
evaluated here offline from the coverage epoch lists — WISE's
6-month fixed-elongation cadence means offsets +/-23..97 d usually
move a window off the visit comb, so the searchable-unit count is
decided by data, not assumed.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
D = REPO / "surveys" / "wise-crossings"
INPUTS = {
    "hypotheses": D / "hypotheses.md",
    "elongation_gate": D / "configs" / "elongation_gate_v1.json",
    "coverage": D / "results" / "coverage_v1_events.ecsv",
    "epochs": D / "results" / "epochs_v1.json",
    "saturation": D / "results" / "saturation_cut_v1.ecsv",
}
SEED, DEV_FRACTION = 20260825, 0.30
BANDS = ("W1", "W2", "W3", "W4")
# Amendment v1.1 (2026-08-24, pre-search): the ZTF/PS1 offset pool
# (+/-23..97, redraws +/-113/127) was designed for 0.6-11.6 d windows
# and cannot clear this rung's ~116 d annual windows by construction
# (an offset must exceed the window length yet stay inside the ~249 d
# inter-window gap). The pool is rescaled to the rung's window length:
# a frozen deterministic scan +/-130, +/-135, ... +/-245 in 5 d steps,
# first 8 offsets satisfying the unchanged validity rule. This
# generalizes the redraw mechanism; the validity rule itself (no
# overlap with any same-rung window + >= 2 epoch-supported windows) is
# untouched. v1.0 (original pool) yielded 0 searchable units and is
# superseded before any pixel was touched.
OFFSET_SCAN = [s * o for o in range(130, 246, 5) for s in (1, -1)]
MIN_OFF_EPOCHS = 8
BASELINE_CLIP = 120


def sha(p: Path) -> str:
    return "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()


def windows_of(cov, tid, band):
    sub = cov[np.asarray([str(x) == tid for x in cov["target_id"]])]
    return [(float(r["t_ca_mjd"] - r["window_days"] / 2),
             float(r["t_ca_mjd"] + r["window_days"] / 2),
             int(r[f"n_{band}"])) for r in sub]


def valid_offsets(wins, ts):
    """Frozen rule: an offset is valid iff no shifted window overlaps
    any real same-rung window AND the shifted set contains >= 2 windows
    with >= 1 primary epoch."""
    ts = np.asarray(ts)
    chosen = []
    for off in OFFSET_SCAN:
        ok = True
        for lo, hi, _ in wins:
            slo, shi = lo + off, hi + off
            for lo2, hi2, _ in wins:
                if slo < hi2 and shi > lo2:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            n_win = sum(1 for lo, hi, _ in wins
                        if ((ts >= lo + off) & (ts <= hi + off)).any())
            ok = n_win >= 2
        if ok:
            chosen.append(off)
        if len(chosen) == 8:
            break
    return chosen


def main():
    cov = Table.read(INPUTS["coverage"])
    cov = cov[np.asarray([s == "queried" for s in cov["status"]])]
    sat = Table.read(INPUTS["saturation"])
    epochs = json.loads(INPUTS["epochs"].read_text())
    searchable = {(str(r["target_id"]), b):
                  r[f"status_{b}"] in ("ok", "marginal")
                  for r in sat for b in BANDS}

    units, singles, blocked = [], [], []
    for tid in sorted(set(str(t) for t in cov["target_id"])):
        for b in BANDS:
            if not searchable.get((tid, b), False):
                continue
            wins = windows_of(cov, tid, b)
            n_win = sum(1 for _, _, n in wins if n > 0)
            if n_win == 0:
                continue
            ts = epochs.get(tid, {}).get(b, [])
            in_any = np.zeros(len(ts), bool)
            ts_arr = np.asarray(ts)
            for lo, hi, _ in wins:
                in_any |= (ts_arr >= lo) & (ts_arr <= hi)
            n_off = int((~in_any).sum())
            offs = valid_offsets(wins, ts)
            u = {"target_id": tid, "band": b, "n_windows": n_win,
                 "n_off_epochs": n_off, "n_valid_offsets": len(offs),
                 "valid_offsets": offs}
            if n_win >= 2 and n_off >= MIN_OFF_EPOCHS and len(offs) >= 8:
                units.append(u)
            elif n_win == 1:
                singles.append(u)
            else:
                blocked.append(u)

    rng = np.random.default_rng(SEED)
    strata = {}
    for u in units:
        strata.setdefault(u["target_id"], set()).add(u["band"])
    strata = {t: "".join(sorted(bs)) for t, bs in strata.items()}
    dev = []
    for cls in sorted(set(strata.values())):
        ts = sorted(t for t, c in strata.items() if c == cls)
        n_dev = max(1, round(DEV_FRACTION * len(ts))) if len(ts) > 1 else 0
        dev += sorted(rng.choice(ts, size=n_dev, replace=False).tolist())
    conf = sorted(set(strata) - set(dev))

    n_searches = len(units)
    freeze = {
        "freeze_version": "wise-crossings-thresholds-v1.1",
        "frozen_at": "2026-08-24",
        "hypothesis_version": "wise-crossings-hypotheses-v1.0",
        "input_hashes": {k: sha(p) for k, p in INPUTS.items()},
        "scope": "channel A, 1.0 AU rung only (elongation gate); "
                 "channels B and A-0.1 are geometrically null and carry "
                 "zero trials",
        "statistics": {
            "substrate": "matched-filter forced photometry on L1b int "
                         "cutouts (v2 total-flux kernel, exact fatal-mask "
                         "template) at the per-epoch propagated star "
                         "position; per-frame magzp flux scale; fluxes "
                         "rescaled to a common Vega ZP",
            "primary": "S_A = weighted stack over covered windows of the "
                       "per-window mean excess vs the same star's "
                       "off-window baseline (robust mean, 3x3sigma clip, "
                       "deterministic <= 120-epoch uniform-in-time "
                       "subsample); one-sided positive; WEIGHT_CAP 20x",
            "variance": "per-epoch variances rescaled by "
                        "k = median(r^2/v)/0.4549 over the clipped "
                        "off-window epochs (floor 1)",
            "gates": ">= 2 covered windows AND >= 8 primary off-window "
                     "epochs AND >= 8 valid pseudo-window offsets "
                     "(evaluated above from the frozen coverage epoch "
                     "lists); any failure -> constraint-only",
            "single_window_units": "constraint-only, never candidates",
        },
        "controls": {"kind": "temporal pseudo-windows", "n": 8,
                     "offset_scan_days": OFFSET_SCAN,
                     "amendment_v1_1": "offset pool rescaled to the "
                                       "rung's ~116 d windows (see "
                                       "script header); v1.0 pool "
                                       "yielded 0 units pre-pixel",
                     "validity": "same-rung-only overlap exclusion "
                                 "(trivially satisfied: one rung); an "
                                 "offset is valid iff its shifted "
                                 "windows contain >= 2 windows with "
                                 ">= 1 primary epoch"},
        "threshold_rule": {
            "T": "max one-sided control statistic over the 8 designated "
                 "offsets",
            "exceedance": "S > max(T, 0); margin reported",
            "per_search_crossing_probability": "1/9",
            "expected_control_crossings": round(n_searches / 9.0, 1),
            "disposition": "every exceedance individually adjudicated; "
                           "no silent drops",
        },
        "veto_ladder": [
            "1. SkyBoT/known-object census per in-window visit "
            "(formality: blended-star channel)",
            "2. rate test: within a 1-10 d visit an ordinary mover "
            "displaces arcmin-scale; the star (and any station within "
            "the 6-arcsec-PSF blend) is static",
            "3. W1:W2 same-visit color test (calibrated where both "
            "bands are searchable): a laser line is single-band, "
            "stellar variability is broadband-correlated",
            "4. recurrence at a second covered window",
            "annotations: phase-lock structure; single-visit windows",
        ],
        "quality": {"primary": "qual_frame > 0, qa_status A, "
                               "moon_sep >= 15 deg + exact fatal-mask "
                               "template at the star",
                    "strict": "adds saa_sep >= 5 deg, qual_frame >= 10; "
                              "re-run reported for any exceedance"},
        "alpha": {"fwer_alpha": 0.05, "ks_alpha": 0.01},
        "split": {"seed": SEED, "dev_fraction": DEV_FRACTION,
                  "unit": "target",
                  "strata": "searchable-band class",
                  "dev": dev, "confirmatory": conf},
        "search_units": {"n_units": len(units),
                         "n_single_window": len(singles),
                         "n_gate_blocked": len(blocked),
                         "units": units, "single_window": singles,
                         "gate_blocked": blocked},
        "completeness": {
            "injections": ">= 100 per event-window cell, image-level "
                          "via the WISE PRFGrid stamp response through "
                          "the identical kernel, chord temporal profile, "
                          "line spectra at band effective wavelengths, "
                          "v2 recovery window [2,1]",
        },
    }
    body = json.dumps(freeze, indent=1, sort_keys=True)
    freeze["freeze_content_hash"] = (
        "sha256:" + hashlib.sha256(body.encode()).hexdigest())
    out = D / "configs" / "threshold_freeze_v1.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(freeze, indent=1, sort_keys=True) + "\n")
    print(json.dumps({
        "units": len(units), "single_window": len(singles),
        "gate_blocked": len(blocked),
        "unit_list": [f"{u['target_id']}/{u['band']}" for u in units],
        "expected_control_crossings":
            freeze["threshold_rule"]["expected_control_crossings"],
        "dev": dev, "confirmatory": conf,
        "hash": freeze["freeze_content_hash"]}, indent=2))


if __name__ == "__main__":
    main()
