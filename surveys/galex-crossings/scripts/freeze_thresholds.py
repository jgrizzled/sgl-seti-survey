"""GALEX crossings threshold freeze v1.0 (hypotheses v1.0 + v1.1).

Freezes the searched-unit/trial family, the numeric statistic
definitions, the control constructions and the exceedance rule, bound
to content hashes of the freeze documents and the coverage products.
Also runs the frozen section-6 nonlinearity cut (MCAT NUV magnitude of
the three channel-A stars — catalog annotation-class contact, the
rule's declared input, snapshot-disciplined). **No photon is touched
and no signal statistic is formed.**

Constructions (hypotheses D4/D7, made numeric here):

- unit = (target, channel-rung, band); statistic = max over the
  unit's eligible covered events (D1a: the gj-1276 B 2010 event is
  ineligible for S_rate — forced_dev, constraint-only lane).
- S_rate: in-aperture photon count / usable live seconds (r_ap = 8
  arcsec; B: max over the deduplicated z-grid aperture set).
- S_burst: over usable live seconds, for each boxcar width w in
  {0.05, 0.5, 5, 50} s (step w/2): z(w) = (C_max(w) - lam*w) /
  sqrt(max(lam*w, 0.5)) with lam = N_tot/t_live from the same
  series; S_burst = max over w.
- S_period: H-test (de Jager; m <= 20) over a geometric period grid
  from 0.02 s to t_span/3 at 5x-oversampled Fourier spacing
  (delta(1/P) = 1/(5*t_span)); S_period = max H. Gate: >= 10
  in-aperture photons, else constraint-only.
- Controls (8 per trial, threshold T = max over controls, exceedance
  S > max(T, 0); expected control crossings = trials/9):
  B: 8 pseudo-positions on the same visit at the locus boresight
  radius, angles locus + k*40 deg (k = 1..8), each carrying the full
  z-family aperture pattern; a position within 30 arcsec of an MCAT
  source brighter than NUV 21 or within 60 arcsec of the locus
  segment rotates +5 deg until valid (frozen resolution rule).
  A: 8 same-duration contiguous usable segments drawn from the
  star's off-window visits (outside every in-era 0.1 AU window of
  that target-channel), designated by the frozen selection rule
  (visits time-ordered; equally spaced non-overlapping segments;
  seed 20260826); < 8 valid segments -> constraint-only.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.adapters.mast_gphoton import GphotonClient
from sglsurvey.snapshots import SnapshotStore

SURVEY = REPO / "surveys" / "galex-crossings"
OUT = SURVEY / "results"
RUN = REPO / "runs" / "galex-crossings"

BOUND_FILES = ("hypotheses.md",
               "results/coverage_v1_events.ecsv",
               "results/coverage_v1_summary.json",
               "results/coverage_v1.md")

#: frozen section-6 rule
NONLIN_EXCLUDED_BRIGHTER_THAN = 14.5
NONLIN_MARGINAL_BRIGHTER_THAN = 15.5
MCAT_MATCH_ARCSEC = 15.0

APERTURE_ARCSEC = 8.0
BURST_WIDTHS_S = (0.05, 0.5, 5.0, 50.0)
PERIOD_MIN_S = 0.02
PERIOD_MAX_FRACTION = 1.0 / 3.0
PERIOD_OVERSAMPLE = 5
HTEST_MAX_HARMONICS = 20
MIN_PHOTONS_PERIOD = 10
N_CONTROLS = 8
CONTROL_SEED = 20260826
B_CONTROL_STEP_DEG = 40.0
B_CONTROL_ROTATE_DEG = 5.0
B_MCAT_EXCLUDE_ARCSEC = 30.0
B_MCAT_EXCLUDE_NUV_MAG = 21.0
B_LOCUS_EXCLUDE_ARCSEC = 60.0

#: the confirmatory event-units from coverage v1 (freeze section 9)
UNIT_EVENTS = {
    # (target, channel): [(t_ca date, nuv_s, fuv_s)]
    ("gj-1276", "B"): [("2007-03-03", 109, 109),
                       ("2010-03-03", 1637, 0)],
    ("gj-1276", "A"): [("2007-09-05", 92, 92)],
    ("wolf-359", "A"): [("2007-03-03", 97, 97)],
    ("ross-128", "A"): [("2007-03-17", 110, 110)],
}
#: D1a: rate-class contact at recon -> the 2010 event is ineligible
#: for S_rate (forced_dev constraint-only lane)
FORCED_DEV_RATE_EVENTS = {("gj-1276", "B", "2010-03-03")}

#: channel-A star positions at the unit event epoch (events table)
A_STARS = {
    "gj-1276": (343.475, -6.7825),
    "wolf-359": (164.1165, 7.012),
    "ross-128": (176.9356, 0.8033),
}


def sha(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def nonlinearity_cut(client, store):
    """Frozen section-6 rule on the MCAT NUV magnitude of each A star."""
    rows = []
    for tid, (ra, dec) in A_STARS.items():
        mcat = client.mcat_box(ra, dec, 0.01, store)
        best = None
        for objid, nuv, fuv, sra, sdec in mcat:
            if nuv is None or float(nuv) < -99:
                continue
            d = np.hypot((float(sra) - ra) * np.cos(np.radians(dec)),
                         float(sdec) - dec) * 3600.0
            if d <= MCAT_MATCH_ARCSEC and (best is None or d < best[0]):
                best = (d, float(nuv),
                        None if fuv is None or float(fuv) < -99
                        else float(fuv), str(objid))
        if best is None:
            cls = "ok"      # no MCAT source at the star: fainter than
            nuv = None      # any nonlinearity concern by construction
            note = (f"no MCAT NUV source within {MCAT_MATCH_ARCSEC}\" "
                    f"(high-PM M dwarf; per-visit positions verified "
                    f"at dev)")
            d = fuv = objid = None
        else:
            d, nuv, fuv, objid = best
            cls = ("excluded" if nuv < NONLIN_EXCLUDED_BRIGHTER_THAN
                   else "marginal" if nuv < NONLIN_MARGINAL_BRIGHTER_THAN
                   else "ok")
            note = ""
        rows.append({"target_id": tid, "ra_deg": ra, "dec_deg": dec,
                     "mcat_objid": objid, "match_arcsec": d,
                     "nuv_mag": nuv, "fuv_mag": fuv, "class": cls,
                     "note": note})
        print(f"[nonlin] {tid}: NUV {nuv} ({cls}) {note}", flush=True)
    return rows


def main():
    hashes = {f: sha(SURVEY / f) for f in BOUND_FILES}

    client = GphotonClient()
    store = SnapshotStore(RUN)
    nonlin = nonlinearity_cut(client, store)
    (OUT / "nonlinearity_cut_v1.json").write_text(
        json.dumps({"rule": {
            "excluded_brighter_than": NONLIN_EXCLUDED_BRIGHTER_THAN,
            "marginal_brighter_than": NONLIN_MARGINAL_BRIGHTER_THAN,
            "match_arcsec": MCAT_MATCH_ARCSEC},
            "stars": nonlin}, indent=2) + "\n")

    excluded = {r["target_id"] for r in nonlin
                if r["class"] == "excluded"}

    units = []
    n_trials = 0
    for (tid, ch), evs in sorted(UNIT_EVENTS.items()):
        for band, col in (("NUV", 1), ("FUV", 2)):
            live = [(d, (n, f)[col - 1]) for d, n, f in evs
                    if (n, f)[col - 1] > 0]
            if not live:
                continue
            if ch == "A" and tid in excluded:
                units.append({"target_id": tid, "channel": ch,
                              "band": band, "events": live,
                              "statistics": [],
                              "status": "nonlinearity-excluded",
                              "trials": 0})
                continue
            rate_events = [d for d, _ in live
                           if (tid, ch, d) not in FORCED_DEV_RATE_EVENTS]
            stats = []
            if rate_events:
                stats.append({"name": "S_rate", "events": rate_events})
            stats.append({"name": "S_burst",
                          "events": [d for d, _ in live]})
            stats.append({"name": "S_period",
                          "events": [d for d, _ in live],
                          "gate": f"min {MIN_PHOTONS_PERIOD} photons"})
            units.append({"target_id": tid, "channel": ch, "band": band,
                          "events": live, "statistics": stats,
                          "status": "confirmatory",
                          "trials": len(stats)})
            n_trials += len(stats)

    expected_crossings = n_trials / (N_CONTROLS + 1)

    config = {
        "survey": "galex-crossings",
        "version": "threshold_freeze_v1.0",
        "date": "2026-08-26",
        "bound_hashes": hashes,
        "statistics": {
            "aperture_arcsec": APERTURE_ARCSEC,
            "S_rate": "in-aperture photons / usable live seconds; "
                      "B: max over deduplicated z-grid apertures; "
                      "unit statistic = max over eligible events",
            "S_burst": {"widths_s": list(BURST_WIDTHS_S),
                        "step": "width/2",
                        "z": "(C_max - lam*w)/sqrt(max(lam*w, 0.5)), "
                             "lam = N_tot/t_live",
                        "statistic": "max over widths, then events"},
            "S_period": {"grid": "geometric, P_min 0.02 s to "
                                 "t_span/3, delta(1/P) = "
                                 "1/(5*t_span)",
                         "test": "H-test (de Jager), m <= 20",
                         "min_photons": MIN_PHOTONS_PERIOD,
                         "orbital_smear": "not corrected; carried by "
                                          "injections (freeze s4)"},
        },
        "controls": {
            "n": N_CONTROLS,
            "rule": "T = max over controls; exceedance S > max(T, 0)",
            "B": {"construction": "pseudo-positions, same visit, locus "
                                  "boresight radius, angles locus + "
                                  "k*40 deg",
                  "rotate_step_deg": B_CONTROL_ROTATE_DEG,
                  "mcat_exclusion": [B_MCAT_EXCLUDE_ARCSEC,
                                     B_MCAT_EXCLUDE_NUV_MAG],
                  "locus_exclusion_arcsec": B_LOCUS_EXCLUDE_ARCSEC,
                  "apertures": "full z-family pattern, translated"},
            "A": {"construction": "8 same-duration contiguous usable "
                                  "segments from the star's off-window "
                                  "visits; time-ordered, equally "
                                  "spaced, non-overlapping",
                  "off_window": "outside every in-era 0.1 AU window "
                                "of the target-channel",
                  "seed": CONTROL_SEED,
                  "insufficient": "constraint-only"},
        },
        "forced_dev_rate_events": sorted(
            list(e) for e in FORCED_DEV_RATE_EVENTS),
        "nonlinearity_cut": "results/nonlinearity_cut_v1.json",
        "units": units,
        "n_trials": n_trials,
        "expected_control_crossings": round(expected_crossings, 3),
        "split": "D8: dev = off-window pseudo-units only; "
                 "confirmatory = all listed units, blind",
    }
    path = SURVEY / "configs" / "threshold_freeze_v1.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config, indent=2) + "\n")
    print(f"\nfrozen -> {path}")
    print(f"config sha256: {hashlib.sha256(path.read_bytes()).hexdigest()}")
    print(f"units: {sum(1 for u in units if u['trials'])} searched "
          f"({len(units)} listed), trials: {n_trials}, expected "
          f"control crossings: {expected_crossings:.2f}")


if __name__ == "__main__":
    main()
