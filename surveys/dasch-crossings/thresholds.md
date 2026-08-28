# DASCH crossings threshold freeze v1.0

Frozen 2026-08-26 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json`
(`sha256:83dae35e…f7cb0a`), bound to content hashes of
`hypotheses.md` v1.0, `results/coverage_v1_windows.ecsv` and
`results/coverage_v1_summary.json`. No photometric data touched.
The ZTF/PS1/PTF crossings constructions are adopted at freeze time
and adapted to the catalogue-level substrate (D3); constants below
are frozen, per-unit numeric thresholds materialize from the 8
designated controls inside the search scripts (house rule: an
exceedance is S > max(T, 0), 1/9 expected control-crossing rate
per trial).

## Channel B — locus-track hit statistic

Per covered window, per in-window exposure: one `platephot` pull
(APASS refcat) at the event locus; the search region is the
z-family track segment {relay_apparent(z, t): z ∈ 550, 1000, 2500,
5500, 10000 AU} evaluated at the exposure midpoint, dilated by
**r_match = 10″** (sensitivity re-runs at 5″/15″ annotated; DASCH
solved-plate astrometry ~1–2″, faint-end and coarse-series scatter
larger; track motion within an exposure ≤ 0.2″). A **hit** = a
platephot row within the region that passes the fatal-AFLAGS
template (defect/wedge/mult-exp/rejected-blend/uncertain-date,
mask 186647040), has archive `reject_flag` = 0, and is either
**uncatalogued** (blank `ref_number`) or refcat-matched but
**≥ 1.0 mag brighter** than its refcat `stdmag`. Strict template
adds the blend-class + plate-quality bits; its completeness cost is
measured (§ completeness), never silent.

- **S_event** = max over the unit's covered windows of the window
  hit count. **S_stack** = total hits over all covered windows.
- **Controls:** the identical statistic at 8 ring loci offset
  (±60″, ±90″, ±120″ in RA; ±90″ in Dec) from the locus — same
  plates, same windows, same platephot subregions (free with the
  primary pull). T = the per-control max (S_event analogue) /
  total (S_stack analogue); threshold = max over the 8.

## Channel A — on-star window excess

- **Detected regime** (star brighter than the plate limits):
  S = robust z of the in-window `magcal_magdep` rows vs the
  unit's off-window baseline (median difference / scaled MAD /
  √n_in; brightening positive); baseline = all usable off-window
  lightcurve rows outside every covered window padded to 2× the
  half-width; gate ≥ 20 off-window epochs. Rows with nonzero
  archive `reject_flag`, `TOO_BRIGHT`, or `SATURATED` excluded
  (D6).
- **Limits-only regime** (star fainter than the p90 in-window
  `limiting_mag_local` — teegarden): the B hit-count construction
  at the star position, identical r_match and candidate classes.
- **Controls:** 8 temporal pseudo-windows at ±23/47/71/97 d
  (re-draws ±113/127 d), mirror gate, same-rung exclusion — the
  PTF A-construction verbatim (the ±5.8 d windows sit 365 d apart,
  so the offsets stay clear of adjacent annual events).

## Units and trials (from the frozen coverage record)

18 searched units × 2 statistics = **36 trials**: dev 6 (expected
control crossings 0.67), confirmatory 30 (expected 3.33). Dev =
van-maanen A 0.1 AU (`forced_dev`, D1), teegarden A 0.1 AU
(limits-only machinery), wolf-359 B 2.5 R☉ (B-chain machinery).
Confirmatory = the remaining 15, including all three B 1.2 R☉
units and both headline century families (van-maanen B, gj-1276 B)
blind. Ledger (coverage-without-statistic): teegarden B 1.2 R☉
(38 events, 0 covered), ross-128 B 1.2 R☉ (no era events).

## Veto ladder (exceedances; hypotheses §5 verbatim)

Calibrated: known-object census (Horizons/SkyBoT at the exceedance
epoch); trail-morphology test (60-min exposures trail ≳ 15″/hr
movers; fwhm/ellipticity from the row); recurrence on the
recomputed track at the unit's other covered windows. Annotations:
same-plate multi-exposure family, cutout defect inspection
(documented DASCH defect classes), chromatic/colorterm structure.
Survivors are `retained-ambiguous`, never dismissed.

## Completeness

Per covered window: platephot-vs-querycat field-star recovery in
the same subregion → detection-efficiency-vs-magnitude curve
(Wilson intervals); fatal- and strict-template costs measured on
the same population; per-window m90 from the curve (C1 analogue:
`limMag*` columns alone never qualify a constraint). Positive
control at dev: a documented high-amplitude variable through the
identical querycat → lightcurve → statistic chain (≤ 0.3 mag; RY
Cnc fallback).
