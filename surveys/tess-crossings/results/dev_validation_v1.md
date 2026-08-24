# TESS crossings dev-stage validation v1 — results

Run 2026-08-24 under threshold freeze **v1.1** (`dev_validation.py`,
`results/dev_validation_v1.json`). The dev stage was reshaped by two
pre-search findings, both made before any signal statistic existed:

## Finding T1 — the science-array usability screen (amendment v1.1)

TESScut serves collateral (non-science) CCD pixels with aperture = 1
and ~zero flux. Three of eight cutouts sat on or across the science
edge (columns 45–2092): **ross-128 B s42 (the v1.0 dev unit,
x = 6–36: entirely virtual)**, teegarden A s44 (x = 4–18), and
gj-1276 A s42 (star at x = 2094, off the trailing edge). The screen
(`usability_screen.py`, margin 4 px + 12 px ring extent for B) is now
a frozen input; the three rows are `off_science_array`
(nominal-covered / unusable — the PS1 nominal-vs-exact lesson in TESS
form, caught here at dev rather than at confirmatory). Re-freeze:
**6 B units (wolf-359 ×3, teegarden ×3 — all confirmatory), 12
trials, 1.3 expected control crossings; machinery validation
reassigned to the zero-trial A reference rows.**

## Machinery validation on the A reference rows (zero trials)

| row | FWHM | star flux → implied T (nominal ZP) | TIC T | k | S_c ref | S_p ref | >5σ cadences in-window |
|---|---|---|---|---|---|---|---|
| teegarden s71 (200 s, 10,024 cadences) | 2.4 px | 9,489 e-/s → 10.50 | 10.565 | 6.7 | +39.5 | 3.1 | 0 |
| van-maanen s43 (600 s, 2,794 cadences) | 3.0 px (fit ceiling) | 3,186 e-/s → 11.68 | 11.99 | 292 | +82.3 | 13.8 | 24 |

- **Per-cadence photometry, time handling, and flux scale validate**:
  teegarden's star recovers its TIC magnitude to **+0.07 mag**
  through the identical kernel chain at 200 s cadence.
- **Flux-scale gate**: the wolf-359 **B field passes cleanly** (10 TIC
  stars, ZP = 20.25, scatter **0.066 mag** ≤ 0.2) — the discovery
  fields calibrate star-by-star as frozen. The tiny A cutouts are
  calibrator-sparse (1–3 stars); the van-maanen field's 3-star check
  is junk-dominated (scatter 0.56) and its target-star
  self-calibration is off by −0.31 mag with the FWHM fit at its
  ceiling (likely blend/crowding at 21″ pixels): **van-maanen A
  reference depths carry a declared ±0.3 mag scale caveat**; A
  fields use target-star self-calibration, B fields use the field-star
  ZP.
- **Finding T2 — blended-A statistics are systematics-dominated**, as
  every prior survey found for its blended channel: sector-scale
  photometric drift gives k = 6.7–292 and reference S_c of +39/+82 on
  ordinary stars; the 24 ">5σ" van-maanen in-window cadences are
  drift/scattered-light residuals, not flares. This validates the
  freeze's structural call: A rows carry **zero trials** and are
  reported as light curves with reference depths only. For channel B
  the equivalent common-mode drift is absorbed differentially by the
  ring controls, which see the same field — the construction the
  freeze relies on.

## Verdict

Machinery cleared for the confirmatory run: per-cadence matched
filter, BTJD→UTC window bookkeeping, baseline/k/WEIGHT_CAP, both
frozen statistics, and the star-calibrated flux scale all exercised
on zero-trial rows; the usability screen is in place; no amendment to
the statistics or thresholds required. Next: confirmatory B (6 units,
teegarden + wolf-359, grazing rungs included), completeness, report.
