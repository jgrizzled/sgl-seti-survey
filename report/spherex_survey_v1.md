---
title: "SPHEREx survey v1 — all 77 universal-list corridors"
date: 2026-08-20
status: "77 corridors / 88 endpoints; 9,808 constraints, 0 candidates"
---

# SPHEREx survey v1

Extends the 3-star pilot (`spherex_pilot_v1.md`) first to the 15
southern corridors ZTF cannot reach (2026-08-20,
`surveys/spherex/results/scaleup_v1_summary.md`, where the estimator
was developed to v3), then to the whole universal list (2026-08-21,
`surveys/spherex/results/allsky_v1_summary.md`). The numbers below are
the all-sky run (`calib_v4`, AnalysisRun `run-ddd97a27faad`); the
southern-only bullets that follow them are retained as the record of
how the rule set was reached.

## All-sky result

- **Coverage.** 77 corridors / 88 endpoints; 60,426 Level-2 detector-
  exposures (2025-05 → 2026-08), 74,762 usable/partial evaluations,
  33,879 slim cutouts from S3 with none missing.
- **Controls.** Flux scale within ±0.19 mag in every detector on
  208,541 star measurements. 1,226 searches × 16 controls: 64
  exceedances (5.2 %, under the 1/16 FAR), all vetoed (56 phase
  split, 4 role coincidence, 3 cross-detector complement, 1 template
  coverage). Layer-1: nothing follows a parallax track. **0
  Candidates.**
- **Constraints.** 9,808 records. Joint six-detector m90 (duty ≥ 0.5,
  |µ| ≤ 1″/yr): median 20.75 AB, 21.5–21.8 on the cleanest corridors
  (GJ 625, 82 Eri, GJ 1061, GJ 66, GJ 338), 16.9–17.9 on the
  Galactic-plane ones (α Cen A/B, GJ 11068). On a median corridor at
  550–720 AU: grey reflectors ≳ 1.3 × 10⁵ km (albedo 0.1), 700 K
  emitters ≳ 17 km, 1000 K emitters ≳ 5 km.
- **Next.** Re-run as each quick release adds a parallax phase (the
  single-season and season-weighted corridors gain a working phase
  test); candidate SED discriminator only if something survives.

## Southern scale-up (record of the rule set)

- **Coverage.** 24,099 public Level-2 detector-exposures (2025-05 →
  2026-08), 31,962 usable/partial precise evaluations, 14,579 slim
  cutouts from S3 with no missing products; every corridor covers
  550–10,000 AU at every visit.
- **Controls.** Flux scale on 54,265 catalogued-star measurements:
  within ±0.19 mag in every detector (0.054 mag offset applied).
  266 searches (228 per-detector + 38 joint) × 8 offset controls: 37
  exceedances (14 %, the empirical FAR of the rule), all vetoed —
  31 by the parallax-phase split, 3 single-season, 1 Rx/Tx role
  coincidence, 2 by the template-coverage rule (GJ 13157 D2 at
  z = 10,000 AU, 10–19″ from a W1 = 7.7 star in a b = −3° field, with
  the static template defined for only half the epochs there).
  **0 Candidates.**
- **Constraints (v2 estimator, `calib_v3`).** 2,128 records: 1,824
  per-detector and 304 joint six-detector (flat-F_ν) stacks. Per-
  detector m90 (duty ≥ 0.5, |µ| ≤ 1″/yr) 18.2–21.2 AB; joint m90
  19.1–21.9 AB — 21.2–21.9 on the clean high-latitude corridors
  (Lalande 21185, σ Dra, GJ 625, GJ 338, Struve 2398, GJ 687),
  19.5–19.8 on the Galactic-plane ones (HD 219134, GJ 13157, 61 Cyg).
  Physically, at 550–720 AU on the clean corridors: grey reflectors
  ≳ 1.0 × 10⁵ km (albedo 0.1), 700 K emitters ≳ 14 km, 1000 K
  emitters ≳ 4 km; ×2–3 larger on the crowded ones.
- **What the v2 estimator bought.** Calibrating epoch variances to the
  template's residual scatter (instead of the confusion-inflated
  per-cutout MAD) gained a median 0.44 mag per detector (up to 1.7 in
  the deep field); the joint stack a further 0.2–0.9 mag. The v1
  record (`calib_v2`: 1,824 constraints, 28 vetoed exceedances) is
  retained.
- **v3 estimator (`calib_v4`).** Sixteen offset controls (FAR < 1/16:
  21 exceedances of 266, 8 %) at zero median depth cost, and a cross-
  detector season-complement phase test for the deep fields where the
  per-detector split is degenerate. All 21 vetoed (17 phase, 2 cross-
  detector, 1 template coverage, 1 role coincidence). **0 Candidates.**
  The 2,128 constraints were re-issued under these thresholds; this is
  the frozen rule set for the northern corridors.
- **What limits depth now.** Confusion, not epoch count: the
  Galactic-plane corridors remain ~1.5 mag shallower, and the deep-
  field GJ 687 corridor is only marginally deeper than clean 60-epoch
  fields.
- **Next.** Re-run when the next quick release adds the missing
  parallax phase (Struve 2398 is single-season; GJ 1221/GJ 251/GJ 687
  are season-weighted); then the 61 remaining universal-list corridors
  shared with ZTF/WISE; remaining estimator items (cross-detector
  phase test, 16 controls, candidate SED discriminator).

AnalysisRuns `run-4f8c9fd41ce9` (v1 estimator, `calib_v2`),
`run-4fbed2e52668` (v2, `calib_v3`) and the `calib_v4` run (v3; re-issued
when the northern corridors complete); records under `runs/spherex/`.
