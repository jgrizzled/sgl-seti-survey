# ZTF scale-up v1 — stage summaries (run 2026-08-20 → 21)

Hypotheses `ztf-hypotheses-v1.0` (unchanged from the pilot); registry
`pilot_wise_2026.yaml` v1.5; membership = ZTF overlay v1 (62 corridors
with antipode Dec > −28°, 69 endpoints, grid grading 29 ok / 18 edge /
15 gap). Driver `scripts/run_scaleup.sh`; wall time ≈ 23 h (one resume
after an uncaught IRSA 503 in screening); `runs/ztf/` = 234 GB,
regenerable.

| stage | run | numbers |
|---|---|---|
| coarse | `coarse_v1` | 147,063 public quadrant-exposures (r 77,635 / g 66,466 / i 2,962; 7,315 bad-quality), 327,390 evaluations, 174,876 hits |
| precise | `precise_v1` | 148,032 usable / 337 partial / 10,322 unusable / 16,185 unknown (product not served or transient 503 — 217 of the latter) |
| screening | `screen_v1` | 1,686,774 psfcat ScreenMatch records; fixed-z point filter: 130 track-consistent-or-mixed, 3 static-background-star, 4 indeterminate — none with a single follower |
| cutouts | `products/cut` | 67,003 sci/diff/msk sets (65,651 with difference images) |
| tensors | `calib_v1/tensors` | 138 endpoint×role tensors, 65,827 flux maps, 147,429 epoch-samples, median 449 epochs per pair-band, difference-image weight at peak median 0.99 |
| calibration | AnalysisRun `run-1dc29a973d0c` | 330 pair-band searches; **2,640 Constraint records** (2,574 recovery curves, 66 not-constrainable — all in gap-graded corridors); T percentiles 3.9 / 5.7 / 10.7; 263 of 330 pair-bands have both parallax phases at the peak cell |
| candidates | `calib_v1/records/candidate*.jsonl` | 24 exceedances → 14 phase-vetoed at calibration; 10 adjudicated (A–E): 2 fail held-out halves, all 10 fail the survey-wide trials test → **0 surviving candidates** |

## Depths (90 % recovery, duty ≥ 0.5, |µ| ≤ 1″/yr, control-corrected AB)

| band | recovery-curve cells | median m90 | 90th pct | best |
|---|---|---|---|---|
| g | 1,050 | 22.54 | 22.95 | 23.57 (GJ 11547) |
| r | 1,075 | 22.40 | 22.84 | 23.44 (GJ 11547) |
| i | 449 | 20.05 | 20.83 | 21.05 |

By grid grade: ok 22.55, edge 22.33, gap 21.28 (median m90); every
not-constrainable cell is in a gap corridor. Lowest-z interval
(550–624 AU) constrainable in 328/330 searches, highest (3,258–10,000
AU) in 315/330. Deepest corridors: GJ 11547, α Cen A/B, 61 Vir, GJ 318
(m90 ≈ 23.2–23.6). The pilot's Ross 128 (edge, 40 % science-image
regime) sits ≈ 1 mag shallower than ok-grade corridors with full
difference-image coverage — the reference-strip effect quantified.

## Look-elsewhere (`calib_v1/look_elsewhere.json`)

Leave-one-out control ratios S/T_loo over all 330 searches: control
exceedance rate 12.5 % (41 expected) vs real 7.3 % (24 observed); real
S/T percentiles (50/90/99/100) 0.84/0.98/1.16/1.35 vs controls
0.85/1.02/1.24/1.86. The real trajectories are, if anything, quieter
than the controls. Per-draw global p of the strongest real exceedance
(gj-229-a tx r, S/T 1.35) is 0.004 → 1.4 expected among 330 searches.

## Adjudication (`calib_v1/adjudication.md`)

Tests A (top-3 epochs > 50 % of stack), B (significant epochs cluster
on the sky while the track spans > 10″), C (odd/even and first/second-
half held-out stacks each ≥ 2σ), D (science-image regime default), E
(survey-wide trials: expected chance exceedances N_trials × p ≥ 1).
gj-876 tx r and gj65-a tx r fail C; all ten fail E (expected 5–35
chance exceedances of their strength). Note that high-z cells (z >
7,000 AU: gj-3112, gj-66-a, gj-783) have track spans of 2–5″, so test
B cannot separate a relay from a faint static residual there — the
difference image and the trials test carry the decision.

## Lessons beyond the pilot

1. With difference-image coverage the optical stack reaches m90 ≈ 22.5
   on ~450 epochs, ≈ 1.5 mag short of the Gaussian-noise limit — T is
   set by subtraction residuals, so a residual-aware variance (or
   ZOGY-style scorr images) is the next sensitivity step.
2. Gap-graded corridors are worth ~1.3 mag less and hold every
   not-constrainable cell; they are the SPHEREx hand-off set together
   with the 15 invisible corridors.
3. Eight-control thresholds at 330 searches guarantee ~40 chance
   exceedances; the leave-one-out look-elsewhere test must be part of
   every calibration run, not an afterthought.
4. IRSA IBE serves gzip bodies (masks ~13 KB) and throws sporadic 503s
   (0.3 %); every fetch in every stage needs the retry-then-record path.
