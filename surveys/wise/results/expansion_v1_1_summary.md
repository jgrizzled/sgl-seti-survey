---
title: "v1.1 expansion — Proxima Cen, Wolf 359, Ross 248, GJ 65 A/B through the full pipeline"
date: 2026-08-19
analysis_run: "run-20040d16f249 (supersedes run-88887a3ac8d9 for the calibration stage)"
---

# v1.1 expansion results

Five endpoints added under registry v1.1 (`sha256:27895a11…`) and
hypotheses v1.1 (endpoint set expanded; all physics parameters
unchanged), then run through the complete scripted chain — coarse
discovery → precise pass → catalog screening → recurrence → cutouts →
sample tensors → full injection-calibration rerun for all **24
endpoint-roles / 12 endpoints / 9 corridors**. Total archive-data
footprint now 4.8 GB (all regenerable).

## Curation (see registry_curation.md §v1.1 for evidence)

Proxima Cen / Wolf 359 / Ross 248: Gaia DR3 (Wolf 359 RV from Fouqué
et al. 2018). GJ 65 A/B: GRAVITY 2024 orbit + constructed mass-weighted
Gaia barycenter, validated against 2MASS (0.34" over 17.4 yr) and the
Gaia A→B relative position (0.013" closure through sglseti).

## New-corridor pipeline outcomes

| Corridor | Coarse rows | Usable W1 epochs (rx) | Screen matches (rx+tx) | Recurrence verdict |
| --- | --- | --- | --- | --- |
| proxima | 2,427 | 416 | 7,256 | chance/static only |
| wolf359 | 1,822 | 334 | 2,297 | chance/static only |
| ross248 | 2,417 | 446 | 6,555 | chance/static only |
| gj65 | 2,232×2 | 402 | 7,351 | chance/static only |

All cryo W3/W4 frames usable in every new corridor. Screening added
23,459 ScreenMatch records (67,374 total); the fixed-z recurrence
filter shows only chance-level and static-star-pattern peaks — no
track-consistent candidate.

## Calibration rerun (all 24 pairs, registry v1.1)

**768 Constraint records** (24 pair-roles × 4 bands × 8 z-intervals,
all `recovery_curve`). New-corridor median 90%-recovery depths
(duty ≥ 0.5, Vega mag): GJ 65 W1 16.3–17.1 and W2 15.3–15.9 (sparse
field — the deepest corridors yet, with W3 ~13.1); Wolf 359 W1
15.1–16.6; Ross 248 W1 15.7–15.9; Proxima W1 14.8 (dense field near
the Galactic-plane-adjacent α Cen antipode region).

**Threshold exceedances: 4 of 96 pair-bands**, all vetoed with
Candidate records:

- `sirius-a/rx/W2` S=871 vs T=451 — phase-imbalanced (941/144; a relay
  is phase-balanced), bright-star-artifact corridor, and a real source
  at this stacked S (~11–12 mag) is excluded by the catalog screening.
- `proxima-cen/tx/W2` S=250 vs T=214 — textbook single-phase static
  star (phase S = −6.7 / +284.8).
- `barnard-star/tx/W3` (8.8 vs 8.7) and `ross-154/tx/W3` (7.1 vs 7.0) —
  1% over an 8-draw null on single-visit cryo cells; placement
  variance.

(The v1.0 calibration run had checked only the two stack_v1-flagged
Lalande cells against thresholds; the full-exceedance census is a
procedure improvement introduced with this rerun, and the two Lalande
cells remain vetoed in it.)

**Search conclusion across all 12 endpoints: no surviving candidate.**

## Records census (cumulative, runs/wise/)

| Record | Count |
| --- | --- |
| QuerySnapshot | 187 |
| Observation | ~18,900 unique frame-bands |
| IntersectionEvaluation (coarse + precise) | ~76,000 |
| ScreenMatch | 67,374 |
| AnalysisRun (stack, calib v1.0, calib v1.1) | 3 |
| Constraint (current, v1.1) | 768 |
| Candidate (all vetoed) | 6 |
