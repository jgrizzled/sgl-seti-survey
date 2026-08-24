---
title: "DECam southern survey (frozen decision rule, blind confirmatory run)"
date: 2026-08-24
status: "complete — confirmatory 14 endpoints / 86 cells: 0 candidates at family-wise α = 0.05; development 5 / 30: 0 candidates; 1,904 injection-calibrated constraints; first survey of the 15 southern corridors unreachable by PS1/ZTF, and the first built on the v2 design with no v1 exploratory phase"
---

# DECam southern survey

**Pins:** registry v1.5 · hypotheses `surveys/decam/hypotheses.md` v1.0
(freeze `sha256:c3c7ece419c20067…`, seed 20260824) · T0 = 58500 ·
engine `sglsurvey/` · adapter `sglsurvey/adapters/noirlab_decam.py`
· AnalysisRuns `run-7565b9349301` (development) / `run-3d3cf9840678`
(confirmatory, blind) in `runs/decam/v2/records/` · numbers from
`surveys/decam/results/report_tables.md` · recon
`surveys/decam/notes/decam_recon_2026-08-24.md`.

## Summary

DECam/NOIRLab instcal exposures (2012–2026, CTIO Blanco, grizY) are
the only 1″-class multi-epoch optical archive covering the 15
universal-list corridors below δ = −30° — corridors whose optical cell
was previously closed by nothing, and which include the
engineering-backbone picks σ Dra, HD 219134 and Lalande 21185. This
survey ran Pipeline A end-to-end (16,283 Observations, 3,061 precise
evaluations, 93,362 ScreenMatch records, 1,223 usable exposure product
sets) and the v2 statistical chain under a rule frozen before any
analysis touched the data — the first survey of the programme with no
v1 exploratory phase.

**No cell of the blind confirmatory set (86 cells over 11 corridors /
14 endpoints) reaches the family-wise threshold** R̃_FWER = 1.564;
11 cells have R > 1 against a ring expectation of 11.9. Development
set (3 pilot corridors forced + struve2398 drawn): 30 cells, R̃_FWER
= 1.247, 4 exceedances vs 4.0 expected, 0 candidates; quality-mask
sensitivity (primary/strict/loose) changes nothing. Median
persistent-source m90 on the confirmatory set is **g 23.2 / r 22.1 /
i 22.7 / z 22.1 / Y 21.0 AB** (threshold completeness; see table for
the four temporal models) — the deepest optical corridor constraints
of the programme where cells are constrainable, ~0.7 mag past the ZTF
survey medians. Coverage is genuinely heterogeneous PI data: 378 of
712 confirmatory threshold constraints are `not_constrainable`
(single-phase, sub-floor epoch counts, or bright limits), and Y is
void in the development set entirely.

## Rule (frozen)

As WISE v2.1 with DECam bindings: 360 × 5 × 5 grid (uniform in 1/z,
~1.0″ spacing; L∞ µ ≤ 1″/yr; T0 = 58500), single-epoch clip
|S_e| ≤ 5, per-frame weight cap, ring 20/30/40″ (48 offsets), R̃ =
R/q95, survey-wide max-R̃ null → FWER α = 0.05; veto = flux-consistent
NSC DR2 catalogued static (≥ 3 detections). Frozen exposure selection:
EXPTIME ≥ 30 s, grizY, obs_type = object, deepest exposure per night
per band (the gj-1221 LMC cap). Fatal dqmask codes: all but
INTERPOLATED; strict/loose masks reported as a systematic. Flux scale:
per-frame NSC-star calibration through the identical matched filter
(median star-ZP MAD 0.045 mag; header MAGZERO measured unreliable —
outliers to +3.3 mag — and never used). Epoch hold-out annotation at
MJD 60400 (archive still accumulating). Multi-CCD loci are mosaicked
per-CCD-filtered onto a TAN canvas (§8.7); 0 cells needed a
cross-track dimension (max σ_xt(99%) = 0.13″ vs the 0.5″ threshold).

## Results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | candidates |
| --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 14 / 86 | 2 | 1.564 | 11 (11.9) | **0** |
| development | 5 / 30 | 0 | 1.247 | 4 (4.0) | 0 |

Top confirmatory cell: wolf-1069/tx/g at R̃ = 1.10 (global p = 0.695).
Mask sensitivity (development, 22 common cells): R̃_FWER 1.25 / 1.14 /
1.25, 0 candidates under each, ΔR median −0.004 (strict).

## Completeness (confirmatory, AB; median m90 threshold / final-candidate)

| band | persistent | flicker | visit | block |
| --- | --- | --- | --- | --- |
| g | 23.18 / 23.28 | 23.01 / 23.00 | 23.63 / 23.65 | 23.64 / 23.64 |
| r | 22.13 / 23.38 | 23.34 / 23.35 | 23.25 / 23.25 | 22.22 / 21.47 |
| i | 22.73 / 22.80 | 22.91 / 22.91 | — | 22.16 / 22.66 |
| z | 22.08 / 22.08 | 21.89 / 21.89 | 21.94 / 21.94 | — |
| Y | 21.00 / 21.00 | 20.99 / 20.99 | 20.59 / 20.59 | — |

Injections are image-level (Moffat β = 3 at each frame's fitted FWHM,
into each CCD image before the matched filter; median PRF-vs-Gaussian
throughput 0.737), ≥ 400 per cell over the frozen magnitude windows;
constraints carry Wilson intervals in the ledger. The static
flux-consistent veto fired on 4.3 % of threshold-recovered
confirmatory injections (its measured selection cost).

## Positive control

(60000) Miminko — the PS1 control, for cross-survey continuity — via
CADC SSOIS: 57 DECam instcal exposures 2013–2019 (g/r/i/z, Horizons
V 19.1–20.9), 32 surviving fetch + in-field star calibration.
Recovered in all four bands under the no-clip scoring at the exact
predicted position (peak offset 0.0″): S_max 216–330, R̃ 3.2–5.0, rank
p = 0.020–0.041; recovered stack magnitudes are 0.25–0.40 mag fainter
than the Horizons + solar-colour prediction, consistent with trailing
loss for a ~0.5″/min mover in 90 s exposures (irrelevant for the
station-kept survey hypothesis). Under the frozen 5σ single-epoch
clip the individually-detected mover epochs are removed — the same
documented interaction as the ZTF v2 control.

## Screening (layer 1)

NSC DR2 via the Data Lab `astro-datalab` client (anonymous sync;
object-cone → meas-by-objectid — direct meas cones time out): 93,362
ScreenMatch records. NSC DR2 is **time-partial** (catalogued epochs
end 2017–2019 in every probed corridor; recorded per corridor in
`runs/decam/screen_v1/catalog_stats.json`), so catalogued-detection
absence after ~2019 is meaningless and the flux-consistent static
veto leans on NSC mean objects + the deep forced photometry.
Recurrence triage on the pilot corridors found only chance
track-crossings (hd219134's corridor sits at b ≈ +2.5°).

## Coverage and data notes

- All 15 southern corridors have usable data (overlay
  `surveys/decam/targets/overlay_v1.md`: 10 ok / 5 crowded / 0 sparse
  under the frozen cut); calendar-month spread 4–12 → both parallax
  phases exist single-archive for most corridors, unlike PS1.
- Heterogeneous PI cadence is the dominant systematic: per-cell epoch
  counts run 5–88, and `not_constrainable` intervals (53 % of
  confirmatory threshold constraints) concentrate in Y/z and
  single-season cells. gj-687 and gj-1221 (LMC-adjacent) carry the
  richest cadence and survived the frozen deepest-per-night cap with
  88 / 69 epochs per role.
- Archive failure modes (~1 % of files; all recovered or recorded):
  `?hdus=` 500s → full-file fallback; sibling HDU-order mismatches →
  per-file order from the `api/header` page; served-bytes md5 drift →
  content-verified acceptance by EXPNUM; one pilot image is
  unrecoverably corrupt server-side (recorded product-corrupt).
- Observer: CTIO/W84 terrestrial site, validated against an
  independent astropy computation to ≤ 1.7 mas at z = 550 AU
  (`surveys/decam/results/observer_validation.json`).

## Interpretation

For the 15 southern corridors this is the first optical constraint of
any kind: the reflected-sunlight and self-luminous optical cells,
closed at g/r ≈ 23.2–23.4 AB (persistent, 90 % threshold
completeness) where the cadence supports the stack, were previously
open. Null results apply to the frozen hypothesis grid only (compact
station-kept relay, 550–10,000 AU, duty ≥ 0.5, |µ_resid| ≤ 1″/yr,
per temporal model); cells and intervals marked `not_constrainable`
in the ledger are exactly that. The archive keeps accumulating
(~yearly refresh extends the epoch hold-out and Y/z coverage), and the
three-archive joint stage can now add a southern optical member.
