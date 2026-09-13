---
title: "PTF/iPTF corridor survey (frozen decision rule, blind confirmatory run)"
date: 2026-09-10
status: "complete — confirmatory 31 endpoints / 66 cells: 0 candidates at family-wise α = 0.05; development 14 / 36: 0 candidates; 1,632 injection-calibrated constraints; the PS1-era optical corridors re-observed from a second site with both parallax phases on half the cells"
---

# PTF/iPTF corridor survey

**Pins:** registry v1.5 · hypotheses `surveys/ptf/hypotheses.md` v1.0
(freeze `sha256:4915a38e50401cec…`, seed 20260910) · T0 = 56000 ·
engine `sglsurvey/` · adapter `sglsurvey/adapters/irsa_ptf.py`
· AnalysisRuns `run-a55afd948de5` (development) / `run-ef554e3f697d`
(confirmatory, blind) in `runs/ptf/v2/records/` · numbers from
`surveys/ptf/results/report_tables.md` · overlay
`surveys/ptf/targets/overlay_v1.md` · recon checks
`surveys/ptf/results/{wcs_linearity_check,flux_scale_check}.json`.

## Summary

The public PTF/iPTF level-1 archive (Palomar P48, 2009-03 → 2015-01,
g + Mould R, 60-s CCD exposures) is the only other 1″-class
multi-epoch optical archive of the PS1 era on the northern sky, from
the same telescope and site ZTF later used but with a campaign
cadence that revisits fields across the year. It was adopted for the
parallax phases PS1's single-phase 3π cadence lacks before 2018 and
for an era-independent re-observation of the PS1 corridors — not for
depth (R ~ 21 per epoch). Pipeline A ran end-to-end (5,933
Observations over 69 endpoints, 7,689 precise evaluations, 3,456
usable exposure product sets, 3,127 star-calibrated frames) and the
v2 statistical chain under a rule frozen before any analysis touched
the data — the second survey built directly on the v2 design with no
v1 exploratory phase (DECam was the first).

**No cell of the blind confirmatory set (66 cells over 27 corridors /
31 endpoints) reaches the family-wise threshold** R̃_FWER = 1.724;
8 cells have R > 1 against a ring expectation of 8.3, and 5 cells are
void (unstable ring null). Development set (four pilot corridors
forced + nine drawn): 36 cells, R̃_FWER = 1.422, 1 exceedance vs 3.9
expected, 0 candidates; quality-mask sensitivity (primary/strict/
loose) changes nothing. Median persistent-source m90 on the
confirmatory set is **g 21.1 / R 21.0 AB** (threshold completeness;
final-candidate identical — the static veto fired on 0.0 % of
threshold-recovered injections), about 2 mag shallower than PS1/ZTF
and 0.4 mag past the median single-epoch clip limit (20.7 AB).
Coverage is campaign-driven: 40 of the 62 Palomar-visible corridors
are searchable at the 5-epoch floor, per-cell epoch counts run
6–192, R carries 58 of the 66 confirmatory cells, and 118 of 528
confirmatory threshold constraints are `not_constrainable`.

## Rule (frozen)

As WISE v2.1 with PTF bindings: 192 × 5 × 5 grid (uniform in 1/z;
L∞ µ ≤ 1″/yr; T0 = 56000), single-epoch clip |S_e| ≤ 5, per-frame
weight cap, ring 20/30/40″ (48 offsets), R̃ = R/q95, survey-wide
max-R̃ null → FWER α = 0.05; veto = flux-consistent PS1 DR2
catalogued static (≥ 3 detections, transformed per band). Search
image = the sky-subtracted science image (the archive serves no
difference images), so the static veto applies at full weight. Fatal
dmask template 65533 (Laher et al. 2014 Table 15) + metadata seeing
≤ 4″ as the primary mask; strict adds seeing ≤ 2.5″ and |moonillf| ≤
0.8. Flux scale: per-frame PS1 DR2 mean-star calibration through the
identical matched filter — g = PS1 g, R = Jordi (2006) Cousins R
+ 0.21 mag to AB, dwarf-locus calibrator rule, gate ≥ 5 stars and
scatter ≤ 0.2 mag (3,127 of 3,456 frames pass; median scatter 0.048 g
/ 0.062 R with 26 / 30 stars; header MAGZPT differs by +3.9 mag and
is never used). Full TAN-SIP sampling (the linear Jacobian errs up to
0.39″ inside the locus span, 1.5″ at stamp corners); 9 of 90
endpoint-roles carry a cross-track dimension (max σ_xt,99 = 11.5″ —
eps-ind-b, ez-aqr, sirius, wise-0855). No epoch hold-out (closed
archive). Dev/confirmatory split stratified by confusion class with
ross128, kapteyn, ltt1445 and vanmaanen forced into development.

## Results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | candidates |
| --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 31 / 66 | 5 | 1.724 | 8 (8.3) | **0** |
| development | 14 / 36 | 0 | 1.422 | 1 (3.9) | 0 |

Top confirmatory cell: gj-783/rx/R at R̃ = 1.29 (global p = 0.336;
S 7.2 / 6.8 at the two parallax phases — a two-phase cell whose
phase-scramble and trajectory annotations both flag it, p_phase 0.010
/ p_traj 0.000, with R̃ well below the family threshold). Mask
sensitivity (development, 30 common cells): R̃_FWER 1.42 / 1.59 /
1.42, 0 candidates under each, ΔR median −0.005 (strict). Void cells:
eps-ind-b/tx/R, gj-526/rx/R, gj-783/tx/g, wolf-359/tx/g,
wolf-359/tx/R.

## Completeness (confirmatory, AB; median m90 threshold / final-candidate)

| band | persistent | flicker | visit | block |
| --- | --- | --- | --- | --- |
| g | 21.08 / 21.08 (15) | 20.68 / 20.68 (14) | 20.67 / 20.67 (6) | 20.67 / 20.67 (13) |
| R | 21.03 / 21.03 (373) | 20.57 / 20.57 (216) | 20.53 / 20.53 (134) | 20.47 / 20.47 (78) |

Development set: g persistent 21.94 (55 intervals; ross128 with 99 g
epochs and vanmaanen with 55 dominate), R 20.85 (208). Injections are
image-level (Moffat β = 3 at each frame's header SEEING, into the
science cutout before the matched filter; median PRF-vs-Gaussian
throughput 0.68), 400 per cell over the frozen 21.5 ± 2 mag windows;
constraints carry Wilson intervals in the ledger. 118 of 528
confirmatory threshold constraints are `not_constrainable` (78
insufficient recovery fit — sub-floor epochs or fits brighter than
the clip limit; 40 in the five void cells).

## Positive controls

The ZTF/PS1/DECam control (60000) Miminko has only 10 PTF frames over
2009–2015 and 4 survive the ZP gate — PTF's coverage of any one
asteroid's track is too sparse — so both controls were found inside
the pilot-corridor exposures by SkyBoT (ross128 and vanmaanen lie on
the ecliptic; `surveys/ptf/scripts/asteroid_control_recon.py`).

- **(798452) 2012 QR36**, V 20.3–20.5, 46 frames / 11 nights
  (catalogue regime: median single-frame S/N 8.5, 39/40 frames above
  5σ): recovered under the frozen rule at R̃ 6.7 (g, 27 epochs) /
  3.1 (R, 13), rank p = 0.020 (0/48 ring controls above); without the
  clip the stack magnitude is 20.50 / 20.22 vs 20.39 / 20.13 predicted
  from Horizons V + solar colours (+0.12 / +0.10 mag), pinning the
  flux scale. Under the clip the recovered magnitude is 1.1–1.3 mag
  fainter — the documented removal of individually detected mover
  epochs, as in the ZTF and DECam controls.
- **(388125) 2005 UP482**, V 21.6, 35 frames / 11 nights, 16 R frames
  calibrated (stack regime: below the single-frame limit): recovered
  at R̃ 4.8, rank p = 0.020, 21.52 vs 21.32 predicted (+0.20 mag).
  Without the clip a designated control lands on a star (T = 12) and
  the ratio drops to R̃ 0.66 — the clip regime is the one the search
  operates in.

## Coverage and data notes

- 62 Palomar-visible corridors queried; 40 searchable (25 with ≥ 20
  epochs in the best cell and ≥ 4 30-day day-of-year bins, 14 thin,
  1 single-phase), 13 below the 5-epoch floor, 9 with no PTF exposure
  on the corridor at all (`targets/overlay_v1.md`). The 22
  non-searchable corridors enter the ledger as coverage-without-
  statistic.
- R dominates (3,123 of 3,456 usable frames); g exists on ~⅓ of the
  corridors and forms 8 confirmatory cells. 33 of 66 confirmatory
  cells are two-phase (both parallax phases present), the cell the
  survey was adopted for; the rest are single-phase and carry no
  phase annotation.
- Archive integrity: 3,550 mask + 3,456 science cutouts served with
  zero 404s; published MD5s on every full product. Gate attrition is
  the only data loss: 321 R frames fail the 0.2-mag ZP-scatter gate
  (bright-sky and poor-seeing nights).
- Observer: Palomar P48 = the ZTF site model (validated 16 mas at the
  ZTF pilot); exposures are 60 s, so mid-exposure evaluation moves
  the locus < 5 mas.

## Interpretation

For the 40 searchable northern corridors this adds a second-site,
second-cadence optical constraint in the PS1 era: persistent relays
brighter than R ≈ 21.0 / g ≈ 21.1 AB (90 % threshold completeness,
median over cells and distance intervals) are excluded at family-wise
α = 0.05 on the frozen hypothesis grid (compact station-kept relay,
550–10,000 AU, duty ≥ 0.5, |µ_resid| ≤ 1″/yr, per temporal model),
and 33 cells now hold both parallax phases before 2015 — the
persistence test PS1 alone cannot run. Null results apply to the
frozen grid only; cells and intervals marked `not_constrainable` in
the ledger are exactly that. A three-archive joint stage (PS1 + PTF +
ZTF) on a common µ reference epoch would be a rebuild, not a relabel
(learnings §10), and is the natural follow-on; the archive is closed,
so there is no yearly refresh.
