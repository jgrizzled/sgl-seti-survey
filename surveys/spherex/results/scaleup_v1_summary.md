# SPHEREx scale-up v1 — all 15 ZTF-inaccessible corridors (2026-08-20)

Hypotheses `spherex-hypotheses-v1.0` (+ role-coincidence amendment);
registry v1.5; corridors = the 15 universal-list systems at Dec < −28°
(19 endpoints: the 3 pilot stars + ross248, cyg61 A/B, struve2398 A/B,
grb34 A/B, gj1221, gj338 A/B, gj625, gj251, hd219134, wolf1069,
gj3512, gj13157). Same chain as the pilot; calibration products under
`runs/spherex/calib_v2` (pilot `calib_v1` untouched); discovery,
precise, control and screen runs are the shared `*_v1` directories
extended in place. Data volume 6.7 GB, 14,579 slim cutouts, 0 missing.

## Coverage

24,099 detector-exposures discovered, 32,526 precise evaluations:
29,880 usable, 2,082 partial, 564 unusable. Every corridor covers
550–10,000 AU at every visit. Struve 2398's corridor has a single
observing season so far (phase split 0 / 63 in D3) and GJ 1221, GJ 251
and GJ 687 are strongly season-weighted (phase-1 counts 2–5 in D3).

## Flux-scale control (all 15 corridors)

54,265 star measurements: median Δm per detector −0.19 (D1), 0.00,
+0.04, +0.11, +0.07, +0.15 (D6); scatter 0.07–0.11 mag; 0.054 mag
all-band offset applied.

## Layer-1 screen

2.1 million single-epoch 5σ peaks, 87,506 within 10″ of a track, 3,529
position clusters, 2,790 static; every one of the 739 unresolved
clusters has ≤ 2 epochs — nothing follows a parallax track.

## Calibration (`calib_v2`, AnalysisRun `run-4f8c9fd41ce9`)

228 (endpoint, role, detector) searches, T = 1.4–7.7. **28
exceedances → 28 Candidate records, all vetoed, 0 retained**
(20 phase-split, 2 single-phase, 6 role-coincidence). 28/228 = 12 %
matches the FAR < 1/8 of the max-of-8-controls rule.

**1,824 Constraint records**, all `recovery_curve`. Median m90 per
detector (duty ≥ 0.5, 90 % recovery, AB; rx–tx range):

| endpoint | D3 epochs | D1 | D2 | D3 | D4 | D5 | D6 |
|---|---|---|---|---|---|---|---|
| 61-cyg-a | 57 | 18.0 | 17.7–18.3 | 18.2 | 19.0–19.1 | 18.7–19.0 | 18.3–18.5 |
| 61-cyg-b | 57 | 18.0–18.3 | 18.1–18.4 | 18.3–18.4 | 19.1–19.4 | 18.8–18.9 | 18.6–18.7 |
| gj-1221 | 193 | 18.7–18.8 | 18.3–18.6 | 19.0–19.1 | 19.1–19.5 | 19.1 | 19.2–19.4 |
| gj-13157 | 75 | 18.0 | 17.7 | 17.5–17.7 | 18.4–18.8 | 18.5–18.6 | 18.2–18.4 |
| gj-251 | 45 | 18.0–18.1 | 17.6–17.8 | 17.6–18.3 | 18.1–18.7 | 18.4–18.8 | 18.3–18.4 |
| gj-338-a | 54 | 19.8–19.9 | 19.7–20.5 | 20.0–20.4 | 20.1–20.5 | 19.2–19.5 | 18.8–19.0 |
| gj-338-b | 54 | 19.9–20.1 | 19.5–20.3 | 19.7–19.8 | 20.3–20.5 | 19.2–19.5 | 19.1 |
| gj-3512 | 58 | 18.5–18.9 | 18.7–18.9 | 19.3–19.5 | 19.1–19.4 | 18.3–18.8 | 18.6–18.7 |
| gj-625 | 171 | 19.9–20.0 | 19.8–20.1 | 20.4–20.9 | 20.5–20.6 | 19.2–19.6 | 19.4–19.5 |
| gj-687 | 2522 | 19.4–20.2 | 19.2–19.7 | 19.7–19.9 | 20.3–20.8 | 19.5 | 20.2–20.7 |
| groombridge-34-a | 82 | 18.9–19.1 | 19.1–19.5 | 19.4 | 19.8 | 19.2 | 18.8–18.9 |
| groombridge-34-b | 91 | 18.8–19.3 | 19.4 | 19.4–19.8 | 19.4–19.7 | 19.3–19.6 | 18.8–19.0 |
| hd-219134 | 68 | 17.8–18.0 | 17.4–17.9 | 17.6–17.8 | 18.0–18.3 | 17.5–18.5 | 17.5–18.2 |
| lalande-21185 | 67 | 19.5 | 19.4–19.8 | 20.0–20.3 | 20.1–20.5 | 19.5–19.6 | 19.0–19.3 |
| ross-248 | 78 | 19.0–19.6 | 18.9–19.2 | 19.4–19.6 | 19.9–20.1 | 19.1–19.3 | 18.8–19.1 |
| sigma-dra | 108 | 19.8–19.9 | 19.5–19.6 | 20.1–20.3 | 19.9 | 19.9–20.0 | 18.9–19.3 |
| struve-2398-a | 74 | 19.8–20.0 | 19.2–20.0 | 20.2–20.3 | 20.2–20.7 | 19.3–19.7 | 19.0–19.1 |
| struve-2398-b | 75 | 19.6–20.0 | 19.7–19.9 | 19.9–20.3 | 20.6–20.7 | 19.4 | 19.0–19.2 |
| wolf-1069 | 74 | 18.9–19.2 | 18.9–19.4 | 19.1–19.3 | 19.8–20.0 | 19.4–19.7 | 18.5–18.9 |

Depth is set by confusion, not epoch count: the four shallowest
corridors (HD 219134, GJ 13157, GJ 251, 61 Cyg; 17.5–18.5 AB) are the
low-Galactic-latitude ones (|b| ≲ 15°), and the deep-field GJ 687
corridor with 2,500 epochs is no deeper than clean high-latitude
corridors with 60.

## Carried-forward caveats

- Struve 2398 (single season) and the season-weighted corridors need
  the next quick release before their phase test has power.
- Component endpoints of the four binaries share a corridor; their
  loci differ by arcseconds, so their constraints are not independent.
- v2 items from the pilot stand: template-residual variance weights,
  cross-detector phase test, 16 controls.

## v2 estimator rerun (`calib_v3`, AnalysisRun `run-4fbed2e52668`, 2026-08-21)

Same cutouts, templates and screen; two estimator changes and two
adjudication rules (hypotheses amendment 2):

1. **Variance calibration.** Epoch variance = per-cutout MAD variance ×
   the template's per-node reduced χ² (clipped 0.05–20). Per-band m90
   improves by a median **+0.44 mag** (range −0.19 to +1.70; GJ 687
   D3 +1.1, D1 +0.7 — the deep field was the most under-weighted).
2. **Joint six-detector stack** (band `ALL`, flat-F_ν injection): a
   further **+0.2 to +0.9 mag** over the best single detector on 37 of
   38 endpoint×role stacks (GJ 687 tx −0.2: its D6 stack alone is
   already the deepest). ALL-band m90 19.1–21.9 AB; per-band 18.2–21.2.
3. **Candidates.** 266 searches (228 per-band + 38 joint), 37
   exceedances, all vetoed: 31 phase-split, 3 single-season, 2
   template-coverage, 1 role-coincidence. The two template-coverage
   vetoes are GJ 13157 rx/tx D2 (S = 4.1 / 3.5 vs T = 2.5 / 3.1) at
   z = 10,000 AU — the grid edge where the track moves 3 px/yr — with
   the template defined for only 48–51 % of epochs there, phase-0
   counts of 3–7, a failed phase split in the joint stack (0.8 / 4.7),
   and a W1 = 7.7 mag star (CatWISE J095137.88−591717.6) plus three
   12–13 AB static sources within 10–19″ in a b = −3° field: PSF-wing
   residual, not a relay. Both rules were added after this case and
   are frozen as amendment 2.
4. **Constraints.** 2,128 records (1,824 per-band + 304 joint). Joint
   inner-interval (550–720 AU) limits on the clean corridors: grey
   reflectors ≳ 1.0 × 10⁵ km (albedo 0.1), 700 K emitters ≳ 14 km.

| endpoint | best single band m90 (rx / tx) | joint ALL m90 (rx / tx) |
|---|---|---|
| lalande-21185 | 20.4 / 20.5 | 21.2 / 21.2 |
| gj-687 | 20.8 / 21.1 | 21.2 / 20.9 |
| sigma-dra | 20.7 / 20.8 | 21.2 / 21.3 |
| gj-625 | 21.2 / 21.0 | 21.9 / 21.6 |
| struve-2398-a / b | 20.9 / 21.1 · 21.1 / 21.0 | 21.5 / 21.2 · 21.7 / 21.4 |
| gj-338-a / b | 20.8 / 20.7 · 20.5 / 20.4 | 21.5 / 21.5 · 21.3 / 21.3 |
| groombridge-34-a / b | 20.2 / 20.3 · 20.3 / 20.4 | 20.7 / 21.0 · 20.9 / 21.1 |
| ross-248 | 20.3 / 20.6 | 20.9 / 20.9 |
| wolf-1069 | 20.5 / 20.1 | 20.8 / 20.7 |
| gj-1221 | 19.8 / 19.8 | 20.1 / 20.1 |
| gj-3512 | 19.7 / 19.9 | 20.0 / 19.9 |
| gj-251 | 19.8 / 19.5 | 20.1 / 20.2 |
| 61-cyg-a / b | 19.4 / 19.6 · 19.6 / 19.6 | 19.7 / 19.6 · 19.8 / 19.8 |
| hd-219134 | 19.2 / 18.8 | 19.5 / 19.1 |
| gj-13157 | 18.9 / 19.2 | 19.6 / 19.6 |

The calibration-stage candidate census now uses five vetoes (phase
split, single season, role coincidence, template coverage, bright
static neighbour); `calib_v2` is retained as the v1-estimator record.

## v3 estimator rerun (`calib_v4`, 2026-08-21): 16 controls + cross-detector phase test

Same cutouts/templates; hypotheses amendment 3.

- **16 offset controls** (eight added: ±50″ diagonals, ±75″ Dec, two
  skew offsets). Exceedances 37 → **21 of 266 searches (8 % ≈ 1/16)**.
  Depth cost: median 0.00 mag per (endpoint, role, band) stack, worst
  −0.69 where a new control raised T; joint-stack m90 range 19.1–21.8.
- **Cross-detector season-complement test**: where the exceeding
  stack's own phase split is undefined (< 3 epochs in one phase), the
  lacking phase is stacked across all six detectors at the same cell.
  It adjudicated the two deep-field exceedances (GJ 687) that had
  previously been "single-season"; both fail (S < 2).
- 21 Candidate records, all vetoed: 17 phase split, 2 cross-detector
  complement, 1 template coverage, 1 role coincidence. **0 retained.**
- 2,128 Constraints re-issued under the 16-control thresholds; this is
  the record the northern run extends (`calib_v4`).
