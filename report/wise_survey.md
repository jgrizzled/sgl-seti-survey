---
title: "WISE/NEOWISE SGL Survey (v2 re-analysis under a frozen decision rule)"
date: 2026-08-22
status: "complete — blind confirmatory run of 61 endpoints / 488 cells: 0 candidates at family-wise α = 0.05; development set 27 endpoints / 216 cells: 0 candidates; 8,448 injection-calibrated constraints with confidence intervals; supersedes the withdrawn v1–v3 reports (removed in the v1 retirement, 2026-08-24; git history)"
---

# WISE/NEOWISE SGL Survey

**Pins:** registry v1.5 `sha256:794d9f90…` · hypotheses **v2.1**
(`surveys/wise/hypotheses.md`, freeze
`configs/v2_0_freeze.json` content hash `sha256:3c83e91e…`) · physics
cell unchanged from v1.0 · sglseti `tusay2022_eq5_7_v1` 1.1.0, observer
`wise-l1b-spacecraft` (`sha256:5afe2681…`) · AnalysisRuns
`run-af14a6b73cac` (confirmatory) and `run-ed152d2ac6cc` (development)
· every number below is generated from the ledger by
`surveys/wise/scripts/report_tables.py` (`surveys/wise/results/report_tables.md`).

## Summary

The v1 analysis (reports v1–v3, git history) claimed a 90 %-complete exclusion with an "empirical FAR
< 1/8". Its scientific review (2026-08-21) showed that the 8-control
threshold was a per-search rank statistic (≈ 78 chance crossings
expected, 70 seen), that completeness never saw the vetoes, and that
the geometry and provenance claims were metadata. Erratum v3.1 (appended to the v3 report, now in git
history) withdrew those claims. This report is the v2 re-analysis of the same
frames, cutouts and registry under **one decision rule frozen before
the confirmatory data were touched**
(`surveys/wise/v2_plan.md`; hypotheses v2.1):

- a 48-offset exchangeable null per cell, a per-cell normalised
  statistic R̃ = R / q95 and a survey-wide threshold R̃_FWER calibrated
  by 10,000 pseudo-experiments at **family-wise α = 0.05**;
- **image-level injections** of the IRSA empirical PRF into the
  calibrated frames (400 per cell, four temporal models, continuous z,
  L∞ motion box, physical spectra with colour corrections), classified
  with the same rule, giving threshold and final-candidate completeness
  with bootstrap intervals on gap-free reciprocal-distance intervals;
- only **calibrated vetoes** with injection-measured selection
  functions; everything else is an annotation;
- full covariance propagation (3 of 176 cells needed a cross-track
  dimension), the spacecraft observer, an asteroid positive control,
  content-hashed products and invariant tests.

**Result.** On the **confirmatory set** (54 corridors, 61 endpoints,
488 cells, run once, blind) **no cell reaches the family-wise
threshold** (R̃_FWER = 2.444; the highest cell, gj-1087/rx/W1, has
R̃ = 1.33 and global p = 0.96). On the development set (216 cells)
likewise 0 of 216 (R̃_FWER = 2.390). The v1-style count of cells with
R > 1 is 41 on the confirmatory set against a ring expectation of
60 — the "exceedances" of v1 were the threshold's own rank statistics.
The calibrated survey is **≈ 1–1.5 mag shallower** than v1 claimed:
median 90 %-completeness depth for a persistent source W1 13.4 / W2
12.3 / W3 9.7 / W4 7.4 (Vega), and for the worst of the four temporal
models (the duty ≥ 0.5 coverage) W1 12.5 / W2 11.1 / W3 7.9 / W4 4.3.
That is the price of a controlled error rate over ~500 heavy-tailed
local nulls, and it is now measured rather than assumed.

## 1. Decision rule (frozen; hypotheses v2.1)

Cell = endpoint × role × band (704). Per cell, S(z, µ) is the capped
inverse-variance stack over the 64 × 5 × 5 grid, S_max its maximum,
T the maximum over the 8 designated offset controls, R = S_max/T. The
null is the 48-offset ring (radii 20/30/40″ × 16 angles; the real
track is one more ring member under the null). R̃ = R / q95(ring);
the family statistic is max R̃ over the set's cells; R̃_FWER is its 95th
percentile under pseudo-experiments drawing one ring R̃ per cell (cells
taken independent — conservative). Cells whose ring is heavy-tailed
(max > 2.5 × q95) or radius-dependent (inner vs outer KS α = 0.01) are
void: 70 of 488 confirmatory cells (53 heavy-tail; 24 W1, 11 W2, 24 W3,
11 W4) and 27 of 216 development cells.

Two further null constructions are reported as annotations, not
pooled: the phase-coherence scramble (degenerate when a single static
star dominates one parallax phase) and trajectory randomisation
(donor tracks sample other parts of the corridor; differ from the ring
at KS α = 0.01 in 55 % of cells). The per-epoch scramble of the plan
destroys the static-sky coherence of the annual parallax return and is
only a diagnostic floor.

Calibrated vetoes: (a) a catalogued CatWISE2020 source (≥ 3 detections)
whose flux, pushed through the band's PRF-vs-Gaussian radial response
at every epoch's separation and stacked with the cell's weights,
accounts for the measured S within a factor 2 overall and in the
dominant phase; (c) the W3/W4 confirmation procedure (W1/W2 forced
photometry at the same node on the cryo frames against the 300 K
W2/W3 prediction at 3σ, post-cryo W1/W2 persistence, catalogue
counterpart). The held-out-epoch prediction test of v2.0 was demoted to
an annotation at step F (§5).

## 2. Confirmatory result (blind)

| quantity | value |
| --- | --- |
| endpoints / cells | 61 / 488 |
| void cells | 70 |
| R̃_FWER (α = 0.05) | 2.444 (pseudo-experiment max: median 1.89, 99 % 2.47) |
| candidates | **0** |
| highest cell | gj-1087/rx/W1: R̃ 1.33, R 1.33, global p 0.960, rank 1/48 (p_cell 0.02), p_phase 0.05 (phase split 787 / 171) |
| cells with R > 1 | 41 (ring expectation 60.2) |

The ten highest cells (R̃ 1.06–1.33) all have global p ≥ 0.96. Several
are the familiar single-phase static-star cells of v1 (gj-625 rx/tx
W1/W2: S_phase 1.6 / 211, 3.4 / 49; gj-588/tx/W1: 537 / −11) — here
they are simply not unusual relative to their own ring null, and no
veto was needed. Adjudication therefore produced no Candidate record
(`runs/wise/v2/records/candidate_confirmatory.jsonl` is empty by
construction; `adjudication_confirmatory.json` records the run).

## 3. Development-set result (exploratory)

27 endpoints / 216 cells; R̃_FWER = 2.390; **0 candidates**; highest
cell lalande-21185/tx/W2 at R̃ 1.03 (global p 0.999); 16 cells with
R > 1 against 26 expected. Quality-mask sensitivity on 214 common
cells: R̃_FWER 2.39 / 2.34 / 2.31 (primary / strict / loose), 0
candidates under each, ΔR vs primary median 0.00 (16–84 %: −0.05 to
+0.02 strict, −0.005 to +0.007 loose) — the mask choice is not a
systematic at the level of the result.

## 4. Completeness

400 image-level injections per cell (100 per temporal model), the IRSA
second-pass PRF of the frame's focal-plane element placed at the exact
sub-pixel position of a continuous (z, µ, cross-track) track, added in
DN before background estimation and the matched filter. Recovery =
R̃_peak ≥ R̃_FWER within 2 z-nodes and 1 µ-node of the injection.
Median 90 %-completeness magnitudes over cells × intervals
(confirmatory set; Vega; development set in parentheses):

| band | persistent | flicker p = 0.5 | visit-scale | long block | worst (duty ≥ 0.5 coverage) |
| --- | --- | --- | --- | --- | --- |
| W1 | 13.44 (13.46) | 12.71 (12.68) | 12.53 (12.50) | 12.69 (12.75) | **12.5** |
| W2 | 12.28 (12.41) | 11.44 (11.59) | 11.12 (11.38) | 11.36 (11.31) | **11.1** |
| W3 | 9.70 (10.04) | 8.76 (9.01) | 7.94 (8.41) | 8.23 (8.65) | **7.9** (threshold sensitivity only) |
| W4 | 7.35 (7.48) | 6.30 (6.46) | 4.28 (3.87) | 5.66 (6.68) | **4.3** (threshold sensitivity only) |

Threshold and final-candidate curves coincide: the flux-consistent
static-source veto fired on 3 of 42,830 threshold-recovered
confirmatory injections (0 of 20,414 on the development set) — an
injected source on the track is not explained by catalogued stars.
Every Constraint carries the per-interval 68 % bootstrap interval and
the injection count; 4,987 of 5,856 confirmatory records are
`recovery_curve`, the rest `not_constrainable` (void null, or a
logistic fit extrapolated outside the injected range). W3/W4 records
carry `completeness_kind = threshold`, `exclusion_claim = false`;
W3/W4 "long block" injections frequently have no epoch on at all (the
cell spans one or two cryo visits), so W3/W4 cannot constrain
long-block intermittency.

Measured throughput of the Gaussian matched filter on the empirical
PRF: W1 0.80, W2 0.75, W3 0.46, W4 0.58 (median over all injections
0.66). The v1 W4 kernel was also twice too wide in pixels (erratum v3.1
item 9, git history); v2 reads the pixel scale from the header.

## 5. Calibrated vetoes, annotations and the held-out test

Selection functions (confirmatory injections, threshold-recovered):
static flux-consistent veto 0.0 %; held-out-epoch prediction test pass
rate persistent 0.989, flicker 0.989, visit 0.891, **block 0.242**;
on null trajectories the same test passes 45 % overall and **63 % of
those above their cell's q95**. A test that rejects three quarters of
in-scope long-block sources while passing most null exceedances is not
a calibrated veto under the plan's own criterion (§6: it may not reject
a present, detectable in-scope signal), so at step F — on the
development set, before the confirmatory set was touched — it was
demoted to an annotation carried with its rates (hypotheses v2.1 §8).
Parallax-phase split, single-epoch share, companion-band significance,
cryo visit count, catalogue proximity and grid-edge fit are
annotations on every Candidate record; none was needed.

## 6. W3/W4

Searched and ranked with the same rule as W1/W2 (no cadence veto).
Void fraction is higher (24 of 122 W3 cells) because a single cryo
visit gives the ring little room. No W3/W4 cell approached R̃_FWER
(highest eps-eri/tx/W4 at R̃ 1.07, global p 1.0), so the confirmation
procedure (§1c) has not been exercised on a real exceedance; W3/W4
results are **raw threshold sensitivity for a 300 K blackbody with the
empirical PRF — no exclusion of warm structures is claimed**. For
orientation only: at the confirmatory median W3 m90 for a persistent
source (9.7, ≈ 4.3 mJy after colour correction) the v3 §6 arithmetic
gives a 300 K radiator radius of ~150 km at 550 AU; at the worst-model
depth (7.9) ~340 km.

## 7. Geometry, observer and positive control

- Covariance propagation (sglseti Monte Carlo, N = 2,000, seed
  20260822, z = 550 / 10,000 AU, epochs 2010 / 2017 / 2024): the 99 %
  cross-track half-width is ≤ 2.9″ for 173 of 176 endpoint-roles; the
  ε Ind Ba/Bb photocentre (rx 4.2″, tx 16.8″) and EZ Aqr tx (7.8″)
  exceed 0.5 × FWHM and were searched with a cross-track dimension
  (statistic and nulls maximised over offsets).
- Independent astropy-only locus (156 linear-astrometry cells):
  median residual 0.013″, max 3.3″; Tx agrees to ≲ 0.05″ once the 2d/c
  light time is included, Rx differs by the model's 2µz/c relay
  light-time term (largest for Barnard's Star at 10,000 AU).
- Observer: spacecraft position from the L1b headers (40,782 frames;
  geocentric radius 6,726–6,917 km; frame-checked to 3 km against
  astropy); locus shift ≤ 17 mas at 550 AU (v1 stated 2.5 mas).
- Positive control: asteroid (14000), Horizons ephemeris with WISE as
  observer, 127 W1 / 128 W2 frames over 2010–2024, driven blind
  through the production chain. Recovered at the Horizons position
  (peak offset ≤ 1″): W2 R̃ = 5.5 (a candidate under the survey rule),
  W1 R̃ = 1.49 (rank 1/48 but below R̃_FWER — at stack magnitude 16.8
  with 127 frames it is ~3 mag fainter than the W1 m90 of a ~800-epoch
  cell, consistent with the injections). Throughput against the NEOWISE
  pipeline's own profile-fit photometry of the same frames: −0.42 mag
  in W2 (MAD 0.09, 71 frames), −0.52 ± 0.40 in W1 — the PRF-vs-Gaussian
  response plus sub-pixel sampling, measured end to end.

## 8. Scope statement

Targeted coverage of the frozen 88-endpoint portfolio (77 corridors)
under the v1.0 physics cell: 550–10,000 AU, component-wise |µ| ≤ 1″/yr,
unresolved source, duty ≥ 0.5 under the four temporal models above.
Nothing here constrains the prevalence of SGL relays or the
Picky-Network hypothesis; population inference would need a separately
pre-registered generative model. Cross-archive confirmation (ZTF /
PS1 / SPHEREx at the predicted position) remains the primary
independent evidence for any future candidate.

## 9. Records and reproducibility

`runs/wise/v2/`: 176 + 7 tensors (9 kept trajectories, per-trajectory
summaries for 99, content hash of their cutouts in the header —
`--only-missing` rebuilds a tensor whose inputs changed), 183
injection files, null ensembles, completeness products, PRFs
(`sha256` per band), observer table, control run with verbatim
Horizons / TAP snapshots, `records/` (2 AnalysisRuns whose
`observation_set_hash` is the hash of the cutout digests actually
used; 8,448 Constraints, each `supersedes` the overlapping v1 record
of `run-1b2e86e9219e`; 0 Candidates). `tests/test_v2_invariants.py`
(15 tests: kernel and PRF normalisation, DN↔mag and Fν↔mag round trips,
zero-flux injection, stamp response = full re-convolution, gap-free
interval tiling, leave-one-out ratios, scramble invariance, FWER /
rank / BH, manifest stale detection, record fields, corridor table,
v1 threshold reproduction) run before every stage. Stages:
`surveys/wise/scripts/run_v2.sh A–I`.

## 10. Review findings — disposition

| finding | disposition |
| --- | --- |
| 1 false-alarm rate | closed: exchangeable ring null, R̃, FWER α = 0.05, rank statements with Wilson intervals, blind confirmatory set |
| 2 completeness without vetoes | closed: image-level injections through the full rule and the calibrated vetoes; threshold vs final-candidate curves with CIs |
| 3 99 % locus not propagated | closed: MC envelopes per cell, cross-track dimension where needed, spacecraft observer, independent check |
| 4 injection model | closed: empirical PRF, sub-pixel, DN before photometry, physical spectra with colour corrections, positive control; the Gaussian-kernel throughput is now measured (0.46–0.80) |
| 5 distance / motion / duty mismatches | closed: log-uniform draws, L∞ declared, gap-free intervals, four temporal models, continuous placement |
| 6 uncalibrated vetoes | closed: calibrated static/halo test and W3/W4 procedure only; held-out test and all other rules are annotations with measured rates |
| 7 quality mask | closed: primary / strict / loose reported (no effect on the result) |
| 8 population inference | closed by scope statement (§8) |
| 9 reproducibility | closed: content hashes, stale detection, invariant tests, ledger-generated tables |
| open | the W3/W4 confirmation procedure has not been exercised on a real exceedance (none occurred); W3/W4 long-block coverage is nil by cadence; three cross-track cells use a 3–5-offset approximation to the envelope |
