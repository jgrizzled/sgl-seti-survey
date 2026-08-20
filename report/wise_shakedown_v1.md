---
title: "WISE/NEOWISE SGL Shakedown — Survey Report v1"
date: 2026-08-19
status: "Complete: plan §4 steps 1–8; registry v1.1; no surviving candidate; 768 constraints"
---

# WISE/NEOWISE SGL shakedown — survey report v1

> **Superseded by [`wise_survey_v2.md`](wise_survey_v2.md)**
> (2026-08-19), which covers the expanded 26-endpoint survey, the
> universal target portfolio, and batches 2–3. This file remains as
> the 12-endpoint shakedown snapshot.

**Summary.** We searched 14.5 years of WISE/NEOWISE single-exposure
imaging (2010–2024) for compact persistent sources on the solar
gravitational lens focal lines of nine nearby stellar systems (twelve
endpoint hypotheses, receive and transmit roles each). The pipeline ran
end-to-end from frozen hypotheses through injection-calibrated
constraints: 19,444 archival frame-band products were discovered and
pinned, 21,385 precise usable-pixel intersections evaluated, 67,374
catalog detections screened against predicted tracks, and joint
(relay-distance × residual-motion) forced-photometry stacks searched
under predeclared empirical thresholds. **No candidate survives.** Six
threshold or marginal excesses were individually adjudicated and vetoed
(static field stars, a bright-star halo, and control-level placement
variance). The survey delivers **768 injection-calibrated constraint
records** — 90%-recovery flux limits over 550–10,000 AU for every
endpoint, role, and WISE band — with full content-addressed provenance.

## 1. Search definition (frozen hypotheses)

Hypothesis freeze `wise-hypotheses-v1.1` (parameters frozen 2026-08-18;
endpoint set expanded 2026-08-18 after the pilot passed):

- **Endpoints (12, component-level for multiples):** Barnard's Star,
  Ross 154, Lalande 21185, Proxima Cen, Wolf 359, Ross 248 (linear
  provider); α Cen A, α Cen B, Sirius A, Sirius B, GJ 65 A, GJ 65 B
  (two-body orbit provider). Registry v1.1
  `sha256:27895a11f4df…` with per-value provenance; curation evidence in
  `surveys/wise/notes/registry_curation.md` (notably: the hip2 Sirius
  solution validated *as* the barycenter against 2MASS to 0.16″; the
  GJ 65 barycenter constructed from mass-weighted Gaia components and
  validated to 0.34″/17.4 yr, orbit convention closed to 0.013″).
- **Hypothesis cell:** station-kept relay on the Sun–star focal line,
  Rx and Tx as separate hypotheses; relay distance 550–10,000 AU
  (log-uniform prior); unresolved point source; duty cycle ≥ 0.5;
  residual motion |µ| ≤ 1″/yr; 99% confidence + 10″ padding.
- **Geometry:** sglseti `tusay2022_eq5_7_v1` v1.1.0 (pinned commit
  `21f6f3d`), Earth-center observer (declared ~2.5 mas approximation),
  `astropy_builtin` ephemeris.
- **Out of scope (separate model IDs, not covered by this report's
  nulls):** barycenter and planetary endpoints, duty cycle < 0.5,
  |µ| > 1″/yr, extended morphologies, relay swarms, inactive relics.

## 2. Data and usable coverage

Archive: IRSA `wise.neowiser_merge_p1bm_frm` (verified union of the
4-band cryo, 3-band cryo, post-cryo, and NEOWISE-R phases; 62.7M
frame-band rows), with products from the unified IBE `merge` tree
(md5-verified). Usable coverage after exact-WCS + mask + quality tests
(counts are W1 epochs for the Rx role; Tx within a few frames; W2
essentially equal to W1):

| Corridor (antipode of) | Usable W1 epochs | Span | W3 / W4 (2010 cryo) |
| --- | --- | --- | --- |
| Barnard's Star | 396 | 2010-03 – 2024-02 | 32 / 16 |
| Ross 154 | 342 | 2010-04 – 2024-04 | 15 / 15 |
| Lalande 21185 | 404 | 2010-05 – 2024-06 | 15 / 15 |
| α Cen A / B | 428 / 426 | 2010-02 – 2024-01 | 38 / 19 |
| Sirius A / B | 422 | 2010-04 – 2024-05 | 18 / 18 |
| Proxima Cen | 416 | 2010-02 – 2024-01 | 40 / 20 |
| Wolf 359 | 334 | 2010-06 – 2024-06 | 14 / 14 |
| Ross 248 | 446 | 2010-01 – 2024-06 | 22 / 22 |
| GJ 65 A / B | 402 | 2010-01 – 2024-06 | 19 / 19 |

Masks restructure coverage substantially: 44% of geometric hits carry
2–6 **disjoint** covered relay-distance intervals (frame edges and
fatal-bit masks splitting the locus) — interval-valued coverage is the
norm, and every intersection record carries its exact intervals.

## 3. Pipeline and records

Five immutable, content-addressed record streams (IDs via sglseti's
provenance hashing; raw archive responses snapshotted verbatim):

| Record | Count | Produced by |
| --- | --- | --- |
| QuerySnapshot | 267 | discovery + screening TAP queries |
| Observation | 19,444 | coarse discovery (`coarse_v1`) |
| IntersectionEvaluation | 52,440 coarse + 21,385 precise | nominal-WCS pass; exact WCS+mask `covered_z_intervals` pass |
| ScreenMatch | 67,374 | single-exposure catalog screening (`screen_v1`) |
| AnalysisRun | 3 | stack (superseded grid), calibration v1.0, calibration v1.1 (`run-20040d16f249`) |
| Constraint | 768 | injection calibration (`calib_v1`) |
| Candidate | 6 (all vetoed) | stack + calibration adjudications |

Search geometry: 64-node relay-distance grid **uniform in 1/z**
(constant ~5.6″ on-sky spacing — an early log-spaced grid left ~17″
gaps at low z and was replaced) × 5×5 residual-motion grid.

## 4. Search results

**Catalog screening:** ~90–200 chance matches per visit per
endpoint-role, as expected under a 10″ screen. A fixed-z point-recurrence
filter found peaks up to 13/27 visits, every one resolved as a **static
background star** by the parallax-phase test: their detections occupy a
single ~5-day day-of-year window each year (sky rms 0.17–0.28″), where
a genuine fixed-z relay is carried to *both* alternating parallax
phases by the locus. The phase test is this survey's key cheap veto.

**Forced-photometry stacks:** matched-filter fluxes on 8,200+ int/unc
cutout pairs, stacked per (z, µ) cell over all usable epochs at a
common zero point. Detection thresholds were predeclared as the maximum
grid significance over **8 offset-control trajectories** per
endpoint-role-band (empirical false-alarm rate < 1/8 per grid search).
The W1/W2 search is field-star-contamination limited (control maxima
S ≈ 20–2,650), not noise limited.

**Adjudications — all six candidates vetoed:**

| Cell | Excess | Verdict |
| --- | --- | --- |
| lalande/rx/W1, lalande/tx/W2 (stack marginals) | +9% / +18% over 1 control | far inside 8-control null (S=30 vs T=77; S=17 vs T=74) |
| sirius-a/rx/W2 | S=871 vs T=451 | phase-imbalanced 941/144; bright-star halo; screening excludes a real ~11–12 mag source |
| proxima-cen/tx/W2 | S=250 vs T=214 | single-phase static star (−6.7 / +284.8) |
| barnard/tx/W3, ross-154/tx/W3 | +1% over threshold | single-visit cryo cells; placement variance |

**Conclusion: no track-consistent source in any of the 24
endpoint-role hypotheses.**

## 5. Completeness-calibrated constraints

Injections: Gaussian-PSF point sources at every grid node via
matched-filter linearity, randomized duty-cycle realizations
(32 repeats, seed 20260818), declared grid/µ mismatch factors;
recovered when the injected cell exceeds the predeclared threshold.
Each of the 768 Constraint records states: *sources brighter than m₉₀
recovered with ≥ 90% probability* for its endpoint × role × band ×
z-interval, duty ∈ [0.5, 1], |µ| ≤ 1″/yr.

Median m₉₀ across 550–10,000 AU (duty ≥ 0.5; Vega mag; Rx / Tx):

| Endpoint | W1 | W2 | W3 | W4 |
| --- | --- | --- | --- | --- |
| Barnard's Star | 15.6 / 14.1 | 15.0 / 14.2 | 11.2 / 12.7 | 10.8 / 10.6 |
| Ross 154 | 15.5 / 15.3 | 14.4 / 13.8 | 12.8 / 12.8 | 10.5 / 10.5 |
| Lalande 21185 | 16.9 / 17.4 | 16.7 / 15.4 | 13.6 / 12.5 | 11.0 / 10.3 |
| α Cen A | 15.0 / 14.0 | 14.4 / 14.1 | 11.7 / 10.8 | 9.8 / 9.9 |
| α Cen B | 15.3 / 14.0 | 14.4 / 14.3 | 11.0 / 10.9 | 9.8 / 9.9 |
| Sirius A | 13.0 / 15.2 | 13.5 / 13.5 | 10.5 / 9.2 | 10.2 / 9.6 |
| Sirius B | 13.6 / 14.5 | 13.7 / 13.2 | 10.6 / 10.5 | 10.1 / 9.8 |
| Proxima Cen | 14.8 / 14.8 | 13.1 / 14.4 | 10.7 / 10.8 | 10.0 / 9.5 |
| Wolf 359 | 16.6 / 15.1 | 16.5 / 15.7 | 11.7 / 11.3 | 9.6 / 8.8 |
| Ross 248 | 15.7 / 15.9 | 14.7 / 14.9 | 12.5 / 13.4 | 9.6 / 10.8 |
| GJ 65 A | 16.5 / 17.1 | 15.3 / 15.9 | 13.1 / 12.0 | 10.6 / 9.3 |
| GJ 65 B | 16.3 / 17.1 | 15.3 / 16.0 | 13.2 / 11.9 | 10.6 / 9.3 |

Depth is set by corridor crowding and bright-star artifacts, not photon
noise: sparse corridors (Lalande, GJ 65, Wolf 359) reach W1 ≈ 16.5–17.4;
the Sirius corridor bottoms at W1 ≈ 13. Full 64-node curves per pair ×
band × duty: `runs/wise/calib_v1/m90_curves.npz`.

## 6. Physical interpretation of the flux limits

Band physics is not interchangeable (hypotheses §6). Converting
representative depths (W1 15.5 / W2 14.4 / W3 12.6 / W4 10.4) to flux
density: ~195 / 300 / 290 / 580 µJy. What the nulls do and do not
constrain, for a source at relay distance z (observer distance ≈ z):

| Emitter model | z = 550 AU | z = 1,000 AU | z = 3,000 AU |
| --- | --- | --- | --- |
| **300 K radiator** (W3 12.6) — radius ≥ | **39 km** | 71 km | 212 km |
| 300 K radiator (W2 14.4) — radius ≥ | 230 km | 418 km | 1,260 km |
| 150 K radiator (W3 12.6) — radius ≥ | 311 km | 566 km | 1,700 km |
| 12 K solar-equilibrium blackbody | unconstrained at any physical size (peaks ~240 µm) |
| Reflected sunlight, aᵍᵉᵒ = 0.1, full phase (W1 15.5) — diameter ≥ | 4.8×10⁵ km | 1.6×10⁶ km | — |

Honest reading: this survey meaningfully constrains **warm, actively
powered structures** — a ~40 km-scale 300 K radiator at the near focal
line would have been recovered at ≥90% probability (W3, 2010 epochs
only; duty ≥ 0.5 within those epochs) — and W1/W2 constrain hot
components, nonthermal/beacon emission at these flux densities over the
full 14.5-yr baseline. It does **not** constrain passive cold
infrastructure (12 K equilibrium) at any plausible size, nor
planet-scale reflectors; those require far-IR/submm archives or
occultation/dynamical methods (plan §5 "dark infrastructure").

## 7. False-alarm and trials accounting

Declared search family: 24 endpoint-roles × 4 bands × 1,600 (z, µ)
cells, maximized per pair-band. The threshold rule (max of 8 matched
control trajectories per pair-band) fixes the empirical per-search
false-alarm rate below 1/8 without Gaussian assumptions — necessary
because per-pixel uncertainties exclude confusion noise, making
nominal S values inflated. Observed: 4 of 96 pair-band searches
exceeded threshold; all were vetoed by phase structure, cross-band
consistency, or catalog-screening exclusion (§4). FAR resolution is the
main statistical limitation (see §8).

## 8. Caveats and gaps

1. Gaussian-PSF model (search and injections): measures pipeline
   throughput, not PSF-wing mismatch; an empirical-PSF upgrade would
   modestly change absolute depths.
2. 8 control trajectories per pair-band → coarse FAR resolution;
   production runs want O(100) sky-rotated/time-scrambled corridors.
3. Constant-background photometry biases crowded corridors (α Cen,
   Sirius); injections absorb this to first order.
4. No known-moving-object positive control was run (SSO association
   snapshots are retained for this); synthetic injections stand in.
5. W3/W4 constraints rest on a single 2010 cryo visit: no recurrence
   or phase power; duty cycle is sampled only within that visit.
6. Catalog screening layer is not itself a null (catalog incompleteness);
   the constraints rest on the forced-photometry stack.
7. Diagonal astrometric covariance by declared assumption (propagated
   99% loci ≤ 0.1″ against 10″ padding — immaterial here).

## 9. Reproducibility

Everything regenerates from tracked code + pinned identities: registry
v1.1 `sha256:27895a11f4df…`; hypotheses v1.1 (content hash on every
record); sglseti `21f6f3d`; geometry model `tusay2022_eq5_7_v1` v1.1.0;
calibration AnalysisRun `run-20040d16f249` (seed 20260818). Raw archive
responses (267 snapshots), product checksums, and per-stage run configs
live under `runs/wise/` (4.8 GB, gitignored, rebuildable via
`surveys/wise/scripts/` in pipeline order). Ledger-ready records:
observation, intersection, analysis-run, constraint, and candidate
JSONL streams per plan goal 5.

## 10. What the shakedown taught us (inputs to the next adapters)

- **The merged L1b inventory + IBE tree is a one-stop discovery
  surface**; frame-inventory WCS supports server-side coarse discovery.
  Mask files carry full WCS — masks alone give exact astrometry and
  valid pixels at 30× less volume than images.
- **Interval-valued z-coverage is the norm** (44% of hits disjoint):
  the record schema's explicit intervals are load-bearing.
- **Parallax-phase structure is a decisive, nearly free veto** for any
  archive with ~6-month revisit cadence; build it into stage-2
  likelihoods from the start.
- **Uniform-in-1/z gridding** is the right relay-distance
  parameterization for any imaging archive (constant on-sky spacing).
- **Crowding, not depth, limits corridor sensitivity**; corridor
  background belongs in target prioritization (plan §5) ahead of raw
  epoch counts.
- Compute/volume: full pilot ≈ 5 GB archive data, ~8 h wall-clock,
  dominated by TAP queries and per-frame FFTs — comfortably scalable
  to ~20 systems before needing archive-side compute.
