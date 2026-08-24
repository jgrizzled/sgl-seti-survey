---
title: "sgl-seti-survey — Project Plan"
date: 2026-08-21
status: "v0.8 — v2 programme steps 1–7 done: all five surveys under the frozen rule (0 candidates in every blind confirmatory run); v1 retired 2026-08-24; next: step 8 (archive-family expansion on the v2 design)"
tags:
  - SETI
  - technosignatures
  - solar-gravitational-lens
  - archival-search
---

# sgl-seti-survey — Project Plan

## 1. Introduction

### 1.1 Background

A technologically advanced interstellar network might place communication
relays on or near the gravitational focal lines of stars, using stellar
lenses to reduce the power and aperture required for links between
neighboring systems. For the Sun, the focal line for rays grazing the
photosphere begins near 547.8 AU, so a relay serving a target star would
sit roughly on the anti-star ray at ~550–10,000 AU heliocentric distance.

This hypothesis is unusually testable for SETI: it predicts _where_ on
the sky a relay associated with a given star should appear, _how_ that
position moves (a huge annual parallax of ~20–375 arcsec depending on
distance, plus secular drift opposite the target star's proper motion),
and _when_ Earth crosses hypothesized beam geometries between the Sun and
nearby stars. The strategy, signal classes, and prior art (Gillon,
Hippke, Tusay et al., Marcy et al.) are laid out in the sglseti repo
under `notes/seti_strategies_for_sgl_technosignatures.md`.

### 1.2 Motivation for this project

Decades of archival astronomy data already exist that incidentally
covers these predicted sky regions and time windows — all-sky infrared
surveys, multi-epoch optical imaging, radio pointings at nearby stars,
event-level high-energy archives. Nobody has systematically asked which
of those observations intersected an SGL corridor or crossing window,
and what they would have seen if a relay were there. Mining archives is
far cheaper than new telescope time, closes real cells of the search
parameter space, and produces the target lists and priors that any
future dedicated observing campaign would need anyway.

**sgl-seti-survey is the archive-facing consumer of sglseti**: it
discovers which archival observations intersect SGL predictions,
cross-matches them precisely, and reports the results. (A full
observation-and-coverage ledger is a separate future project; this
project may produce inputs for it but does not build it.)

### 1.3 The sglseti package

`sglseti` (assumed available in an adjacent folder, `../sglseti`; not yet
on PyPI) is the stateless, offline, reproducibility-focused geometry engine
for this project.

The resulting package boundary provides:

- **Historical and interval target generation** — exact archival point
  epochs or observation intervals, keyed by caller-supplied IDs, with
  adaptive Rx/Tx loci over relay distance and conservative swept regions.
- **Versioned target-state models** — linear astrometry, acceleration,
  orbital solutions, and externally sampled ephemerides, with full
  endpoint identity and provenance.
- **Uncertainty products** — propagated role- and epoch-dependent
  uncertainty kept distinct from assumed search padding.
- **Observer-state models** — Earth center, terrestrial sites, and
  time-dependent spacecraft states from pinned kernels or state tables.
- **Beam-crossing geometry** — closest approach, impact parameter,
  transverse speed, side of axis, and parameterized windows.
- **Reproducibility** — pinned resources, checksummed scientific inputs,
  versioned models, and manifests.

`sglseti` still deliberately ends at generic target products. It does not
query archives, interpret detector footprints, run signal searches, manage
candidates, or claim coverage. Those archive- and instrument-specific
responsibilities begin here.

## 2. Goals

1. Determine which existing public archival observations intersect
   uncertainty-aware SGL relay corridors in space and time, or pass close
   to a hypothesized beam axis, for a curated set of nearby endpoints.
2. Search the usable intersecting data for SGL-consistent signatures,
   especially sources that recur along the predicted parallax,
   anti-target-motion, and role-dependent track.
3. Measure the end-to-end sensitivity and false-positive behavior of each
   search through synthetic-source injection, positive controls, and
   matched control corridors.
4. Publish reproducible per-target, per-archive reports with
   completeness-qualified constraints, for example: "persistent point
   sources following model M were recovered with at least 90% probability
   above flux F over these relay-distance intervals and epochs."
5. Emit stable observation, intersection, analysis-run, constraint, and
   candidate records suitable for a future coverage ledger without
   building the ledger service in this project.

Non-goals for now: new telescope observations, a public
coverage-ledger service, candidate follow-up campaigns beyond basic
vetting.

## 3. Scientific scope and methodology

### 3.1 Baseline hypothesis and search grid

The baseline local-artifact hypothesis is a compact, actively station-kept
relay on the Sun's target-star focal line at 550–10,000 AU. Rx and Tx are
separate hypotheses. The survey does not assume that every relay is
continuously visible or that every signal class follows the same emission
model.

Every analysis must freeze and identify at least:

- target endpoint and target-state model;
- Rx or Tx role;
- relay-distance prior and sampling or adaptive tolerance;
- geometry-model and observer-state versions;
- uncertainty confidence level and any additional search padding;
- source morphology and spectral model;
- persistence or duty-cycle assumption;
- allowed stationkeeping residuals or other motion parameters; and
- detection pipeline and decision threshold.

Broader hypotheses — planetary endpoints, relay swarms, off-axis
infrastructure, or inactive relics — receive separate model IDs and
constraints. A null result for one cell of this grid must not be described
as coverage of the others.

### 3.2 Coverage states and records

The pipeline distinguishes discovery, geometry, analysis, and scientific
constraint:

| Record/state            | Meaning                                                                   |
| ----------------------- | ------------------------------------------------------------------------- |
| Observation             | Immutable archive metadata and data-product identity                      |
| Intersection evaluation | Coarse candidate and exact hit/miss against one SGL hypothesis            |
| Usable intersection     | The locus crosses valid detector response after masks and quality cuts    |
| Analysis run            | A versioned pipeline actually searched specified data and parameter space |
| Constraint              | Injection-calibrated sensitivity or a qualified null result               |
| Candidate               | A retained event or track with its competing-model tests                  |

Coarse-pass misses remain useful audit records but are not scientific
coverage. A footprint hit is not coverage until the relevant data are
usable, analyzed, and assigned a measured detection efficiency.

### 3.3 Pipeline A — corridor discovery and precise intersection

For each target, role, archive time span, relay-distance range, and
uncertainty model, use `sglseti` to generate a conservative discovery
envelope. This replaces a fixed 20–30 arcmin cone with a computed region
that has a documented angular error bound.

Two-pass join:

1. **Coarse discovery.** Represent the discovery envelope as a spatial MOC
   or another conservative index and query archive metadata. Use a generic
   ObsCore/SIA/TAP adapter where the archive implements those semantics
   adequately, and archive-specific adapters otherwise. Preserve the raw
   response alongside the normalized observation record.
2. **Precise intersection.** Evaluate the adaptive or swept SGL locus over
   each observation's start-to-stop interval and intersect it with the
   detector's exact spherical WCS footprint. Account for distortion, chip
   gaps, dithers, masks, and spatially varying usable response. Return all
   intersected relay-distance intervals, which may be disjoint.

A midpoint is sufficient for short WISE frames; it is not the general
observation model. MOCs and nominal field shapes accelerate discovery but
do not replace exact WCS and valid-pixel tests.

### 3.4 Detection and motion-model screening

Use a layered search:

1. query source, reject, variability, and moving-object tables for cheap
   candidate screening;
2. perform forced photometry or source extraction at the predicted
   positions;
3. use difference imaging and trajectory-aware shift-and-stack when the
   source may lie below single-exposure thresholds;
4. fit the continuous relay distance and bounded residual motion rather
   than testing only a fixed grid;
5. compare the SGL hypothesis with inertial-background, ordinary stellar
   parallax/proper-motion, Keplerian Solar-System orbit, and instrumental
   models; and
6. require successful prediction of held-out epochs for a persistent
   candidate.

A real source may be absent from a static catalog, split into multiple
catalog entries, present as a high-motion or poor-fit source, or retained
only in a reject table. Catalog absence alone is neither a detection nor a
null result.

Every pipeline is calibrated by blind synthetic injections spanning relay
distance, flux, source morphology, duty cycle, background, and detector
position. Positive controls use known moving objects where applicable.
Negative controls use sky-rotated or time-scrambled corridors matched in
ecliptic latitude, Galactic latitude, source density, and observing
conditions. Empirical false-alarm rates must account for all target, role,
distance, residual-motion, and signal-shape trials.

### 3.5 Pipeline B — beam-axis proximity

Execution status and sequence: §11.

Run `sglseti` crossing searches over each archive's actual temporal
coverage rather than imposing a universal 1990 start. For every potentially
relevant observation, record the minimum impact parameter during the
observation, closest-approach time, uncertainty, transverse speed, and side
of the axis.

Apply effective beam radius, annular illumination, scan pattern, wavelength,
receiver response, and duty cycle afterward as explicit hypotheses.
Ingress/egress windows are convenient derived views, not universal binary
crossings.

Query both the remote target and local Rx/Tx corridors when the archive and
signal hypothesis justify them. Include Sun avoidance, daylight, field
visibility, and detector response; temporal coincidence alone does not make
an observation constraining.

### 3.6 Archives and data products of interest

| Priority | Archive or product family                                                                                                          | Rationale                                                                                                                                            |
| -------- | ---------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1        | WISE merged L1b: cryogenic WISE + NEOWISE single exposures; source/reject tables; CatWISE, unWISE, VarWISE, moving-object products | **Done (§4).** Mature all-sky multi-epoch shakedown; W1/W2 motion and hot-source sensitivity; W3/W4 carry the warm waste-heat test                   |
| 2        | SPHEREx public spectral images                                                                                                     | **Pilot active (§6).** Repeated all-sky coverage in 102 near-IR channels; spectral discrimination and a second IRSA adapter                          |
| 4        | Rubin/LSST, Pan-STARRS, DECam/NOIRLab, ATLAS, Catalina, SkyMapper, HSC/CFHT                                                        | **PS1 pilot complete (§7).** Multi-epoch optical detections and image-level shift-and-stack; Rubin is now a current rather than purely future source |
| 5        | IRAS, AKARI, Herschel, Planck, Spitzer; 2MASS and photographic plates such as POSS/DASCH                                           | Longer-wavelength thermal tests and long time baselines; some require spacecraft observers or older target solutions                                 |
| 6        | Breakthrough Listen; VLA/VLASS/COSMIC; RACS/EMU, LoTSS, GLEAM, NVSS/SUMSS; suitable ALMA/MeerKAT data                              | Local traffic, narrowband or broadband signals, and wide-field radio coverage of antipodal corridors                                                 |
| 7        | TESS FFIs, Kepler/K2, and other high-cadence imaging                                                                               | Time-domain searches of corridors and crossing events, with interval-aware spacecraft geometry                                                       |
| 8        | MPC isolated observations and tracklets, JPL/MPC known-object ephemerides, survey reject tables                                    | Candidate and veto inputs; not a complete search because 550–10,000 AU reflex rates may fall below ordinary intra-night tracklet thresholds          |
| 9        | GALEX; Chandra/XMM/eROSITA/Swift; Fermi-LAT event products                                                                         | Opportunistic UV and high-energy coincidence and persistent-source tests                                                                             |
| 10       | ESO/Keck and other spectral archives                                                                                               | Continuous or pulsed laser-line searches at the target and local corridor                                                                            |

## 4. WISE/NEOWISE — shakedown and first survey (v1 — retired 2026-08-24, see §10 step 7; reports and scripts in git history)

The first milestone — a completeness-calibrated end-to-end run of
Pipeline A against the WISE merged L1b products — was executed
2026-08-18 → 21. All eight shakedown steps (hypothesis freeze, pilot
registry, snapshotted discovery, precise WCS+mask pass, catalog
screening, forced-photometry shift-and-stack, injection calibration,
report) passed their success criteria, and the pipeline was then run
over the full universal target list. Canonical documents:

- Hypothesis freeze v1.0 + addenda: `surveys/wise/hypotheses.md`
- Per-step status and result summaries: `surveys/wise/README.md`,
  `surveys/wise/results/*.md`
- Registry curation and corridor notes: `surveys/wise/notes/`
- Reports: `report/wise_shakedown_v1.md` (pilot, 12 endpoints),
  `report/wise_survey_v2.md` (47 endpoints / 38 corridors, 8 pc),
  `report/wise_survey_v3.md` (88 endpoints / 77 corridors, universal
  v2 / picky-network portfolio, corrected depth scale — current)
- Target selection: `targets/universal_v2.md`,
  `notes/sgl_seti_star_ranking_methods.md`,
  `notes/picky_network_hypothesis_sgl_seti.md`

**Outcome.** Registry v1.5 (88 endpoints / 77 corridors; batch 5 —
the universal v2 picky-network expansion — complete 2026-08-21).
5,632 ledger-ready Constraint records under calibration v0.2.1
(`run-1b2e86e9219e`, corrected total-flux scale); every threshold
exceedance in every batch individually vetoed by the tracked census
(`adjudicate_exceedances.py`) and no surviving candidate. W1
90%-recovery depths span ~11–16 Vega mag (corrected scale; median
13.9), set by field crowding rather than noise. Physically, the W3
nulls exclude warm (~300 K) structures ≳ ~100 km at 550 AU;
reflected-light and 12 K-equilibrium hypotheses are essentially
unconstrained by WISE.

**Learnings that shape later adapters** (details in the linked docs):

1. _Masks carry the WCS._ ~130 KB `-msk` products give exact astrometry
   and usable pixels without the image; precise passes are cheap.
2. _Parallax-phase test is the decisive veto._ A static background
   source recurs at one day-of-year window; a real relay must appear at
   both alternating parallax phases. Built into screening and stacking;
   loses power in dense corridors where both phase positions are
   occupied.
3. _Confusion, not noise, sets depth_ at 6" resolution. Control
   trajectories (8 offsets) define thresholds; per-pixel uncertainties
   understate the real floor.
4. _Effective-epoch floor._ One frame can dominate a stack
   (lacaille-8760 S=398/399); calibration v0.2.0 caps per-epoch weight
   at 20× median and reports N_eff.
5. _z-grid must be uniform in 1/z_, not log z — a log grid left 17"
   gaps at small relay distance.
6. _Endpoint curation dominates effort._ Published Hipparcos/hip2
   solutions for Sirius and Procyon are already barycentric; Gaia
   component solutions for tight binaries (GJ 65, RUWE ~11) are
   orbit-corrupted; propagate CNS5 epochs before cross-matching;
   always Kepler-check extracted orbits.
7. _Registry is survey-agnostic._ Endpoint solutions, corridor grouping
   (`wise_corridors.py`) and the universal list are reused verbatim by
   later archives; only per-survey overlays change ordering.
8. _Erratum (2026-08-20)._ The shared matched filter returned PSF peak
   amplitude, not total flux; WISE depth labels are optimistic by
   1.9 mag (W1–W3) / 3.3 mag (W4). Detection statistics are unaffected.
   See `surveys/wise/results/erratum_flux_scale_2026-08-20.md`.

**DONE (2026-08-21):** after WISE batch 5 (universal v2 portfolio:
88 endpoints / 77 corridors, `surveys/wise/results/batch5_summary.md`),
all 176 sample tensors were rebuilt under the total-flux kernel and
recalibrated as v0.2.1 `run-1b2e86e9219e`; the exceedance census was
rerun uniformly over all batches by the new
`adjudicate_exceedances.py`; `report/wise_survey_v3.md` (supersedes v2)
carries corrected depths (survey-wide Rx medians W1 13.9 / W2 13.0 /
W3 10.6 / W4 7.1; best W1 15.9), the regenerated figure, and rescaled
physical limits (300 K radiator ≥ ~99 km at 550 AU via W3); the
published artifact is updated. Lesson 9: the AnalysisRun config must
include the photometry convention — the first recalibration reproduced
the old run id because the config hash ignored the kernel.

**Scientific review (2026-08-21):** `surveys/wise/scientific_review.md`
recommends major revision before the survey is presented as a
calibrated exclusion experiment. Findings, in severity order: the
8-control threshold gives a 1/9 per-search crossing probability (≈78
expected crossings vs 70 seen — no survey-wide false-alarm control);
completeness is threshold-crossing only and never passes through the
vetoes; 326/352 W3/W4 cells are auto-vetoed by the cryo-cadence rule,
so the warm-structure exclusion is unsupported; the 99 % locus is
metadata, not propagated geometry; injections are analytic Gaussians on
sampled tensors, not image-level with the empirical PRF; grid/interval
gaps (9.2 % of prior mass), L2-vs-L∞ motion bound, single temporal
model; vetoes uncalibrated; quality masks untested; provenance claims
exceed the implementation. Response: `surveys/wise/v2_plan.md`; the
four overreaching claims are withdrawn by a v3.1 erratum before any v2
computation. **These findings apply, with the substitutions listed in
`v2_plan.md` §12, to every v1 survey** — which is why §10 redoes them
all.

## 5. Second adapter: ZTF (v1 — retired 2026-08-24, see §10 step 7)

**Why ZTF second.** The plan's second adapter exists to stress the
interfaces with a genuinely different archive and to open a physical
cell WISE left closed. ZTF does both: ground-based observer state,
CCD-quadrant footprints, three filters, seeing/airmass/moon quality
dimensions, and reference-subtracted difference images as a first-class
product; scientifically, 1" pixels remove the confusion floor, g/r/i is
the reflected-sunlight band, and 8 years × ~1,000 epochs per corridor
give ~16 parallax cycles for the phase test and enough cadence to fit
continuous relay distance + residual motion. It is also ~90% of a
Rubin adapter. Probed 2026-08-20: IRSA IBE `ztf/products/sci` returns
1,132 frames (g 483 / r 494 / i 155, 2018-06 → 2026-06) at the
Ross 128 corridor. Coverage: corridors at Dec ≳ −30°, i.e. **52 of 76
universal-list systems** (the southern stars).

**Pilot (3 corridors) — complete 2026-08-20.** Report
`report/ztf_pilot_v1.md`, stage detail `surveys/ztf/results/pilot_v1_summary.md`.
Outcome: 80 Constraint records, 0 candidates, m90 ≈ 21.2–21.5 AB;
lessons — fixed-grid CCD gaps (all three antipodes affected; ≈13 % of
sky), reference-image edge strips (hybrid search image), per-cell epoch
floor + single-epoch clip, and the asteroid control's flux-scale catch.
The steps as executed:

1. Freeze `surveys/ztf/hypotheses.md`: same prior, roles, duty cycle
   and |µ_resid| bound as WISE v1.0; add filter set, Palomar observer
   site, reflected-light and hot-source interpretations, per-exposure
   quality cuts, and an explicit statement of what a g/r limit means
   physically.
2. Choose three corridors from the in-footprint set spanning low/mid
   Galactic latitude and both hemispheres of the ecliptic; reuse the
   existing registry entries unchanged.
3. Adapter `sglsurvey/adapters/irsa_ztf.py`: discovery via IBE
   `products/sci` metadata (snapshotted), nominal footprint from the
   quadrant CD matrix, exact footprint from the science-image WCS +
   mask product, cutouts via IBE `?center=&size=`.
4. Precise pass, then screening against `ztf_objects_dr24` (and the
   per-frame PSF source lists) with the parallax-phase test.
5. Forced photometry on difference-image cutouts along the track;
   joint z × µ stack with phase-split and offset controls; injection
   calibration on a 1/z-uniform grid with the effective-epoch cap.
6. Positive control: recover a known main-belt asteroid through the
   same shift-and-stack machinery (§3.4 has lacked a positive control).
7. Report `report/ztf_pilot_v1.md`; then scale to all 52 in-footprint
   systems via a ZTF overlay of the universal list.

Success criteria are those of §4 (reproducible joins, declared
recovery probability over a nontrivial flux × distance region, controls
inside thresholds, constraints traceable to usable pixels), plus: the
adapter interface in `sglsurvey/adapters/base.py` survives without
WISE-specific leakage, and the ground-based observer state validates
against sglseti's terrestrial-site model at the arcsecond level. All
met (topocentric shift ≤ 0.016″; interface unchanged, two estimator
generalisations in the shared calibration).

**Scale-up v1 — complete 2026-08-21.** Overlay
(`surveys/ztf/targets/overlay_v1.md`: 62 visible corridors, 29 ok /
18 edge / 15 gap) → full chain over 69 endpoints: 147,063 exposures,
2,640 Constraint records, 24 exceedances, 0 surviving candidates after
stage-7 adjudication and the survey-wide leave-one-out look-elsewhere
test (24 observed vs 41 expected chance exceedances). m90 ≈ 22.5 AB
median (g/r), 23.6 best. Report `report/ztf_survey_v1.md`; summary
`surveys/ztf/results/scaleup_v1_summary.md`. Deferred to v2: residual-
aware variance / ZOGY score images (≈1.5 mag headroom), per-corridor
references for edge strips, fainter second control. Gap-graded (15) +
invisible (15) corridors hand off to SPHEREx.

## 6. Third adapter: SPHEREx (v1 — retired 2026-08-24, see §10 step 7)

SPHEREx quick-release data are served at IRSA (`spherex.obscore` TAP
view, collections `spherex_qr2` and `spherex_qr2_deep`; QR1 rows are
superseded by the QR2 reprocessing). Probed 2026-08-20: ~300–650
Level-2 spectral images per sky position for 2025-05-24 → 2026-08-11
outside the deep fields, ~7,000 inside them. It is the right third
adapter rather than second: its 0.75–5 µm, 6.15"-pixel regime largely
overlaps W1/W2 physically, and 14 months yield only ~2–3 parallax
phases, but it adds (a) per-pixel wavelength (linear-variable filter,
R ≈ 40–130) and hence spectral discrimination of any candidate against
a stellar SED, (b) coverage of the 15 southern-corridor systems ZTF
cannot reach (Dec < −28°), and (c) a third, spectral-image-shaped data
product (multi-extension IMAGE/FLAGS/VARIANCE/ZODI/PSF-cube/WAVE-table
files with SIP WCS, surface-brightness units and a space-based
observer) to test the adapter interface. Running the SPHEREx pilot in
parallel with the ZTF scale-up and the last WISE batch is deliberate:
both of those are limited by IRSA download throughput, and SPHEREx
products are additionally mirrored in a public S3 bucket
(`nasa-irsa-spherex`) that supports byte-range reads at ~20 MB/s, so
the adapter fetches cutouts there and leaves IRSA's bandwidth to the
other two surveys.

**Pilot (3 stars) — complete 2026-08-20.** Report
`report/spherex_pilot_v1.md`, stage detail
`surveys/spherex/results/pilot_v1_summary.md`. Outcome: 288 Constraint
records, 0 candidates (8 exceedances, all vetoed), m90 ≈ 19.3–20.8 AB
in all six detectors; flux scale verified on 13,657 star measurements
to ±0.2 mag. Lessons — sub-pixel-phase matched filter for the
under-sampled PSF (0.4 mag), static-sky template in place of
difference images (~3 mag), detector-dependent deep-field seasons,
Rx/Tx role-coincidence veto. Sub-project `surveys/spherex/`,
adapter `sglsurvey/adapters/irsa_spherex.py`, hypotheses
`surveys/spherex/hypotheses.md` v1.0 (WISE v1.0 physics, Earth-centre
observer with the ≤ 0.02" LEO offset carried as a budget term, L2
FLAGS fatal template, per-detector bands D1–D6 with the per-epoch
wavelength carried on every sample). Corridors chosen from the
ZTF-inaccessible set, reusing registry entries unchanged:

| endpoint        | why                                                                                                                                                                                              |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `lalande-21185` | nearest star ZTF cannot reach (2.55 pc); corridor at Dec −36°, b −68° (clean field); 349 exposures                                                                                               |
| `gj-687`        | corridor 3° from the south ecliptic pole, inside the SPHEREx deep field: 7,308 exposures, continuous cadence — the parallax-phase and spectral-sampling stress case; engineering-backbone basket |
| `sigma-dra`     | engineering-backbone top pick (5.8 pc, quiet G9V); corridor at Dec −70° outside the deep field; 643 exposures                                                                                    |

Steps mirror §5 with SPHEREx specifics: (1) freeze hypotheses;
(2) TAP ObsCore discovery snapshotted, nominal footprint from the
ObsCore polygon; (3) precise pass on S3 range-read cutouts (IMAGE,
FLAGS, VARIANCE, ZODI, nearest PSF plane, wavelength table) against the
SIP WCS + FLAGS usable-pixel footprint; (4) layer-1 screening of bright
single-epoch sources along the track against CatWISE2020 (SPHEREx QR2
publishes no source catalogs); (5) matched-filter forced photometry in
surface-brightness units converted to µJy via the per-exposure pixel
solid angle, joint z × µ stack with phase-split and offset controls,
injection calibration per detector; (6) positive control: flux-scale
recovery of catalogued 2MASS/CatWISE stars in the same cutouts through
the same estimator (the lesson of the ZTF asteroid control);
(7) report `report/spherex_pilot_v1.md`. **Scale-up complete
2026-08-20:** all 15 ZTF-inaccessible corridors / 19 endpoints,
0 candidates (`report/spherex_survey_v1.md`). **v2 estimator
2026-08-21** (`calib_v3`): epoch variances calibrated to the static
template's residual scatter (+0.44 mag median) and a joint six-
detector flat-spectrum stack (+0.2–0.9 mag): 2,128 Constraints, joint
m90 21.2–21.9 AB on clean corridors, 19.5–19.8 AB in the Galactic
plane; 37 exceedances all vetoed under five rules (phase split, single
season, role coincidence, template coverage, bright static neighbour).
**v3 estimator 2026-08-21** (`calib_v4`): 16 offset controls (FAR
< 1/16 at no median depth cost) and a cross-detector season-complement
phase test; 21 exceedances, all vetoed. **All-sky run complete 2026-08-21**
(`report/spherex_survey_v1.md`, `surveys/spherex/results/allsky_v1_summary.md`):
77 corridors / 88 endpoints, 60,426 exposures, **9,808 Constraints,
0 Candidates** (64 exceedances of 1,226 searches, all vetoed); joint
m90 median 20.75 AB, 21.5–21.8 on clean corridors, 16.9–17.9 in the
Galactic plane. SPHEREx is now the first archive to cover the full
universal list. Next: re-run per quick release; candidate SED
discriminator only if something survives.

Success criteria: those of §4/§5, plus: the adapter interface absorbs
a multi-extension spectral product and a wavelength-per-sample axis
without changes to `base.py`; the flux scale is verified against
catalogued stars to ≤ 0.2 mag before any depth is quoted.

## 7. Fourth adapter: Pan-STARRS1 (v1 — retired 2026-08-24, see §10 step 7; first non-IRSA archive)

Started 2026-08-20 while IRSA throughput was saturated by the ZTF
scale-up and the SPHEREx pilot: PS1 DR2 warps (2009–2014, grizy,
δ > −30°) are served entirely by MAST (`ps1filenames.py` listing,
`fitscut.cgi` cutouts, catalogs API), so the pilot ran in parallel
without touching IRSA. Deliberately the **same three corridors as the
ZTF pilot** (Ross 128, ε Ind A/B, Proxima) so PS1 extends each ZTF
corridor by a 5–10-year baseline. Sub-project `surveys/panstarrs/`,
adapter `sglsurvey/adapters/mast_ps1.py`, hypotheses
`surveys/panstarrs/hypotheses.md` v1.0 (WISE/ZTF v1.0 physics,
Haleakalā observer, IPP mask template 16255, star-calibrated flux
scale), report `report/ps1_pilot_v1.md`.

Result: 1,351 warps → 1,355 usable precise evaluations (617 warps,
~75 % of geometric hits usable — no ZTF-style grid-gap losses) →
6,040 DR2 detection matches, nothing track-following → 320 Constraints
(m90 ≈ 21.0–21.4 AB g/r/i, 19.9 z, 19.0 y), 6 exceedances all vetoed
single-phase, **0 surviving candidates**; asteroid (60000) recovered at
0″ with magnitudes matching Horizons + solar colours to ≤ 0.1 mag.

Lessons that feed the next stages:

1. _Per-frame star-calibrated zero point_ (DR2 `mean` stars through the
   identical matched filter) replaces header ZP + a-posteriori
   throughput correction: it measured FPA.ZP + 2.5 log EXPTIME − 0.5 mag
   per frame to 0.07 mag and is verified by the asteroid control. Adopt
   as the default flux-scale method wherever a per-epoch star catalog
   exists.
2. _Parallax phase is cross-archive._ 3π revisits each field at the
   same season (97:3 phase split on every corridor), so PS1 alone cannot
   run the static-background veto; the stage-2 likelihood must combine
   archives at both phases (PS1 + ZTF on these corridors).
3. Full skycell masks (3.3 MB fpack) give exact WCS + usability, as the
   WISE `-msk` trick; `CONV.BAD` marks resampled OTA-gap bands with
   finite image values — masks, not NaNs, define usability.
4. `detection.obsTime` trails warp `MJD-OBS` by ~50–60 s; ~15 % of warp
   epochs have no catalogued detections.

**Scale-up (2026-08-21, complete):** overlay
`surveys/panstarrs/targets/overlay_v1.md` (62 of 77 corridors δ > −30°;
minor parallax phase ≤ 12 % everywhere; 30 calibrator-sparse), batched
driver `run_scaleup.sh` with per-batch purge (≈ 150 GB transient, 18 GB
retained), lazy calibration cutouts (97.6 % of 5,516 flux maps
star-calibrated). AnalysisRun `run-eaa6d89a1ec9`: 24,852 warps →
16,527 usable evaluations → 134,874 DR2 matches → **5,520 Constraints**
(median m90 g 21.1 / r 20.9 / i 20.7 / z 19.9 / y 18.9 AB) → 78
exceedances (11 %, the chance rate by construction), 72 phase-vetoed, 3
adjudicated-vetoed (two catalogued stars, one non-persistent), **3
marginal** (82 Eri rx y, Fomalhaut rx y, GJ 526 tx z: at-threshold /
grid-edge / single-filter, no catalogued counterpart) carried as
qualified nulls for the ZTF cross-archive test. Report
`report/ps1_survey_v1.md`, summary
`surveys/panstarrs/results/scaleup_v1_summary.md`. New stage-7 tool
`adjudicate_candidates.py` (split-half, static-star, other-band,
grid-edge tests) is reusable by every adapter.

**Stage 2, joint PS1 + ZTF (2026-08-21, complete):** sub-project
`surveys/joint/` (`joint_ps1_ztf.py`, `marginal_ztf_test.py`), report
`report/joint_ps1_ztf_v1.md`, AnalysisRun `run-fd75b2c982c2`. Joint
weighted stack of the per-archive tensors on the exact common
trajectory (station-kept relay, µ_resid = 0 — the tensors' µ reference
epochs differ, so µ ≠ 0 stays per-archive), ZTF on the PS1
star-calibrated scale, 8 shared controls, common phase reference: 138
endpoint-roles / 414 cells, **321 with both parallax phases** (PS1
alone ≈ 0), 3,312 Constraints with m90 ≈ 23.3 (g, r) / 21.5 (i) AB over
2009–2026, 35 exceedances (8.5 %, below chance) all vetoed (14
single-phase, 16 phase-split, 3 non-persistent, 2 other-band, 2
catalogued stars on i-band tracks at stage 7). The three PS1 marginal
cells vetoed by direct ZTF forced photometry along the PS1-fitted
(z, µ) tracks (S = 0.2 / 1.4 / −0.1 over 900–1,200 frames). **PS1
survey and joint stage close with 0 candidates.** Lessons: ZTF i is
too sparse/clustered to carry the phase test; the catalogued-star test
belongs in the automatic rules; a common µ reference epoch is needed
for any joint µ-grid analysis.

**Stage 2 v2 (2026-08-21, complete):** PS1 tensors rebuilt on the
common µ reference epoch T0 = 59800 (`run_common_t0.sh`,
`runs/panstarrs/calib_t0_59800`), joint stack over the full 5 × 5 µ
grid (`joint_ps1_ztf_mugrid.py`, AnalysisRun `run-ac08543c5b29`):
3,312 Constraints, on-grid m90 ≈ 23.0 (g) / 22.9 (r) / 21.1 (i), 43
exceedances (10 %, below chance), 4 retained then vetoed at stage 7
(field-wide systematic; faint stars along the track; i-only pair
absent in ZTF g+r) — **still 0 candidates**. Lesson: at T0 = 59800 the
PS1 epochs are 7–13 yr from the reference, so the 0.5″/yr µ step
under-samples the family there (off-grid-marginalised i depth
collapses; g/r lose 0.45 mag) — v3 needs a mid-baseline T0 with
≈ 0.1″/yr µ sampling or analytic interpolation. The
catalogued-static-source test is now automatic in every adapter's
census (`sglsurvey/vetting.py`: PS1 DR2 / ZTF DR24 snapshots offline,
CatWISE via VizieR for WISE; ≥ 3 catalog detections required).

**TODO — three-archive and southern extension of the joint stage:**

- _WISE into the joint stage._ PS1 + ZTF close only the reflected /
  self-luminous optical cell; the thermal cell is WISE's, and the WISE
  batches already cover these corridors at both parallax phases
  (NEOWISE 6-month cadence). Adding WISE to `surveys/joint` gives a
  three-archive test with a colour axis (a real relay must be
  consistent from 0.5 to 22 µm) and lets the phase veto and the
  W1:W2 static-source veto act on the same records. Needs: WISE
  tensors on a common µ reference epoch (see above), Vega→AB and
  surface-brightness conventions reconciled in `vetting`/`joint`.
- _Southern corridors._ 15 of 77 universal-list corridors lie below
  δ = −30° (Lalande 21185, Ross 248, 61 Cyg, Struve 2398, Groombridge
  34, GJ 1221, GJ 338, GJ 625, GJ 687, GJ 251, σ Dra, HD 219134,
  Wolf 1069, GJ 3512, GJ 13157), unreachable by PS1 and ZTF, so the
  optical cell there is currently closed by nothing — and several
  top engineering-basket targets (σ Dra, HD 219134, Lalande 21185)
  are among them. Options: SPHEREx (piloted, all-sky, 0.75–5 µm) as
  the in-hand near-IR route; DECam/NOIRLab Astro Data Lab (reachable,
  unprobed) as the only 1″ multi-epoch optical equivalent of PS1/ZTF.
  First step: a DECam recon + 3-corridor pilot on the ZTF/PS1 pattern
  (Lalande 21185, σ Dra, HD 219134), then a southern overlay.

## 8. Open questions and research directions

- **Endpoint hypotheses.** A network may aim at a stellar component,
  system barycenter, planet, or orbital acquisition region. Determine which
  hypotheses are physically motivated for each nearby system and never use
  component astrometry as an unlabeled barycenter substitute.
- **Emission models.** Translate flux limits into constraints separately
  for reflected sunlight, equilibrium and actively heated thermal
  emission, beacons, leakage, and pulsed or intermittent sources. Choose
  observing bands from those models rather than labeling all infrared data
  as equivalent waste-heat coverage.
- **Stationkeeping and inactive relays.** Quantify plausible deviations
  from an ideal focal line, secular drift, and failure-state orbits.
  Determine residual-motion priors that are broad enough to be physical
  without making the search statistically unconstrained.
- **Model accuracy budgets.** Allocate angular error among target-state
  propagation, observer ephemerides, relativistic geometry, timing, and
  numerical locus sampling. Validate that budget at the oldest epochs and
  smallest relay distances before assigning coverage.
- **Footprint representation.** Establish a common discovery
  representation, probably MOCs, while retaining exact spherical WCS,
  detector masks, chip gaps, dithers, and spatial sensitivity for the
  precise pass.
- **Detection access at scale.** Coarse metadata queries are cheap;
  image-level searches and epoch photometry are not. Benchmark bulk
  downloads, TAP/SIA/IBE services, cloud-hosted collections, and archive
  compute before scaling the target list.
- **Motion-model comparison.** Build calibrated likelihoods for the SGL
  track, a static or ordinary stellar source, a Keplerian Solar-System
  object, and detector artifacts. Test whether continuous distance and
  residual-motion fits are identifiable at each archive's cadence.
- **Statistical trials.** Define the global search family before candidate
  selection and propagate trials over targets, roles, relay distance,
  motion residuals, epochs, bands, source shapes, and persistence models.
- **Target prioritization.** Score systems by endpoint quality, distance,
  multiplicity, physical interest, astrometric history, corridor
  background, and actual calibrated archive coverage. Do not let proximity
  alone determine the order.
- **Crossing physics.** Replace generic ingress/egress language with
  explicit transmitter and receiver geometries, beam profiles, scan
  strategies, wavelength dependence, and duty cycles. Determine when
  remote-target versus local-corridor data can constrain each case.
- **Radio scope.** Decide whether this repository performs voltage or
  spectrogram searches, delegates them to archive-specific pipelines, or
  publishes only geometry and coverage products. In every case, preserve
  frequency, drift-rate, polarization, time, and sensitivity dimensions.
- **Dark infrastructure.** Explore constraints from occultations,
  microlensing, reflected-light phase behavior, thermal emission at longer
  wavelengths, and gravitational or dynamical effects for objects that do
  not transmit.
- **Population inference.** Develop a hierarchical framework that combines
  heterogeneous, model-specific completeness curves. A count of archive
  intersections is not a population constraint.

## 9. Practical notes

- **Repo layout:** retain `notes/` for research and planning,
  `sglsurvey/` (or similar) for the pipeline package, `registries/` for
  curated target registries, and `runs/` or `build/` for generated products
  that are ignored except for compact, reviewable summaries.
- **Core records:** model immutable `Observation`,
  `IntersectionEvaluation`, `AnalysisRun`, `Constraint`, and `Candidate`
  records separately. Do not overload one table with discovery, analysis,
  and interpretation states.
- **Archive adapters:** normalize identity, timing, footprint, band,
  calibration, data-quality, and product-location fields, but preserve raw
  archive metadata. Keep discovery, download, and instrument-specific
  valid-pixel logic behind explicit adapter interfaces.
- **sglseti dependency:** during development, pin the adjacent checkout
  (`../sglseti`) to an exact commit and record whether it was dirty. Every
  run also pins target registry, geometry model, observer state, kernels,
  and the model IDs introduced by the improvements roadmap.
- **Reproducibility:** snapshot archive queries and responses, checksums,
  calibration files, software environments, configuration, random seeds,
  and manifests. Archive services and catalogs evolve even when survey code
  does not.
- **Data retention:** keep lightweight metadata and derived measurements
  permanently; record content hashes and durable archive identifiers for
  large images or event files; retain local cutouts and injected products
  according to a documented regeneration policy.
- **Implementation sequence:** ~~complete the required `sglseti` gates,
  implement a WISE adapter and the five core records, run the 3–5-endpoint
  pilot through injection recovery and reporting~~ (done, §4); ~~ZTF as the
  second adapter via a 3-corridor pilot~~ (done, §5; scale-up running);
  ~~SPHEREx third via a 3-star pilot~~ (done, §6; all 15
  ZTF-inaccessible corridors surveyed); ~~Pan-STARRS1 fourth via
  the ZTF pilot corridors, off MAST~~ (done, §7; started early because
  IRSA was saturated); ~~expand archive families only after the interfaces
  survive all four~~ → superseded by the v2 programme (§10): finish the
  in-flight v1 runs, rebuild WISE as v2, carry the fixes to the other
  surveys, retire v1, and only then expand archive families.

## 10. v2 programme — sequence and cross-survey transfer (2026-08-21)

The WISE scientific review (§4) found that the v1 surveys are sound as
exploratory searches but not as calibrated exclusion experiments, and
the same constructions — 8/16-offset thresholds, tensor-level analytic
injections, hand-tuned vetoes, nominal-locus geometry, filename-based
provenance — were copied from WISE into ZTF, SPHEREx, PS1, and the
joint stage. v2 repairs the statistical experiment once, on WISE, and
then applies the repaired design to every survey. The full WISE design
is `surveys/wise/v2_plan.md`; its §12 is the transfer table.

### 10.1 Sequence

| Step | What                                                                                                                                                                                                                                                                                                                                                                                                                               | Exit condition                                                                                                                               |
| ---- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| 1    | **Let the in-flight v1 runs finish**: SPHEREx northern run (62 corridors, frozen v3 rules, §6) and the joint PS1+ZTF common-T0 µ-grid calibration (§7). No new v1 rules or batches start.                                                                                                                                                                                                                                          | both runs summarised under `surveys/{spherex,joint}/results/` and the v1 reports closed with their own "exploratory, pending v2" status note |
| 2    | **Errata for the v1 reports.** WISE v3.1 (withdraws the four claims in `v2_plan.md` §0); the same relabelling — "FAR < 1/N" → rank statement, "90 %-complete exclusion" → threshold sensitivity, vetoes → heuristic review — for the ZTF, SPHEREx, PS1, and joint reports.                                                                                                                                                         | errata committed; no numbers recomputed                                                                                                      |
| 3    | **Shared-code modifications** (`v2_plan.md` §8.5): `records` fields, `photometry` injection seam, `vetting` flux-consistent static test + parallax-phase test, `geometry` per-epoch observer, new `inject`, `nulls` (promoted from ZTF `look_elsewhere.py`), `manifest`; promote `wise_corridors.py` out of `surveys/wise/` (ZTF/PS1/SPHEREx import it). All additive; v1 scripts keep running. `tests/` with the invariant tests. | tests green against v1 products; v1 pipelines unchanged in output                                                                            |
| 4    | **WISE v2 in `surveys/wise/`** per `v2_plan.md` §10 (steps A–I): hypothesis v2.0 freeze with hold-out split → covariance + positive control → null ensemble → image-level injections → frozen rule on the development set → blind confirmatory run → `report/wise_survey_v4.md`.                                                                                                                                                   | v4 report built from the ledger; review findings 1–9 each closed or explicitly deferred                                                      |
| 5    | **Retrospective**: update `v2_plan.md` §12 with what actually transferred, what cost more than planned, and any shared-code change forced by the WISE run.                                                                                                                                                                                                                                                                         | §12 revised                                                                                                                                  |

**Status 2026-08-22:** steps 1–5 done. Errata appended to all five v1 reports (WISE v3.1); shared code in `sglsurvey/` (`nulls`, `inject`, `manifest`, `corridors`, additive changes to `records`, `photometry`, `vetting`, `geometry`) with `tests/test_v2_invariants.py` green; WISE v2 in `surveys/wise/` ran A–I under hypotheses v2.1 — 0 candidates on both the development (216 cells) and the blind confirmatory (488 cells) sets, `report/wise_survey_v4.md`; retrospective in `v2_plan.md` §12.1 (ring-only null, per-cell normalisation with heavy-tail exclusion, held-out test demoted to an annotation, ≈ 1.5 mag depth cost of family-wise error control, W4 pixel-scale erratum). Step 6 done 2026-08-23 via the survey-agnostic engine `sglsurvey/v2/` (profiles in `surveys/{ztf,panstarrs,spherex}-v2/`, joint in `surveys/joint/`): blind confirmatory runs ZTF 228 cells / PS1 450 / joint 230 / SPHEREx (see `report/spherex_survey_v2.md`) — 0 candidates everywhere; reports `report/{ztf,ps1,joint_ps1_ztf,spherex}_survey_v2.md`; lessons in `v2_plan.md` §12.2 (static-source double counting in difference/template images, the single-epoch clip as a bright completeness limit, full-WCS sampling for SPHEREx). Step 7 next.
| 6 | **Other surveys v2**, in `surveys/{ztf,panstarrs,spherex,joint}-v2/`, in the order ZTF → PS1 → joint → SPHEREx (ZTF has the richest cadence and an existing positive control, so it validates the transfer fastest; the joint stage needs both optical v2s; SPHEREx last because its next quick release adds a parallax phase anyway). Each is a scaled-down `v2_plan`: reuse v1 discovery/precise/screen and cutouts, rebuild only null, injection, vetting, geometry check, manifests, report. | each v2 report built from its ledger |
| 7 | **Retire v1.** Delete `surveys/{wise,ztf,panstarrs,spherex,joint}/` and `runs/<survey>/` products no longer referenced; rename each `*-v2` to the plain name; move v1 reports to `report/archive/`; remove stale references (README, this plan, `targets/` overlay builders, memory notes) and dead code paths in `sglsurvey/` kept only for v1 compatibility. v1 stays in git history. | `grep -r "wise-v2\|_v1\|v1_" ` clean except history notes; tests green |

**Step 7 done 2026-08-24** (user amendment: v1 reports deleted rather than archived, and the final reports renamed to plain names — `wise_survey_v4.md` → `wise_survey.md`, etc.). What was done, and the deliberate deviations from the row above:

- _Reports:_ only the five canonical reports remain (`report/README.md`); all pilots/v1–v3 reports and `report/figures/` deleted (git history).
- _Surveys:_ each `*-v2` renamed to the plain name. Kept inside the renamed dirs because they are still current: the discovery/precise/screen/fetch scripts (Pipeline A stages that v2 reuses and step 8 needs), corridor/targets overlay builders, per-archive `notes/`, `surveys/wise/{v2_plan.md,scientific_review.md}`, PS1's `purge_products.py` and SPHEREx's `static_template.py` (generators of kept products). Deleted: the standalone v1 hypothesis files after their still-active physics cells were incorporated directly into the v2 hypothesis documents; all v1 analysis scripts (sample_tensor, injection/forced-stack calibration, adjudicators, look_elsewhere, figures); and v1 `results/` summaries.
- _Runs:_ v1 tensor products and superseded calibration runs deleted (~73 GB: all `calib_*/tensors`, `stack_v1`, SPHEREx `calib_v1–v3`); kept because still referenced — coarse/precise records, cutout products (rebuild inputs), screening records/snapshots (catalogue loaders + the ScreenMatch ledger), `m90_curves.npz` + `records/` of the final v1 calibrations (injection windows + supersession targets), PS1 `zeropoints.jsonl`, SPHEREx `calib_v4/templates`, ZTF `control_v1`. v2 products moved to `runs/<survey>/v2/`. SPHEREx supersession links re-derived against `calib_v4` (they had pointed at the deleted pilot-era `calib_v1` ledger; `supersedes` is outside the constraint identity hash).
- _Code:_ v1-only vetting retired (`static_source_test` proximity veto and helpers); all path constants updated; every survey's ledger-driven report tables rebuild cleanly; tests green (the one v1-tensor regression test now skips, its input being deleted).
- _Exit-criterion note:_ the stage directories keep their `coarse_v1` / `precise_v1` / `screen_v1` names — they are the still-current first versions of live stages that v2 reuses by design, not superseded analyses.
- _Provenance note:_ the freeze files (`surveys/*/configs/v2_*freeze.json`) are kept byte-identical to the pre-registered versions — every AnalysisRun pins their hashes. The path-rename sweep and the 2026-08-24 documentation restoration edited the hypothesis documents after those freezes, so the `hypotheses_hash` values recorded in the freezes refer to the documents as of freeze time (git history), not the current self-contained files; future freezes hash the complete active hypothesis document.
  | 8 | Resume archive-family expansion (DECam south, three-archive joint) on the v2 design. | — |

**Step 8 started 2026-08-24 — DECam south.** Recon
(`notes/decam_recon_2026-08-24.md`): all 15 southern corridors have
instcal coverage (~15.6k exposures 2012–2026, 4–12 calendar months →
both parallax phases single-archive, unlike PS1); Astro Archive API +
`?hdus=` single-CCD fetches verified; `astro-datalab` (pinned 2.22.1;
2.24.0 has undeclared imports) queryClient works anonymously for
NSC DR2 (async/MyDB need an account — deferred to scale-up); Data Lab
TAP unusable both probe days; Data Lab cutout service strips TPV → not
used for astrometry. Built: adapter
`sglsurvey/adapters/noirlab_decam.py` (EXPNUM-joined triplets, dqmask
md5-verified, per-exposure EXTNAME→HDU map with fetch-time assert),
`GeometryContext.decam_default()` (CTIO W84), static focal-plane
layout `surveys/decam/configs/decam_focal_plane_v1.json` (validated vs
2012/2024 exposures, worst corner 10.2"), `surveys/decam/` with draft
hypotheses v0.1 (freeze pending — first survey with no v1 exploratory
phase, v2 discipline from day one). Pilot Pipeline A run (lalande /
sigmadra / hd219134): coarse 306 exposures / 612 evaluations; precise
222 hit exposures, 412 evaluations, ~all usable (5 i-band 404s);
screen 11,430 ScreenMatch (NSC DR2 confirmed time-partial: catalogued
epochs end 2017–2019 per corridor — recorded in
`runs/decam/screen_v1/catalog_stats.json`); recurrence triage clean
(dense hd219134 corridor is Galactic-plane, b ≈ +2.5°); single-CCD
image+wtmap products fetched for stage 2. Next: hypothesis freeze
v1.0 (user), positive-control asteroid, star flux-scale check, then
the v2 chain and the southern overlay.

Rules during the programme: no v1 result is cited as an exclusion;
nothing in a `-v2` directory may import from a v1 survey directory
(shared code goes through `sglsurvey/`); every v2 hypothesis freeze is
a new version with a declared hold-out before any v2 script touches
the confirmatory set.

### 10.2 What each v2 survey must change (from the WISE learnings)

Common to all (the review's findings 1, 2, 4, 5, 6, 9):

- **Null ensemble and global error rate** replace "T = max of N
  controls": extended spatial offsets + time scrambling + trajectory
  randomisation per cell; leave-one-out exceedance ratio R; survey-wide
  max-R null → FWER ≤ 0.05 candidate threshold plus BH q-values; per-cell
  rank statements with binomial intervals. The joint stage calibrates
  one FWER across archives.
- **Image-level injections** into the retained cutouts before the
  matched filter, with the archive's empirical PSF and a physical
  spectrum; continuous z (log-uniform), µ (L∞ box), cross-track (from
  the propagated envelope), four temporal models; ≥ 400 per cell;
  threshold and final-candidate completeness reported separately with
  Wilson intervals; gap-free reciprocal-distance intervals.
- **Calibrated vetoes vs annotations**: only tests with independent
  evidence and an injection-measured selection function may reject
  (flux-consistent catalogued-static / halo test; held-out-epoch
  prediction; a physical cross-band consistency test where a second band
  exists). Phase split, season, role coincidence, significance-ratio
  rules become annotations. Candidates that survive are reported as
  `retained-ambiguous`, not vetoed by discretion.
- **Covariance propagation** with an empirical envelope-vs-PSF check; a
  cross-track tensor dimension where the 99 % envelope exceeds half a
  PSF. Optical PSFs (1–2″) make this test _stricter_ than for WISE —
  expect more PS1/ZTF endpoints to need the extra dimension.
- **Hypothesis freeze with hold-out**: a random endpoint split
  stratified by confusion class (development vs confirmatory), seed in
  the freeze; epoch hold-out pre-registered where new epochs keep
  arriving (ZTF, SPHEREx) and calibrated post-hoc where the mission is
  over (WISE, PS1 — see `v2_plan.md` §1.8); FWER α = 0.05 as the
  candidate rule with BH q-values informational; the L∞ motion bound
  stated; the distance prior distinguished from the grid.
- **Positive control** through the full chain (asteroid; ZTF has one,
  PS1 and SPHEREx need one).
- **Quality-mask sensitivity** (primary / strict / loose) reported as a
  systematic.
- **Content-hashed manifests**, mandatory stale-product detection,
  invariant tests, report tables generated from the ledger.

Survey-specific:

- _ZTF_: `look_elsewhere.py` becomes the shared `nulls` module; nightly
  cadence supports all four temporal models; hold out the last
  observing year for the prediction test; per-quadrant PSF from the
  sci header; diff-image injections must go into the _science_ image
  before differencing, not into the diff.
- _PS1_: cutouts are purged after tensoring (`purge_products.py`), so
  image-level injection re-fetches per batch as `run_common_t0.sh`
  does; single-phase cadence means visit-scale ≈ persistent and the
  phase test cannot carry weight — the joint stage supplies phases;
  per-skycell PSF; star-calibrated ZP already in place.
- _SPHEREx_: 16 controls → same 1/17 problem; PSF cube already in the
  product; inject an SED, not a magnitude (per-pixel wavelength);
  template-coverage and bright-static-neighbour rules become
  annotations / flux-consistent vetoes; northern-vs-southern run is the
  natural hold-out split; re-run when QR3 adds a phase.
- _Joint_: inject the same physical source into all archives; colour
  consistency 0.5–22 µm becomes a calibrated veto; WISE v2 tensors on
  the common T0 make the three-archive stage possible.

## 11. Pipeline B — crossings surveys (started 2026-08-23)

Execution of §3.5. The survey-independent geometry is computed once per
observer class and each archive intersects it; hypotheses (beam radius,
wavelength, duty cycle) are frozen per survey before any data are
touched, with the v2 statistical discipline (frozen thresholds,
dev/confirmatory splits, control-based exceedance budgets, injection
completeness) inherited from §10.

### 11.1 Completed

| Step | What | Outcome |
| ---- | ---- | ------- |
| 1 | **Universal crossing list** (`sglsurvey/crossings.py` → `crossings/universal_v1/`, `xng-a09e2db7681d`). Every b(t) minimum for the 88-endpoint registry × both link directions, Earth-center observer, 1980→2028; no beam radii, no b cut; boundary minima flagged. | 16,586 events, 0 invalid; 960 with b < 0.01 AU; ~semiannual minima (annual, b ≈ 1 AU for the 5 ecliptic-pole targets) |
| 2 | **ZTF crossings survey** (`surveys/ztf-crossings/`, report `report/ztf_crossings.md`). Channels: A = uplink interception at the star near opposition; B = downlink pre-lens interception at the antipode (z-track 550–10,000 AU); sunward combinations out of scope. Coverage intersect, saturation cut, threshold freeze v1.0 + dev-driven amendments v1.1 (parallax-factor systematics template, same-rung pseudo-windows) and v1.2 (empirical variance rescale), confirmatory run, stamp-response injection completeness. | **0 candidates.** 36 searchable units, 3 exceedances vs 4.0 expected control crossings, all adjudicated. Depths: B median m90 ≈ 21.8 (relay power ≳ ~130 W through the 2.5 R☉ solar cone during covered windows); A contrast-limited m90 ≈ 15.8–17.6 (10-m uplink ≳ 56–400 kW), wise-0855 empty field ≥ 22 |

Lessons carried forward (details in the report): the crossing geometry
phase-locks to the sidereal year, so wide-beam rungs have window ≈
observing season and are constraint-only by construction — design
narrow-rung-first for optical archives; blended on-star channels are
systematics-dominated (PM dipole vs multi-year references, k factors
10²–10⁴) and the template + rescale recipe is reusable; the antipode
channel is the workhorse; the 1/9 exceedance budget behaved exactly as
designed.

### 11.2 Next steps

| Step | What | Exit condition |
| ---- | ---- | -------------- |
| 3 | **PS1 crossings.** Intersect the universal list with the PS1 era (2010–2014: windows disjoint from ZTF's, extending the covered-window record backwards); reuse the MAST adapter and the ztf-crossings channel constructions; expect no i-band-only complications since PS1 grizy spans the same hypothesis range. _In progress 2026-08-24 (`surveys/ps1-crossings/`): hypotheses v1.0, coverage intersect (B wide 14/47 events on all 7 targets; grazing 1/27+1/34 but the covered event is van-maanen at b = 0.28 R☉ with a TTI i-band pair — deepest graze in the programme), grizy saturation cut, threshold freeze v1.0 (19 units, 2.1 expected control crossings, warp-direct substrate, ZTF v1.1/v1.2 rules adopted at freeze, A 1.0 AU rung constraint-only by declaration, TTI-pair mover veto; B-wide dev = teegarden + wolf-359) all frozen. Dev search complete: 0 exceedances, no amendment (5 wolf-359 units clean; teegarden unit `track_masked` — correlated CONV.BAD chip-gap attrition, finding P1; saturation gate PASS, frozen levels stand). Next: confirmatory run._ | Frozen hypotheses + thresholds, confirmatory run, completeness, report |
| 4 | **Spacecraft-observer crossing lists.** Derivative `crossings/` products with per-epoch spacecraft observers (WISE, SPHEREx, TESS; observer code already in `sglsurvey/geometry.py`) — the Earth-center list is invalid at b ≲ R☉ for L2/Earth-trailing observers. | Per-observer events tables with manifest parity to `universal_v1` |
| 5 | **WISE crossings survey** on the derivative list: IR bands dodge the 532 nm-only constraint; NEOWISE cadence (~6-month visits) vs window durations decides which rungs are viable — check before freezing. | Same chain as step 2 |
| 6 | **Joint crossings stage.** Cross-archive coincidence on shared covered windows (ZTF × PS1 × WISE), colour consistency as a calibrated veto, one ledger of covered windows per target per rung across all archives. | Unified covered-window ledger + joint report |
| 7 | **Radio-scope decision** (§2 goals): whether corridor/crossing geometry products feed Breakthrough Listen / VLASS-class searches directly or stay geometry-only. | Decision recorded here |

Universal-list caveats that gate later steps: impact-parameter
uncertainty is not propagated (`sglseti.crossing_uncertainty` exists
for per-event follow-up); the 2028 window end means yearly refresh runs
as archives extend.
