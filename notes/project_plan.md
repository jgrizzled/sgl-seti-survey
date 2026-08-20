---
title: "sgl-seti-survey — Project Plan"
date: 2026-08-20
status: "v0.3 — WISE complete, ZTF pilot starting"
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

| Priority | Archive or product family                                                                                                          | Rationale                                                                                                                                   |
| -------- | ---------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| 1        | WISE merged L1b: cryogenic WISE + NEOWISE single exposures; source/reject tables; CatWISE, unWISE, VarWISE, moving-object products | **Done (§4).** Mature all-sky multi-epoch shakedown; W1/W2 motion and hot-source sensitivity; W3/W4 carry the warm waste-heat test        |
| 2        | SPHEREx public spectral images                                                                                                     | Repeated all-sky coverage in 102 near-IR channels; spectral discrimination and a second IRSA adapter                                        |
| 4        | Rubin/LSST, Pan-STARRS, DECam/NOIRLab, ATLAS, Catalina, SkyMapper, HSC/CFHT                                                   | Multi-epoch optical detections and image-level shift-and-stack; Rubin is now a current rather than purely future source                     |
| 5        | IRAS, AKARI, Herschel, Planck, Spitzer; 2MASS and photographic plates such as POSS/DASCH                                           | Longer-wavelength thermal tests and long time baselines; some require spacecraft observers or older target solutions                        |
| 6        | Breakthrough Listen; VLA/VLASS/COSMIC; RACS/EMU, LoTSS, GLEAM, NVSS/SUMSS; suitable ALMA/MeerKAT data                              | Local traffic, narrowband or broadband signals, and wide-field radio coverage of antipodal corridors                                        |
| 7        | TESS FFIs, Kepler/K2, and other high-cadence imaging                                                                               | Time-domain searches of corridors and crossing events, with interval-aware spacecraft geometry                                              |
| 8        | MPC isolated observations and tracklets, JPL/MPC known-object ephemerides, survey reject tables                                    | Candidate and veto inputs; not a complete search because 550–10,000 AU reflex rates may fall below ordinary intra-night tracklet thresholds |
| 9        | GALEX; Chandra/XMM/eROSITA/Swift; Fermi-LAT event products                                                                         | Opportunistic UV and high-energy coincidence and persistent-source tests                                                                    |
| 10       | ESO/Keck and other spectral archives                                                                                               | Continuous or pulsed laser-line searches at the target and local corridor                                                                   |

## 4. WISE/NEOWISE — shakedown and first survey (complete)

The first milestone — a completeness-calibrated end-to-end run of
Pipeline A against the WISE merged L1b products — was executed
2026-08-18 → 20. All eight shakedown steps (hypothesis freeze, pilot
registry, snapshotted discovery, precise WCS+mask pass, catalog
screening, forced-photometry shift-and-stack, injection calibration,
report) passed their success criteria, and the pipeline was then run
over the full universal target list. Canonical documents:

- Hypothesis freeze v1.0 + addenda: `surveys/wise/hypotheses.md`
- Per-step status and result summaries: `surveys/wise/README.md`,
  `surveys/wise/results/*.md`
- Registry curation and corridor notes: `surveys/wise/notes/`
- Reports: `report/wise_shakedown_v1.md` (pilot, 12 endpoints),
  `report/wise_survey_v2.md` (47 endpoints / 38 corridors)
- Target selection: `targets/universal_v2.md`,
  `notes/sgl_seti_star_ranking_methods.md`,
  `notes/picky_network_hypothesis_sgl_seti.md`

**Outcome.** Registry v1.5 (89 endpoints / 77 corridors, batch 5 in
progress at time of writing). Through v1.4: 3,008 ledger-ready
Constraint records, 34 Candidate records — every one vetoed — and no
surviving candidate. W1 90%-recovery depths span ~13–17.5 Vega mag,
set by field crowding rather than noise. Physically, the W3 nulls
exclude warm (~300 K) structures ≳ 39 km at 550 AU; reflected-light and
12 K-equilibrium hypotheses are essentially unconstrained by WISE.

**Learnings that shape later adapters** (details in the linked docs):

1. *Masks carry the WCS.* ~130 KB `-msk` products give exact astrometry
   and usable pixels without the image; precise passes are cheap.
2. *Parallax-phase test is the decisive veto.* A static background
   source recurs at one day-of-year window; a real relay must appear at
   both alternating parallax phases. Built into screening and stacking;
   loses power in dense corridors where both phase positions are
   occupied.
3. *Confusion, not noise, sets depth* at 6" resolution. Control
   trajectories (8 offsets) define thresholds; per-pixel uncertainties
   understate the real floor.
4. *Effective-epoch floor.* One frame can dominate a stack
   (lacaille-8760 S=398/399); calibration v0.2.0 caps per-epoch weight
   at 20× median and reports N_eff.
5. *z-grid must be uniform in 1/z*, not log z — a log grid left 17"
   gaps at small relay distance.
6. *Endpoint curation dominates effort.* Published Hipparcos/hip2
   solutions for Sirius and Procyon are already barycentric; Gaia
   component solutions for tight binaries (GJ 65, RUWE ~11) are
   orbit-corrupted; propagate CNS5 epochs before cross-matching;
   always Kepler-check extracted orbits.
7. *Registry is survey-agnostic.* Endpoint solutions, corridor grouping
   (`wise_corridors.py`) and the universal list are reused verbatim by
   later archives; only per-survey overlays change ordering.
8. *Erratum (2026-08-20).* The shared matched filter returned PSF peak
   amplitude, not total flux; WISE depth labels are optimistic by
   1.9 mag (W1–W3) / 3.3 mag (W4). Detection statistics are unaffected.
   See `surveys/wise/results/erratum_flux_scale_2026-08-20.md`.

## 5. Second adapter: ZTF (active)

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

**Scale-up (next):** ZTF overlay of the universal list with primary-
grid gap position and secondary-field epoch counts per corridor; custom
per-corridor reference images for the edge strips; a fainter second
asteroid control; then the 52 in-footprint systems.

## 6. Third adapter: SPHEREx (next)

SPHEREx QR1/QR2 are served at IRSA (`spherex.obscore`, SIA collections
`spherex_qr1/qr2`; ~330 exposures per sky position for 2025-06 →
2026-06). It is the right third adapter rather than second: its
0.75–5 µm, 6.2"-pixel regime largely overlaps W1/W2 physically, and one
year yields only ~2 parallax phases, but it adds (a) 102-channel
spectral discrimination of any candidate against a stellar SED, (b)
coverage of the 24 northern-star / southern-corridor systems ZTF cannot
reach, and (c) a third, spectral-cube-shaped data product to test the
adapter interface. Plan: after the ZTF pilot report, freeze a SPHEREx
hypothesis document reusing the WISE band-specific interpretations,
start with the 24 ZTF-inaccessible systems, and grow the baseline as
quick-release data accrue.

## 7. Open questions and research directions

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

## 8. Practical notes

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
  pilot through injection recovery and reporting~~ (done, §4); ZTF as the
  second adapter via a 3-corridor pilot (§5), SPHEREx third (§6); expand
  archive families only after the interfaces survive all three.
