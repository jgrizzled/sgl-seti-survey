---
title: "sgl-seti-survey — Project Plan"
date: 2026-08-21
status: "v1.2 — Pipeline A (§4): all six surveys complete on the v2 design (WISE, ZTF, SPHEREx, PS1, joint, DECam), 0 candidates everywhere; Pipeline B (§5): programme closed 0 candidates incl. TESS (§5.6, complete 2026-08-24); expansion queue §5.8 next item: ATLAS/ASAS-SN. Execution record: notes/project_history.md; lessons: notes/learnings.md"
tags:
  - SETI
  - technosignatures
  - solar-gravitational-lens
  - archival-search
---

# sgl-seti-survey — Project Plan

Companion documents: `notes/project_history.md` (chronological
execution record — what ran, when, with what outcome) and
`notes/learnings.md` (consolidated design lessons). This plan keeps
goals, methodology, sequencing, and current status.

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

Execution status and sequence: §5.

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

| Priority | Archive or product family                                                                                                          | Rationale                                                                                                                                      |
| -------- | ---------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| 1        | WISE merged L1b: cryogenic WISE + NEOWISE single exposures; source/reject tables; CatWISE, unWISE, VarWISE, moving-object products | **Done (§4.1).** Mature all-sky multi-epoch shakedown; W1/W2 motion and hot-source sensitivity; W3/W4 carry the warm waste-heat test           |
| 2        | SPHEREx public spectral images                                                                                                     | **Done — all-sky (§4.3).** Repeated all-sky coverage in 102 near-IR channels; spectral discrimination and a second IRSA adapter                |
| 4        | Rubin/LSST, Pan-STARRS, DECam/NOIRLab, ATLAS, Catalina, SkyMapper, HSC/CFHT                                                        | **PS1 (§4.4) and DECam (§4.6) done.** Multi-epoch optical detections and image-level shift-and-stack; Rubin is now a current source (§6)       |
| 5        | IRAS, AKARI, Herschel, Planck, Spitzer; 2MASS and photographic plates such as POSS/DASCH                                           | Longer-wavelength thermal tests and long time baselines; some require spacecraft observers or older target solutions                           |
| 6        | Breakthrough Listen; VLA/VLASS/COSMIC; RACS/EMU, LoTSS, GLEAM, NVSS/SUMSS; suitable ALMA/MeerKAT data                              | **Geometry-only decision (§5.5).** Local traffic, narrowband or broadband signals, and wide-field radio coverage of antipodal corridors |
| 7        | TESS FFIs, Kepler/K2, and other high-cadence imaging                                                                               | **TESS crossings done (§5.6).** Time-domain searches of corridors and crossing events, with interval-aware spacecraft geometry                 |
| 8        | MPC isolated observations and tracklets, JPL/MPC known-object ephemerides, survey reject tables                                    | Candidate and veto inputs; not a complete search because 550–10,000 AU reflex rates may fall below ordinary intra-night tracklet thresholds    |
| 9        | GALEX; Chandra/XMM/eROSITA/Swift; Fermi-LAT event products                                                                         | Opportunistic UV and high-energy coincidence and persistent-source tests                                                                       |
| 10       | ESO/Keck and other spectral archives                                                                                               | Continuous or pulsed laser-line searches at the target and local corridor                                                                      |

2026-08-24 addendum: the crossings expansion queue (§5.8) and future
projects (§6) add archives not itemised above — ASAS-SN, PTF/iPTF,
Gattini-IR/WINTER, gPhoton photon-level GALEX, heliospheric imagers
(LASCO/STEREO-HI/WISPR), radio epoch-metadata (VAST/RACS, LoTSS),
Rubin alert stream, and Gaia DR4 epoch products.

## 4. Pipeline A surveys

Execution of §3.3–3.4 over the universal registry (v1.5: 88
endpoints / 77 corridors), one subproject per archive (joint stages
combine archives). Six surveys complete; **every blind confirmatory
run returned 0 candidates**.

The four 2026-08 v1 surveys were exploratory-grade: the WISE scientific
review (findings condensed in `notes/learnings.md` §1) found their
shared statistical constructions could not support exclusion claims,
and the **v2 program** (2026-08-21 → 24) rebuilt the experiment once
on WISE, applied it to every survey through the survey-agnostic engine
`sglsurvey/v2/` (since flattened to `sglsurvey/`), and retired v1 (reports and scripts to git history).
Sequence, retirement inventory, and execution record:
`notes/project_history.md` §5; design rules and transfer amendments:
`notes/learnings.md` §§1–3. Standing rules: no v1 result is cited as
an exclusion; every hypothesis freeze is a new version with a declared
hold-out before any script touches its confirmatory set; each survey's
frozen decision rule lives in its own `hypotheses.md`.

**Open items:** ~~the three-archive joint stage~~ — **done 2026-08-25
as joint v3.0** (§4.5: 0 candidates blind; the stack-regime asteroid
control also closed, `notes/learnings.md` §10 item 1);
SPHEREx re-run per quick release (QR3 adds a true parallax-phase
hold-out); the stack-regime positive controls and template refits
listed in `notes/learnings.md` §10; Rubin and Gaia DR4 (§6).

### 4.1 WISE/NEOWISE

**Description.** All-sky mid-infrared (W1–W4) merged L1b single
exposures, 2010–2024. First adapter and end-to-end shakedown of
Pipeline A; W1/W2 carry the reflected-light / hot-source cell, W3/W4
the warm waste-heat sensitivity (threshold statements only — no
exclusion, per the review). The registry, corridor grouping and
universal target list built here are reused verbatim by every later
survey.

**Status.** Complete — v1 2026-08-18 → 21; rebuilt as v2 2026-08-22
after the scientific review. Report `report/wise_survey.md`; narrative
`notes/project_history.md` §§1, 5.

**Result.** **0 candidates** (development 216 + blind confirmatory 488
cells, family-wise α = 0.05); W1 depths confusion-limited.

### 4.2 ZTF

**Description.** Northern 1″ optical time domain (g/r/i, 2018–,
Dec ≳ −30°: 52 of 76 universal-list systems) with difference images;
opens the reflected-sunlight cell with ~16 parallax cycles and
stresses the adapter interface with a ground-based observer and
CCD-quadrant footprints; ~90 % of a Rubin adapter. Adapter
`sglsurvey/adapters/irsa_ztf.py`.

**Status.** Complete — pilot + scale-up 2026-08-20 → 21, v2 rerun
2026-08-23. Report `report/ztf_survey.md`; narrative
`notes/project_history.md` §2.

**Result.** **0 candidates** (v2 confirmatory 228 cells); m90 ≈ 23.2
AB (zg, persistent).

### 4.3 SPHEREx

**Description.** All-sky near-infrared spectral survey (0.75–5 µm,
per-pixel wavelength — spectral discrimination of any candidate);
covers the 15 southern corridors ZTF cannot reach; third data-product
shape (multi-extension spectral images, SIP WCS, S3 byte-range
cutouts). Adapter `sglsurvey/adapters/irsa_spherex.py`.

**Status.** Complete — pilot, southern scale-up and all-sky run
2026-08-20 → 21, v2 rerun 2026-08-23; standing item: re-run per quick
release (QR3 adds a true parallax-phase hold-out). Report
`report/spherex_survey.md`; narrative `notes/project_history.md` §3.

**Result.** **0 candidates** (all-sky, 77 corridors / 88 endpoints —
the first archive to cover the full universal list; v2 rerun blind).

### 4.4 Pan-STARRS1

**Description.** PS1 DR2 warps (2009–2014, grizy, δ > −30°) via MAST —
the first non-IRSA archive — extending each ZTF corridor back by a
5–10-year baseline; single-phase 3π cadence means the parallax-phase
test needs the joint stage (§4.5). Adapter
`sglsurvey/adapters/mast_ps1.py`.

**Status.** Complete — pilot + scale-up 2026-08-20 → 21, v2 rerun
2026-08-23. Report `report/ps1_survey.md`; narrative
`notes/project_history.md` §4.

**Result.** **0 candidates** (v2 confirmatory 450 cells); asteroid
(60000) positive control recovered to ≤ 0.1 mag.

### 4.5 PS1 + ZTF (joint)

**Description.** Joint weighted stack of the PS1 and ZTF tensors on
the common trajectory and µ reference epoch (`surveys/joint/`),
supplying the parallax phases PS1 alone lacks (321 of 414 cells
both-phase) and a cross-archive persistence test over 2009–2026.

**Status.** Complete — v1 2026-08-21, v2 rerun 2026-08-23, **v3
(three archives, adding the WISE W1/W2 colour axis) 2026-08-24 → 25**:
freeze v3.0 with the optical statistic unchanged; WISE
joint-conventions tensor rebuild (T0 59800, AB ZP 25,
`surveys/wise/profile.py` → `runs/wise/v3/`); three-archive
same-source injections; the colour-consistency test as a calibrated
veto; plus the stack-regime asteroid positive control (220000) through
the per-archive and joint chains. Report `report/joint_ps1_ztf.md`;
narrative `notes/project_history.md` §§4, 8.

**Result (v3).** **0 candidates** (blind confirmatory 45 endpoints /
230 cells, R̃_FWER 1.778) — optical numbers bit-identical to v2; all
230 cells W1/W2-annotated, 0 would-fire colour vetoes, measured
false-veto rate 0 in the fitted population (5/46,302 raw, all above
the optical bright limit); control recovered in all three joint
families (R̃ 2.21/3.04/1.77 vs threshold 1.54).

**Result.** **0 candidates** (v2 confirmatory 230 cells), including
the PS1 marginal cells vetoed by direct ZTF forced photometry.

### 4.6 DECam

**Description.** NOIRLab DECam instcal exposures (2012–2026, CTIO,
grizY) — the only 1″-class multi-epoch optical archive covering the
15 corridors below δ = −30° (previously closed by nothing; includes
σ Dra, HD 219134, Lalande 21185). First survey with no v1 exploratory
phase — v2 discipline from day one. Adapter
`sglsurvey/adapters/noirlab_decam.py`.

**Status.** Complete — recon through blind confirmatory run,
2026-08-24. Report `report/decam_survey.md`; narrative
`notes/project_history.md` §5; recon
`surveys/decam/notes/decam_recon_2026-08-24.md`.

**Result.** **0 candidates** (blind confirmatory 14 endpoints /
86 cells; development 30 cells); the program's deepest optical
corridor constraints (median persistent m90 g 23.2 / r 22.1 / i 22.7 /
z 22.1 / Y 21.0 AB) where cells are constrainable — 378/712
confirmatory threshold constraints `not_constrainable` in the
heterogeneous PI coverage.

## 5. Pipeline B — crossings surveys

Execution of §3.5, one subproject per archive. Shared geometry
products feed every survey: the universal crossing list
(`crossings/universal_v1/`: 16,586 b(t)-minimum events, 88 endpoints ×
both link directions, 1980→2028) and the spacecraft-observer lists
(`crossings/{wise,tess,spherex}_v1`: 5,170 / 3,390 / 1,130 events;
the TESS HEO correction validated at |Δb| ≤ 0.33 R☉). Channels:
A = uplink interception at the star near opposition; B = downlink
pre-lens interception at the antipode; sunward combinations out of
scope (revisited in §5.8 item 9). Hypotheses (beam radius, wavelength,
duty cycle) are frozen per survey before any data are touched, with
the v2 statistical discipline inherited from §4.

The 2026-08-23 → 24 programme closed with **0 candidates** across all
archives (tally: 47 searched trials, 5 exceedances vs ~6.1 expected);
execution record `notes/project_history.md` §6; lessons
(sidereal-year phase locking, the two structural gate theorems,
on-star systematics templates, the antipode workhorse, correlated
exact-mask attrition, the warp-direct substrate) in
`notes/learnings.md` §8.

### 5.1 ZTF

**Description.** First crossings survey: channels A and B on the
difference-image substrate, 2018–2026 — the design (coverage
intersect, saturation cut, frozen thresholds, pseudo-window controls,
stamp-response injection completeness) that later archives adopted.

**Status.** Complete 2026-08-23. Report `report/ztf_crossings.md`.

**Result.** **0 candidates**; 36 searchable units, 3 exceedances vs
4.0 expected, all adjudicated. Downlink relay power ≳ ~130 W through
the 2.5 R☉ cone in covered windows; channel A contrast-limited
(10-m uplink ≳ 56–400 kW).

### 5.2 Pan-STARRS1

**Description.** The ZTF channel constructions on the warp-direct
star-calibrated substrate (PS1 has no difference images), grizy,
2009–2015 — extends the covered-window record back a decade, disjoint
from the ZTF era.

**Status.** Complete 2026-08-24. Report `report/ps1_crossings.md`.

**Result.** **0 candidates**; 11 searched B units, 2 exceedances vs
2.1 expected (ross-128 r vetoed flux-consistent-static; gj-1276 i
retained-ambiguous, handed to the joint stage). The grazing rung ended
unconstrained — correlated chip-gap mask attrition took the
van-maanen b = 0.28 R☉ event.

### 5.3 WISE

**Description.** Crossings over the WISE spacecraft-observer list
(W1/W2, 2010–2024) — the elongation-90° surveyor case.

**Status.** Complete 2026-08-24 as a structural null: no pixel was
searched. Report `report/wise_crossings.md`.

**Result.** **0 searchable units** by two pre-pixel gates — Gate I
(elongation theorem: channels B and A-0.1 invisible in principle) and
Gate II (control-geometry theorem: annual ~116 d windows defeat the
temporal pseudo-window family). The 91 covered windows enter the joint
ledger as coverage-without-statistic.

### 5.4 ZTF + PS1 + WISE (joint)

**Description.** Unified covered-window ledger across the three
archives (593 rows, 42 multi-archive target-channels, full disposition
census) plus the programme's one designated follow-up: the gj-1276
recurrence test over 11 ZTF event-bands.

**Status.** Complete 2026-08-24 under pre-registration
(`surveys/joint-crossings/plan.md`). Report
`report/joint_crossings.md`.

**Result.** **Programme closed with 0 candidates.** The flat-SED
persistent-relay interpretation of the PS1 gj-1276 anomaly refuted at
90 % (surviving interpretations — red/line SEDs or non-persistence —
not promotable); structurally open cells enumerated in the report.

### 5.5 Breakthrough Listen + VLASS (radio)

**Description.** Metadata-only intersection of BL Open Data pointings
(93 position cones at stars + antipodes) and VLASS per-tile epoch
dates with the universal crossing windows, to decide the repo's radio
scope.

**Status.** Complete 2026-08-24. **Decision: geometry-only — no radio
signal search in this repo.** Report `report/radio_crossings.md`.

**Result.** Narrow-rung on-star radio coverage: zero; antipodes:
zero in-window BL pointings — the workhorse channel is archivally
virgin (a future-observation recommendation); the wide rung is already
answered by BL's published nulls. Micro-follow-ups handed off (VLASS
six-cutout quick-look; 3 teegarden APF spectra to the spectral-archive
family, §3.6 row 10).

### 5.6 TESS

**Description.** FFI cutouts of the crossing windows from the
Ecliptic sectors — the first archive able to resolve complete
ingress→egress crossing light curves and test the pulse-period cell
(two grazing-rung events: wolf-359 b = 0.55 R☉, teegarden
b = 1.04 R☉; two complete near-grazing channel-A windows). The era's
two deepest grazes (van-maanen 0.078 R☉, gj-1276 0.49 R☉) fell
between sectors and stay open.

**Status.** Complete — coverage gate through blind confirmatory,
completeness (SPOC PRF, both temporal models) and report, 2026-08-24.
Report `report/tess_crossings.md`; execution log
`notes/project_history.md` §7.

**Result.** **0 candidates** (6 confirmatory B units × 2 statistics).
The pulse-period cell closes clean on both resolved grazing crossings
— the first constraint of its kind (≥1-cadence pulses ≳ 2 kW through
the 1.2 R☉ cone); persistent chord limits ~0.6–4 kW through the
grazing cones. Chord statistic systematics-dominated (design lesson:
v2 needs a detrending layer); 1 vetoed exceedance, 1 budget-absorbed,
1 retained-ambiguous (teegarden 0.1 AU — recurrence at a future
ecliptic sector is the designated follow-up).

### 5.7 Standing maintenance

Yearly crossing-list refresh (the 2028 window end) as archives extend;
fold the TESS rows into the covered-window ledger at that refresh
(report `tess_crossings.md` §6 enumerates them; the joint-stage v1
ledger predates TESS and stays as pre-registered) and re-run the
teegarden 0.1 AU recurrence test when a future ecliptic sector covers
the target;
revisit SPHEREx crossings when 3+ annual windows exist per target
(decision 2026-08-24: Gates I+II apply and the 1.5-yr archive is too
thin); re-run the radio intersection if BL's MeerKAT holdings reach
the public archive; the v2-style null-ensemble redesign is the
designated route should the d = 1 wide-beam rung ever be searched
with calibrated error rates. Universal-list caveat: impact-parameter
uncertainty is not propagated (`sglseti.crossing_uncertainty` exists
for per-event follow-up).

### 5.8 Archive-expansion queue (adopted 2026-08-24)

Datasets adopted from the 2026-08-24 archive brainstorm, ordered by
which structurally open cell (`report/joint_crossings.md` §3) each
opens and by adapter cost. TESS (§5.6) — the queue's first item — is
complete; item 2 (ATLAS + ASAS-SN) is next. Access details below
marked _unprobed_ are from general
knowledge, not verified endpoints — each item starts with a
DECam-style reachability recon before any freeze.

| #   | Dataset                                                           | Open cell / rationale                                                                                                                                                                                                                                                                          | Access route (status)                                                                          |
| --- | ----------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| 1   | **TESS FFIs**                                                     | Pulse periods between exposure and window length; first resolved crossing light curves                                                                                                                                                                                                         | **Done (§5.6)** — 0 candidates, pulse cell closed clean                                        |
| 2   | **ATLAS forced photometry + ASAS-SN Sky Patrol**                  | Window _coverage fraction_ — nightly/quad-nightly all-sky cadence (ATLAS δ > −50°, o/c ~19.5, 2015–; ASAS-SN all-sky, g/V ~18, 2012–) incl. southern targets ZTF never sees; duty-cycle/flare-beacon rung, not depth                                                                           | Public forced-photometry APIs at arbitrary positions — no image handling (_unprobed_)          |
| 3   | **PTF/iPTF (2009–2016)**                                          | Fills the PS1-era → ZTF-era temporal gap in the north with difference-imaging-era data (R ~20.5–21)                                                                                                                                                                                            | IRSA IBE `ptf/products/` — same machinery as `irsa_ztf.py`, adapter mostly config (_unprobed_) |
| 4   | **GALEX time-tagged photons (gPhoton, 2003–2013)**                | Pulse-period cell in the UV, in an era predating everything but PS1; 5 ms photon time-stamps, light curves at arbitrary positions                                                                                                                                                              | gPhoton photon database via MAST (_unprobed_)                                                  |
| 5   | **DASCH scanned plates (1885–1992, DR7)**                         | A century of annually-recurring windows before the 1980 list start (B ~15–17). Prerequisite: extend `sglsurvey.crossings` window backward — cheap, and exercises the §7 model-accuracy budget at old epochs                                                                                    | DASCH lightcurve + cutout services at Harvard (_unprobed_)                                     |
| 6   | **Palomar Gattini-IR (J, 2018–, J ~16) + WINTER (Y/J/Hs, 2023–)** | The >900 nm time-domain cell (1064 nm falls in Y); shallow — framed as opening the cell, not deep exclusion. The 1550 nm line hypotheses stay with the spectral-archive family (§3.6 row 10: APOGEE H-band, ESO NIRPS, alongside the APF hand-off)                                             | Archive access uncertain — recon first (_unprobed_)                                            |
| 7   | **Radio metadata extensions**                                     | Repeat the VLASS geometry-only construction (radio decision unchanged, §5.5) on archives that revisit fields: ASKAP VAST + RACS per-field epoch dates, LoTSS pointing dates; the virgin antipode channel is cheap to keep testing                                                              | CASDA / LoTSS DR services (_unprobed_); BL-MeerKAT re-check stays a §5.7 maintenance item      |
| 8   | **Kepler/K2 superstamps + FFIs**                                  | One-afternoon footprint intersect: K2 ecliptic campaign fields (2014–2018) + Kepler prime field vs `universal_v1` antipodes; proceed only on a hit (30 min continuous cadence in the pre-ZTF era)                                                                                              | MAST (well-known); intersect needs no adapter                                                  |
| 9   | **Heliospheric imagers (SOHO/LASCO, STEREO/HI, PSP/WISPR)**       | First a one-page geometry study: are any of the sunward channel combinations (declared out of scope at the ZTF freeze) _only_ visible at small elongation? If yes, these are the unique archival substrate (decades of continuous coverage, V ~8–13); if no, close the question for the record | Geometry study first; archive recon only on a positive answer                                  |

## 6. Future projects (both pipelines)

Adopted 2026-08-24; blocked on external releases/approvals, so
tracked here rather than queued:

- **Rubin alert stream.** Alerts have streamed to brokers since
  2026-02; Early DP2 (2026-07-27) has coadds + catalogs, with visit
  and difference images targeted late 2026; the six-month DR1 was
  cancelled in favour of a full Year-1 release. Near-term project:
  broker-based forced-photometry-style screening of the southern
  corridors (Pipeline A) and crossing windows (Pipeline B); the
  bulk-image adapter waits for visit/difference images. **A request
  for Rubin data-platform access is pending — if approved, a
  dedicated Rubin project (full adapter on the v2 design) follows.**
- **Gaia DR4 (not before 2026-12).** All per-transit epoch photometry
  _and astrometry_ for every source: a qualitatively new test —
  astrometric SGL-track fitting against the 20–375″ relay parallax —
  plus an L2-observer crossing list via the now-generic
  `register_spacecraft_table_observer` machinery. Candidate for the
  2027 crossing-list refresh cycle.

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
- **Radio scope.** ~~Decide whether this repository performs voltage or
  spectrogram searches, delegates them to archive-specific pipelines, or
  publishes only geometry and coverage products.~~ Decided 2026-08-24
  (§5.5): geometry and coverage products only. In every case,
  preserve frequency, drift-rate, polarization, time, and sensitivity
  dimensions.
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
