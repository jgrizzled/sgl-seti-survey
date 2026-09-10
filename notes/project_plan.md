---

## title: "sgl-seti-survey — Project Plan"
date: 2026-09-06
tags:
  - SETI
  - technosignatures
  - solar-gravitational-lens
  - archival-search

# sgl-seti-survey — Project Plan

Companion documents: `notes/project_history.md` (chronological
execution record — what ran, when, with what outcome) and
`notes/learnings.md` (consolidated design lessons). This plan keeps
goals, methodology, sequencing, and current status. The programme-wide
index of every constraint and open cell is `report/programme_ledger.md`
(§6).

## 1. Introduction

### 1.1 Background

A technologically advanced interstellar network might place communication
relays on or near the gravitational focal lines of stars, using stellar
lenses to reduce the power and aperture required for links between
neighboring systems. For the Sun, the focal line for rays grazing the
photosphere begins near 547.8 AU, so a relay serving a target star would
sit roughly on the anti-star ray at ~550–10,000 AU heliocentric distance.

This hypothesis is unusually testable for SETI: it predicts *where* on
the sky a relay associated with a given star should appear, *how* that
position moves (a huge annual parallax of ~20–375 arcsec depending on
distance, plus secular drift opposite the target star's proper motion),
and *when* Earth crosses hypothesized beam geometries between the Sun and
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

## 4. Pipeline A surveys



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
cutouts). Adapter `sglsurvey/adapters/irsa_spherex.py`. The deferred
controls (former ledger item O3) were two items the v2 rerun left
open: (i) the six-detector joint cell from the stored per-tensor
accumulators, and (ii) a template refit with injections, because
slow-source absorption by the static template was unmodelled. QR3 was
not out when they ran (IRSA `spherex.obscore` served only
`spherex_qr2`, checked 2026-09-04) and neither depended on it.

**Status.** Complete — pilot, southern scale-up and all-sky run
2026-08-20 → 21, v2 rerun 2026-08-23; deferred controls complete
2026-09-04 (`report/spherex_joint6.md`; hypotheses v2.1 for the joint
family; `notes/learnings.md` §10, §13). The v3 templates were found
lost in the v1 retirement (a dangling calib_v4 link) and regenerated
from the cutouts, verified against a rebuilt tensor. One dev-driven
amendment: q95 ≤ 0 cells are void (an engine gap, carried into §4.7).
Report `report/spherex_survey.md`; narrative `notes/project_history.md`
§3 and §20. The QR3 re-run is its own future survey, §4.7.

**Result.** **0 candidates** (all-sky, 77 corridors / 88 endpoints —
the first archive to cover the full universal list; v2 rerun blind;
joint cell 0 candidates in 122 confirmatory cells, R̃_FWER 2.249, 11
R > 1 vs 12.7 expected; dev 54 cells, 0). Joint-cell persistent m90
20.93 AB, worst-of-four 19.15. Depths: the v2 persistent m90 (D1–D4
20.1–20.4) are overstated for slow sources because the static
template absorbs a median 49 % of a source present in a node's only
2–3 visits (1,056 cells × a 28-source ladder: Δm 0.73 median, p90
1.3–1.9; cadence-driven — deep field 5–10 %, single-visit corridors
80–90 %; validated end to end to < 0.01); corrected limits D1–D4
19.2–19.4, D5 18.8, D6 18.2, joint 20.02, published report-level with
the records flagged.

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

### 4.7 SPHEREx re-run at QR3 (future; former ledger item O8)

**Description.** Repeat of the §4.3 survey on the next SPHEREx quick
release. QR3 is the first release with a true parallax-phase hold-out
per corridor, so it upgrades the persistence test from the QR2
phase-poor cadence to a genuine held-out-epoch prediction. Carry into
the freeze, from the §4.3 deferred controls: a source-excluded
template (per node, drop epochs with the hypothesised track within
~2 FWHM — removes the 49 % absorption by construction); the q95 ≤ 0
degenerate-null guard in the per-detector engine
(`sglsurvey/nulls_stage.py`); the six-detector joint cell (hypotheses
v2.1) as a standard family; and per-exposure-time, not per-detector,
visit/block on-patterns in the joint injections.

**Status.** Not started; waiting on the release.

**Trigger.** IRSA `spherex.obscore` serving a `spherex_qr3` collection
(QR2 only as of 2026-09-04). Runs as a new hypothesis version with its
own declared hold-out, per the standing rules.

### 4.8 Rubin corridor survey (future; former §6 / ledger item O11)

**Description.** Pipeline A over Rubin/LSST on the southern corridors,
the archive the ZTF adapter was ~90 % of a design for. Rubin alerts
have streamed to brokers since 2026-02; Early DP2 (2026-07-27) has
coadds + catalogs, with visit and difference images targeted late
2026; the six-month DR1 was cancelled in favour of a full Year-1
release. Access granted 2026-08-26 (`RSP_API_TOKEN`) and the recon
completed the same day
(`surveys/rubin/notes/rubin_recon_2026-08-26.md`): TAP / SIA /
DataLink / SODA all verified; DP2 (Early) = 28,698 LSSTCam visits
2025-03 → 2025-12 with DiaSource + per-visit forced difference
photometry at r ~ 24 single-epoch depth; the geometry intersect shows
19/76 corridors covered, phase-limited. The image-level v2 chain wants
the Exposure-cutout substrate (image + mask + variance + PSF in one
call — the simplest v2 substrate probed so far). A broker-based,
forced-photometry-style screening of the corridors remains the
near-term option should one be wanted before the images.

**Status.** Tabled 2026-08-26: Early DP2's 1–3 calendar-month corridor
coverage is parallax-phase-poor, so the staging decision ran the
Pipeline B catalog-level survey first (§5.13) and holds Pipeline A for
the images.

**Trigger.** The visit/difference-image release (late 2026), revisited
at the Year-1 release alongside Gaia DR4 (§4.9). The Pipeline B
image-level re-run (§5.21) shares the adapter.

### 4.9 Gaia DR4 astrometric SGL-track test (future; former §6 / ledger item O12)

**Description.** Gaia DR4 publishes all per-transit epoch photometry
*and astrometry* for every source — a qualitatively new Pipeline A
test: astrometric SGL-track fitting against the 20–375″ relay parallax
and the anti-target-motion drift, rather than a photometric stack. The
matching Pipeline B product (an L2-observer crossing list) is §5.22.

**Status.** Not started; adopted 2026-08-24.

**Trigger.** The DR4 release (not before 2026-12); candidate for the
2027 crossing-list refresh cycle (§5.27).

### 4.10 Long-wavelength thermal archives (not adopted; former §3.6 row 5)

**Description.** IRAS, AKARI, Herschel, Planck and Spitzer: the
longer-wavelength thermal tests of the waste-heat cell that W3/W4
(§4.1) only reach as threshold statements. Some of these require
spacecraft observers (Herschel and Planck at L2, Spitzer on its
Earth-trailing orbit) or older target solutions; the generic
spacecraft-observer machinery built for TESS/SOHO/STEREO/PSP (§5)
covers the first need.

**Status.** Not adopted; no recon. Priority 5 on the original list.

**Trigger.** A DECam-style reachability recon, which every new archive
starts with before any freeze; none is scheduled.

### 4.11 2MASS and photographic-plate corridors (not adopted; former §3.6 row 5)

**Description.** Long time baselines for the persistent-source track
on the near-infrared and photographic substrates — 2MASS (1997–2001)
and the POSS plates. DASCH already served the photographic century for
Pipeline B (§5.14); its Starglass substrate (exposure lists, century
lightcurves, cutouts) is the obvious plate substrate should a corridor
survey be adopted.

**Status.** Not adopted; no recon.

**Trigger.** None scheduled; would start with a reachability recon
like every new archive.

### 4.12 MPC observations, tracklets and known-object products (not adopted; former §3.6 row 8)

**Description.** MPC isolated observations and tracklets, JPL/MPC
known-object ephemerides and survey reject tables as candidate and
veto inputs to the Pipeline A screens (§3.4 step 1 and 5). Not a
complete search on its own: 550–10,000 AU reflex rates may fall below
ordinary intra-night tracklet thresholds, so an unlinked relay would
sit in the isolated-observation pool rather than in a tracklet.

**Status.** Not adopted as a survey. Known-object ephemerides already
serve as positive controls inside the completed surveys (the asteroid
controls of §4.4, §4.5, §5.9 and §5.13).

**Trigger.** A candidate that survives a confirmatory run — none has.

## 5. Pipeline B — crossings surveys



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
family, §5.17). Follow-ons: the CASDA/LoTSS extension (§5.7) and the
cutout quick-look (§5.8).

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

### 5.7 Radio metadata extensions — ASKAP + LoTSS

**Description.** The §5.5 geometry-only construction repeated on
archives that revisit fields: CASDA ObsCore (RACS / VAST / EMU / FLASH
and every ASKAP continuum collection, per-SBID on-sky intervals;
VAST-pilot dates via `casda.observation`, anonymous TAP) and ASTRON VO
`lotss_dr3.pointings` (`dateallobs`), intersected with the universal
crossing windows (`scripts/coverage_intersect_v2.py`). Adopted because
the antipode channel was archivally virgin and cheap to keep testing.

**Status.** Complete 2026-09-04. Report `report/radio_crossings_ext.md`;
execution log `notes/project_history.md` §16. Decision unchanged:
geometry-only. Re-runs at each yearly refresh (§5.27).

**Result.** Grazing rungs: 0 observations in 65 ASKAP + 30 LoTSS
in-span events. **The antipode channel is no longer archivally virgin
at the 0.1 AU rung**: 6 validated wide-field ASKAP epochs on 5 targets
sit inside antipode windows (wolf-359 2023-09 and 2024-09 at 0.8° from
the VAST_2257-06 centre; ross-154 RACS-high 2021-12 at 0.63°;
ross-128, van-maanen, teegarden nearer the footprint edge); on-star A
0.1 AU: gj-1276 VAST 2023, teegarden RACS-mid 2024 + LoTSS 2023 (in
beam), van-maanen VAST 2024 — the first archival radio coverage of the
workhorse channel. The "nobody has ever looked" statement now holds
for targeted radio and the grazing rungs only.

### 5.8 Radio quick-look

**Description.** A broadband transient check at the relay position —
a look, not a search, the geometry-only decision unchanged: RACS/VAST
catalogue cone-search first, image cutout second, on the 6 validated
ASKAP B 0.1 AU antipode epochs and the 4 on-star narrow-rung epochs of
§5.7, plus the 4 strict + 2 tolerance-edge VLASS on-star epochs.

**Status.** Done 2026-09-04. Report `report/radio_quicklook.md`;
execution log `notes/project_history.md` §18. The ASKAP stage ran with
the OPAL login (`.env`); only the REJECTED FLASH SB83234 is closed.

**Result.** 17 epochs, 0 coincident catalogue components (10
catalogues, 7 epoch-resolved), 0 pixels above 3σ in 12 ASKAP + 26
VLASS cutouts; both wolf-359 antipode epochs empty at 0.2 mJy rms,
the teegarden 2026 antipode empty at 0.07 mJy (2-h FLASH); exact tile
times put the two VLASS tolerance-edge rows outside the window. The
gj-1276 2.7σ in-window pixel is the maximum of a 16-epoch noise
series. No residual.

### 5.9 ATLAS + ASAS-SN

**Description.** The window coverage-fraction / duty-cycle-beacon
rung: nightly/quad-nightly all-sky forced photometry (ATLAS o/c
~19.5, all-sky since 2022, δ > −50° since 2015; ASAS-SN g/V ~18,
2013–) over the universal crossing windows. The survey's new cell is
**recurrence** — 22 in-era windows per grazing-rung target, incl. the
van-maanen b = 0.24 R☉ family lost by PS1 masks and TESS sector
gaps; recurrence-stacked statistic primary. ATLAS is the workhorse
for both channels (authenticated arbitrary-position difference-flux
forced photometry; built-in MPC positive control); ASAS-SN's role is
narrowed by its access model (Sky Patrol v2 anonymous but
catalogued-sources-only) to on-star channel A plus a
coverage-fraction ledger; v1 is reCaptcha-gated — a manual
instrument for pre-registered follow-ups only.

**Status.** COMPLETE 2026-09-04 (`report/atlas_asassn_crossings.md`).
Resumed 2026-09-03 on the new server. Reachability recon
(`surveys/atlas-asassn-crossings/notes/atlas_asassn_recon_2026-08-25.md`)
and hypothesis freeze v1.0 (D1–D8, incl. the wolf-359 `forced_dev`
remedy for pre-freeze data contact) done 2026-08-25; the resume
checklist (`notes/server_migration_resume.md`) cleared 2026-09-03.
Pre-data amendments (hypotheses §12): **A1** — the freeze's "22
semiannual windows" counted the sunward axis crossing of each year,
which the universal declaration already excludes; the searchable
population is **11 annual windows per target per channel** (B 1.2 R☉
44 / 2.5 R☉ 55 / 0.1 AU 77; A 0.1 AU 77), stack depth ~0.4 mag
shallower than pre-declared; **A2** — ATLAS era end MJD 61286,
ASAS-SN v2 ceiling re-measured unchanged at 2025-06-16 (reprocessing
lag); **A3** — no radec-list batching exists; a windowed task runs in
~2.2 min, so the 308-task fleet (231 B incl. the D3 mini-track ends,
77 A) drained in 7.4 h with 0 failures (`scripts/atlas_drain.py`).
**Coverage ledgers** (`results/coverage_v1_`*, `asassn_ledger_v1_*`):
the grazing rungs are attrition-dominated — B 1.2 R☉ 8/44 windows
covered (van-maanen 4/11), 2.5 R☉ 25/55, against 0.1 AU 74/77 (B)
and 68/77 (A) — because 0.5–1.3 d windows sit below ATLAS's nightly-
at-best revisit; the pre-declared ~100–180-exposure grazing stacks
are ~18–30 exposures in reality (≈ 20.6 stacked), the 0.1 AU stacks
reach ~21.5. ASAS-SN v2 field-level epoch lists show the same shape.
Threshold freeze v1.0 (hypotheses §13; seed 20260904; dev = wolf-359

- gj-908; per-statistic gates, mirror-gated control validity, ≥ 4
valid controls per trial): 14 searched units / 37 trials; S_stack
constraint-only at every grazing rung. Dev 2 exceedances (one
single-exposure epoch, adjudicated); blind confirmatory 5 vs 3.2
expected on 3 events — a retained-ambiguous 2.8σ quad-consistent
van-maanen 1.2/2.5 R☉ chord (2017-04, no recurrence), a
difference-image baseline systematic (gj-1276 0.1 AU), and the
channel-A PM-dipole systematic (teegarden). D6 control (15000) CCD:
scale gate pass, ±0.1 mag.

**Result.** **0 candidates.** First per-window constraints through
the van-maanen b = 0.24 R☉ photosphere-grazing cone: ≳ 78 W (S_event,
o ≈ 20.5) in 4 of 11 windows; 2.5 R☉ cones 280–590 W in 4–5 windows
(teegarden, van-maanen, wolf-359); 0.1 AU cones 26–55 kW per window,
26–41 kW recurrence-stacked, all seven targets; pulse cell 140 W–890 W
(grazing) / 36–79 kW (0.1 AU) per 30-s exposure. Channel A delivers
no constraint (PM dipole against the template). The recurrence cell
the survey was adopted for is not calibratable at nightly cadence
for 0.3–1.3 d windows (≤ 3 valid pseudo-stacks); coverage fractions
18 % / 45 % / 95 % by rung are the structural finding.

### 5.10 PTF/iPTF

**Description.** IRSA IBE level-1 epochal images, 2009-03 → 2015-01
(2.8M CCD exposures; late iPTF was never publicly released), g +
Mould R, on a scie-direct substrate (no public difference images;
level-2 reference coadds are annotation inputs only). PS1-era-parallel
from a different site, cadence and filter set: the covered-window
record's first **independent re-observations of PS1-era windows**,
including the van-maanen April deep-graze family (b 0.27–0.37 R☉)
that PS1 lost to chip-gap masks. Adapter
`sglsurvey/adapters/irsa_ptf.py`; the archive publishes MD5s and
served 710/710 cutout fetches with zero 404s — the
highest-integrity interface in the programme so far.

**Status.** Complete — recon, hypothesis freeze v1.0 (D1–D8),
coverage, threshold freeze, dev and blind confirmatory + completeness
all 2026-08-26, entirely on the dev machine. One amendment (v1.1:
384-px cutout calibrator support — a surface-density fix, statistics
untouched). Report `report/ptf_crossings.md`; execution log
`notes/project_history.md` §10; recon
`surveys/ptf-crossings/notes/ptf_recon_2026-08-26.md`.

**Result.** **0 candidates** (blind confirmatory 6 units / 10
trials, 0 exceedances vs 1.11 expected control crossings). First
pre-2015 constraint on the photosphere-grazing cell: relay power
≳ ≈190 W (van-maanen 2011 b = 0.28 R☉, g) / ≲110 W (ross-128, R)
through the 2.5 R☉ cone in the covered windows; kW-class (8–30 kW)
0.1 AU-cone limits across the 2009–2013 recurrences; the
recurrence-stack cell closes clean on 4 units (m90 20.6–≥22.0 AB);
the B 1.2 R☉ rung is structurally uncovered (coverage-without-
statistic ledger). Dev also resolved gj-1276 A constraint-only at
the frozen offset gate — campaign cadence defeats on-star temporal
controls (lesson recorded, report §4).

### 5.11 SOHO/LASCO sunward survey

**Description.** The sunward channel combinations, excluded from
every prior survey and revisited by the 2026-08-25 geometry study
(answer **YES**): S1 downlink post-lens (apparent source = the
solar-limb graze point at 19′–40′) and S2 uplink past the Sun (star
at ε ≤ 5.7° for the 0.1 AU rung) are confined to ε ≲ 6° for every
frozen rung except the 1 AU uplink outer skirt — coronagraphs and
heliospheric imagers are the **unique substrate**; MW-class power
cell. LASCO C2/C3 first (SDAC/NRL anonymous trees, 1996–2026) with a
mandatory SOHO observer list — the L1 halo offset is 0.93 R☉, so
`crossings/soho_v1` was built (10,714 events; grazing-family
|Δb| ≤ 0.21 R☉ vs Earth-center, membership unchanged; the van-maanen
October family deepens to b ≈ 0.13 R☉); 128 LASCO-era grazing
windows on the van-maanen/wolf-359/teegarden/gj-1276 family.

**Status.** **COMPLETE 2026-08-27** — report
`report/lasco_crossings.md`. Chain:
geometry study + LASCO recon + `soho_v1` + hypotheses v1.0 +
coverage 2026-08-25 (S1 1.2 R☉ `not_constrainable` behind the C2
occulter; S1 2.5 R☉ 148/153 windows covered+visible in the
[2.2, 2.5] R☉ wings; 0.1 AU rungs 96–98 % covered at ~1,000 C3
frames/window); threshold freeze v1.0 + amendments v1.1/v1.2; dev
closed and blind confirmatory run once, both 2026-08-26; completeness

- controls + report 2026-08-27. Later: STEREO HI-1 (own observer
list, §5.16, complete) and WISPR (PSP geometry pass, §5.18, complete;
survey §5.19). Survey docs:
`surveys/heliospheric-crossings/`.

**Result.** **0 candidates blind** (`results/confirmatory_v1.md`):
30 searched trials over 13 units, 2 exceedances vs 3.3 expected —
under budget, every recurrence statistic null. Adjudications:
gj-1276 S_pulse = one 12-min frame (cosmic-ray class, vetoed by the
frozen persistence rule); van-maanen S_event = one 2020-10-05 window
whose ~3 h all-PA burst follows a CDAW-catalogued C2 CME by 48 min
(adjudicated CME-period systematic, retained, non-promotable).
Structural finding: ross-128's S1 antipode is **permanently blended**
(fixed sky point 121″ from a VT 7.2 star — the S1 source is a fixed
ICRS point the Sun sweeps past) → constraint-only,
nominal-covered/resolution-blended. Dev (11 trials, 1 exc vs 1.2):
ross-154's Venus 8-yr-synodic stray-light family vetoed; gj-908
(V 8.98) detected window-locked at the predicted S2 position and
nulled by the stellar template — the in-situ positive control.
Completeness (injections into the real null series, coronal
systematics-limited): recurrence-stack m90 V 4.5–8.4 → **the first
sunward-cell constraints: pulse ≥ 6.7 MW per 12-min through the
2.5 R☉ cone; 40–530 MW recurrence-stacked downlink; 3–26 GW through
the 0.1 AU cone; 10-m uplink 0.7–1.1 GW** (±0.3 mag band-conversion
systematic declared; Uranus recovered at 20σ stacked with its
methane-spectrum deficit documented; flux gate carried by 217
unbiased star measurements + era stability 0.064 mag).

### 5.12 GALEX/gPhoton

**Description.** The UV pulse-period cell: GALEX time-tagged
photons (5 ms stamps) from the gPhoton photon database at MAST,
2003-06 → 2013-05, queried at arbitrary positions over the
universal crossing windows. Recon (2026-08-26) found the grazing
rungs (1.2/2.5 R☉, ±0.35 d windows) **structurally uncovered** by
GALEX's sparse visit cadence — they enter the ledger as
coverage-without-statistic — while the **0.1 AU rung** (±5.8 d
windows) has 5 usable in-window units at 4 target-channels
(gj-1276 B 2007+2010 incl. a 1.6 ks NUV visit; wolf-359 A,
gj-1276 A, ross-128 A), 4 with simultaneous FUV — the programme's
first two-band photon-level units, in the pre-ZTF 2007–2010 era.
Substrate = direct photon-event counting from the DB (no images);
photon pull verified in-window (1,695 NUV photons, ms-tick stamps).
Earth-center `universal_v1` valid (LEO, standing 0.010 R☉ budget).

**Status.** **Complete** — recon, hypothesis freeze v1.0 (D1–D8),
amendments v1.1 (aspect `flag % 2` gate; restores the 1.6 ks
gj-1276 B 2010 visit from the pervasive flag-64 state, astrometry
verified) and v1.2 (S_period tick-jitter — the 5 ms/200 Hz
quantization artifact), coverage (47 snapshots; 5 units / 2,045
NUV s confirmed; grazing rungs 0 covered), threshold freeze v1.0,
dev stage (ZP_eff 19.592 ± 0.088; background 0.209 cts/s/8″;
segment gate → gj-1276 A both bands + wolf-359 A FUV
constraint-only, family 5 band-units / 15 trials; wolf-359 flares =
in-situ positive control + FRED templates), blind confirmatory,
completeness and report — all 2026-08-26, entirely on the dev
machine. Report `report/galex_crossings.md`; execution log
`notes/project_history.md` §11; survey docs
`surveys/galex-crossings/`.

**Result.** **0 candidates** (15 frozen trials, 14 effective;
**1 exceedance vs 1.56 expected control crossings** — the wolf-359
A NUV 5-photon 0.5 s burst, adjudicated retained-ambiguous: flare
veto cannot positively fire at that fluence, no in-archive
recurrence window; micro-flare leading, recorded for any future UV
mission able to re-cover a wolf-359 window). First UV constraints
of the programme: persistent relays ≳ 20–33 kW through the 0.1 AU
cone in the covered 2007–2010 windows (NUV m90 21.3–21.9, FUV
≈ 21.2–21.8; wolf-359 variability-limited at ≈ 415 kW); single UV
pulses ≳ 1.7 MW (0.5 s) / 17 MW (0.05 s) peak; **coherent pulse
trains P = 0.5–50 s ≳ ≈ 6 kW time-averaged** on the 1.6 ks gj-1276
B 2010 window — the programme's first sub-second-cadence crossing
constraint; 266 nm (quadrupled Nd:YAG) the first in-band harmonic
of the 1064 nm family covered anywhere in the programme.

### 5.13 Rubin DP2 crossings (catalog-level)

**Description.** The first Rubin survey: Pipeline B over DP2 (Early)
— 28,698 LSSTCam visits, 2025-03-26 → 2025-12-08, on the
**catalog-level substrate** (DiaSource difference-image detections +
per-visit forced difference photometry; no visit/difference images
until late 2026). The recon
(`surveys/rubin/notes/rubin_recon_2026-08-26.md`) found one
grazing-family event with an in-window epoch — **ross-128 downlink,
b = 1.86 R☉, one r visit at ~24 AB single-epoch depth** in the
±5.8 d 0.1 AU window (the ±0.35 d grazing windows themselves are
uncovered) — potentially an order-of-magnitude-deeper power
constraint on the workhorse antipode channel than the ZTF ~130 W
record. ross-154 A is expected saturation-excluded (V 10.4 vs the
~16 mag 30 s limit); A 1.0 AU keeps its programme-wide deferred
status (§5.27). New structural cell: the first survey whose
detection substrate is another pipeline's difference-image catalog —
the freeze must handle association-not-forced-photometry,
archive-side injection flags, and real/bogus scores without absolute
tuning.

**Status.** Complete — recon through blind confirmatory, all
2026-08-26. Report `report/rubin_crossings.md`; execution log
`notes/project_history.md` §13; survey docs
`surveys/rubin-crossings/` (freeze v1.0 + amendment v1.1: control
locus-avoidance 10″ → 2.5″, found at dev).

**Result.** **0 candidates** (1 blind trial — ross-128 B 0.1 AU r —
0 exceedances vs 0.11 expected; zero DiaSources in the discovery
cone). Deepest wide-rung single-epoch flux threshold in the
covered-window record: broadband-r relay power ≳ ~1.3 kW through the
0.1 AU cone at the b(t) = 20.9 R☉ rim sample, z 1000–10,000 AU —
a magLim-referenced threshold statement, **not injection-calibrated**
(no archive injections at the field; no pixels to inject into), so
no exclusion is claimed. SkyBoT asteroid positive control passed
through the frozen chain (0.308″, S_det 31.3). Ledger: grazing rungs
structurally uncovered; ross-154 A saturation-limited (16 visits);
designated follow-up — image-level re-run at the visit/difference-
image release upgrades the unit to an injection-calibrated exclusion.

### 5.14 DASCH

**Description.** The pre-1980 century: Harvard scanned plates
(~430,000, DR7, ~1885–1990; B ~12–16 typical, ~18 best), accessed
through the anonymous Starglass REST API (exposure lists,
century-long lightcurves with per-epoch non-detection limits, FITS
cutouts, per-plate subregion photometry incl. uncatalogued
detections — a catalogue-level-first substrate, new to the
programme alongside §5.13's). The prerequisite backward extension
`crossings/universal_1885_v1` (Earth center, 1885→1993, 37,096
events; shared-era rows match `universal_v1` to 23 s / 3×10⁻⁷ R☉)
revealed the survey's cell: **van-maanen, gj-1276 and wolf-359 are
photosphere-grazing (b ≲ 1.2 R☉) on essentially every annual
crossing of the DASCH century**, with ~15–28 calibrated-exposure
windows per target-channel at ±1 d (~120 across the family +
teegarden/ross-128 at 2.5 R☉) — an order of magnitude more
grazing-rung windows than all previous surveys combined, ending
where the 1980 list begins. §7 old-epoch model budget measured at
recon: ephemeris ≤ 1×10⁻⁵ R☉ vs DE440S; astrometric sensitivity
≤ 2×10⁻⁴ R☉ — the pre-1941 `long_propagation_span` annotations are
conservative bookkeeping. Known-issue set (plate defects, blends,
source splitting, missing lightcurve points) drives the
single-detection vetting design; APASS-B is the science refcat
(ATLAS-g has documented false trends); exposure-interval window
overlap (60-min median exposures) and a `time_accuracy_days` gate
(logbook-dated plates are date-only) belong to the coverage stage.

**Status.** Complete — recon + backward extension, hypothesis
freeze v1.0 (D1–D8) + amendments v1.1–v1.3 (all dev-stage,
pre-confirmatory), coverage (18 searched units), threshold freeze,
dev, blind confirmatory, completeness and report, all 2026-08-26 on
the dev machine. Report `report/dasch_crossings.md`; execution log
`notes/project_history.md` §12; recon
`surveys/dasch-crossings/notes/dasch_recon_2026-08-26.md`.

**Result.** **0 candidates** (blind confirmatory 15 units / 30
trials; 2 exceedances vs 3.33 expected control crossings, both
adjudicated defect-class in the plate pixels). The programme's
first pre-1980 power constraints — grazing cones ≳ 20 kW on the
deep-plate windows (van-maanen/wolf-359 B 1.2 R☉), sub-MW across
the shallow-plate century bulk back to the 1890s, MW-class through
the 0.1 AU cones, over ~10–60 independent annual recurrences per
unit; the century recurrence-stack cell closes clean on every
searched unit. Completeness measured from pooled field-star
recovery: 90 % depth sits 1.5–2.5 mag above the archive's limMag
columns (C1 at catalogue level). Both
positive controls passed (RY Cnc eclipse chain; (7) Iris 1911
recovered at 0.2″). Amendments discovered by dev: APASS refcat
PM-less for high-µ stars (v1.1 routing), DASCH's defect classifier
eats real single-plate transients (v1.2 — `SUSPECTED_DEFECT`
demoted to annotation), sub-limit stars are their own hit
background (v1.3 quiescent gate).

### 5.15 Kepler/K2 footprint intersect

**Description.** A one-afternoon footprint intersect: the K2 ecliptic
campaign fields (2014–2018) and the Kepler prime field against the
crossing list, to proceed to a survey only on a hit (30-min continuous
cadence in the pre-ZTF era). Observer correction: Kepler trailed Earth
by 0.04–1.14 AU (|Δt_ca| median 33 d), so Earth-center is invalid at
every rung and the intersect ran on a spacecraft list
(`crossings/kepler_v1`, Horizons −227), not `universal_v1`.

**Status.** Closed 2026-09-04 with no survey — no hit possible. Note
`surveys/kepler-crossings/notes/kepler_k2_footprint_intersect_2026-09-04.md`;
execution log `notes/project_history.md` §17.

**Result.** K2fov silicon polygons × MAST CAOM actual campaign ranges
→ 0 of 1,435 field-active events on silicon at any b. K2 boresights
stayed 61°–158° from the Sun — never within 22° (37° with the true C14
start) of the spacecraft's anti-sun point — while ≤ 0.1 AU channel
positions sit within 3.5° of it; the prime field is 66°+ off. The
elongation gate, fourth confirmation. Five registry stars (wolf-359
C14, ross-128 C1, ross-154 C7, van-maanen C8, gj-876 C3) have K2 light
curves at b ≈ 0.9–1.0 AU — wide rung, recorded only.

### 5.16 STEREO-A HI-1 sunward survey

**Description.** Second heliospheric substrate for the two sunward
combinations (S1 downlink post-lens, S2 uplink past the Sun), after
LASCO §5.11. STEREO-A HI-1A level-2 frames 2007–2026 (RAL/UKSSDC
processing, NASA SSC mirror; 72″/px, 40-min sums, 630–730 nm) on the
spacecraft-observer list `crossings/stereoa_v1` (Horizons −234). The
camera's 4° inner edge excludes the grazing cones (`not_constrainable`);
the 0.1 AU rungs give a 44-h / ~66-frame arc per window at ε 4.4°–6°,
21 windows per unit, plus a 19-d star-fixed baseline transit — and the
mission's pointing eras (east 2007–2014, rolled **west** 2015-11 →
2023-08, east again) put the arc before or after t_ca respectively.

**Status.** **COMPLETE 2026-09-04** — `report/stereo_hi_crossings.md`.
Chain in one day: observer list + census → recon (header WCS sub-pixel,
0.05-mag 2-D colour-corrected ZP, arc depth V 10–10.5/frame) → coverage
(272/295 events; the 2014–15 gap) → freeze v1.0 (D1–D9, user-approved)
→ dev with two pre-confirmatory amendments (v1.1 level-2 substrate +
24-px margin after a ridge/edge background bias on level 1; v1.2 2°
bright-body proximity veto replacing the frame-wide veto) → blind
confirmatory once → completeness + Pallas positive control → report.
Survey docs `surveys/stereo-hi-crossings/`.

**Result.** **0 candidates blind** (`results/confirmatory_v1.md`): 27
trials over 9 units, 3 exceedances vs 3.0 expected — two extended
coronal-transient fronts in the pixels (gj-1276 S2 2009-01, teegarden
S1 2010-08) and one non-point-source background-offset plateau
(gj-1276 S2 2012-10), all non-recurrent; every S_stack null. ross-128
S1 constraint-only (antipode Tycho-blended, as in LASCO). Dev: 12
trials, 2 adjudicated (gj-908 S1 plateau; ross-154 single-frame
flare-class). Completeness (injections into the real null series):
recurrence-stack m90 V_eq 11.9–12.9 → **downlink 12–29 MW through the
0.1 AU cone, 10-m uplink 3–61 MW** (100–1000× deeper than LASCO on the
same cells); pulse cell open but GW-class (control patches carry
bright-star excursions). Pallas recovered at +0.01 / −0.02 mag on two
arc passages. Honest framing: the "sub-MW downlink" hope for this
substrate was a grazing-cone figure; HI-1 cannot see the grazing
cones. Hand-offs: v2 baseline/offset design notes, WISPR (geometry pass
§5.18, done: the 0.1 AU cells at 14°–54° from PSP; survey §5.19),
yearly refresh of `stereoa_v1` (§5.27).

### 5.17 Spectral-archive family — laser lines in the crossing windows

**Description.** The first non-imaging survey: archived
high-resolution spectra of the deep-family stars taken inside their
channel-A windows (Earth on the Sun–star axis, star at opposition),
searched for an unresolved emission line at rest in the star's frame
(cells 532 / 1064 / 1550 nm and generic). Substrate: ESO phase-3
HARPS, ESPRESSO, NIRPS; CFHT SPIRou APERO products via CADC; CARMENES
GTO DR1 VIS — all anonymous. Unit = target × instrument × event; the
null ensemble is 60 seeded same-star same-instrument out-of-window
spectra (leave-one-out template; PSF-consistent matched-filter
statistic; T = ensemble maximum). Channel B has no spectral substrate
(sky fibres). Docs `surveys/spectral-archives/` (recon note,
hypotheses D1–D9, thresholds v1.0 + amendment v1.1, dev machinery log,
results, adjudication, completeness).

**Status.** **COMPLETE 2026-09-05** — `report/spectral_archives.md`.
Recon → freeze (user-approved) → fetch (961 files, 12 GB) → dev
(ross-154, gj-908; four passes that produced the v1.1 gates:
wavelength-solution, spectrum/pixel SNR, PSF-consistent statistic,
own-error σ floor) → blind confirmatory once (wolf-359, ross-128,
teegarden) → adjudication → injections → report. One machinery defect
surfaced in adjudication (ESPRESSO phase-3 wavelengths are vacuum, not
air; the two ESPRESSO combos were re-run, results unchanged in
substance).

**Result.** **0 candidates.** 24 units, 77 trials, 18 exceedances vs
4.7 expected — 12 sky/telluric ([O I] 5577, OH airglow incl. the
CARMENES 882.95/775.28 nm features of two independent years, O₂ A band,
1.35/1.8 µm water bands, a CO₂-line correction residual at 1569.50 nm
present in all eight wolf-359 NIRPS spectra), 4 from one ESPRESSO
spectrum carrying the Hg-lamp pentad, 2 from a 2-px cosmic hit; the
excess over expectation is the barycentric ensemble's blindness to
observer-frame emission (v2 design item). Completeness (injections,
measured NIRPS SED for wolf-359): **12 W at 1550 nm and 22 W at
1064 nm through the 1.2 R☉ cone** (wolf-359 NIRPS 2025-03-03), 190–220 W
(SPIRou 2021), **37 W at 532 nm through the 2.5 R☉ cone** (ross-128
HARPS 2021-03-17), 53–95 W (wolf-359 NIR) and 400 W (ross-128 SPIRou
1550 nm) through 2.5 R☉, 0.1–7 kW through 0.1 AU on 13 units —
the programme's first line-SED constraints and the first at
1064/1550 nm. CARMENES-VIS and every generic cell beyond 200 % of the
continuum are unconstrained at v1. Hand-offs: raw-frame grazing hits
(HIRES 2010, MAROON-X 2021/2024, CRIRES 2009, APF 2016), SOPHIE 2022
(date-stripped header), SPIRou teegarden 2025 (proprietary to
~2026-11; pre-registered), HPF teegarden 2019–20 and CARMENES NIR
(data asks); van-maanen future-observation recommendation (one
HARPS/ESPRESSO/NIRPS exposure at t_ca ± 7 h on Oct 6–7); yearly
refresh of the NIRPS/SPIRou/ESPRESSO epochs (§5.27).

### 5.18 PSP/WISPR observer-geometry pass

**Description.** The prerequisite named by the sunward geometry study
(§5.11, caveat 1): Parker Solar Probe is an inner-heliosphere observer
(0.046–0.73 AU), so the sunward apparent sources — the star (S2) or
its antipode (S1), at elongation ε = atan(b_e / r_along) — subtend far
larger angles than from 1 AU, and WISPR (two ram-side telescopes,
13.5°–108.5° from Sun center, observing whenever r < 0.25 AU) might see
them. Question: per encounter, do the grazing-cone sources ever cross
the 13.5° inner edge, and what does the 0.1 AU rung look like?

**Status.** **COMPLETE 2026-09-05** —
`surveys/wispr-crossings/notes/psp_wispr_geometry_pass_2026-09-05.md`,
`results/psp_census_v1.json`. Observer list `crossings/psp_v1`
(`xng-c339df63a608`, 7,448 events, 2018-08-15 → 2026-12-01) on the
generic spacecraft machinery with two orbit-forced changes: 10-min
Horizons sampling fetched in yearly chunks (PSP sweeps ~34°/6 h at
perihelion) and 0.5-d coarse bracketing (the star-side and
antipode-side minima of one axis are ~2 d apart around perihelion).
WISPR field edges verified from the NRL Data Users Guide v5; ram side
and half-heights are literature values flagged for header recon.

**Result.** PSP crosses every axis twice per orbit at a per-target
repeating heliocentric radius, so each target's geometry recurs every
encounter (29 encounters E1–E29, q 35.7 → 9.85 R☉). The post-t_ca half
of every window is on the ram side, the pre-t_ca half never is.
**Grazing cones:** `not_constrainable`**.** 1.2 R☉ ≤ 6.9° at any encounter
(mission floor); 2.5 R☉ reaches the edge only at the E22+ perihelia
for gj-1111 (13.7°, 2–3 min per window at b_e 2.42–2.49 R☉) and falls
0.2° short for wolf-359 (13.3°) — a knife-edge on the baffle-limited
F-corona strip; LASCO C2 remains the unique grazing-cone substrate.
**0.1 AU rungs: open, every encounter.** 533 S1 / 561 S2 event-rungs
on 21 targets per channel have a WISPR-I arc (median 23 h per event,
ε 14°–54°, b_e 2.5–18 R☉ on perihelion-side crossings; WISPR-O adds
50°–108° from E6 on); ~700–1,100 WISPR-I hours per encounter; of the
seven-system family only wolf-359-S2, gj-1276-S1 and ross-154-S1
never appear. Public L2 through E27. Hand-offs: the WISPR 0.1 AU
survey (§5.19, recon-first: in-field point-source depth vs r, header
verification, cadence, encounter-local baselines, Mercury/Venus
control); `psp_v1` refresh with the yearly cycle (SPK end 2026-12-01);
Solar Orbiter SoloHI (5°–45°, 0.28–0.9 AU) as the natural fourth
substrate on the same machinery (`--observer solo`, Horizons −144).

### 5.19 PSP/WISPR sunward survey

**Description.** The survey the geometry pass (§5.18) enabled: the two
sunward combinations on the 0.1 AU rung from inside 0.25 AU, where
the rung subtends 14°–65° and its in-beam arc runs through WISPR-I
for ~23 h per event on 21 targets per channel, every encounter.
Approved at the freeze (D0 run) for the cells it adds — the inner
cone from a new vantage, first sunward constraints on new targets,
the first in-band 532 nm sunward substrate, a 5-min pulse cell — at
LASCO-class depth, not as a deepening of the HI-1 cells.

**Status.** **COMPLETE 2026-09-06** — `report/wispr_crossings.md`;
survey docs `surveys/wispr-crossings/` (recon, hypotheses D0–D10,
thresholds v1.0 + amendments v1.1–v1.5, dev/confirmatory tables,
completeness, controls); execution record `notes/project_history.md`
§22. Chain: recon → freeze (user-approved) → dev with five
pre-confirmatory amendments (template completion, proximity veto,
bright-star class, template response, ZP-uncertainty gate + two-frame
pulse persistence) → blind confirmatory once → adjudication with a
Gaia DR3 static-content check and a latitude-curvature test →
completeness → report.

**Result.** **0 candidates blind**: 30 searched units, 588 events,
90 trials, 15 exceedances vs 10.0 expected, all adjudicated
non-promotable — 9 static-content / latitude-curvature systematics
(Tycho-incomplete stars and bright neighbours in the source
apertures, matched by Gaia to ~20 %; near-plane sources on the
brightness ridge), 3 bright-star class (gj-783 S2), 2 arc-exit ramps,
1 extended coronal front (wolf-1061 S1, pixels flat). The LASCO/HI-1
family stacks are null. Completeness: recurrence-stack m90 V_eq
7.8–12.0 (median 9.3) → **downlink 0.08–3.9 GW through the 0.1 AU
cone, 10-m uplink 0.07–4.4 GW** (±0.3 mag flux-scale systematic);
pulse cell opened but systematics-limited (control thresholds 4–70σ).
Vesta E22 recovered at −0.19 ± 0.32 mag. 5 units not covered (sources
above the field's +15° edge), fomalhaut S2 unsearchable. Honest
framing: proximity to the beam buys visibility, not flux — the survey
confirms the 0.1 AU sunward cells at LASCO class from a third
instrument and extends them to 16 new targets; the machinery lessons
(Gaia template with proper-motion propagation, quadratic-in-latitude
controls, measured field extent) are the hand-off for any v2, along
with SoloHI on the same machinery and the E28+ refresh.

### 5.20 Palomar Gattini-IR + WINTER

**Description.** The > 900 nm time-domain cell: Palomar Gattini-IR
(J, 2018–, per-epoch J ≈ 14.5 Vega; ~1.17–1.33 µm, so J only) and
WINTER (Y/J/Hs, 2023–) — **1064 nm needs WINTER Y**. Shallow
(MW-class at 0.1 AU), so framed as opening the cell rather than as a
deep exclusion. The 1550 nm line hypotheses went to the
spectral-archive family (§5.17) instead.

**Status.** **Blocked — recon 2026-09-04**
(`surveys/gattini-crossings/notes/gattini_winter_recon_2026-09-04.md`;
execution log `notes/project_history.md` §15). PGIR DR1 (Data Lab,
2018-10 → 2022-10) is a 2MASS-source light-curve catalog, not images:
forced at 2MASS-epoch positions (every in-era target reads as empty
sky; no antipode has a source in the PSF), float32 times quantised to
0.25 d; the epochal stacks are not public. WINTER has no public
release (login-only portal). Visit lists recorded: 14 units, 80–224
epochs each; wolf-359 B in-window at 2.5 R☉ twice — PGIR imaged every
in-era target-channel.

**Blocker.** A data ask to the PGIR team (user action) for epochal
images or PM-propagated forced photometry, or a future PGIR
epochal-stack release / DR2; for WINTER, a public release (re-checked
at each yearly refresh: IRSA holdings, Data Lab, `winter.caltech.edu`,
the instrument-paper series).

### 5.21 Rubin image-level crossings

**Description.** Image-level re-run of the §5.13 catalog-level unit
(ross-128 B 0.1 AU r, plus whatever later visits add) on Rubin
visit/difference images, with pixel injections — upgrading the
magLim-referenced ~1.3 kW threshold to an injection-calibrated
exclusion. The grazing rungs, structurally uncovered in DP2, re-enter
the coverage stage with each later release. Shares the adapter with
§4.8.

**Status.** Designated follow-up of §5.13; not started.

**Trigger.** The visit/difference-image release (late 2026), then the
Year-1 release.

### 5.22 Gaia DR4 L2-observer crossing list

**Description.** An L2-observer crossing list for Gaia via the
now-generic `register_spacecraft_table_observer` machinery, so DR4's
per-transit epoch photometry can be searched inside the crossing
windows on the same footing as the spacecraft surveys. Companion to the
Pipeline A track test, §4.9.

**Status.** Blocked, not started.

**Trigger.** The DR4 release (not before 2026-12); candidate for the
2027 crossing-list refresh (§5.27).

### 5.23 Solar Orbiter SoloHI sunward survey

**Description.** The natural fourth heliospheric substrate after LASCO
(§5.11), STEREO-A HI-1 (§5.16) and WISPR (§5.19): SoloHI's 5°–45°
field from an observer at 0.28–0.9 AU, on the same spacecraft-observer
machinery (`--observer solo`, Horizons −144) with the PSP lessons —
10-min Horizons chunks and 0.5-d bracketing for a fast observer, a
Gaia template with proper-motion propagation, quadratic-in-latitude
controls, measured field extent.

**Status.** **COMPLETE 2026-09-06** — `report/solohi_crossings.md`;
survey docs `surveys/solohi-crossings/` (geometry pass + recon note,
hypotheses D0–D9, thresholds v1.0 + amendments v1.1–v1.2, dev
machinery log, dev/confirmatory tables, completeness, controls);
execution record `notes/project_history.md` §23. Observer list
`crossings/solo_v1` (`xng-dea7100a725c`, 4,624 events). Chain in one
day: observer list → geometry pass + recon (field model measured from
headers: anti-ram side, orbit-plane frame, 0.49° seams; grazing cones
`not_constrainable`; 0.1 AU rungs open every orbit, pre-t_ca half,
with an off-beam same-tile baseline before every arc) → freeze
(user-approved) → dev with two pre-confirmatory amendments
(saturation mask at 0.85 DSATVAL; structure-noise epoch gate
err ≤ 1.5 units + colour calibrators V ≥ 4.5) → blind confirmatory
once → adjudication with pixel stamps, the Horizons `@-144` census,
aligned stacks and a post-blind median-interpolation diagnostic →
completeness → report.

**Result.** **0 candidates blind**: 18 searched units (of 27), 75
events, 54 trials, 8 exceedances vs 6.0 expected, all adjudicated
non-promotable — 4 control-interpolation artifacts (a bright control
star losing its core to the saturation mask, extrapolated into the
source by the frozen quadratic on one-sided ladders: gj-1111 S2 ×3,
gj-251 S2), 2 uncatalogued moving objects through the patch
(teegarden S2), 2 extended level offsets of alternating sign (gj-876
S1). The LASCO/HI-1/WISPR family stacks are null (van-maanen S1/S2,
wolf-359 S2, ross-128 S2, gj-1276 S1, teegarden S2). Completeness:
recurrence-stack m90 V_eq 9.1–13.9 (median 10.5) → **downlink
14 MW – 1.2 GW through the 0.1 AU cone (median 340 MW), 10-m uplink
33–350 MW (median 89 MW)** — LASCO-class, ~1 mag deeper than WISPR's
stacks (±0.3 mag flux-scale systematic); pulse cell opened but
systematics-limited (V_eq 4.8–7.9). Vesta recovered at −0.05 ± 0.14
mag. 9 units `constraint_only` (arcs inside the sunward-edge band:
saturation plateau in the ≥ 45-s regime, unresolved corona structure
inside ε ≈ 8°); 61-vir S2 saturated. Honest framing: the survey
confirms the 0.1 AU sunward cells at LASCO class from a fourth
instrument, at the other window phase, and extends them to 13 new
targets; the machinery lessons (robust control interpolation, a
brighter control-star mask under saturation, the edge band) are the
hand-off for any v2, with P13+ at the yearly refresh.

### 5.24 SPHEREx crossings

**Description.** The crossings survey on the SPHEREx spacecraft
observer list (`crossings/spherex_v1`, 1,130 events). SPHEREx is an
elongation-90° surveyor like WISE, so Gates I and II of §5.3 apply
(channels B and A-0.1 invisible in principle; annual windows defeat
the temporal pseudo-window family).

**Status.** Deferred by decision 2026-08-24: the 1.5-yr archive is too
thin.

**Trigger.** 3+ annual windows per target in the public archive.

### 5.25 S2 1 AU outer-skirt blended search

**Description.** The one sunward cell the 2026-08-25 geometry study
left outside the coronagraph substrate: the S2 uplink (star → relay
beam intercepted by Earth on the anti-target side of the Sun) at the
1 AU rung, whose outer skirt sits at ε ≳ 30–40° and is reachable by a
night-sky archive (ATLAS, §5.9) as a blended search on the star. Recon
(2026-09-07, `surveys/skirt-crossings/notes/skirt_recon_2026-09-07.md`)
reshaped the cell: Earth's axis distance is exactly r_E sin ε, so the
1 AU rung has no temporal signature (ledger only) and the testable
quantity is the **elongation-locked step of sub-1 AU uplink beams**
(0.90 / 0.95 AU, edges 64° / 72°) — on in the skirt below the edge and
again past opposition, off in the quadrature band — above the
measured ATLAS solar-elongation floor (~50°), tested with a null
ensemble of 8 field stars per target in the same exposures.

**Status.** **COMPLETE 2026-09-10** — report `report/skirt_crossings.md`.
Chain: recon → freeze v1.0 (D1–D8, approved 2026-09-07) → 162
reduced-mode full-history ATLAS tasks with server-side proper motion
(~30 h drain) → dev with amendments A1–A4 (airmass-layer gate
rewording; **parallax response correction with a measured PSF-fit
response exponent k = 2.42** — the parallactic displacement peaks at
quadrature and is itself a ~1 % elongation-locked step; the
colour-unmatched rule for the five reddest targets) → cut (16 of 18
targets unsaturated) → threshold freeze (28 searched units, 81 trials,
9.0 expected crossings) → blind confirmatory → adjudication →
completeness (A5: no depth for exceedance trials) → report.

**Result.** **0 candidates.** 24 exceedances vs 9.0 expected — 13
asymmetric (S2-side-only or opposite-sign, which a beam cannot do),
2 single-cycle, 1 ensemble-shared, 8 `retained_ambiguous_colour_unmatched`
(teegarden, gj-3512: red targets with no colour-matched control); the
9 colour-matched targets ran to budget (6 vs 5.0). 90 % depths on the
57 non-exceedance trials: **0.6–24 MW through the 0.90–0.95 AU beams
(S_sym median 6–8 MW; 0.6–1.5 MW on gj-1276 / eps-ind-b)**, a
1.5–6.5 mag contrast below each star — the first constraint on the
sub-1 AU uplink rungs. gj-13157 crowded (G 9.4 star at 14″) and
systematics-limited. 1 AU ledger: every conjunction cycle has 7–31
S2-side nights at ε ≥ 50°.

### 5.26 High-energy archives

**Description.** Chandra / XMM / eROSITA / Swift and Fermi-LAT event
products for opportunistic high-energy coincidence with the crossing
windows and persistent-source tests on the corridors. The UV half of
the original row ran as GALEX/gPhoton (§5.12). Recon (2026-09-07,
`surveys/highenergy-crossings/notes/highenergy_recon_2026-09-07.md`)
reshaped the row: **XMM-Newton (solar aspect 70–110°) and eROSITA
(scan at 90° elongation) cannot observe any crossing window** — both
channels sit at elongation 180° ± 6° — and join WISE/SPHEREx under
Gate I (§5.3); Chandra has no in-window observation in 27 years; the
pointed archives contribute one unit (Swift XRT 4.2 ks + UVOT UV grism
on wolf-359, 2018-03-06, A 0.1 AU). The window substrate is the two
all-sky monitors: **Fermi-LAT** (≥ 100 MeV; every grazing window
2008 → 2018-03 has 1–19 ks of in-FoV livetime, measured from the
weekly spacecraft files; after the 2018-03 solar-array anomaly one
grazing window in three has none) and **Swift/BAT** (half of the
grazing windows and 95 % of the 0.1 AU windows carry ≥ 1 ks within 30°
of the boresight) — a high-energy pulse/burst coincidence cell,
photon-starved and near-background-free at the LAT (7 photons > 100
MeV within 1° in 3 d at the van-maanen antipode; arbitrary
position/time pulls verified through the LAT data server). Corridor
half: catalogue-level and ready — eROSITA-DE upper limits at
arbitrary positions by API (DR1 eRASS1; DR2 eRASS:3 released
2026-07-31, catalogue-only) for the 36 corridors in the western
Galactic hemisphere (UL median 3.4 × 10⁻¹⁴ erg cm⁻² s⁻¹), plus
5XMM-DR15 / CSC 2.1 / 2SXPS-LSXPS / BAT-157m / 4FGL-DR4 cones; sources
within 7′ are common (34/88 corridors at eRASS:3 depth, LMC-direction
corridors crowded), so the screen is an SGL-track test, not a
presence test — first case: a persistent X-ray source 1.0′ from the
teegarden antipode (LSXPS + eRASS:3).

**Status.** **COMPLETE 2026-09-07** — report `report/highenergy_crossings.md`.
Chain in one day: recon → freeze v1.0 (D1–D11) → amendments v1.1
(Jeffreys rate floor; A gj-1276 1.2 R☉ demoted, 145 trials, T = 3.451)
and v1.2 (per-trial ensemble exceedance gate: 9 B van-maanen trials
constraint-only, 3C 279 at 1.8° from the antipode) → coverage (580
windows, 518 searchable) → threshold freeze (55,474 pseudo-windows) →
dev (controls GRB 130427A and the 3C 454.3 flare PASS) → blind
confirmatory → completeness → Swift/BAT constraint-only arm (HEASoft
`batsurvey` in a container, 238 pointings) → XRT look → corridor
screen. Execution log `notes/project_history.md` §26.

**Result.** **0 candidates blind**: 136 calibrated trials, 0
exceedances vs 0.05 expected (lane L 10,778 photons vs 10,517
expected); the 9 constraint-only B van-maanen trials all exceed
through the 3C 279 April-2014 flare, vetoed `known_source`. 90 %
limits (Γ = 2): recurrence-stacked 100 MeV–1 GeV flux 1.4–3.6 × 10⁻⁷
ph cm⁻² s⁻¹ through the 1.2 R☉ cone (**0.12–0.33 MW**), 0.2–0.5 MW at
2.5 R☉, 2–9 MW at 0.1 AU; single-window 0.5–1.2 MW (grazing); 1-ks
pulses Φ50 ≈ 3 × 10⁻³ ph cm⁻². BAT 14–195 keV: 175 grazing windows at
median ~30 mCrab (1.6 / 5.9 MW), no threshold. XRT look: wolf-359 at
3.6× its quiescent rate, flat, nothing pulse-like. Corridor screen: 46
corridors empty within 7′, 212 of 281 sources `fixed_sky` by the
SGL-track test, 69 unresolved (LMC-direction and CSC-only); eRASS1 /
eRASS:3 upper limits at all 36 western corridors. Hand-offs: the UVOT
UV grism to §5.17; the yearly refresh re-pulls the LAT era (data
server, ~90 s per position) and the post-2018 livetime table.

### 5.27 Standing maintenance

The yearly crossing-list refresh (the 2028 window end) as archives
extend, and the items that ride on it:

- ~~Fold the TESS and GALEX rows into the covered-window ledger~~ —
**done 2026-09-10** as part of the programme ledger v2 (§6), which
folds every completed survey; the joint-stage v1 ledger predates TESS
and stays as pre-registered. At the refresh, re-run the ledger build
(§6) after the per-survey re-runs.
- Re-run the teegarden 0.1 AU recurrence test when a future ecliptic
sector covers the target (§5.6); the wolf-359 A NUV
retained-ambiguous burst (§5.12) awaits a future UV mission for its
recurrence test.
- Radio: re-run the v2 CASDA/LoTSS arm
(`surveys/radio-crossings/scripts/coverage_intersect_v2.py`, ~12 min)
plus the v1 VLASS arm (epochs 3.2 / 4.x) — VAST's weeks-cadence over
the wolf-359/gj-1276 antipode field makes the 0.1 AU antipode windows
a recurring test; re-run the intersection if BL's MeerKAT holdings
reach the public archive. The v2-style null-ensemble redesign is the
designated route should the d = 1 wide-beam rung ever be searched
with calibrated error rates.
- Heliospheric: refresh `stereoa_v1`, `psp_v1` and `solo_v1` (SPK ends 2026-12-01 / 2030-11-20;
WISPR E28+ public L2; SoloHI P13+ — perihelion 2026-08-19 — public L2 after 2026-04-10).
- Spectral archives (§5.17): pull the newly public
NIRPS/SPIRou/ESPRESSO/HARPS epochs of the 2026–27 windows (ensembles
and thresholds cached in `runs/spectral-archives/v1/`) and re-check
the SPIRou teegarden 2025-11-08 release.
- Re-probe the near-IR time-domain substrates blocked in §5.20 (WINTER
release; any PGIR epochal-stack release or DR2).
- ~~Propagate impact-parameter uncertainty in the universal crossing
list~~ — **done 2026-09-06** (`sglsurvey/crossings_uncertainty.py`,
`crossings/universal_v1/uncertainty*.{ecsv,npz,json}`; history §24):
128-draw target-state Monte Carlo per event, no rung membership
changes at 95 % (half-width ≤ 0.012 R☉ overall, ≤ 1e-4 R☉ on the
grazing family; `t_ca` σ ≤ 264 s; side 100 % consistent), i.e. the
target-state term sits at or below the declared 0.010 R☉ LEO
observer floor (only eps-ind-b's worst event, 0.012 R☉, exceeds it). Observer state, ephemeris and the model floor
remain `not_propagated`. Re-run per list at the refresh
(`--product <name>`); the 84 boundary rows still at the window edge
are the 2028-end events the refresh will complete.
- Channel A at the 1.0 AU rung keeps its programme-wide deferred
status.



## 6. Programme ledger (goal 5)

**Description.** The cross-survey index of every constraint and open
cell: one harvest file per report transcribing what the report states
(per target where it gives per-target numbers), merged into a single
table keyed by target × channel × rung × band with the programme's
coverage states (`searched`, `constraint_only`, `ledger_only`,
`not_constrainable`, `structurally_open`, `vetoed_known_source`,
`retained_ambiguous`, `no_survey`), native units plus the derived
transmitter power where a report gives one, and the report-section or
results-file citation on every row. Emits the shared Constraint /
AnalysisRun / Candidate records for the Pipeline B surveys, which had
only survey-specific results files. An index, not a re-analysis: no
limit is re-derived. Survey docs `surveys/programme-ledger/`
(`harvest/SCHEMA.md`, `scripts/`, `results/`); records and the full
table under `runs/programme-ledger/v2/`.

**Status.** **COMPLETE 2026-09-10** — report `report/programme_ledger.md`.
Decisions 1–5 approved 2026-09-10 (scope, states, native units,
provenance, ECSV + report); six parallel harvest passes over the 27
reports; 22 of the 28 harvest files regenerate byte-identically from
generators kept in `scripts/harvest/`, six are hand transcriptions.

**Result.** 28,588 rows from 28 sources, **0 candidates**: 20,541
searched, 191 constraint-only, 344 ledger-only, 7,301 structurally
open, 87 not constrainable, 97 no-survey, 16 retained-ambiguous and 11
vetoed rows; 28 AnalysisRun, 21,163 Constraint, 27 Candidate records.
Pipeline A: all 88 endpoints searched, 61 in six surveys. Pipeline B:
102 distinct (target, channel, rung) cells searched over 37 targets, 35
with a published power limit; **14 registry targets have no per-target
crossings row in any report** (82-eri, alpha-cen-a/b, gj-1061, gj-367,
gj-66-a/b, gj-687, kapteyn-star, lp-145-141, proxima-cen, sigma-dra,
struve-2398-a/b — programme-wide statements only); the A 1.0 AU rung
is ledger-only bar one ZTF unit; the S1 grazing rungs are
not-constrainable in every substrate but LASCO C2 at 2.5 R☉. The twelve
retained-ambiguous cells are enumerated (report Table 7). Refresh: re-run
the generators, then `build_ledger_v2.py` and `render_tables.py`.

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

