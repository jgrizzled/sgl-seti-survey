---
title: "sgl-seti-survey — Project history"
date: 2026-08-24
status: "execution record, moved out of notes/project_plan.md on 2026-08-24; append-only"
---

# Project history

Chronological execution record of the surveys and programmes. The plan
(`notes/project_plan.md`) keeps goals, methodology, sequencing, and
current status; this file keeps what was actually run, when, and with
what outcome. Cross-cutting lessons are consolidated in
`notes/learnings.md`. Final reports are in `report/`; the v1 reports
and analysis scripts named below were deleted at v1 retirement (§5)
and live in git history.

Section references of the form `§N` inside the moved text refer to
`notes/project_plan.md` unless otherwise noted.

## 1. WISE/NEOWISE v1 — shakedown and first survey (2026-08-18 → 21; retired 2026-08-24)

The first milestone — a completeness-calibrated end-to-end run of
Pipeline A against the WISE merged L1b products — was executed
2026-08-18 → 21. All eight shakedown steps (hypothesis freeze, pilot
registry, snapshotted discovery, precise WCS+mask pass, catalog
screening, forced-photometry shift-and-stack, injection calibration,
report) passed their success criteria, and the pipeline was then run
over the full universal target list. Canonical documents (v1 reports
and the standalone v1 hypothesis file now in git history):

- Hypothesis freeze v1.0 + addenda: `surveys/wise/hypotheses.md`
- Per-step status and result summaries: `surveys/wise/README.md`,
  `surveys/wise/results/*.md`
- Registry curation and corridor notes: `surveys/wise/notes/`
- Reports: `report/wise_shakedown_v1.md` (pilot, 12 endpoints),
  `report/wise_survey_v2.md` (47 endpoints / 38 corridors, 8 pc),
  `report/wise_survey_v3.md` (88 endpoints / 77 corridors, universal
  v2 / picky-network portfolio, corrected depth scale)
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

**Recalibration (2026-08-21):** after WISE batch 5 (universal v2
portfolio: 88 endpoints / 77 corridors,
`surveys/wise/results/batch5_summary.md`), all 176 sample tensors were
rebuilt under the total-flux kernel and recalibrated as v0.2.1
`run-1b2e86e9219e`; the exceedance census was rerun uniformly over all
batches by the new `adjudicate_exceedances.py`; `report/wise_survey_v3.md`
(superseded v2) carried corrected depths (survey-wide Rx medians
W1 13.9 / W2 13.0 / W3 10.6 / W4 7.1; best W1 15.9), the regenerated
figure, and rescaled physical limits (300 K radiator ≥ ~99 km at
550 AU via W3); the published artifact was updated.

**Scientific review (2026-08-21):** recommended major revision before
the survey is presented as a calibrated exclusion experiment; the nine
findings are condensed in `notes/learnings.md` §1 (full review and the
`v2_plan.md` response in git history:
`git show 8fb226f:surveys/wise/scientific_review.md`). The four
overreaching claims were withdrawn by a v3.1 erratum before any v2
computation. The findings applied to every v1 survey — which is why
the v2 programme (§5 below) redid them all.

## 2. ZTF v1 — second adapter (2026-08-20 → 21; retired 2026-08-24)

**Why ZTF second.** The plan's second adapter existed to stress the
interfaces with a genuinely different archive and to open a physical
cell WISE left closed. ZTF did both: ground-based observer state,
CCD-quadrant footprints, three filters, seeing/airmass/moon quality
dimensions, and reference-subtracted difference images as a first-class
product; scientifically, 1" pixels remove the confusion floor, g/r/i is
the reflected-sunlight band, and 8 years × ~1,000 epochs per corridor
give ~16 parallax cycles for the phase test and enough cadence to fit
continuous relay distance + residual motion. It is also ~90% of a
Rubin adapter. Probed 2026-08-20: IRSA IBE `ztf/products/sci` returns
1,132 frames (g 483 / r 494 / i 155, 2018-06 → 2026-06) at the
Ross 128 corridor. Coverage: corridors at Dec ≳ −30°, i.e. 52 of 76
universal-list systems.

**Pilot (3 corridors) — complete 2026-08-20.** Report
`report/ztf_pilot_v1.md`, stage detail
`surveys/ztf/results/pilot_v1_summary.md` (git history). Outcome: 80
Constraint records, 0 candidates, m90 ≈ 21.2–21.5 AB; lessons in
`notes/learnings.md`. The steps as executed:

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
   same shift-and-stack machinery (the pipeline's first positive
   control).
7. Report; then scale to all 52 in-footprint systems via a ZTF overlay
   of the universal list.

Success criteria were those of the WISE shakedown (reproducible joins,
declared recovery probability over a nontrivial flux × distance region,
controls inside thresholds, constraints traceable to usable pixels),
plus: the adapter interface in `sglsurvey/adapters/base.py` survives
without WISE-specific leakage, and the ground-based observer state
validates against sglseti's terrestrial-site model at the arcsecond
level. All met (topocentric shift ≤ 0.016″; interface unchanged, two
estimator generalisations in the shared calibration).

**Scale-up v1 — complete 2026-08-21.** Overlay
(`surveys/ztf/targets/overlay_v1.md`: 62 visible corridors, 29 ok /
18 edge / 15 gap) → full chain over 69 endpoints: 147,063 exposures,
2,640 Constraint records, 24 exceedances, 0 surviving candidates after
stage-7 adjudication and the survey-wide leave-one-out look-elsewhere
test (24 observed vs 41 expected chance exceedances). m90 ≈ 22.5 AB
median (g/r), 23.6 best. Report `report/ztf_survey_v1.md` (git
history). Deferred to v2: residual-aware variance / ZOGY score images
(≈1.5 mag headroom), per-corridor references for edge strips, fainter
second control. Gap-graded (15) + invisible (15) corridors handed off
to SPHEREx.

## 3. SPHEREx v1 — third adapter (2026-08-20 → 21; retired 2026-08-24)

SPHEREx quick-release data are served at IRSA (`spherex.obscore` TAP
view, collections `spherex_qr2` and `spherex_qr2_deep`; QR1 rows
superseded by the QR2 reprocessing). Probed 2026-08-20: ~300–650
Level-2 spectral images per sky position for 2025-05-24 → 2026-08-11
outside the deep fields, ~7,000 inside them. It was the right third
adapter rather than second: its 0.75–5 µm, 6.15"-pixel regime largely
overlaps W1/W2 physically, and 14 months yield only ~2–3 parallax
phases, but it added (a) per-pixel wavelength (linear-variable filter,
R ≈ 40–130) and hence spectral discrimination of any candidate against
a stellar SED, (b) coverage of the 15 southern-corridor systems ZTF
cannot reach (Dec < −28°), and (c) a third, spectral-image-shaped data
product (multi-extension IMAGE/FLAGS/VARIANCE/ZODI/PSF-cube/WAVE-table
files with SIP WCS, surface-brightness units and a space-based
observer) to test the adapter interface. Running the SPHEREx pilot in
parallel with the ZTF scale-up and the last WISE batch was deliberate:
both of those were limited by IRSA download throughput, and SPHEREx
products are additionally mirrored in a public S3 bucket
(`nasa-irsa-spherex`) that supports byte-range reads at ~20 MB/s, so
the adapter fetches cutouts there and leaves IRSA's bandwidth to the
other two surveys.

**Pilot (3 stars) — complete 2026-08-20.** Report
`report/spherex_pilot_v1.md` (git history). Outcome: 288 Constraint
records, 0 candidates (8 exceedances, all vetoed), m90 ≈ 19.3–20.8 AB
in all six detectors; flux scale verified on 13,657 star measurements
to ±0.2 mag. Lessons in `notes/learnings.md`. Sub-project
`surveys/spherex/`, adapter `sglsurvey/adapters/irsa_spherex.py`,
hypotheses `surveys/spherex/hypotheses.md` v1.0 (WISE v1.0 physics,
Earth-centre observer with the ≤ 0.02" LEO offset carried as a budget
term, L2 FLAGS fatal template, per-detector bands D1–D6 with the
per-epoch wavelength carried on every sample). Corridors chosen from
the ZTF-inaccessible set, reusing registry entries unchanged:

| endpoint        | why                                                                                                                                                                                              |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `lalande-21185` | nearest star ZTF cannot reach (2.55 pc); corridor at Dec −36°, b −68° (clean field); 349 exposures                                                                                               |
| `gj-687`        | corridor 3° from the south ecliptic pole, inside the SPHEREx deep field: 7,308 exposures, continuous cadence — the parallax-phase and spectral-sampling stress case; engineering-backbone basket |
| `sigma-dra`     | engineering-backbone top pick (5.8 pc, quiet G9V); corridor at Dec −70° outside the deep field; 643 exposures                                                                                    |

Steps mirrored the ZTF pilot with SPHEREx specifics: (1) freeze
hypotheses; (2) TAP ObsCore discovery snapshotted, nominal footprint
from the ObsCore polygon; (3) precise pass on S3 range-read cutouts
(IMAGE, FLAGS, VARIANCE, ZODI, nearest PSF plane, wavelength table)
against the SIP WCS + FLAGS usable-pixel footprint; (4) layer-1
screening of bright single-epoch sources along the track against
CatWISE2020 (SPHEREx QR2 publishes no source catalogs);
(5) matched-filter forced photometry in surface-brightness units
converted to µJy via the per-exposure pixel solid angle, joint z × µ
stack with phase-split and offset controls, injection calibration per
detector; (6) positive control: flux-scale recovery of catalogued
2MASS/CatWISE stars in the same cutouts through the same estimator
(the lesson of the ZTF asteroid control); (7) report.

**Scale-up complete 2026-08-20:** all 15 ZTF-inaccessible corridors /
19 endpoints, 0 candidates (`report/spherex_survey_v1.md`, git
history). **v2 estimator 2026-08-21** (`calib_v3`): epoch variances
calibrated to the static template's residual scatter (+0.44 mag
median) and a joint six-detector flat-spectrum stack (+0.2–0.9 mag):
2,128 Constraints, joint m90 21.2–21.9 AB on clean corridors,
19.5–19.8 AB in the Galactic plane; 37 exceedances all vetoed under
five rules (phase split, single season, role coincidence, template
coverage, bright static neighbour). **v3 estimator 2026-08-21**
(`calib_v4`): 16 offset controls (FAR < 1/16 at no median depth cost)
and a cross-detector season-complement phase test; 21 exceedances, all
vetoed. **All-sky run complete 2026-08-21:** 77 corridors / 88
endpoints, 60,426 exposures, 9,808 Constraints, 0 Candidates
(64 exceedances of 1,226 searches, all vetoed); joint m90 median
20.75 AB, 21.5–21.8 on clean corridors, 16.9–17.9 in the Galactic
plane. SPHEREx became the first archive to cover the full universal
list.

Success criteria: those of the earlier pilots, plus: the adapter
interface absorbs a multi-extension spectral product and a
wavelength-per-sample axis without changes to `base.py`; the flux
scale is verified against catalogued stars to ≤ 0.2 mag before any
depth is quoted. All met.

## 4. Pan-STARRS1 v1 and the joint PS1+ZTF stage (2026-08-20 → 21; retired 2026-08-24)

First non-IRSA archive. Started 2026-08-20 while IRSA throughput was
saturated by the ZTF scale-up and the SPHEREx pilot: PS1 DR2 warps
(2009–2014, grizy, δ > −30°) are served entirely by MAST
(`ps1filenames.py` listing, `fitscut.cgi` cutouts, catalogs API), so
the pilot ran in parallel without touching IRSA. Deliberately the
same three corridors as the ZTF pilot (Ross 128, ε Ind A/B, Proxima)
so PS1 extends each ZTF corridor by a 5–10-year baseline. Sub-project
`surveys/panstarrs/`, adapter `sglsurvey/adapters/mast_ps1.py`,
hypotheses `surveys/panstarrs/hypotheses.md` v1.0 (WISE/ZTF v1.0
physics, Haleakalā observer, IPP mask template 16255, star-calibrated
flux scale), report `report/ps1_pilot_v1.md` (git history).

Pilot result: 1,351 warps → 1,355 usable precise evaluations (617
warps, ~75 % of geometric hits usable — no ZTF-style grid-gap losses)
→ 6,040 DR2 detection matches, nothing track-following → 320
Constraints (m90 ≈ 21.0–21.4 AB g/r/i, 19.9 z, 19.0 y), 6 exceedances
all vetoed single-phase, 0 surviving candidates; asteroid (60000)
recovered at 0″ with magnitudes matching Horizons + solar colours to
≤ 0.1 mag. Lessons (star-calibrated zero points, cross-archive
parallax phase, skycell masks / `CONV.BAD`, `detection.obsTime` lag):
`notes/learnings.md`.

**Scale-up (2026-08-21, complete):** overlay
`surveys/panstarrs/targets/overlay_v1.md` (62 of 77 corridors
δ > −30°; minor parallax phase ≤ 12 % everywhere; 30
calibrator-sparse), batched driver `run_scaleup.sh` with per-batch
purge (≈ 150 GB transient, 18 GB retained), lazy calibration cutouts
(97.6 % of 5,516 flux maps star-calibrated). AnalysisRun
`run-eaa6d89a1ec9`: 24,852 warps → 16,527 usable evaluations →
134,874 DR2 matches → 5,520 Constraints (median m90 g 21.1 / r 20.9 /
i 20.7 / z 19.9 / y 18.9 AB) → 78 exceedances (11 %, the chance rate
by construction), 72 phase-vetoed, 3 adjudicated-vetoed (two
catalogued stars, one non-persistent), 3 marginal (82 Eri rx y,
Fomalhaut rx y, GJ 526 tx z: at-threshold / grid-edge /
single-filter, no catalogued counterpart) carried as qualified nulls
for the ZTF cross-archive test. Report `report/ps1_survey_v1.md` (git
history). New stage-7 tool `adjudicate_candidates.py` (split-half,
static-star, other-band, grid-edge tests) reusable by every adapter.

**Stage 2, joint PS1 + ZTF (2026-08-21, complete):** sub-project
`surveys/joint/` (`joint_ps1_ztf.py`, `marginal_ztf_test.py`), report
`report/joint_ps1_ztf_v1.md` (git history), AnalysisRun
`run-fd75b2c982c2`. Joint weighted stack of the per-archive tensors on
the exact common trajectory (station-kept relay, µ_resid = 0 — the
tensors' µ reference epochs differ, so µ ≠ 0 stays per-archive), ZTF
on the PS1 star-calibrated scale, 8 shared controls, common phase
reference: 138 endpoint-roles / 414 cells, 321 with both parallax
phases (PS1 alone ≈ 0), 3,312 Constraints with m90 ≈ 23.3 (g, r) /
21.5 (i) AB over 2009–2026, 35 exceedances (8.5 %, below chance) all
vetoed (14 single-phase, 16 phase-split, 3 non-persistent, 2
other-band, 2 catalogued stars on i-band tracks at stage 7). The three
PS1 marginal cells were vetoed by direct ZTF forced photometry along
the PS1-fitted (z, µ) tracks (S = 0.2 / 1.4 / −0.1 over 900–1,200
frames). PS1 survey and joint stage closed with 0 candidates.

**Stage 2 v2 (2026-08-21, complete):** PS1 tensors rebuilt on the
common µ reference epoch T0 = 59800 (`run_common_t0.sh`,
`runs/panstarrs/calib_t0_59800`), joint stack over the full 5 × 5 µ
grid (`joint_ps1_ztf_mugrid.py`, AnalysisRun `run-ac08543c5b29`):
3,312 Constraints, on-grid m90 ≈ 23.0 (g) / 22.9 (r) / 21.1 (i), 43
exceedances (10 %, below chance), 4 retained then vetoed at stage 7
(field-wide systematic; faint stars along the track; i-only pair
absent in ZTF g+r) — still 0 candidates. Lesson (common µ reference
epoch / T0 under-sampling): `notes/learnings.md`. The
catalogued-static-source test became automatic in every adapter's
census (`sglsurvey/vetting.py`: PS1 DR2 / ZTF DR24 snapshots offline,
CatWISE via VizieR for WISE; ≥ 3 catalog detections required).

**The TODOs this stage left**, and how they resolved:

- _WISE into the joint stage_ (three-archive test with a 0.5–22 µm
  colour axis; needs WISE tensors on a common µ reference epoch and
  Vega→AB / surface-brightness conventions reconciled) — still open;
  tracked in plan §4 (open items) and `notes/learnings.md` §10.
- _Southern corridors._ 15 of 77 universal-list corridors lie below
  δ = −30° (Lalande 21185, Ross 248, 61 Cyg, Struve 2398, Groombridge
  34, GJ 1221, GJ 338, GJ 625, GJ 687, GJ 251, σ Dra, HD 219134,
  Wolf 1069, GJ 3512, GJ 13157), unreachable by PS1 and ZTF; several
  top engineering-basket targets (σ Dra, HD 219134, Lalande 21185)
  among them. SPHEREx covered them in the near-IR (§3 above); the
  optical cell → the DECam survey (complete 2026-08-24,
  `report/decam_survey.md`; recon, pilot and runs in §5 below).

## 5. v2 programme (2026-08-21 → 24)

**Rationale.** The WISE scientific review (§1) found that the v1
surveys were sound as exploratory searches but not as calibrated
exclusion experiments, and the same constructions — 8/16-offset
thresholds, tensor-level analytic injections, hand-tuned vetoes,
nominal-locus geometry, filename-based provenance — were copied from
WISE into ZTF, SPHEREx, PS1, and the joint stage. v2 repaired the
statistical experiment once, on WISE, and then applied the repaired
design to every survey. The full WISE v2 design and its retrospectives
(`v2_plan.md`, `scientific_review.md`) were deleted at reorganization
and live in git history (`git show 8fb226f:surveys/wise/v2_plan.md`);
the durable design rules, transfer amendments, and lessons are
consolidated in `notes/learnings.md`, and each survey's frozen
decision rule is live in its own `hypotheses.md`.

**Sequence as planned** (each step's execution recorded below):

| Step | What | Exit condition |
| ---- | ---- | -------------- |
| 1 | **Let the in-flight v1 runs finish**: SPHEREx northern run (62 corridors, frozen v3 rules, §3 above) and the joint PS1+ZTF common-T0 µ-grid calibration (§4 above). No new v1 rules or batches start. | both runs summarised and the v1 reports closed with their own "exploratory, pending v2" status note |
| 2 | **Errata for the v1 reports.** WISE v3.1 (withdraws the four claims in `v2_plan.md` §0); the same relabelling — "FAR < 1/N" → rank statement, "90 %-complete exclusion" → threshold sensitivity, vetoes → heuristic review — for the ZTF, SPHEREx, PS1, and joint reports. | errata committed; no numbers recomputed |
| 3 | **Shared-code modifications** (`v2_plan.md` §8.5): `records` fields, `photometry` injection seam, `vetting` flux-consistent static test + parallax-phase test, `geometry` per-epoch observer, new `inject`, `nulls` (promoted from ZTF `look_elsewhere.py`), `manifest`; promote `wise_corridors.py` out of `surveys/wise/`. All additive; v1 scripts keep running. `tests/` with the invariant tests. | tests green against v1 products; v1 pipelines unchanged in output |
| 4 | **WISE v2** per `v2_plan.md` §10 (steps A–I): hypothesis v2.0 freeze with hold-out split → covariance + positive control → null ensemble → image-level injections → frozen rule on the development set → blind confirmatory run → report. | v4 report built from the ledger; review findings 1–9 each closed or explicitly deferred |
| 5 | **Retrospective**: update `v2_plan.md` §12 with what actually transferred and what it cost. | §12 revised |
| 6 | **Other surveys v2**, in the order ZTF → PS1 → joint → SPHEREx (ZTF has the richest cadence and an existing positive control, so it validates the transfer fastest; the joint stage needs both optical v2s; SPHEREx last because its next quick release adds a parallax phase anyway). Each reuses v1 discovery/precise/screen and cutouts; rebuild only null, injection, vetting, geometry check, manifests, report. | each v2 report built from its ledger |
| 7 | **Retire v1**: delete superseded v1 surveys/runs/reports, rename `*-v2` to plain names, remove stale references. v1 stays in git history. | stale-reference grep clean; tests green |
| 8 | Resume archive-family expansion (DECam south, three-archive joint) on the v2 design. | — |

**Rules during the programme:** no v1 result is cited as an exclusion;
nothing in a `-v2` directory may import from a v1 survey directory
(shared code goes through `sglsurvey/`); every v2 hypothesis freeze is
a new version with a declared hold-out before any v2 script touches
the confirmatory set.

**Steps 1–5 done 2026-08-22.** Errata appended to all five v1 reports
(WISE v3.1); shared code in `sglsurvey/` (`nulls`, `inject`,
`manifest`, `corridors`, additive changes to `records`, `photometry`,
`vetting`, `geometry`) with `tests/test_v2_invariants.py` green; WISE
v2 in `surveys/wise/` ran stages A–I under hypotheses v2.1 —
0 candidates on both the development (216 cells) and the blind
confirmatory (488 cells) sets, `report/wise_survey_v4.md` (now
`report/wise_survey.md`); retrospective lessons in
`notes/learnings.md` §2.

**Step 6 done 2026-08-23** via the survey-agnostic engine
`sglsurvey/v2/` (now `sglsurvey/`; profiles in `surveys/{ztf,panstarrs,spherex}-v2/`,
joint in `surveys/joint/`): blind confirmatory runs ZTF 228 cells /
PS1 450 / joint 230 / SPHEREx — 0 candidates everywhere; reports
`report/{ztf,ps1,joint_ps1_ztf,spherex}_survey_v2.md` (now the plain
names); lessons in `notes/learnings.md` §2.

**Step 7 (retire v1) done 2026-08-24** (user amendment: v1 reports
deleted rather than archived, and the final reports renamed to plain
names — `wise_survey_v4.md` → `wise_survey.md`, etc.). What was done,
and the deliberate deviations from the step-7 row:

- _Reports:_ only the five canonical reports remain; all pilots/v1–v3
  reports and `report/figures/` deleted (git history).
- _Surveys:_ each `*-v2` renamed to the plain name. Kept inside the
  renamed dirs because they are still current: the
  discovery/precise/screen/fetch scripts (Pipeline A stages that v2
  reuses and step 8 needs), corridor/targets overlay builders,
  per-archive `notes/`, `surveys/wise/{v2_plan.md,scientific_review.md}`
  (later deleted in the `ff0e46a` reorganization; recover via
  `git show 8fb226f:surveys/wise/v2_plan.md`), PS1's
  `purge_products.py` and SPHEREx's `static_template.py` (generators
  of kept products). Deleted: the standalone v1 hypothesis files after
  their still-active physics cells were incorporated directly into the
  v2 hypothesis documents; all v1 analysis scripts (sample_tensor,
  injection/forced-stack calibration, adjudicators, look_elsewhere,
  figures); and v1 `results/` summaries.
- _Runs:_ v1 tensor products and superseded calibration runs deleted
  (~73 GB: all `calib_*/tensors`, `stack_v1`, SPHEREx `calib_v1–v3`);
  kept because still referenced — coarse/precise records, cutout
  products (rebuild inputs), screening records/snapshots (catalogue
  loaders + the ScreenMatch ledger), `m90_curves.npz` + `records/` of
  the final v1 calibrations (injection windows + supersession
  targets), PS1 `zeropoints.jsonl`, SPHEREx `calib_v4/templates`, ZTF
  `control_v1`. v2 products moved to `runs/<survey>/v2/`. SPHEREx
  supersession links re-derived against `calib_v4` (they had pointed
  at the deleted pilot-era `calib_v1` ledger; `supersedes` is outside
  the constraint identity hash).
- _Code:_ v1-only vetting retired (`static_source_test` proximity veto
  and helpers); all path constants updated; every survey's
  ledger-driven report tables rebuild cleanly; tests green (the one
  v1-tensor regression test now skips, its input being deleted).
- _Exit-criterion note:_ the stage directories keep their `coarse_v1` /
  `precise_v1` / `screen_v1` names — they are the still-current first
  versions of live stages that v2 reuses by design, not superseded
  analyses.
- _Provenance note:_ the freeze files (`surveys/*/configs/v2_*freeze.json`)
  are kept byte-identical to the pre-registered versions — every
  AnalysisRun pins their hashes. The path-rename sweep and the
  2026-08-24 documentation restoration edited the hypothesis documents
  after those freezes, so the `hypotheses_hash` values recorded in the
  freezes refer to the documents as of freeze time (git history), not
  the current self-contained files; future freezes hash the complete
  active hypothesis document (`notes/learnings.md` §9).

**Step 8 started 2026-08-24 — DECam south.** Recon
(`surveys/decam/notes/decam_recon_2026-08-24.md`): all 15 southern corridors have
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
hypotheses v0.1 (first survey with no v1 exploratory phase, v2
discipline from day one). Pilot Pipeline A run (lalande / sigmadra /
hd219134): coarse 306 exposures / 612 evaluations; precise 222 hit
exposures, 412 evaluations, ~all usable (5 i-band 404s); screen 11,430
ScreenMatch (NSC DR2 confirmed time-partial: catalogued epochs end
2017–2019 per corridor — recorded in
`runs/decam/screen_v1/catalog_stats.json`); recurrence triage clean
(dense hd219134 corridor is Galactic-plane, b ≈ +2.5°); single-CCD
image+wtmap products fetched for stage 2 (CCDs selected by the covered
locus — 51/213 arc centres sit in chip gaps; 38 exposures span 2–3
CCDs; 1 exposure lost to a server-side-corrupt image). Pre-freeze
validations (`surveys/decam/results/`): CTIO site model ≤ 1.7 mas vs
independent astropy (PASS); flux scale — MAGZERO is counts convention
and per-frame unreliable (outliers +3.3 mag), per-frame NSC-star
calibration mandatory (star-ZP MAD 0.014–0.10 with 35–6,300
stars/CCD); southern overlay v1 grades all 15 corridors 10 ok /
5 crowded / 0 sparse (`surveys/decam/targets/overlay_v1.md`).

**DECam freeze + dev set (2026-08-24, later).** User approved the
hypotheses §8 decisions; v1.0 frozen
(`surveys/decam/configs/v2_freeze.json`, seed 20260824, new additive
engine mechanism `forced_dev` records the pilots as forced development
corridors): dev = lalande + sigmadra + hd219134 + struve2398 (drawn);
confirmatory 11 corridors / 14 endpoints, untouched. Positive control
chosen and run: (60000) Miminko (SSOIS-discovered DECam coverage, 57
exposures 2013–2019 g/r/i/z, Horizons V 19.1–20.6;
`configs/asteroid_control_v1.json`) — recovered in all four bands
under the no-clip scoring (S_max 216–330, rank p = 0.020–0.041,
throughput deficit +0.25–0.39 mag consistent with trailing loss; the
frozen 5σ single-epoch clip removes individually-detected mover
epochs, same as the ZTF v2 control). v2 geometry stage: 38 pairs, max
σ_xt(99%) 0.13" < 0.5" threshold → 0 cross-track cells. Per-frame
star ZPs (`runs/decam/zeropoints.jsonl`, matched-filter NSC
calibration): 373/391 dev frames, median MAD 0.045 mag; lalande's
sparse field needed the relaxed star window (16.0–21.0,
class_star > 0.7, ≥ 6 stars, matching the control script). Dev
tensors: 10 endpoint-role pairs (lalande 11 epochs/role, hd219134
~29, sigmadra ~54, struve2398 50). Dev chain complete (AnalysisRun
`run-7565b9349301`): 30 cells / 480 Constraints / 314 recovery
curves; null ensemble stable (0 heavy-tail; real exceedances at the
ring expectation, R̃_FWER = 1.25); median m90 (threshold) g 23.28 /
r 23.38 / i 22.81 / z 22.23 AB (Y below the epoch floor) — deeper
than ZTF's 22.5 median as expected; 98 static-flux-consistent vetoes
fired; 0 candidates, 0 retained-ambiguous.

**DECam blind confirmatory run (2026-08-24, complete —
`report/decam_survey.md`).** Rule confirmed from dev with no
amendment; Pipeline A over the 11 confirmatory corridors (16,283
Observations, 3,061 precise evaluations, 93,362 ScreenMatch records,
1,223 usable exposure product sets survey-wide), then the blind run
(AnalysisRun `run-3d3cf9840678`): **confirmatory 14 endpoints /
86 cells, 0 candidates** at family-wise α = 0.05 (R̃_FWER = 1.564;
11 cells R > 1 vs 11.9 expected; top cell wolf-1069/tx/g at R̃ = 1.10,
global p = 0.695); 1,904 injection-calibrated constraints; median
persistent m90 g 23.2 / r 22.1 / i 22.7 / z 22.1 / Y 21.0 AB — the
programme's deepest optical corridor constraints where cells are
constrainable (378/712 confirmatory threshold constraints
`not_constrainable` in the heterogeneous PI coverage); quality-mask
sensitivity (primary/strict/loose) changed nothing. First survey of
the 15 southern corridors, closing the step-8 southern optical cell.

## 6. Crossings programme, steps 1–7 (2026-08-23 → 24; closed)

Execution of plan §3.5 (Pipeline B is now plan §5). Reports are
live in `report/`.

| Step | What | Outcome |
| ---- | ---- | ------- |
| 1 | **Universal crossing list** (`sglsurvey/crossings.py` → `crossings/universal_v1/`, `xng-a09e2db7681d`). Every b(t) minimum for the 88-endpoint registry × both link directions, Earth-center observer, 1980→2028; no beam radii, no b cut; boundary minima flagged. | 16,586 events, 0 invalid; 960 with b < 0.01 AU; ~semiannual minima (annual, b ≈ 1 AU for the 5 ecliptic-pole targets) |
| 2 | **ZTF crossings survey** (`surveys/ztf-crossings/`, report `report/ztf_crossings.md`). Channels: A = uplink interception at the star near opposition; B = downlink pre-lens interception at the antipode (z-track 550–10,000 AU); sunward combinations out of scope. Coverage intersect, saturation cut, threshold freeze v1.0 + dev-driven amendments v1.1 (parallax-factor systematics template, same-rung pseudo-windows) and v1.2 (empirical variance rescale), confirmatory run, stamp-response injection completeness. | **0 candidates.** 36 searchable units, 3 exceedances vs 4.0 expected control crossings, all adjudicated. Depths: B median m90 ≈ 21.8 (relay power ≳ ~130 W through the 2.5 R☉ solar cone during covered windows); A contrast-limited m90 ≈ 15.8–17.6 (10-m uplink ≳ 56–400 kW), wise-0855 empty field ≥ 22 |
| 3 | **PS1 crossings survey** (`surveys/ps1-crossings/`, report `report/ps1_crossings.md`, executed 2026-08-24). ZTF channel constructions on the warp-direct star-calibrated substrate (PS1 has no difference images), grizy, era 2009–2015 (disjoint from ZTF, extends the covered-window record back); A 1.0 AU rung constraint-only by freeze; ZTF v1.1/v1.2 rules adopted at freeze — dev run (teegarden + wolf-359) validated them with **no amendment**; saturation levels verified at dev (true saturation ≥ 1 mag brighter than frozen E). | **0 candidates.** 11 searched B units (channel A entirely constraint-only via the frozen gates — single-phase cadence starves the pseudo-window controls); 2 exceedances vs 2.1 expected: ross-128 r vetoed (flux-consistent static, DR2 *stack* catalog), gj-1276 i retained-ambiguous (night-consistent i ≈ 22.8, no stack counterpart, single covered window — recurrence handed to the joint stage). m90 median 21.6 AB (4 units grid-censored ≥ 22); downlink ≳ ~280 W (g median) to ≲ 37 W (censored i) through the 2.5 R☉ cone in 2010–2014 windows. ~40 % of covered unit-rows mask-lost to correlated CONV.BAD chip-gap attrition, **including the van-maanen b = 0.28 R☉ grazing event** — the grazing rung ends unconstrained (nominal-covered / mask-unusable in the ledger) |
| 4 | **Spacecraft-observer crossing lists** (executed 2026-08-24). `sglsurvey/crossings.py` generalized to per-observer runs (`--observer wise\|tess\|spherex`, generic `register_spacecraft_table_observer` in `sglsurvey/geometry.py`); TESS SSB vectors fetched from JPL Horizons (−95, 6 h sampling, raw response + table sha-pinned under `crossings/observers/`); WISE reuses the v2 L1b-header spacecraft table; SPHEREx runs Earth-center with the LEO offset declared as a budget (a coarse table of a 95-min orbit would be false precision). | **Done — manifest parity with `universal_v1`.** `crossings/{wise,tess,spherex}_v1`: 5,170 / 3,390 / 1,130 events, 0 invalid. TESS validation (van-maanen): matched events shift \|Δb\| ≤ 0.33 R☉ (median 0.18) and \|Δt_ca\| ≤ 3.2 h vs Earth-center — the correction the step existed for; the 13.7-d HEO wobble adds shallow extra minima for the ecliptic-pole b ≈ 1 AU family (physical, every-minimum rule) |
| 5 | **WISE crossings survey** (`surveys/wise-crossings/`, report `report/wise_crossings.md`, executed 2026-08-24). Chain ran hypotheses freeze → elongation gate → coverage intersect (45 TAP cones, 91 covered events / 43 targets, rich 211–343-epoch baselines) → saturation cut (W1/W2 searchable = 11 faint targets; every ordinary 10 pc M dwarf saturates) → threshold freeze v1.0 → v1.1 (control pool rescaled to the rung's ~116 d windows, pre-pixel). | **Structural null: 0 searchable units, 0 trials, no pixel searched.** Gate I (elongation theorem): channels B and A-0.1 are invisible in principle to an elongation-90° surveyor — 0/220 and 0/101 events; only A-1.0 survives (102 events). Gate II (control-geometry theorem): the rung's annual ~116 d windows defeat the 8-offset pseudo-window family even after pool rescaling (2–6 achievable; merged 232 d windows admit zero) — the third independent confirmation that the d=1 wide-beam rung cannot be temporally controlled under this design. The 91 covered-window rows enter the step-6 ledger as coverage-without-statistic. Two freeze defects (nonexistent qa_status value; mis-scaled pool) were caught and amended pre-pixel — the ordering discipline worked |
| 6 | **Joint crossings stage** (`surveys/joint-crossings/`, report `report/joint_crossings.md`, executed 2026-08-24 under pre-registration `plan.md` v1.0 + amendment v1.1). Unified covered-window ledger (593 rows: ZTF 86 / PS1 316 / WISE 191; 42 multi-archive target-channels; full disposition census) + the programme's one designated follow-up: the gj-1276 recurrence test — frozen ZTF channel-B stack at z fixed to 550 AU over 11 ZTF event-bands (2019–2026, retained cutouts), single pre-registered joint trial, bad-epoch anomaly removed by the established single-epoch clip (record-based amendment, both versions reported). | **Programme closes with 0 candidates.** J = 0.15 vs T_J = 1.12 (no exceedance); joint recovery at 22.8 AB = 0.94, joint m90 = 22.93 → **flat-SED persistent-relay interpretation of the PS1 anomaly refuted at 90 %**; surviving interpretations narrowed to red/line SEDs (untestable in current optical archives) or non-persistence — not promotable. Programme tally: 47 searched trials, 5 exceedances vs ~6.1 expected, structurally open cells enumerated in the report |
| 7 | **Radio-scope decision** (`surveys/radio-crossings/`, report `report/radio_crossings.md`, executed 2026-08-24). Metadata-only intersection: BL Open Data API (93 position cones at stars + antipodes, beam-classified, snapshotted) and VLASS per-tile epoch dates × the universal crossing windows. | **Decision: geometry-only — no radio signal search in this repo.** Narrow-rung on-star radio coverage: zero (the one cell where re-analysis would beat BL's blind search has no data); antipodes: zero BL pointings within 0.5° in-window (the workhorse channel is archivally virgin — a future-observation recommendation); wide rung: 143 in-beam BL obs / 22 events / 18 targets, already answered by BL's published nulls (cited). Micro-follow-ups handed off: 4 strict VLASS narrow-rung on-star epochs (+2 tolerance-edge incl. one antipode) as a six-cutout quick-look check; 3 teegarden APF optical spectra 5.3 d inside a narrow window (spectral-archive family, plan §3.6 row 10). **Crossings programme closed, steps 1–7** |

Lessons carried forward: `notes/learnings.md` §8.

## 7. TESS crossings survey — execution log (2026-08-24; complete)

Reopened 2026-08-24 as the first item of the crossings
archive-expansion queue (plan §5.8). The coverage gate
(`surveys/tess-crossings/results/coverage_gate_v1.md`) PASSED
(corrected same day for a HEASARC sector-table year typo that had
falsely doubled the grazing count): the Ecliptic sectors (42–44, 71,
91; 600 s/200 s FFIs) hold complete ingress→egress light curves of two
grazing-rung events (wolf-359 b = 0.55 R☉ s42; teegarden b = 1.04 R☉
s91 at 200 s) plus two complete near-grazing channel-A windows
(gj-1276 b = 0.97 R☉, van-maanen b = 0.52 R☉) — the first archive able
to resolve crossing light curves and test the pulse-period cell. The
era's two deepest grazes (van-maanen 0.078 R☉, gj-1276 0.49 R☉, 2022
spring) fell between sectors and stay open.

Executed through the threshold freeze 2026-08-24
(`surveys/tess-crossings/`): TESScut recon + hypotheses v1.0 (channel
A constraint-only — ≥9 d windows in 27 d sectors defeat the temporal
8-control family, 4th appearance; channel B = discovery via spatial
ring controls), per-cadence coverage refinement (8 cubes; grazing
windows hold 81–405 clean cadences), TIC cut (gj-908 excluded
T=7.10), threshold freeze v1.0 (7 B units × 2 statistics — chord
amplitude + pulse max — = 14 trials, 1.6 expected control crossings;
dev = ross-128, all grazing units fully confirmatory; A = 4
reference-only light-curve rows). Dev stage done → freeze v1.1: the
science-array usability screen found TESScut serving collateral CCD
pixels — 3/8 cutouts off-science incl. the ross-128 dev unit → 6 B
units (all confirmatory), 12 trials, 1.3 expected crossings;
machinery validated on the zero-trial A rows (teegarden star
recovered to +0.07 mag; wolf-359 B field ZP scatter 0.066; blended-A
systematics k up to 292 confirm the zero-trial A design).

Confirmatory run 2026-08-24
(`surveys/tess-crossings/results/confirmatory_v1.{json,md}` +
`adjudication_v1.json`): 6 B units × 2 statistics, nested z family
explicit. **No candidates.** The pulse statistic — the survey's new
cell — is null in all six units (first pulse constraint, 200 s →
window, on fully resolved grazing crossings). The chord statistic is
systematics-dominated (sector-scale scattered-light drift; S, T in
the tens; one star-contaminated control at 927): 3 exceedances vs 1.3
expected — teegarden 2.5 R☉ vetoed by the frozen chord-shape test
(split-half ratio 0.31, background anti-correlation −0.785);
wolf-359 2.5 R☉ adjudicated control-crossing (margin +0.46, nested
windows contradict a real source: 0.1 AU superset S = −19); teegarden
0.1 AU retained-ambiguous, non-promotable (truncated window weakens
the shape test; same sector-end cadences as its vetoed sibling;
recurrence requires a future ecliptic sector — the designated
follow-up). Design lesson recorded: an undetrended chord statistic
saturates its error budget with drift — a v2 needs a differential/
detrending layer.

Completeness 2026-08-24 (`results/completeness_v1.json`,
`scripts/completeness_v1.py`): exact stamp-response injections with
the SPOC per-camera/CCD PRF (4 files fetched from
archive.stsci.edu/missions/tess/models/prf_fitsfiles, nearest grid
point to each cutout's CCD position, sha-recorded under
`runs/tess-crossings/products/prf/`) through the identical kernel
chain; 200 draws per unit per temporal model (d = 1 chord and
d = 0.1 boxcar on the frozen per-unit period grids), magnitude grid
T 10–18, recovery against each unit's actual frozen thresholds, full
2 × 2 model × statistic recovery matrices. Finding C1 (flux scale
only, no statistic touched): the dev-stage zp_check used a hardcoded
1.5 px kernel where the freeze requires the identical per-cube
kernel — re-measured through the unit kernel with injection
responses normalized at the calibration reference. Wolf-359
unit-kernel ZP 20.80 (scatter 0.112, pass); teegarden s91 20.03
(scatter 0.493, FAILs the 0.2 gate → declared ±0.5 mag scale caveat).
Depths (natural pairing): wolf-359 chord ≥ 18.0 grid-censored at
1.2/2.5 R☉ (0.1 AU 15.2, star-contaminated control), pulse
16.5/16.3/12.1; teegarden chord 17.5/17.6/≥ 18.0, pulse
14.8/14.7/14.9; A reference rows (S = 5, threshold-free) chord
≥ 18.0, pulse 14.1 (teegarden) / 13.2 (van-maanen, ±0.3). Pulse
recovery plateaus 0.93–0.98 at all brightnesses — real duty-cycle
miss probability (a d = 0.1 train can land entirely in
straylight-masked cadences), included in the curves. Physics:
persistent downlink through the grazing cones excluded above
~0.6–4 kW during the resolved crossings; first pulse-cell limits
~2 kW / 1.3 MJ per ≥ 1-cadence pulse (wolf-359 1.2 R☉).

Survey report `report/tess_crossings.md` written 2026-08-24 —
**TESS crossings survey COMPLETE, 0 candidates**; the teegarden
0.1 AU retained-ambiguous row awaits future-ecliptic-sector
recurrence (plan §5.7 maintenance, alongside folding the TESS rows
into the covered-window ledger at the next refresh). Queue item 2
(ATLAS + ASAS-SN) is next.

## 8. Joint Pipeline A v3 — WISE colour axis (2026-08-24; frozen, in flight)

The plan §4 leading open item executed up to its freeze, in place in
`surveys/joint/`. (a) **WISE engine port**: `surveys/wise/profile.py`
+ `run.py` bind the flattened engine to the v1 WISE inputs under the
joint conventions — T0 = 59800 (v2 tensors at 57800 cannot be
relabelled), AB on ZP 25 via frame Vega MAGZP + Vega→AB offsets (W1
+2.699 / W2 +3.339; a flat-Fν source now has equal tensor flux in
every band of every archive), 9-node µ grid (0.25″/yr; PS1 grid =
[::2] subgrid — the PS1 T0 under-sampling lesson applied), W1/W2
only; v2 observer table, PRF grids and covariance envelopes reused.
Full rebuild: 146 tensors / 69 endpoints / 62 corridors in ~1.3 h on
7 workers, zero failures. (b) **`joint.py` generalised in place** (v3,
`runs/joint/v3/`): optical statistic and family byte-unchanged from
v2; added per-cell W1/W2 annotations at the joint node, the
colour-consistency veto — deficit D_b = (f_pred − f_w)/√(σ_w²+σ_pred²)
≥ 5 in every usable W band, σ_w = max(1/√B, ring MAD) so the empirical
confusion floor caps the significance — charged to the final-candidate
curve, and a `wise-inject` stage running the WISE injection chain once
per optical band family with the family's exact PS1/ZTF union window.
Verified on wolf-359: all 400 W1 draws pair bit-exactly (z, mag, µ,
model) with the frozen ZTF zg injection products — the j-th injection
is one physical source in all three archives with no change to the
frozen optical products. (c) **Freezes before any confirmatory
product**: `surveys/joint/hypotheses.md` v3.0 (self-contained) +
`configs/v3_freeze.json` (sha256:bb68869c…, superseding v2
9bb5791f…); WISE operational freeze `surveys/wise/configs/
v3_freeze.json` (joint split pinned; `hypotheses.md` §9 records the
rebuild as operational-only). Next: three-family injections (running),
then dev nulls (false-veto-rate check) → confirmatory once.

**Stack-regime positive control (2026-08-25, §8 continued).** The
learnings §10 item-1 control executed and passed: asteroid (220000) on
nights with predicted V 21.9–23.05 — below the single-frame limits, so
the frozen |S_e| ≤ 5 clip removes (almost) nothing and only the
ephemeris-weighted stack can recover it. ZTF (35 frames, all sub-clip):
zr S 9.6 / R̃ 2.90, zg S 5.3 / R̃ 1.95, throughput +0.23/+0.29 mag vs
Horizons + solar colours. PS1 (25 warps; the one 5.6σ frame clipped by
the rule): i S 5.8 / R̃ 1.77 at −0.14 mag; z an honest non-detection
(median single-frame S/N 0.34 × √8); g/r/y too few epochs standalone.
Joint g/r/i families through the exact v3 statistic: R̃ 2.21 / 3.04 /
1.77 — every family above the dev threshold 1.54, i.e. the control
would be promoted by the family rule. Scripts
`surveys/{ztf,panstarrs}/scripts/asteroid_stack_fetch.py` +
`surveys/joint/scripts/asteroid_control_stack.py`; summaries under
`runs/{ztf/v2,panstarrs/v2,joint/v3}/control/220000_stack/`.

**v3 development pass (2026-08-25, §8 continued).** Full dev chain run
(nulls all three masks → completeness → adjudicate): **0 candidates**;
the optical machinery reproduces v2 bit-identically (102 cells, 1 void,
R̃_FWER 1.5401, 1,632 constraints, median m90 g 22.11 / r 22.05 /
i 21.33 threshold — equal to the v2 dev files to the digit). Colour
axis validated: all 102 cells W1+W2-usable; on real (null) peaks the
deficit statistic spans median −0.49 to max 2.53 (no would-fire vetoes,
ν = 5); on injections the measured false-veto rate is 1/23,117
(4×10⁻⁵) — the one firing a mag 16.8 long-block source whose WISE
on-window under-overlaps (× PRF throughput), a regime already above the
optical single-epoch clip limit and hence excluded from the fitted
constraints (0/23,116 in the fitted population). Ring confusion floor
dominates σ_w in bright corridors as intended (alpha-cen σ_ring ≈ 164
vs σ_pix ≈ 5 flux units). Records under runs/joint/v3/. Next: the
confirmatory set, once.

**v3 confirmatory run and report (2026-08-25, §8 closed).** The blind
confirmatory set analysed once (nulls → completeness → adjudicate):
**0 candidates** — 45 endpoints / 230 cells, 5 void, R̃_FWER 1.778,
14 R > 1 vs 25.2 expected, 3,648 constraints; every optical number
bit-identical to the v2 confirmatory run (medians to the digit), as
the freeze declared. Colour axis on the confirmatory set: 230/230
cells W1+W2-usable, real-peak deficit D median −0.33 / max 3.73, zero
would-fire vetoes; injections 5/46,302 vetoed, all at mag 16.8–18.8 —
above each cell's optical bright limit, so **0 false vetoes in the
fitted population** (mechanism: PRF throughput × on-window overlap,
confined to the catalogue-layer regime). Ledger-generated tables
`surveys/joint/results/report_tables.md` (report stage extended with
guarded colour lines); survey report `report/joint_ps1_ztf.md`
rewritten in place as the v3 three-archive report, superseding the v2
text (git history). Plan §4 updated; the §4 three-archive open item
and learnings §10 items 1 and 3 are closed.

## 9. Heliospheric sunward-channel geometry study (2026-08-25; plan §5.8 item 9 — answer YES)

One-page geometry study answering the question frozen at the ZTF
crossings freeze: are the two sunward (direction × side) combinations
*only* visible at small solar elongation? Executed with ATLAS (§5.8
item 2) paused for the server migration. Script
`surveys/heliospheric-crossings/scripts/sunward_geometry.py` over
`crossings/universal_v1` (16,586 events); products
`results/sunward_geometry_v1.json` + study note
`notes/sunward_geometry_study_2026-08-25.md` (same subproject).

**Answer: yes — recon justified.** The four combinations partition
cleanly: searched channels A/B have apparent sources at elongation
≥ 91.7° (median ~150°); the sunward pair — S1 = downlink post-lens
(outbound, target side; apparent source = the solar-limb graze point
at 19′–40′ for the 1.2–2.5 R☉ rungs, independent of Earth's in-beam
offset) and S2 = uplink past the Sun (inbound, anti-target side;
source = the star at ε ≈ arcsin(b/1 AU), occulted at b < 1 R☉:
144 events / 3 targets) — is confined to ε ≲ 6° for every frozen rung
except the 1 AU uplink outer skirt (ε ≳ 30° twilight-visible, recorded
as a low-priority channel-A extension). Analytic forms validated
against the frozen list to ≤ 0.032° on the b ≤ 0.1 AU population
(673 events, all ε ≤ 3.29° at t_ca). Census: LASCO-era (1996–) S1
grazing windows 128 (≤ 1.2 R☉) / 160 (≤ 2.5 R☉) on the
van-maanen/wolf-359/teegarden/gj-1276 (+ross-128) family — ~4× the
whole ZTF-era channel-B count, minutes cadence, on exactly the family
lost to PS1 chip-gap masks and TESS sector gaps. Sensitivity framing:
opens the sunward cell at MW-class power (V ~ 8–13 substrates), not
deep exclusion. Caveats recorded: L1/SOHO-observer crossing list
(halo cross-track ~R☉) precedes any grazing freeze; STEREO-A needs its
own observer list; WISPR needs a PSP-observer geometry pass; S1
post-lens annulus beam profile is a hypothesis-freeze question.
Recommended recon order LASCO → STEREO HI-1 → WISPR (deferred). Plan
§5.8 row 9 updated; recon does not block on the server migration.

**LASCO reachability recon (2026-08-25, same day; recon note
`surveys/heliospheric-crossings/notes/lasco_recon_2026-08-25.md`,
star check `scripts/recon_star_check.py` →
`results/recon_star_check_v1.json`; sample frames
`runs/heliospheric-crossings/recon/`).** All services anonymous and
reachable: SDAC `pub/lasco_level05/` (1996-02 → 2025-02, YYMMDD/c2|c3
plain-HTTP daily trees), NRL `lz/level_05/` mirror (current through
2026-06; bare top-level listing times out — use dated paths), NRL
`lz/level_1/` calibrated tree (1996 → 2017-08-31 only: BUNIT MSB,
derolled solar-north-up, time-corrected), JPL Horizons −21 (SOHO SSB
vectors at both era ends). Product facts from real frames: level 0.5
is helioprojective-only WCS, raw DN; C3 Clear 56″/px (contains
532 nm), C2 synoptic Orange 11.9″/px (532 nm just outside); ~12 min
modern cadence both cameras; one live silent-truncation fetch failure
(size+FITS-verify mandatory). Go/no-go validation on a real 19 s L1
C3 frame: synthetic celestial WCS (Sun RA/Dec + Meeus P-angle +
CROTA; orientation convention resolved empirically = East left, solar
north up) pattern-matched 3/3 bright Taurus stars at 0.0/0.4/0.7 px —
det S/N 1124 (V 1.65) / 759 (V 3.0) / 411 (V 4.3), crude single-frame
depth ~V 8–9 — and the common translation (−15.6, −2.6) px = 14.6′ ≈
0.93 R☉ projected is the SOHO halo transverse offset measured
directly, confirming the observer-list prerequisite. Pre-freeze open
items in the note (build `crossings/soho_v1` + census re-run; C2
occulter fate of the 1.2 R☉ rung; coronal-background control design;
0.1 AU-rung volume policy ~0.5 TB at full cadence → frozen
subsampling; 2017-09+ star-ZP chain validated against the L1 overlap
era). Plan §5.8 row 9 updated.

**SOHO-observer crossing list `crossings/soho_v1` (2026-08-25, same
day; `xng-298d55d0ce6b`).** `sglsurvey/crossings.py` gained
`--observer soho` / `--fetch-observer soho` on the TESS pattern:
Horizons −21 SSB vectors, 6 h, 1996-01-01 → 2026-10-01 (stop bounded
by the Horizons SPK end 2026-10-05, not the universal 2028 window —
extend at the yearly refresh; table 44,925 rows sha-pinned under
`crossings/observers/`, SOHO–Earth distance sanity-checked at
0.008–0.011 AU). 10,714 events, 0 invalid, manifest parity. Census +
Earth-center validation
(`surveys/heliospheric-crossings/scripts/soho_census.py` →
`results/soho_census_v1.json`): sunward grazing pairs shift |Δb|
median 0.124 / max 0.211 R☉ (TESS-scale, as the recon's measured
0.93 R☉ halo offset predicted at axis-aligned geometry); large-b
events shift up to 2.7 R☉ (the 0.01 AU L1 radial offset projected) —
so wide-rung units must also be defined on this list. Grazing-rung
membership unchanged (S1 ≤ 1.2 R☉: 117 ev / same 4 targets; ≤ 2.5 R☉:
153 / same 5); S1 0.1 AU 216 ev / 7 targets. Notable: the van-maanen
annual October S1 family deepens in the SOHO frame — 2009–2015 events
drop from b ≈ 0.27 to 0.13–0.23 R☉, strengthening the headline cell.
Next: hypothesis freeze (remaining pre-freeze items in the recon
note's addendum).

**Hypothesis freeze v1.0 DRAFTED (2026-08-25, same day;
`surveys/heliospheric-crossings/hypotheses.md`) — pending user
approval of decisions D1–D9.** Input pinned to `crossings/soho_v1`;
channels S1/S2 supersede the universal sunward out-of-scope
declaration for this survey only. Key frozen-at-draft structure:
level-0.5 uniform search substrate with level 1 as validation+masks
(D1); S1 filled-cone windows with the annulus double-passage as an
adjudication discriminator (D2); radial+temporal detrend layer (D3);
8 same-radius PA-ring controls primary / temporal pseudo-windows
secondary (D4); Uranus conjunction passage as positive control (D5);
full-cadence grazing + 3-hourly 0.1 AU subsampling ≈ 200 GB transient
(D6); S2 1.0 AU out of scope + S1 1.2 R☉ gated on the measured C2
inner radius (D7); dev = ross-154/gj-908/ross-128, confirmatory = the
four grazing-family targets (D8). S2 is an empty-field excess search
(all seven targets V ≥ 9 > the C3 single-frame limit — no saturation
cut exists). Pre-freeze data contact verified clean: the three recon
frames' date lies in zero sunward windows. Declared sensitivity:
opens the sunward cell at MW-class (0.3–1 MW recurrence-stacked
downlink through the grazing cones incl. van-maanen 0.13 R☉).

**Hypotheses FROZEN + coverage intersect DONE (2026-08-25, same day;
freeze v1.0 approved as recommended D1–D9; coverage
`surveys/heliospheric-crossings/results/coverage_v1.{json,md}`,
occulter gate `results/occulter_radii_v1.json`, 4,320 listing
snapshots under `runs/heliospheric-crossings/coverage/`).** Occulter
annuli measured from level-1 frames across eras: C2 usable
[2.2, 6.3] R☉ (bright diffraction ring at 1.9–2.0 R☉; occulted core
NOT zero-filled — the automated valid-pixel metric fails, brightness-
profile inspection rules), C3 usable [4.4, 29] R☉ (post-1997 core =
constant fill plateau, MAD 0). Coverage (703 event-rung rows,
day-granularity from SDAC/NRL dated listings): **S1 1.2 R☉
not_constrainable (0/117 visible — D7 confirmed, ledgered as
nominal-covered/occulter-unusable); S1 2.5 R☉ 148/153 covered+visible
in its [2.2, 2.5] wings (visfrac 0.12–0.25, ~28–60 C2 frames/window,
29–31 windows/target × 5 targets); 0.1 AU rungs the workhorse —
S1 208/216, S2 212/217, ~1,000+ C3 frames/window full-cadence,
visfrac 0.80 (grazing-depth targets, 4.4 R☉ core cut) to 1.00
(gj-908).** All 23 losses explained: the 1998 attitude-loss gap (one
full annual generation, family-wide), the 2026-07+ tail past NRL
currency (refresh-recoverable), 2 isolated day gaps. Next: threshold
freeze (searchable units 5 C2 + 14 C3; split per D8, seed declared
there).

**Threshold freeze v1.0 (2026-08-25, same day;
`surveys/heliospheric-crossings/thresholds.md` +
`configs/threshold_freeze_v1.json` sha256:919af49c…, seed 20260825).**
19 searchable units (5 × S1 2.5 R☉ C2 wing-regime + 7 × S1 0.1 AU +
7 × S2 0.1 AU C3), 43 trials (C2 units 3 statistics, C3 units 2 —
S_pulse full-cadence-only per D6), expected control crossings 4.8
(dev 11/1.2, confirmatory 32/3.6). No z-grid exists — both channels'
predicted positions are z-independent (axis geometry only), z enters
via the beam-radius rungs alone. Statistics = S_stack (recurrence
stack, primary) / S_event / S_pulse on detrended series (azimuthal-
median radial profile subtraction 0.1 R☉ bins → 30 d temporal median
→ variance rescale k); rule S > max(T,0); T = max over 8 same-radius
PA-ring controls (ΔPA ±25/50/75/100°), temporal pseudo-windows
secondary-reported. Gates: ≥ 5-star WCS fit + ZP scatter ≤ 0.2;
aperture ≥ 90 % valid (catches pylon/edge/telemetry blocks);
planet-proximity 10 px; ≥ 10 visible epochs/event, ≥ 8 events/stack,
≥ 200 baseline epochs. Split per D8: dev = ross-154 + gj-908 (S1+S2
0.1 AU) + ross-128 (S1 2.5 R☉); confirmatory = the four deep-graze
targets' 12 units + ross-128's two C3 units. Next: dev search.

**Dev stage complete — amendments v1.1 + v1.2 frozen, 0 candidates
(2026-08-26; `surveys/heliospheric-crossings/results/dev_search_v11.json`,
engineering log `notes/dev_machinery_log_2026-08-25.md`).** Dev data:
32,843 frames / 62 GB (keep-alive fetch fix after a 1-connection
serialization), full measure 0 worker errors, astrometry C2 80.3 %
(2-star L1 rule) / C3 99.99 %. The frozen v1.0 statistic taken
literally was numerically invalid (occulted-core fill-plateau epochs
→ zero formal error); the amended chain (v1.1: astrometry-gated epoch
validity, ring-differential statistic with night medians +
baseline-median centring + empirical night σ, measured colour system
c_C3 = 0.429 mag/(B−V) with per-frame scatter 0.53 → 0.38, S2
stellar template for V ≤ 12.5 targets, Tycho-2 VT ≤ 11 mask,
high-background class; v1.2: S2 source-star mask exemption,
bright-planet in-FOV veto) closed dev at **11 trials, 1 exceedance vs
1.2 expected**. Adjudications: ross-154 S1's z = 20–32 events lay on
the exact 8-year Venus synodic cycle (2000/2008/2016/2024, elongation
5.9–8.0°, ephemeris-verified) — Venus stray light, vetoed; ross-128
S1 S_event 1.45 vs T 1.20 = the 2015-03-17 St. Patrick's Day CME
storm window — adjudicated control-crossing, retained. Positive
control achieved in-situ: gj-908 (V 8.98) detected window-locked at
the predicted S2 position, modelled by the stellar template, statistic
null after subtraction. Next: blind confirmatory run (14 units) under
the frozen v1.0+v1.1+v1.2 chain — fetch, measure, reduce once.

## 10. PTF/iPTF reachability recon (2026-08-26; plan §5.8 item 3)

Live probes from the dev machine, all anonymous
(`surveys/ptf-crossings/notes/ptf_recon_2026-08-26.md`). The plan
row's `ptf/products/` path does not exist — the archive is IRSA IBE
`ptf/images/level1` (epochal CCD exposures; 175 metadata columns
incl. `pfilename`, MD5 `pchecksum`, full PV WCS + corners, `dmask` +
`sexcat` ancillaries) and `level2` (reference coadds + uncert +
depth-of-coverage + PSF per field/filter/chip). Cutout service
verified (SIP WCS preserved, shifted CRPIX); gid 100 and 101 both
fetch anonymously; TAP `ptf.ptf_procimg` / `ptf_objects` /
`ptf_lightcurves` live (per-object variability stats,
`transient_flag`, sources to R ≈ 21).

Two plan corrections: the public archive spans MJD 54891–57051
(**2009-03-01 → 2015-01-28**, 2,818,277 CCD exposures — late iPTF
never released), and there are **no public difference images** →
warp-direct substrate on the PS1-crossings pattern. Coverage probed
at all five grazing-family targets (both channels, 0.01° boxes):
patchy and campaign-driven — teegarden antipode has zero PTF
coverage; wolf-359/gj-1276 antipodes have epochs but none within
±3 d of a grazing t_ca. The survey case rests on **van-maanen**
(346 antipode epochs; the April b 0.27–0.37 R☉ family — deepest of
the era, destroyed by PS1 chip-gap masks — covered in 2009/2011/2013
with 9 epochs at closest 0.36 d in 2011) and **ross-128** (254
antipode epochs, 4 grazing events covered ±3 d), plus 2 on-star
events each at van-maanen and teegarden. The full in-era window scan
(24 grazing events per target-channel) closed at recon: 6/10/14
family events covered at ±1/3/5 d → ~5–7 searched target-channels.
Timing convention verified from a header (`obsmjd` = shutter-open
UTC; `t_mid = obsmjd + aexptime/2`). Adapter confirmed mostly-config
on `irsa_ztf.py`; no server dependency (cutout-scale pulls). Next:
hypothesis freeze (remaining open items: dmask bit definitions,
in-window LIMITMAG distribution, photometry chain).

**Freeze prep + hypotheses draft (2026-08-26, later).** Pre-freeze
homework closed: dmask bits pinned from Laher et al. 2014 Table 15
(fatal template 65533 — every bit except 2¹ object-detected;
ancillaries must be selected by `anciltype`, slot order varies);
off-window depth sampled at 3 fields (30 offset header stamps —
R LIMITMAG median 21.1, range 19.8–21.7, keyword present in only
~1/3 of frames → the survey computes its own 5σ depths);
`photcalflag` measured ≈ 0 at 6/8 scoped positions → field-star ZP
primary for all frames. Definitive flat-chord era scope
(`surveys/ptf-crossings/results/era_scope_v0.json`): **the 1.2 R☉
rung has zero in-window epochs** (closest 0.36 d vs ±0.3 d edge);
the van-maanen deep-graze family is covered at the 2.5 R☉ rung
(2011 b = 0.278 R☉, 2 epochs); 9 units with coverage, led by
van-maanen B 0.1 AU (4 ev / 25 ep) and ross-128 B (2 ev / 67 ep +
2 ev / 9 ep at 2.5 R☉); teegarden/gj-908/ross-154 antipodes have
zero PTF epochs. `surveys/ptf-crossings/hypotheses.md` v1.0 DRAFTED
(PS1-crossings construction, scie-direct substrate, decisions D1–D8
incl. dev = wolf-359 B + gj-1276 A + ross-128 A, confirmatory = both
headline B families + van-maanen/teegarden A) — **FROZEN 2026-08-26,
all D1–D8 approved as recommended**. Next: coverage stage (fresh
snapshot-disciplined pulls; adapter build on the irsa_ztf pattern).

**Coverage stage complete (2026-08-26, later).** Adapter
`sglsurvey/adapters/irsa_ptf.py` built and validated (explicit-column
IBE discovery — the default column set omits fid/WCS/checksums;
MD5 verification incl. the archive's truncated-achecksum quirk,
prefix-verified; dmask int16→uint16 reinterpretation; exact-footprint
smoke test usable at the survey position). Coverage run
`surveys/ptf-crossings/scripts/coverage_intersect.py` (93 snapshot-
disciplined discovery boxes, ~12 min): **covered events match the era
scope exactly** under the full frozen channel definition (the scoping
scan had omitted the side-of-axis cut — denominators halve, covered
sets unchanged). B: 1.2 R☉ 0/24 (structurally uncovered, as frozen);
2.5 R☉ 3/30 — **van-maanen 2011 b = 0.28 R☉ with a same-night g pair
in-window** + ross-128 ×2; 0.1 AU 7/42 — the full van-maanen April
recurrence family 2009/2010/2011/2013 incl. `evt-4c4ea2c397` (2010),
the exact event PS1's chip-gap masks destroyed, now with 2
independent PTF epochs; ross-128 2012-09 carries 65 epochs / 39
same-night pairs. A: 0.1 AU 4/42; 1.0 AU 76/513 constraint-only.
photcalflag≈0 throughout confirms field-star calibration as the only
chain. Results `surveys/ptf-crossings/results/coverage_v1*`;
snapshots `runs/ptf-crossings/`. Next: threshold freeze (searched-
unit population, saturation cut, statistic thresholds), then dev.

**Threshold freeze v1.0 (2026-08-26, later).** Saturation cut
(`saturation_cut_v1.*`, 86 era targets, frozen rule E_R 14.5 /
E_g 15.0, Mould-R via Jordi 2006): of the four covered A 0.1 AU
units **only gj-1276 R survives** — teegarden R 13.73, van-maanen g
12.40 and ross-128 R 9.86 all excluded-class, as the frozen rule
anticipated. Threshold freeze
(`surveys/ptf-crossings/thresholds.md`,
`configs/threshold_freeze_v1.json` sha256:26e235ee…ef53ae, bound to
hypotheses + coverage + saturation hashes): **7 searched units, 11
trials (dev 1 / confirmatory 10), expected control crossings 1.22**.
Units: van-maanen B 2.5 g (the b = 0.28 R☉ same-night pair, S_event)
+ B 0.1 g/R (2 and 3 events, S_event+S_stack); ross-128 B 2.5 R +
B 0.1 R (2 events each, S_event+S_stack) + B 0.1 g (45-epoch single
event); gj-1276 A R (dev). wolf-359 B = single-epoch class (dev);
substrate = scie-direct with PS1-DR2 field-star ZP gate (≥ 5
calibrators, scatter ≤ 0.2); controls = 8 ring trajectories (B) / 8
pseudo-windows (A); veto ladder incl. the same-night repeat test.
Both headline B families fully confirmatory-blind. Next: dev stage
(machinery on wolf-359 B + gj-1276 A + ross-128 A, bright-star
saturation gate, asteroid positive control), then blind confirmatory.

**Dev stage complete (2026-08-26, later;
`surveys/ptf-crossings/results/dev_search_v1.md`).** Machinery
validated end-to-end on real data (110/110 cutouts, zero 404s;
`build_flux_map_ptf` added to `sglsurvey/photometry.py`). Dev-fixed
realizations: dwarf-locus calibrator restriction (unrestricted PS1
set inflated ZP scatter to ~0.25; restricted 0.047 on the test
frame), 256-px cutouts everywhere (128 px held ~2 calibrators at
these latitudes), asteroid-local calibrator cones. Results: wolf-359
B single-epoch class clean (S 0.101 vs T 2.099, 8/8 ring controls,
k chain live); **gj-1276 A resolved constraint-only under the frozen
offset-validity gate** (1/12 valid offsets — campaign cadence) →
0 searched dev trials, confirmatory = 10 trials / 1.11 expected
control crossings, all channel B. Saturation gate: R boundary
bracketed 12.98(sat)–14.06(clean) vs frozen 14.0 — no amendment;
g one-sided (no searchable consequence); ross-128 exclusion
confirmed (bits 8+6 fire at core). Positive control: (8971)
Leucocephala position-locked in 3/3 exposures, internal RMS
0.025 mag; +0.24 absolute vs Horizons V + assumed V−R within
prediction systematics (Horizons response snapshotted); flux scale
pinned by the PS1-DR2 calibration. Next: blind confirmatory run
(6 B units / 10 trials) — fetch, measure, reduce once, with
completeness injections.

**Confirmatory run + completeness + report — SURVEY COMPLETE
(2026-08-26, later; `report/ptf_crossings.md`).** Blind confirmatory
(600/600 cutout pairs, zero 404s): **10 trials, 0 exceedances** vs
1.11 expected control crossings — **0 candidates**; every
coverage-stage event retained its epochs at measure, k 1.00–1.10.
One amendment (v1.1, frozen before re-reduction): confirmatory
cutouts 256 → 384 px — the 256-px support was inoperable at the
ross-128 antipode (142/144 frames < 5 in-frame calibrators, a
surface-density failure; calibrator rule and statistics untouched;
superseded pass kept, van-maanen statistics moved ≤ 0.13). One
mechanical fix mid-measure: blank in-header MAGZPT on
non-photometric frames crashed `build_flux_map_ptf` (guarded).
Completeness (Moffat-β3 stamp-response, 200 draws × 15 mags, no dead
z): m90 = 21.37 (van-maanen B 2.5 g — the b = 0.28 R☉ unit),
≥ 22.0 ×3 (censored), 21.61 / 20.58 elsewhere. Physical: the
photosphere-grazing cell's first pre-2015 constraint — ≈ 190 W
(van-maanen 2011) / ≲ 110 W (ross-128) through the 2.5 R☉ cone;
kW-class (8–30 kW) through the 0.1 AU cone across the 2009–2013
recurrences; recurrence stacks close clean on 4 units. Lessons
(report §4): calibrator support is a field property; dwarf-locus
restriction load-bearing; PTF = highest-integrity archive interface
so far (710/710 fetches, published MD5s); campaign cadence defeats
on-star temporal controls; side-of-axis matters at scoping.

## 11. GALEX/gPhoton reachability recon (2026-08-26; plan §5.8 item 4)

Queue item 4 started as `surveys/galex-crossings/`: the UV
pulse-period cell via MAST's gPhoton photon database (recon note
`surveys/galex-crossings/notes/galex_recon_2026-08-26.md`, era scope
`results/era_scope_v0.json`). All routes probed live and anonymous
from the dev machine: the Mashup SQL service
(`mastcomp.stsci.edu/.../GalexPhotonListQueryTest` — the `mast.`
host 404s; `fGetTimeRanges` from the old docs does not exist),
`fGetNearbyAspectEq` (per-second aspect with band + boresight
distance; `band` is `'FUV/NUV'` when both detectors are on —
substring-match), the `NUVPhotonsV`/`FUVPhotonsV` photon views
(ra/dec/time box queries; 5 ms-tick stamps verified), the GR6+7
MCAT, and MAST CAOM. Photon-DB era measured 2003-06-07 →
2013-05-01. Era scope: the grazing rungs (1.2/2.5 R☉, ±0.35 d
windows) are **structurally uncovered** at all 10 grazing-family
positions — coverage-without-statistic, incl. the van-maanen
deep-graze family (visits in 2004/2008, none in-window); the
**0.1 AU rung carries the survey**: 5 usable in-window units at 4
target-channels (gj-1276 B 2007 + 2010 — the latter a 1,637 s NUV
visit, photon pull verified: 1,695 photons / 30″ box ≈ 1.0 ct/s
background; wolf-359 A, gj-1276 A, ross-128 A 2007, each ~100 s
AIS-length with simultaneous FUV), plus one rim edge case
(gj-908 A at 36.3′ boresight vs the adopted 33′ cut). Earth-center
`universal_v1` valid (LEO, the standing 0.010 R☉ budget). Next:
hypothesis freeze — pre-freeze items are photon-flag semantics, the
boresight cut, the pulse-statistic calibration route (relative
Poisson counting vs gAperture/gPhoton2 absolute chain), and an
M-dwarf flare veto for the three on-star units.

**Hypothesis freeze v1.0 (2026-08-26, same day).** Pre-freeze probes
(off-window only): `aspect` and `imgrun` schemas pinned; photon flag
census (100 % flag 0 in 8,621 off-window photons at the gj-1276
antipode); aspect-correction sanity on an MCAT field star (20,243
photons, centroid 0.63″ from catalog, RMS 4.2″ ≈ PSF) — with a ×2.6
raw-rate vs naive-ZP discrepancy recorded as a dev-stage
live-time/aperture gate. `surveys/galex-crossings/hypotheses.md`
**FROZEN** (D1–D8 approved as recommended): photon-event substrate
(flag 0, boresight ≤ 33′, MCAT-star in-visit calibration at the
≤ 0.2 mag gate); three statistics per unit — S_rate (d = 1),
S_burst ({0.05, 0.5, 5, 50 s}), S_period (H-test, 20 ms → T_visit/3,
LEO orbital phase smear carried by the injections, not corrected);
FUV a search band on the 4 FUV-live units (≤ 27 trials, FWER
α = 0.05); controls = ≥ 8 same-visit pseudo-positions (z-family
mirrored) + off-window same-position visits; A-channel M-dwarf
flare veto (FRED morphology + two-band discriminator; NUV contains
266 nm = quadrupled Nd:YAG — the programme's first in-band harmonic
of the 1064 nm family). Pre-freeze contact declared in full,
remedy D1: the recon's in-window reachability pull at gj-1276 B
2010 demotes **S_rate only** on that unit to `forced_dev`
(burst/period stay blind); ross-128 A retained with declaration.
D8: dev = off-window pseudo-units only; **confirmatory = all 5
units, blind**. gj-908 A rim-excluded (36.3′ > 33′); A 1.0 AU
deferred programme-wide. Next: coverage stage under snapshot
discipline (fresh pulls, exact per-event per-epoch positions).

**Coverage stage + amendment v1.1 (2026-08-26, same day).** Adapter
`sglsurvey/adapters/mast_gphoton.py` (snapshot-disciplined Mashup SQL
client) + `surveys/galex-crossings/scripts/coverage_intersect.py`:
one era-wide aspect discovery query per (channel, target), then exact
per-event windows with per-second boresight distances to the
per-event predicted positions (B: apparent relay per z-grid point at
the visit epoch). 47 snapshots under `runs/galex-crossings/`.
Population after the side-of-axis filter: A 70 / B 69 events, 7
targets each. First pass (v1.0 gates) zeroed gj-1276 B 2010: all
1,637 in-window seconds carry aspect flag 64 — pervasive at that
field in 2009–2010, while the recon's flag-0 observation was a 2007
visit. gPhoton's own convention (`PhotonPipe.py` L534–536: usable
aspect ⇔ `flag % 2 == 0`; the `flagDiv2` column; the recon's 0.63″
astrometry check ran on flag-64 seconds) established bit 0 as the
only bad-aspect bit → **amendment v1.1** (user-approved):
aspect gate `flag % 2 == 0`; v1.0-gate outputs preserved as
`coverage_v1_flag0_superseded.*`; formal flag-64 astrometric
verification queued for dev. Result: **the 5 frozen units confirmed
exactly** (gj-1276 B 2010 1,637 NUV s dmin 23.7′; gj-1276 B 2007
109 s NUV+FUV at 32.4′; gj-1276 A 92 s; wolf-359 A 97 s; ross-128 A
110 s — totals 2,045 NUV / 408 FUV s); live seconds z-independent at
coverage (2–30″ spread ≪ 33′ gate); grazing rungs 0/39 and 0/49
covered; gj-908 A rim-limited (110 s at 36.3′); teegarden antipode
archive-empty. Record `results/coverage_v1.md`. Next: threshold
freeze (trials tally, control-ensemble sizes, dev pseudo-unit
thresholds), then the blind confirmatory run.

**Threshold freeze v1.0 (2026-08-26, same day).**
`scripts/freeze_thresholds.py` → `configs/threshold_freeze_v1.json`
(`sha256:18346431…017e47`), bound to hypotheses v1.0+v1.1 and the
coverage products; no photon touched. The frozen section-6
nonlinearity cut ran as its declared MCAT input (snapshotted):
wolf-359 NUV 19.28 ok, ross-128 21.18 ok, gj-1276 no MCAT source
within 15″ → ok by construction (per-visit high-PM identification a
dev item) — no unit lost. Family: **8 searched units (4 NUV + 4
FUV) × 3 statistics = 24 trials**; unit statistic = max over
eligible events (S_rate on gj-1276 B 2010 stays forced_dev per
D1a); statistics made numeric (S_rate cts/s in 8″ apertures with
the B z-family max; S_burst boxcar z over {0.05,0.5,5,50 s};
S_period H-test m ≤ 20 on the 5×-oversampled geometric grid
0.02 s → t_span/3, ≥ 10-photon gate); controls 8 per trial —
B pseudo-positions at the locus boresight radius (40° spacing, +5°
MCAT/locus-exclusion rotation rule), A same-duration off-window
segments (seed 20260826) — threshold T = max over controls,
exceedance S > max(T, 0), **expected control crossings 2.67**.
Veto ladder order frozen (flare veto with two-band discriminator
first; SkyBoT at burst exceedances; detector-fixed clustering;
recurrence). Next: dev stage on off-window pseudo-units (flag-64
astrometry, live-time/×2.6 closure + ≤ 0.2 mag MCAT gate, control
census, flare census, statistic machinery + injections), then the
blind confirmatory run.

**Dev stage — CLOSED (2026-08-26, same day; `results/dev_v1.md`).**
Pseudo-units only, off-window guarded (`scripts/dev_stage.py`,
`galexlib.py`). (i) Flag-64 astrometry PASS — 16 source-blocks,
median centroid 0.37″, RMS ≈ PSF: the v1.1 premise formally
verified. (ii) Calibration PASS — ZP_eff 19.592 ± 0.088 over 15
deduped MCAT calibrators (the catalog carries per-visit duplicate
rows; dedup at 3″); the recon ×2.6 fully explained (duplicates +
crude live-time + 0.49 mag aperture/dead-time term); blank sky
0.209 cts/s per 8″ aperture. (iii) B pseudo-position rule valid
(0 rotations × 3 blocks); A segment gate assessed → gj-1276 A
NUV/FUV (4/3 segments; its long visit rim-only) and wolf-359 A FUV
(0; long visits NUV-only) constraint-only — confirmatory family
**5 band-units / 15 trials / 1.67 expected control crossings**.
(iv) Flare census: wolf-359 flared in all 3 long visits (S_burst
105.7/62.0/14.4; FRED templates rise ≤ 10–20 s, decay 30–40 s,
contrast to 21.6×; one slow-rise non-FRED event → the two-band
prong is load-bearing); gj-1276/ross-128 quiet. (v) B pseudo-unit
(1,381 s flag-64 segment) null-clean on all three statistics; the
A pseudo-unit drew the 2009-03-25 flare segment and **all three
statistics fired through the full chain — the conditional positive
control satisfied in-situ**; injections: persistent m90 ≈ NUV 22.3
(1.4 ks class), 0.5 s pulse fluence ~5 photons, trains to P = 10 s.
Dev's structural find: the 5 ms tick (200 Hz) corrupts the H-test
grid at commensurate harmonics (null H 2,700–3,800 vs ~32) →
**amendment v1.2** (user-approved, validated first): S_period on
U(0,5 ms)-jittered times, seed 20260826 — null restored to H 29–38,
P = 10 s recovery 6/6 at the previously-dead amplitude
(`results/jitter_validation_v0.json`). Cleared for the blind
confirmatory run under v1.0+v1.1+v1.2.

**Blind confirmatory + completeness + report — SURVEY COMPLETE
(2026-08-26, later; `report/galex_crossings.md`).** The first
sanctioned in-window contact (`scripts/confirmatory_search.py`;
statistics locked to `results/confirmatory_v1.json` before any
injection): 5 band-units / 15 frozen trials, 14 effective (ross-128
A FUV S_period degraded at the frozen 10-photon control gate; its
S_rate was an exact tie with the control max — strict inequality
holds). **1 exceedance vs 1.56 expected control crossings — 0
candidates.** The exceedance (wolf-359 A NUV S_burst 4.68 vs 4.13)
adjudicated retained-ambiguous: a 5-photon 0.5 s cluster (0.80
expected) — FRED test unresolved at that fluence, two-band
discriminator powerless (0.067 FUV expected), strict ≤ 25′ re-run
unchanged, asteroid excluded by timescale, no in-archive recurrence
window; leading interpretation an ordinary CN Leo micro-flare
(within the star's measured off-window burst activity). D1a lane
reported: 0.235 cts/s at the 2010 locus ≈ background. Completeness
(photon-level injections vs locked thresholds; supplement for the
1,637 s 2010 event + wolf-359's variability-inflated grid):
persistent m90 NUV 21.3–21.9 (wolf-359 18.6, variability-limited),
FUV ≈ 21.2–21.8 (nominal-ZP budget) → **20–33 kW through the 0.1 AU
cone** (415 kW at wolf-359); pulses 90 % at 4–5 photons →
≳ 1.7 MW (0.5 s) / 17 MW (0.05 s) peak; **coherent trains
P = 0.5–50 s at 60 photons / 0.037 cts/s ≈ 6 kW time-averaged on
the 2010 window — the programme's first sub-second-cadence
crossing constraint** (120 photons at P = 0.05 s, drift smear
carried). Lessons (report §4): read flag conventions from pipeline
source; jitter discretized time stamps; dedupe the MCAT; control
supply is per-band; on-star flare-star units want a design-time
two-band fluence gate. The wolf-359 burst is recorded for any
future UV mission able to re-cover a window; ledger rows queued for
the §5.7 refresh.

## 12. DASCH recon + crossings backward extension (2026-08-26; plan §5.8 item 5)

Session goal: open queue item 5 — the pre-1980 century of crossing
windows in the Harvard plates (DR7) — by executing the plan row's
prerequisite (extend the universal crossing list backward from 1980)
and a DECam-style reachability recon of the DASCH services (marked
_unprobed_ in the plan). Both done same-day on the dev machine; full
facts in `surveys/dasch-crossings/notes/dasch_recon_2026-08-26.md`.

**Backward extension.** `crossings/universal_1885_v1`
(`xng-e2d1063af9d0`): Earth center, 1885-01-01 → 1993-01-01, same
registry/model/every-minimum construction — 37,096 events, 0
invalid, 91 MB, ~85 min. Deliberately overlaps `universal_v1` over
1980–92 as a regression check: **all 4,104 shared-era events match
1:1, max |Δt_ca| 23 s, max |Δb| 3.4×10⁻⁷ R☉**. 13,053 events
`degraded`, dominated by `long_propagation_span` (sglseti's declared
75-yr linear-motion bound → everything before ~1941). The §7
old-epoch accuracy budget was then *measured*: astropy builtin vs
JPL DE440S heliocentric Earth ≤ 6 km ≈ 1×10⁻⁵ R☉ at 1885;
grazing-family astrometric sensitivity (5 km/s RV + 0.1 mas/yr µ
perturbations through the rigorous `apply_space_motion` propagation)
≤ 2×10⁻⁴ R☉ and ≲ 5 s in t_ca — 50× under the standing 0.010 R☉
LEO budget, so the degraded flags are conservative bookkeeping, not
physical limits. Unmodelled century-baseline orbital motion stays a
per-target freeze check.

**Recon.** DASCH DR7 rides the Starglass REST API
(`api.starglass.cfa.harvard.edu/public/`, anonymous, JSON POST;
WAF rejects python-urllib UAs — send a curl-like UA; registered
`x-api-key` tier exists if rate limits ever bite). All five science
endpoints probed live end-to-end: queryexps (9,000–15,000 exposures
per probed position, 27 columns incl. per-exposure `limMagApass`,
`exptime` **in minutes**, `wcssource`), querycat, lightcurve
(van Maanen: 3,544 rows 1890–1989, calibrated mags + per-epoch
`limiting_mag_local` on every non-detection row, `time_accuracy_days`
mode 60 s but absent on non-detection rows, AFLAGS/BFLAGS bits
unpinned), cutout (FITS verified: 20′, 1.44″/px, TAN, target
centred; minimal header — CD-diagonal-only, resampling status
unpinned), platephot (subregion detections **including uncatalogued
sources** — catalogue-level channel-B search possible; exactly-50-row
probe response unexplained), mosaic_package (pre-signed S3, binning
1/16). DR7 facts: ~430k plates, 97 %/89 % astrometric/photometric
success, APASS DR8 = B-band science refcat ("excellent long-term
stability") vs ATLAS-refcat2 with documented false long-term trends
(astrometry only); typical depth B 12–16, probed best 18.3; Menzel
gap 1953–68 confirmed at every position; known-issues page (defects,
blends, splitting, missing points) drives the single-detection
vetting design.

**Science yield of the recon.** Coverage probed at the grazing
family's 10 positions (stars + antipodes): ~2.3–3.1k calibrated
exposures each, 150–340 at B ≥ 15, both hemispheres. Intersected
with `universal_1885_v1`: **van-maanen, gj-1276 and wolf-359 are
photosphere-grazing (≤ 1.2 R☉) on essentially every annual crossing
of the DASCH century** (van-maanen b ≈ 0.59–0.69 R☉ throughout —
the modern deep-graze family is the tail of a secular deepening);
teegarden B holds 38 grazing events and the 2.5 R☉ rung adds
ross-128 everywhere. ~15–28 covered windows per target-channel at
±1 d (~120 grazing windows across the family; 0.1 AU rung tracks
±3 d at ~35–48/channel) — roughly an order of magnitude more
grazing-rung windows than all prior surveys combined, ending where
the 1980 list begins. Next: hypothesis freeze (open items listed in
the recon note: flag-bit pinning, platephot semantics, cutout
resampling, series policy, flat-chord scan, timing gate, positive
control, snapshot set).

## 13. Rubin DP2 crossings survey (2026-08-26; plan §5.13 — complete)

Same-day chain following the morning's RSP recon (§ `surveys/rubin/`):
plan §5.13 opened (Pipeline A tabled until the image release, §6);
hypothesis freeze v1.0 (D1–D8; pre-freeze recon contact at the
ross-128 antipode declared, unit kept blind; reliability barred from
absolute cuts; ross-154 saturation expectation frozen as a rule, not
a decision) → coverage (21 snapshots; era re-measured identical;
ross-128 B 0.1 AU the only searchable unit — all five z positions on
detector 102, magLim r 23.366, Δt −5.66 d, b(t) = 20.9 R☉ rim
sample; ross-154 A 16 on-detector visits; grazing rungs 0 visits) →
threshold freeze (1 unit × S_det = 1 trial; Ross 154 G = 9.13 →
excluded per §6) → dev on 62 off-window pseudo-units (blind guard on
the unit visit; null 0/53 associations; **amendment v1.1**: control
locus-avoidance 10″ → 2.5″ after the 10″ rule proved geometrically
impossible at inner-z offsets; SkyBoT positive control PASS —
2006 SE393 at 0.308″, S_det 31.3, diaSourceId = the catalog's own
ssObjectId link; injection census 0 → completeness route ii) →
**blind confirmatory: S_det = 0, zero DiaSources in the discovery
cone, 0 exceedances vs 0.11 expected — 0 candidates.** Report
`report/rubin_crossings.md`. New record classes: first catalog-level
substrate; deepest wide-rung single-epoch flux threshold
(not-injection-calibrated statement, ~1.3 kW through the 0.1 AU
cone); designated follow-up = image-level re-run at the late-2026
visit/difference-image release upgrading the unit to an
injection-calibrated exclusion.

### §12 continuation — freeze through blind confirmatory (same day)

**Freeze chain.** Hypotheses v1.0 drafted post-recon and **frozen on
user approval of D1–D8** (5-target scope; catalogue-level substrate;
532 nm declared unconstrained — out of the blue plates' band — with
355 nm conditionally in-band; no pulse statistic, the century
recurrence stack as the new cell; van-maanen A `forced_dev` after the
recon lightcurve pull). Coverage stage (fresh snapshot pulls,
interval-based overlap, timing gate): **18 searched units** —
including all three B 1.2 R☉ grazing units (van-maanen 12, gj-1276 7,
wolf-359 4 covered windows) and 0.1 AU rungs at 56–63 windows each;
ledger rows teegarden B 1.2 R☉ and ross-128 B 1.2 R☉. Threshold
freeze v1.0: B locus-track hit statistic (r_match 10″, z-family
polyline, 8 same-pull ring controls at ±60–120″), A robust-z with
8 temporal pseudo-windows; 36 trials (dev 6 / confirmatory 30).

**Dev stage (0 exceedances / 6 trials, final chain).** Positive
controls: RY Cnc through the A-chain (85 faint excursions 75 %
phase-locked, depth 0.79 mag, scatter 0.18) and **(7) Iris 1911**
through the B-chain (found via a Horizons coincidence scan of the
teegarden star field after JPL sb_ident 500'd on pre-1950 epochs;
recovered 0.2″ from the timing-extended trail, correctly
trail-elongated). Three dev findings → amendments, all
pre-confirmatory: **v1.1** the DR7 APASS refcat carries pm = 0 dummy
entries for van-maanen and wolf-359 (the APASS "lightcurve" of van
Maanen = 9 spurious rows vs 1,552 real ATLAS ones) → PM-match
routing, era-local ±5 yr baselines, Stouffer S_stack; **v1.2** the
Iris control was killed by the frozen fatal template —
`SUSPECTED_DEFECT` stamps real single-plate transients — demoted to
annotation feeding cutout inspection; **v1.3** the teegarden dev
exceedance adjudicated to (a) the quiescent star marginally
extracted at the plate limit (off-window same-locus measurement:
2/12 deep plates, median B 17.11 ≈ expected 17.3) and (b) a
defect-class row killed by cutout pixels (≤ 3.4σ diffuse vs a 21σ
field star) → the limits-only quiescent-star gate (hit ≥ 1 mag
brighter than measured quiescence; rings never gated).

**Blind confirmatory (15 units / 30 trials): 0 candidates — 2
exceedances vs 3.33 expected, both adjudicated.** Both sat in
wolf-359 A 0.1 AU, whose quiescent measurement worked perfectly
(12/12 deep off-window plates detect the PM-less-refcat star,
median B 15.66): the 1943/ac37753 hit (B 14.58, 0.35 mag above the
plate limit, 7.4″, fwhm 10 px) and the 1982/dny00425 hit (B 13.02,
0.5 mag above limit, 7.1″) both show ≤ 4σ diffuse pixel structure
where real sources on the same cutouts reach 16–20σ — vetoed
defect-class (`results/adjudication_wolf359A_v1.json`); no
recurrence across the unit's 56 covered windows; genuine star
extractions sit at 1.4–2.4″. Annotations: one van-maanen B 0.1 AU
locus hit (1905, S = T = 1, not an exceedance); ross-128 A clean
with heavy-tailed controls (T 7.09). Completeness (field-star
recovery per covered window, Wilson intervals) and report
`report/dasch_crossings.md` close the survey — the programme's
first pre-1980 constraints, kW-class through the grazing cones back
to the 1890s, the century recurrence-stack cell clean on every
searched unit.

**Completeness + constraints (closing the survey).** Field-star
recovery pooled per plate-limit stratum (`completeness_pooled_v1`):
90 % depth sits **1.5–2.5 mag above the archive's limMag columns**
(the C1 lesson at catalogue level); deep-plate strata reach m90
13.5–14.5. Constraint headlines: van-maanen and wolf-359 B 1.2 R☉
≳ 20 kW on deep-plate windows, sub-MW (m90 10.5) over the
shallow-plate century bulk to the 1890s, MW-class 0.1 AU cones;
template costs measured on 8,714 matched rows (fatal 1.2 %, strict
22 %, SUSPECTED_DEFECT on 13.9 % of genuine stars). One flagged
anomaly: the ross-128 antipode field never reaches 90 % recovery
(60–85 % over B 9–14) — constraint rows honest but shallow;
diagnosis queued. Report `report/dasch_crossings.md`.

**Blind confirmatory run — 0 CANDIDATES (2026-08-26;
`surveys/heliospheric-crossings/results/confirmatory_v1.{md,json → confirmatory_search_v1.json}`).**
Fetch 68,330/68,522 frames (92.7 GB new; 21,477 shared from the dev
cache), full measure (0 worker errors), ONE reduce under the frozen
v1.0+v1.1+v1.2 chain: **30 searched trials, 2 exceedances vs 3.3
expected — under budget; every recurrence statistic null.**
Adjudications (frozen ladder): gj-1276 S_pulse 18.7/17.4 = one 12-min
frame (2014-09-04 07:12, +2078 with negative neighbours) → vetoed,
single-frame/cosmic-ray rule; van-maanen S_event 4.51/1.19 = one
window (2020-10-05/06), a ~3 h all-PA annulus disturbance whose onset
follows a CDAW-catalogued C2 CME by 48 min (20:48) — both ±25° rings
swing hundreds in both signs; catalogued CPAs (263°/274°/85°) differ
from the source PA (122°): the azimuthal-median detrend couples all
PAs during a corona-wide transient → adjudicated CME-period
systematic, retained, non-promotable (no recurrence: 27 siblings
≤ 0.17). New structural finding: **ross-128's S1 antipode is
permanently blended** — the fixed S1 sky position sits 121″ (2.2 px)
from a VT 7.2 star, inside the frozen 3 px mask at every epoch
forever → constraint-only, ledger nominal-covered/resolution-blended
(the S1 apparent source is a fixed ICRS point; the Sun sweeps past
it). Bookkeeping blemish recorded (baseline day-picks inside windows,
~5 % conservative contamination; no re-run — blind preserved).
Remaining: completeness (stamp injections + Uranus control + L1-era
ZP validation) → report/lasco_crossings.md.

**Completeness + report — LASCO SURVEY COMPLETE, 0 candidates
(2026-08-27; `report/lasco_crossings.md`; completeness
`surveys/heliospheric-crossings/results/{completeness_v1,power_limits_v1,uranus_control_v1,zp_validation_v1}.json`).**
Stamp-measured response + injections into the real confirmatory null
series vs the fixed confirmatory thresholds (100 draws, seed
20260825): **the coronal night-to-night systematics set the floor —
recurrence-stack m90 V 4.5–8.4**, far shallower than naive stack
scaling (the freeze's "measured, not assumed" clause doing its job).
Physics: downlink ≳ 40–530 MW through the 2.5 R☉ cone
(recurrence-stacked, C2 wings), 3–26 GW through the 0.1 AU cone;
**pulse cell 6.7 MW per ≥ 12-min pulse (gj-1276)**; S2 10-m-class
uplink 0.7–1.1 GW. Controls: 217 single-frame S/N>8 star measurements
unbiased (−0.02 ± 0.39); era stability 0.064 mag year-to-year,
−0.066 across the L1 boundary; Uranus recovered at the predicted
moving position at ~20σ stacked but 1.65 mag faint — its
methane-absorbed red spectrum in the red-weighted Clear band
(diagnosed via same-field stars: chain unbiased; a ±0.3 band-
conversion systematic declared). The flux gate is carried by the
stellar validation + the gj-908 red-dwarf in-situ control. The
sunward cell — archivally virgin before this programme — is closed
at MW-class power; substrate chain validated for STEREO/WISPR
follow-ons. v2 design notes recorded (baseline in-window picks,
CME-robust per-sector detrend, inner-field scale).

## 14. ATLAS + ASAS-SN crossings — coverage stage resumed on the server (2026-09-03; plan §5.9)

**Migration + checklist.** Execution moved to the new Linux server
(`runs/` → 8 TB SSD symlink). Resume checklist cleared the same day:
sglseti 19c8167 clean, `.env` token, ATLAS queue 200, ASAS-SN v2 port
9006 egress open, Python 3.12.2 env. API docs + OpenAPI schema
snapshotted (`runs/atlas-asassn-crossings/docs_snapshot/`).

**Pre-data amendments (hypotheses §12, A1–A3).** Building the task
list exposed a freeze counting defect: the "22 semiannual windows per
target" counted both yearly Earth crossings of the axis, but one of
each pair puts the channel's sky position at solar elongation
0.1–0.5° — the sunward combination the freeze itself excludes
(LASCO owns it). Verified by elongation at every task position
(kept 174–180°). Searchable population = **11 annual windows per
target per channel** (B 44/55/77 by rung, A 0.1 AU 77); stack depth
~0.4 mag shallower than pre-declared; pseudo-window family unaffected
(offsets at elongation ≥ 80°); no archive data at survey positions
had been searched. Eras: ATLAS end MJD 61286; ASAS-SN v2 ceiling
re-measured over 14 positions at **2025-06-16, unchanged after 15
months** (reprocessing lag — v2 is a 2013→2025-06 substrate). Schema
has no radec-list parameter (recon open item closed: one position per
task, serial).

**Coverage drain running.** `scripts/build_task_list.py` → 308
windowed (±110 d) difference-flux tasks: 231 B (t_ca position + D3
mini-track ingress/egress at z = 550, offsets 37.5″) + 77 A (mini-
track collapses by geometric identity — recorded). One windowed task
= 2.2 min server time (vs 25–70 min full-history) → ~12 h serial.
`scripts/atlas_drain.py` submits with `comment` = task key, keeps
≤ 8 in flight, fetch-then-DELETE with the atlas_api sidecar snapshot;
sidecars are the done-list, `inflight.json` + server-side comment
adoption make it crash-resumable. First series (van-maanen B
2016-04 event): 214 rows, 82 % FAQ-mask keep, 5 rows/night median,
mag5sig 19.06 — but 0 exposures inside its 0.34 d grazing window
(10 in the 0.1 AU window): the coverage-fraction attrition the stage
exists to measure. Next: coverage ledger (ATLAS in-window epochs per
event × rung × band + ASAS-SN v2 epoch lists), then the threshold
freeze.

**Coverage stage COMPLETE (2026-09-04).** The 308-task drain
finished at 06:06 UTC (7.4 h; 0 failures; median server time 54 s,
max 499 s; 16 MB of series + sidecars under
`runs/atlas-asassn-crossings/lc/`). ATLAS ledger
(`results/coverage_v1_{events.ecsv,summary.json}`, FAQ mask keep
median 0.93): **the grazing rungs are attrition-dominated** — B
1.2 R☉ 8/44 windows covered (van-maanen 4/11, gj-1276 2, teegarden 1,
wolf-359 1; 4 exposures median per covered window), B 2.5 R☉ 25/55
(van-maanen 7, teegarden 7, gj-1276 4, wolf-359 4, ross-128 3), while
the 0.1 AU rungs are near-complete (B 74/77 with 17 exposures median,
mini-track union identical; A 68/77 with 14). Cause: 0.5–1.3 d
windows against nightly-at-best sampling with weather and sun-gap
losses (ATLAS is a ~1–2 d cadence survey; a window shorter than the
revisit interval is covered only by luck). The §9 pre-declared
recurrence-stack depth (~100–180 exposures) is therefore not
delivered for the grazing rungs: the van-maanen 1.2 R☉ stack holds
18 exposures across 4 windows (≈ 20.6 stacked), 2.5 R☉ 30 across 7;
the 0.1 AU stacks (150–300 exposures) do reach ~21.5. The
photosphere-grazing 0.24 R☉ family is nonetheless **searchable for
the first time** (4 covered windows). Pseudo-window availability
scales the same way (mean 1.4 / 2.7 / 7.0 of 8 offsets with data at
1.2 R☉ / 2.5 R☉ / 0.1 AU, o band) — the 1/9 control budget will be
computed from actual availability at the freeze.
ASAS-SN v2 ledger (`results/asassn_ledger_v1_*`, field-level union
of quality-G image_ids of catalogued sources within 3′; coverage
records only): B 0.1 AU 62/70 in-era windows, 1.2 R☉ 12/40, 2.5 R☉
23/50 (median 1 image), A 0.1 AU 59/69 — the same attrition shape
at ~2 d cadence. Channel-A nearest v2 sources sit 4–8″ from
gj-908/ross-128/ross-154 but 15–36″ from van-maanen/teegarden/
wolf-359 — the high-PM stars' master_list positions are epoch-offset
(the DASCH PM-match lesson); the D2 saturation table must PM-match,
not nearest-match. Next: threshold freeze (seed, D8 split, control
budget from measured availability), then dev → confirmatory.

**Threshold freeze → dev → blind confirmatory (2026-09-04, same
day).** User-approved freeze recommendations (seed 20260904; dev =
wolf-359 + gj-908; per-statistic gates — S_pulse/S_event at ≥ 1
covered window, S_stack at ≥ 3, and ≥ 4 valid mirror-gated controls
for a trial; o searched, c annotation; ASAS-SN ledger-only; D6 MPC
control (15000) CCD). Saturation table on 20-d reduced-mode series
(`scripts/saturation_cut.py`): gj-1276 ok, teegarden marginal
(o 12.59), the other five excluded — channel A searches two stars.
**Freeze** (`configs/threshold_freeze_v1.json`, hypotheses §13): 14
searched units / 37 trials (dev 10, conf 27), 4.55 expected control
crossings; **S_stack is constraint-only at every grazing rung** (≤ 3
valid pseudo-stacks) — the recurrence cell's calibrated test is not
deliverable at nightly cadence; grazing S_pulse/S_event searched for
van-maanen (1.2 + 2.5 R☉), teegarden and wolf-359 (2.5 R☉).
Engine `scripts/search_core.py` (D5 running-median detrend, k
rescale, response-weighted chord, best-R mini-track position per
epoch, nested-z max, mirror-gated pseudo-window controls, max-rule
T); `run_search.py --split`.
**Dev (10 trials, 1.3 expected): 2 exceedances**, one epoch — wolf-359
2.5 R☉ S_event 3.55 / S_pulse 3.18: a single 3.5σ Sutherland exposure
in a degraded-sky quad whose other members read zero (rung 2 fail;
SkyBoT clean; no recurrence) → retained, non-promotable. Control
finding: a 445 µJy chi/N-45 single-frame artefact in a pseudo-window
sets the wolf-359 0.1 AU S_pulse T at 15.8 (FAQ mask has no chi/N
term; frozen max rule keeps the trial valid but insensitive). KS of
standardized baselines p 0.24–0.97; k 1.0–1.5; engine validity
counts match the freeze. No amendments (`results/dev_adjudication_v1.json`).
**Blind confirmatory (27 trials, 3.2 expected): 5 exceedances on 3
events, 0 candidates** (`results/confirmatory_adjudication_v1.json`):
X1 van-maanen 1.2 + 2.5 R☉ S_event 2.78 (evt-fc686a, 2017-04) — a
2.8σ three-exposure quad-consistent chord (+40 µJy) in one of four
covered photosphere-grazing windows, SkyBoT clean, no recurrence →
**retained-ambiguous, non-promotable**; X2 gj-1276 B 0.1 AU S_event
5.09 / S_pulse 3.04 (evt-9c6a08, 2025-03, z = 10000) — in-window
nightly means consistent with zero, the excess manufactured by a
−10 µJy off-window baseline depression the 30-d median cannot
remove; S_pulse quad-inconsistent → adjudicated systematic; X3
teegarden A S_stack 11.8 — the difference flux at the star is the
star's full 34 mJy (5″/yr PM dipole against the template), k 10²–10⁷,
controls −14…+7 → adjudicated channel-A blended systematic. Poisson
P(≥ 5 | 3.2) ≈ 0.23. v2 notes: chi/N mask term; window-mean-vs-zero
test in the ladder; PM/parallax systematics template for channel A.
Completeness (D6 injections, `scripts/completeness_v1.py`) and the
MPC positive control (`scripts/mpc_control.py`) running.

**Completeness + report — ATLAS/ASAS-SN SURVEY COMPLETE, 0 candidates
(2026-09-04; `report/atlas_asassn_crossings.md`).** D6 injections
(`scripts/completeness_v1.py`, 100 draws/mag, seed 20260904, z = 550
geometry, into the real baseline series against the frozen T): B
units recover to o 19.7–20.75 per window, 19.8–20.5 stacked at
0.1 AU; five trials artefact-set (T 15.8–412 from single-frame
pseudo-window outliers → insensitive); both A units systematics-
limited. Positive control (15000) CCD (MPC mode, 2024–2026, V
18.4–20.4): scale gate PASS on bias (apparitions −0.05/−0.10 mag after
one colour term, slope −0.03), scatter photon+rotation-dominated at
S/N 5–8; ±0.1 mag declared. Physics: **van-maanen 0.24 R☉
photosphere-grazing cone ≳ 78 W per window (4 of 11 windows) — first
constraint ever on that family**; 2.5 R☉ cones 280–590 W (4–5 windows
each, three targets); 0.1 AU 26–55 kW / 26–41 kW stacked, seven
targets; pulse cell 140–890 W (grazing). Structural finding: a window
shorter than the revisit interval is covered by luck — 18/45/95 %
coverage by rung — so the recurrence cell is not calibratable at the
grazing rungs on a nightly survey. v2 notes: chi/N mask term; ≥ 2
epochs for S_event; window-mean-vs-zero ladder rung; PM/parallax
template (or reduced-mode PM-propagated photometry) for channel A;
PM-matched ASAS-SN sources. Nothing committed to git during the
session (user's call).

## 15. Gattini-IR + WINTER reachability recon (2026-09-04; plan §5.8 item 6 — BLOCKED)

Live probes from the dev machine, all anonymous
(`surveys/gattini-crossings/notes/gattini_winter_recon_2026-09-04.md`;
scan `scripts/recon_scan.py` → `results/recon_scan_v0.json`;
snapshots `runs/gattini-crossings/recon/`). **Neither archive can
carry a crossings search today; no freeze.**

**PGIR.** The public DR1 (NOIRLab Data Lab `pgir_dr1`: `exposures`
4.13 M quadrant stacks, `photometry` 47.6 G rows, `sources` 148.9 M;
TAP + `queryClient` q3c cones, sub-second) is a **light-curve catalog
of 2MASS point sources**, J < 15.5, era 2018-10-15 → 2022-10-15
(1,130 nights), δ > −28.5°, per-epoch 5σ depth J ≈ 14.5 Vega at the
targets (plan said "J ~16" — that was the nominal AB number),
saturation J 8.5 → 6.0. No image service (SIA 404 where `nsc_dr2`
answers; no storage tree; no IRSA holding). Two defects kill both
channels: (1) photometry is forced at **2MASS-epoch positions** —
every in-era target is a high-PM star 15–120″ from its 2MASS entry,
which reads empty sky (van-maanen 16.6 vs J 11.7, wolf-359 17.1 vs
7.1, teegarden 16.6 vs 8.4) while the flux lands in drifting
neighbour entries; the bright ones (gj-908, ross-154, ross-128
pre-2020) are saturated besides; no antipode has a catalogued source
within 33″ → channel B has no photometry at all; (2) `obsjd` is
**float32 at the source** (paper table types) → 0.25 d bins (2,596
distinct values over 1,130 nights), ~40 % of a grazing window.
Non-detection rows are kept, so the catalog does give the **visit
list**: all 14 in-era target-channels (7 targets × A/B, ≤ 0.1 AU,
4 events each 2019–2022) imaged 80–224 times; the van-maanen
0.25 R☉ April family has zero in-window epochs (closest 1.06 d);
wolf-359 B (b 0.72 R☉) is in-window at 2.5 R☉ in 2020 and 2021 (2
epochs each) and marginally at 1.2 R☉; 0.1 AU rungs 3–26 epochs per
unit. Wavelength correction: PGIR is J only — **1064 nm is not in J**;
that cell needs WINTER Y.

**WINTER.** Y/J/Hs, 1.1″ px, single-visit J 17.6–19.3, operating
since 2023-06 — the right substrate — but closed: no data release, no
proprietary-period statement in any paper (incl. the Dec-2025
instrument paper), no IRSA holding, `winter.caltech.edu` is a
login-only portal.

**Outcome.** Item 6 marked BLOCKED in plan §5.8 (substrate, not
adapter cost); the route is a data ask to the PGIR team (stack
cutouts, or a PM-propagated forced-photometry run over the 56
windows — the visit list is the attachment) or a future image/DR2
release; WINTER re-check added to §5.7 standing maintenance. No
adapter written, no hypotheses drafted, nothing committed to git.
Cone-power framing if a route opens: ~5 mag shallower than ATLAS →
MW-class floors at 0.1 AU (opening the cell, as the plan framed it).

## 16. Radio metadata extensions — ASKAP + LoTSS (2026-09-04; plan §5.8 item 7 — complete)

Geometry-only repeat of the 2026-08-24 VLASS construction (§6 step 7)
on two revisiting archives; decision unchanged, no pixel touched.
Report `report/radio_crossings_ext.md`; script
`surveys/radio-crossings/scripts/coverage_intersect_v2.py` (+
`census_v2.py`); snapshots `runs/radio-crossings/v2/`.

**Recon (same session, all anonymous).** CASDA TAP
(`casda.csiro.au/casda_vo_tools/tap`) `ivoa.obscore` holds every ASKAP
data product with `s_region` polygons and per-SBID `t_min/t_max`;
collections RACS (168 k cubes, 2019-04 →), VAST pilot (45.8 k cubes,
**no `t_min` on the cubes** — dated via `casda.observation` by SBID),
VAST full (264 k, 2022-11 →), EMU, FLASH, WALLABY, DINGO, commissioning.
`dataproduct_type='image'` is a red herring (WALLABY/RACS DR1
mosaics only); the per-epoch products are `cube` +
`cont.restored.t0`. Several cubes per SBID (v1/v2 re-runs, lowres/raw)
→ dedupe per (SBID, field centre), best `quality_level`. `s_fov` is the
image bound (43–99°!), so a declared footprint (2.25° + HPBW/2) sits
beside the archive's own containment. ASTRON VO TAP
(`vo.astron.nl/__system__/tap/run/tap`) has `lotss_dr3.pointings`
(2,551 pointings, `dateallobs` = mid-MJD of every 8-h run; 5,458
observations 2014-05 → 2024-08) — exactly the per-pointing dates the
plan row wanted; DR2 tables carry only `dateobs`.

**Run.** 91 CASDA cones (84 A + 7 B; ~7 s each) + 1 LoTSS table pull;
0 failures; 12 min wall. Universal events 2014–2028, frozen ladder,
actual on-sky interval vs flat-chord window (no tolerance).

**Result.** Grazing rungs (B 1.2 / 2.5 R☉): 0 observations in 65
in-span ASKAP + 30 LoTSS events — ±0.3–0.7 d windows vs 15-min visits
at weeks-to-years cadence. **B 0.1 AU (antipode): 7 of 50 in-span
ASKAP events covered in footprint, 6 validated, 5 targets** —
wolf-359 2023-09-03 (VAST SB52549, UNCERTAIN) and 2024-09-10
(SB65727, GOOD), both 0.80° from the VAST_2257-06 centre; ross-154
2021-12-29 (RACS-high SB34957, 0.63°); ross-128 2021-09-21 (VAST
pilot SB32330, 2.45°); van-maanen 2022-03-29 (RACS-low SB38682,
footprint edge 3.12°); teegarden 2026-05-03 (FLASH SB84179, 2.18°,
UNCERTAIN); van-maanen 2026-03-29 FLASH SB83234 REJECTED. LoTSS: 0
antipode hits (3 antipodes in its footprint). A 0.1 AU: ASKAP 3
(gj-1276 VAST 2023-09-03 at 1.10°; teegarden RACS-mid 2024-11-11 at
2.09°; van-maanen VAST 2024-10-07 at the edge) + LoTSS 1 (teegarden
P043+19 2023-11-04, 1.89°, in beam). A 1.0 AU: ASKAP 211/606 in-span
events (64 targets), LoTSS 13/300 — the surveys' own transient
searches are the constraint, as with BL. Antipode visit lists
(any date): 13–89 ASKAP observations per target.

**Outcome.** The 2026-08-24 "antipode channel archivally virgin"
statement is now qualified: true for targeted radio (BL) and the
grazing rungs; false at the 0.1 AU rung, where wide-field 0.9–1.7 GHz
continuum epochs exist (~0.2 mJy/beam, broadband only — continuum
imaging dilutes CW lines ~10⁹, §5.5). Hand-offs: (a) antipode
catalogue-cone + cutout quick-look on the 6 validated epochs alongside
the VLASS six; (b) the 4 on-star narrow-rung epochs; (c) VLASS arm
refresh + v2 re-run at the yearly refresh (§5.7). Plan §5.5/§5.7/§5.8
updated; nothing committed to git.

## 17. Kepler/K2 footprint intersect (2026-09-04; plan §5.8 item 8 — closed, no hit possible)

Geometry-only, one session. Note
`surveys/kepler-crossings/notes/kepler_k2_footprint_intersect_2026-09-04.md`;
scripts `surveys/kepler-crossings/scripts/{footprint_intersect,onstar_k2_check}.py`;
results `surveys/kepler-crossings/results/`; snapshots
`runs/kepler-crossings/recon/` (26 MAST responses).

**Observer correction first.** The plan row said "vs `universal_v1`
antipodes", but Kepler flew an Earth-trailing heliocentric orbit and
was 0.04 AU (2009) → 1.14 AU (2018) from Earth. Added
`--observer kepler` to `sglsurvey/crossings.py` (Horizons −227, same
pattern as TESS/SOHO) and built `crossings/kepler_v1`
(`xng-07012bdf124d`, 2009-05-01 → 2018-11-01, 3,352 events, 0
invalid, 11 min). Matched Earth-center events shift by |Δt_ca| median
33 d / max 117 d and |Δb| up to 1 AU — invalid at every rung, not only
grazing ones.

**Footprints and dates.** `K2fov` 8.0.1 (added to `pyproject.toml`):
silicon-level channel polygons for C0–C19 (dead modules 3/7, module 4
after C10) and the prime field as campaign 1000 (roll 20° + 90° ×
season, four-season envelope). Actual campaign data ranges from MAST
CAOM TAP (`dbo.obspointing` grouped by `sequence_number`; the sync
endpoint 504s past 60 s, the async UWS endpoint works; one product
mis-tagged C12 → sequence 14, so ranges are clipped to the planned
dates ± 15 d).

**Result.** 0 of 1,435 field-active events on silicon at any b (list
max 1.03 AU). Structural: K2 boresights stayed 61°–158° from the Sun
as seen from the spacecraft (backward-facing campaigns sweep 142° →
61°, forward-facing 61° → 143°; mid-campaign at quadrature), never
within 22° (37° with the true C14 start) of the anti-sun point; the
prime field (ecliptic latitude +65°) never within 66°. Every ≤ 0.1 AU
channel position sits within 3.5° of anti-sun (grazing rungs ≤ 0.7°);
the nearest active boresight to one was 17.3° (29.6° with the true
C14 start) against a 7.6° FOV half-diagonal. Fourth independent
confirmation of the elongation gate (learnings §8), first for a
spacecraft. Sky-only: K2 did image 5 of the 14 ≤ 0.1 AU unit positions
in other campaigns, 1–9 months off the crossings, and five registry
stars have K2 light curves at b ≈ 0.9–1.0 AU (wolf-359 C14 LC+SC,
ross-128 C1, ross-154 C7, van-maanen C8 LC+SC, gj-876 C3; EZ Aqr on
C3 silicon but not targeted) — the wide rung where the archive's own
flare/transient searches are the constraint; recorded, not queued.

**Outcome.** Item 8 closed with no survey and no freeze; the 2009–2015
30-min pulse cell stays open (no continuous-cadence 1 AU mission
pointed at opposition; CoRoT is expected to fail the same gate, not
probed). `crossings/README.md` gained the `kepler_v1` row; plan status
line, §5.8 intro/row 8 and the priority-7 row updated; nothing
committed to git.

## 18. Outstanding-items ledger + radio quick-look (2026-09-04; plan §5.15, item O1 — done)

**Ledger.** With every §5.8 queue item complete, closed or blocked,
the plan gained §5.15: a three-tier table of everything that remains
(Tier 1 ready now: radio quick-look O1, STEREO-A HI-1 O2, SPHEREx
QR2 controls O3, spectral-archive family O4, WISPR geometry pass O5;
Tier 2 small standing items O6–O8; Tier 3 externally blocked O9–O13).
SPHEREx QR3 was checked live: IRSA `spherex.obscore` still serves only
`spherex_qr2`, latest epoch MJD 61241 (≈ 2026-07-20), so the deferred
six-detector joint cell and template refit no longer wait on it.

**Quick-look (O1).** `surveys/radio-crossings/scripts/quicklook_v1.py`
→ `results/quicklook_v1.json`, `report/radio_quicklook.md`, raw
responses + 26 VLASS cutouts under `runs/radio-crossings/quicklook/`.
17 epochs from the v1/v2 in-window tables (10 ASKAP in-footprint at
0.1 AU, 1 LoTSS in-beam, 6 VLASS). Catalogue stage: the VAST
full-survey SBIDs have per-SBID component catalogues in CASDA TAP
(`AS207.vast_extragal_dr1_*_sb<sbid>_components_v01`) — the only truly
epoch-resolved public catalogues — and both wolf-359 antipode epochs
plus gj-1276 A and van-maanen A are empty within 20″ (nearest
components 219–292″; 5σ ≈ 0.9–1.3 mJy). RACS epochs checked against
the release catalogues (low2/high single-epoch, consistent SBIDs;
mid's SBID column shows the 2024-11 teegarden observation is not in
the release); LoTSS DR3 mosaic empty at 71″. VAST-pilot and FLASH
SBIDs have no TAP catalogue. VLASS: CADC CAOM2 gives exact tile
times and SODA cutouts anonymously; nothing above 2.3σ at any
position in 26 tiles, and the two v1 tolerance-edge rows (gj-1276 A,
wolf-359 B, VLASS2.1) are 2.5–3 d outside the strict window.
Closest anything came: a 3.1 mJy RACS-low2 component 29″ from the
van-maanen B 2022 antipode (a fixed-sky deep-graze point; background).

**Authenticated stage (same day).** The user added OPAL credentials to
`.env`; CASDA DataLink accepts them as HTTP basic auth and hands back
a per-product cutout token, after which the SODA async job needs no
auth (the astroquery.casda flow). Stage 1b downloads the level-5
Selavy catalogues (VAST pilot v2, FLASH SB84179, RACS-mid SB67840)
and cone-searches them locally; stage 2 cut out all 12 accessible
ASKAP images (the VAST-pilot v1 products and the REJECTED FLASH
SB83234 are "no permission"). Every position empty: wolf-359
antipode pixels +0.03 / −0.12 mJy at 0.2 mJy rms; teegarden 2026
antipode +0.03 mJy at 0.07 mJy rms in the 2-h FLASH image. One
2.7σ pixel at the gj-1276 star (VAST 2023-09-03, 0.64 mJy) was
followed across all 16 VAST_2257-06 epochs
(`quicklook_gj1276_epochs.py`): mean +0.09 ± 0.06 mJy, no catalogue
component within 60″ in any epoch — the maximum of a noise series.
O1 closed; no residual. Decision §5.5 unchanged. Plan status line,
§5.5, §5.15, both radio reports updated; nothing committed to git.
