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
`sglsurvey/v2/` (profiles in `surveys/{ztf,panstarrs,spherex}-v2/`,
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

## 7. TESS crossings survey — execution log (2026-08-24; in flight)

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
