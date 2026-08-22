---
title: "WISE survey v2 — remediation plan"
status: "draft v0.3 — 2026-08-21; responds to scientific_review.md (2026-08-21); decisions 1–3 of §11 resolved and folded in"
date: 2026-08-21
---

# WISE survey v2 — plan

Response to [`scientific_review.md`](scientific_review.md). The review's
verdict: the v1 infrastructure is sound and transparent, but the
_statistical experiment_ is not — the 8-control threshold has a 1/9
per-search crossing probability (≈78 expected crossings vs 70 seen), the
completeness curves stop at threshold crossing and never see the vetoes,
326/352 W3/W4 cells are auto-vetoed by design, and the 99% locus is
metadata rather than geometry. v2 rebuilds the experiment on the same
data, frames, and registry; it does **not** change the target portfolio
or the physics cell.

Principle for v2: **one frozen decision rule, applied blind to one
end-to-end injection set, one exchangeable null ensemble, and the real
data — in that order.** Everything v1 did by hand after looking at
candidates becomes either a measured selection function or a reported
annotation.

## 0. Scope and what is kept

v2 is built in **`surveys/wise-v2/`** as a sibling of `surveys/wise/`,
so v1 scripts, results, and reports stay runnable and citable while v2
is developed; v1 is deleted and `wise-v2` renamed to `wise` only after
the v4 report ships (see `notes/project_plan.md` §10 for the
cross-survey sequence). Bulk products go to `runs/wise-v2/`; the v1
cutouts under `runs/wise/products/` are read, never modified.

Kept verbatim (no re-run):

- Registry v1.5 (88 endpoints / 77 corridors), `hypotheses.md` §§1–8
  physics (550–10,000 AU, |µ| bound, duty ≥ 0.5, unresolved source).
- Frame discovery (`coarse_v1`), precise pass (`precise_v1`), catalog
  screen (`screen_v1`), and the 27 GB of `products/cut` + `products/msk`
  cutouts. v2 re-samples from these; nothing is re-fetched unless the
  quality-mask sensitivity study (§7) needs frames v1 discarded.
- `sglsurvey/photometry.py` matched filter, `sample_tensor.py` tensor
  layout, `sglsurvey/vetting.py` static-source test.

Withdrawn now, before any new computation, as a v3.1 erratum to
`report/wise_survey_v3.md` (§9 of this plan):

- "90%-complete exclusion" across the stated cell.
- Any statement that the 8-control threshold is an empirical false-alarm
  rate.
- End-to-end W3/W4 ~100 km / 300 K exclusions (relabel as raw threshold
  sensitivity).
- "99% confidence locus" and "every stage content-addressed".

## 1. Hypothesis freeze v2.0 (new version, not an edit)

`hypotheses.md` gets a v2.0 block. Physics unchanged; what changes is
the _decision rule_ and the geometry spec, which v1 left to "the run
config once the pipeline exists" (§9). v2.0 must fix, before any v2
result is looked at:

1. **Motion bound.** Declare the bound **component-wise L∞**
   (|µ*α|, |µ*δ| ≤ 1″/yr), matching the 5×5 grid actually searched.
   Re-declaring as L2 would force a new grid; L∞ is the cheaper honest
   fix and the review says either is acceptable if stated. Record it in
   the Constraint field `residual_motion_bound` with `norm: "linf"`.
2. **Distance prior vs grid.** The log-uniform prior is what _weights_
   completeness statements; the 1/z-uniform grid is the numerical
   tabulation. v2 injections are drawn continuously in z from the
   log-uniform prior (§4), so the prior finally governs placement as
   §3 always claimed.
3. **Temporal models** (four, each its own injection family):
   persistent; independent flicker p=0.5 (v1's model); visit-scale
   on/off (one Bernoulli draw per NEOWISE visit, p=0.5); long-block
   (on for a contiguous ≥50% of the mission span, random start). A cell
   is "covered at duty ≥ 0.5" only where the _worst_ of the four
   reaches the completeness target.
4. **Global detection statistic and error target.** Per cell
   (endpoint × role × band), the statistic is the leave-one-out
   exceedance ratio R = S_max / T_loo, exactly as in
   `surveys/ztf/scripts/look_elsewhere.py`. The survey-wide statistic is
   the maximum R over all cells in the confirmatory set. **Target:
   family-wise P(any false survey detection) ≤ 0.05**, calibrated on the
   null ensemble of §3. No per-cell "exceeds" verdict is reported as a
   detection claim; cells are ranked by R with a global p-value.
5. **Complete candidate rule** (§6): a cell is a _candidate_ iff R
   exceeds the FWER-calibrated global threshold. A candidate may be
   rejected only by a **calibrated veto** frozen here — the
   flux-consistent catalogued-static-source / halo test, the held-out-
   epoch prediction test, or the W3/W4 confirmation procedure (§5) —
   each with an injection-measured selection function. Phase balance,
   W1:W2 significance ratio, cryo epoch count, and bare proximity are
   **annotations** carried on the Candidate record, never grounds for
   rejection.
6. **Observer.** Topocentric WISE position from the L1b frame metadata
   (the `-int` header carries the spacecraft state; if not, use the
   IRSA `scan`/`frame` ephemeris table). Correct the accuracy-budget
   arithmetic: 6,900 km geocentric radius ↔ ~17 mas at 550 AU, not
   2.5 mas.
7. **Confirmatory hold-out.** A **random endpoint split stratified by
   confusion class** (the low / mid / high-confusion corridor classes
   already used to order the v1 batches), drawn with a recorded seed and
   hashed into the v2.0 freeze — not a split by batch, since batches 4–5
   are systematically more crowded than 1–3:
   - _Development set:_ ~30 % of endpoints (≈ 26), every confusion
     class represented in proportion — all rule tuning, threshold
     studies, veto selection-function work.
   - _Confirmatory set:_ the remaining ~70 % (≈ 62) — run once, blind,
     under the frozen rule. If the rule changes after touching this set,
     the confirmatory claim is void and the report must say so.
   Both roles of an endpoint, and both components of a binary, go to
   the same side of the split (they share frames).
8. **Epoch hold-out: calibrated post-hoc, not pre-registered.** NEOWISE
   ended in 2024-08 and no further epochs will exist, so a true
   pre-registered epoch hold-out would cost ~0.15 mag of W1/W2 depth
   (≈ 25 % of epochs) unconditionally and shorten every constraint's
   epoch range, for a test that only matters if a candidate survives.
   Instead: **all epochs are searched.** For any retained candidate, the
   _held-out-epoch prediction test_ refits (z, µ) on epochs ≤ 2021-12-31
   and performs forced photometry at the predicted 2022–24 positions.
   Because the candidate was _selected_ using those epochs, the test is
   optimistically biased; the bias is measured, not assumed: the identical
   refit-and-predict procedure is run on every null-ensemble exceedance
   (§3) and on recovered injections (§4), and its false-pass and
   true-pass rates are frozen with the rule and quoted alongside any
   result. On the **development set only**, the test is additionally run
   as a true hold-out (2022–24 excluded from the search) to measure its
   power at full independence — free there, since the development set is
   exploratory. Cross-archive confirmation (ZTF / PS1 / SPHEREx on the
   same corridor at the predicted position) remains the primary
   independent evidence for any candidate, being the only data never
   part of the WISE search.

Freeze gate: v2.0 hash recorded in `configs/v2_0_freeze.json`; no
script in §3–§8 runs against the confirmatory set before that commit.

## 2. Geometry: propagate the covariance for real

(Review §3.) Stage: `scripts/precise_pass.py` + `sample_tensor.py`.

1. Call the sglseti seeded Monte Carlo locus propagation (the routine
   the hypothesis names and production never calls) per endpoint × role
   with N = 2,000 draws of the full registry covariance (including
   correlations and reference epoch), seed and N in the run config.
2. Emit the 99% cross-track envelope half-width per epoch
   (`sigma_xt_99`) into the precise-stage IntersectionEvaluation, and a
   per-endpoint summary table: max envelope vs W1 PSF FWHM (6.1″).
3. **Decision rule for sampling:** where `sigma_xt_99` ≤ 0.5 × FWHM for
   every epoch (expected for nearly all Gaia/HIPPARCOS endpoints), the
   nominal locus is an adequate search position and the envelope is
   reported as a verified approximation. Where it is not (candidates:
   orbit-derived components — GJ 65 A/B, Sirius B, ε Ind Ba/Bb
   photocenter, Luhman 16), the tensor gains a cross-track dimension
   (3 or 5 offsets spanning the 99% envelope) and the stack statistic
   maximises over it, with that extra trial counted in the null.
4. Empirical coverage check: recompute the locus for 50 endpoints using
   an independent ephemeris (Gaia DR3 catalogue positions propagated
   with `astropy` only) and report residuals against the sglseti locus.
5. Positive control (review §4.5 / priority 4.4): port
   `surveys/ztf/scripts/asteroid_control.py` to WISE — a numbered
   main-belt asteroid near the W1 single-exposure limit, JPL Horizons
   ephemeris as the "trajectory model", driven blind through cutout →
   flux map → tensor → stack → v2 decision rule. Must be recovered as a
   candidate; its recovered flux vs Horizons-predicted H, G gives the
   first empirical throughput check of the whole chain.

Deliverable: `results/v2_geometry_summary.md`; tensors rebuilt only for
endpoints that need the cross-track dimension.

## 3. Null ensemble and the survey-wide false-alarm model

(Review §1 — critical.)

The 8 spatial offsets stay as the _threshold_ definition per cell (they
are fine as a local normaliser), but the _error rate_ comes from a much
larger, exchangeable ensemble built from three independent null
constructions, all preserving local coverage and contamination:

| Null construction          | How                                                                                                                                                              | Trials/cell | What it catches                                 |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | ----------------------------------------------- |
| Spatial offsets (extended) | 8 → 48 offsets on a ring 15–45″ from the locus, both axes, same cutouts                                                                                          | 48          | local background, static confusion              |
| Time scrambling            | permute epochs _within_ a phase and band across the (z, µ) grid so the trajectory coherence is broken but per-epoch noise, cadence, and phase structure are kept | 200         | cadence/phase artefacts, single-epoch dominance |
| Trajectory randomisation   | evaluate the real cutouts on the locus of a _different_ endpoint's (z, µ) sky motion, re-centred on this corridor (motion pattern exchangeable, sky fixed)       | 50          | anything tied to the motion model itself        |

Implementation: `scripts/null_ensemble.py` reading the v2 tensors
(the spatial set needs `sample_tensor.py --offsets ring48`, which
re-samples cutouts only — no re-fetch). Output per cell: the empirical
distribution of R under each construction, and the pooled one.

Statistical products:

1. **Per-cell rank statement** (review §1 rec. 4): the real R's rank among
   N exchangeable controls, with a binomial interval — replaces "T = max
   of 8 controls" language everywhere.
2. **Survey-wide max-R null**: on the development set, draw 10,000
   survey-wide pseudo-experiments (one null R per cell per draw, from
   the pooled per-cell distributions), take max over cells → the FWER
   threshold R_FWER at **α = 0.05**, and the Benjamini–Hochberg q-values
   for the per-cell view. Report both; **FWER is the candidate rule, BH
   is informational** (decided 2026-08-21).
3. **Check of exchangeability**: the three null constructions must give
   compatible R distributions per cell (KS test); cells where they do
   not are flagged `null_unstable` and cannot carry a constraint.
4. **Trial accounting**: the z × µ × (cross-track) grid is inside R
   already (max over grid); endpoints × roles × bands (704, or more
   with cross-track cells) are the family.

Expected outcome: ~70 cells crossing T*loo is exactly what this
construction predicts; the question v2 answers is whether any cell's R
is unusual \_survey-wide*. Under the null none should be.

## 4. End-to-end completeness: inject into images

(Review §2 and §4 — critical/high.)

Injection moves from the tensor to the cutout. New module
`sglsurvey/inject.py` + `scripts/inject_pipeline.py`:

1. **Where:** add the source to each `-int` cutout in DN _before_
   `build_flux_map` (background estimation, masking, matched filter,
   calibration all run unchanged). Masks are untouched — an injected
   source landing on a fatal-bit pixel is lost, as it would be in
   reality.
2. **PSF:** replace the circular Gaussian with the WISE empirical PRF
   (IRSA per-band PRF grid, 9 focal-plane positions; ~dozens of MB,
   fetched once and content-hashed). Place by frame (x, y) with subpixel
   phase drawn uniformly; the matched filter stays Gaussian, so the
   Gaussian-vs-true-PRF throughput loss is _measured_, not assumed.
3. **Spectrum:** inject in physical units. Two spectral families per
   band set: a 300 K blackbody (W3/W4 waste-heat cell) and a flat-Fν /
   G2V-reflected spectrum (W1/W2). Apply the WISE colour corrections
   from the Explanatory Supplement (§IV.4.h) when converting to DN via
   the frame `MAGZP`. Constraint records carry `spectrum_model`.
4. **Placement:** continuous — z from the log-uniform prior on
   [550, 10,000], (µ*α, µ*δ) uniform on the L∞ box, cross-track offset
   from the propagated 99% envelope, duty model from the four families
   of §1.3, flux from a log-uniform band around the v1 m90 (±2 mag).
   No grid nodes, no analytic mismatch factor.
5. **Volume:** per cell (endpoint × role × band) **400** injections
   (vs 32 × 64 nodes analytic). 704 cells × 400 = ~280k injections;
   each touches the ~50–300 cutouts of its track. Cost is dominated by
   matched-filter convolution on modified cutouts — cache the
   _unmodified_ flux map and only re-convolve a stamp around the
   injection (linear operation, so flux map + stamp response is exact).
   Budget: ~1 CPU-day on the dev machine; run in corridor batches like
   `run_scaleup.sh`.
6. **Two completeness curves, always reported separately:**
   - _threshold completeness_ — injection's R ≥ R_FWER;
   - _final-candidate completeness_ — and the injection survives
     every independent rejection test of §5–§6 (run blind on it).
     Each as a function of magnitude, binned in log z with edges at
     reciprocal-distance **midpoints** so the intervals tile 550–10,000 AU
     with zero gap (fixes the 9.2% hole), with Wilson 68/95% binomial
     intervals from the injection counts.
7. **Constraint ledger:** m90 (and m50) per interval from the
   final-candidate curve, worst temporal model, with CI. The Constraint
   schema gains `completeness_kind ∈ {threshold, final_candidate}`,
   `ci_68`, `n_injections`, `prf_model`, `spectrum_model`. v1's 5,632
   records are superseded by link, never deleted.
8. **Invariant tests** (review §9): Gaussian normalisation, PRF
   normalisation, DN↔mag round-trip, interval tiling, "injection of
   zero flux changes nothing", stamp-response = full-reconvolution on a
   random sample. `tests/test_wise_inject.py`, run in CI before any
   batch.

## 5. W3/W4: sensitivity, not exclusion — plus a confirmation path

(Review §2 — critical; priority 3.)

1. Remove the `single-visit cryo cell` veto from the candidate rule. A
   cryo cell has the cadence it has; ~12 exposures in one visit is the
   _expected_ signature, not evidence against.
2. W3/W4 cells are searched and ranked with the same R statistic and
   global threshold as W1/W2. Their null ensemble is dominated by the
   time-scramble and trajectory-randomisation constructions because the
   spatial-offset set is small in a 12-frame visit — report which.
3. **Independent confirmation procedure** for any W3/W4 cell exceeding
   R*FWER (frozen in v2.0):
   a. W1/W2 forced photometry at the \_same* (z, µ) and epoch on the
   cryo frames — a 300 K, ~100 km body predicts a W2/W3 flux ratio;
   absence of W2 at the predicted level at ≥ 3σ rejects the thermal
   interpretation, presence supports it.
   b. Post-cryo NEOWISE W1/W2 along the propagated track (10+ years of
   epochs) — tests persistence, the duty-cycle hypothesis's own
   prediction.
   c. AllWISE / CatWISE / unWISE catalogue search at the W3/W4 peak
   position for a static counterpart (`sglsurvey/vetting.py`).
   d. Only if (a)–(c) are all inconclusive: external archival imaging
   (Spitzer IRAC/MIPS, AKARI) and a follow-up request.
4. **Reporting rule:** W3/W4 Constraint records carry
   `completeness_kind = threshold` only, `exclusion_claim = false`, and
   the physical-radius table is labelled "raw threshold sensitivity,
   300 K blackbody, empirical PRF" — no end-to-end exclusion is written
   until (3) has been exercised on at least one real exceedance and
   on injections (so its true-positive loss is known).

## 6. Vetoes: calibrated rejection tests vs annotations

(Review §6 — high.) The review's objection is to **uncalibrated**
vetoes, not to vetoes as such. The distinction v2 draws:

- A **calibrated veto** is a rule that (i) rests on independent evidence
  that the excess has a non-SGL origin, (ii) cannot reject an in-scope
  signal that is _actually present and detectable_ — it can only lose a
  signal that was genuinely swamped, which is a completeness cost — and
  (iii) has its selection function measured by the §4 injections and
  charged to the final-candidate completeness curve.
- An **annotation** is anything that merely raises suspicion, or that
  _could_ reject a present in-scope signal (an intermittent source really
  can be single-phase). Annotations are carried on the Candidate record
  and reported, never used to reject.

Why the catalogued-static-source test qualifies as a calibrated veto:
an SGL source at 550–10,000 AU has a parallax amplitude of ~20″–375″,
so it cannot _be_ a multi-detection fixed-position catalogue entry
(CatWISE is an all-epoch coadd; a source moving many PSF widths per
year never forms a compact entry). A static star within a PSF width of
the track at the stack-dominant epochs therefore either produced the
excess (veto correct) or swamped a real source at the same place
(completeness loss, measured). The test is only sharp enough to count
as independent evidence with two additions to
`sglsurvey.vetting.static_source_test` (§8.5): a **flux-consistency**
check — the catalogue magnitude pushed through the same matched filter
at the test epochs must account for the measured S within a factor 2
(a distant W1 = 16 star cannot explain S = 40) — and the
**parallax-phase signature** (a static star near the phase-0 cluster
gives high S0 / low S1; a track has both). Proximity alone is an
annotation.

| v1 rule                                            | v2 status                                                                                                                                                                                                          | Selection function measured by                                                                                              |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------- |
| static CatWISE source within 6″ at the major phase | **calibrated veto** when flux-consistent _and_ phase-imbalanced; annotation otherwise                                                                                                                              | injections landing within 1 PSF of catalogued stars (completeness loss charged to the final-candidate curve); null ensemble |
| bright-star halo                                   | **calibrated veto** in its flux-consistent form (PSF-wing flux of the catalogued star at the track position explains S); annotation otherwise                                                                      | same; the 48-offset null samples halos                                                                                      |
| single-phase (S0/S1 < 0.3)                         | annotation — an intermittent in-scope source can be single-phase                                                                                                                                                   | injections under the four temporal models (§4)                                                                              |
| single-epoch dominance                             | annotation; the weight cap already limits it                                                                                                                                                                       | injections with one bright artefact epoch added                                                                             |
| W1:W2 significance ratio "star-like"               | **dropped** — a ratio of significances is not a colour; replaced by the physical W1/W2 and W2/W3 _flux_ ratio test with errors (§5.3a), which is itself a calibrated veto only for the thermal-cell interpretation | injected spectra                                                                                                            |
| single-visit cryo cell                             | **removed** (§5)                                                                                                                                                                                                   | —                                                                                                                           |
| "+<10 % marginal excess"                           | **removed**; superseded by the global p-value                                                                                                                                                                      | —                                                                                                                           |

A candidate (R ≥ R*FWER) may therefore be \_rejected* only by: the
calibrated static/halo test above; the calibrated post-hoc held-out-epoch
prediction test (§1.8); or the W3/W4 confirmation procedure of
§5.3. Anything else leaves it `retained-ambiguous`, reported as such
with its annotations. `adjudicate_exceedances.py` is replaced by
`adjudicate_v2.py` implementing exactly this and nothing discretionary.
Every rejection records which test fired and the injection-measured
loss of that test for the cell's band and z interval.

## 7. Exposure-quality mask and sensitivity analysis

(Review §7 — moderate.)

1. Primary mask (pre-specified): `qual_frame > 0`, `qual_scan ≥ 5`
   (where present), `saa_sep > 0`, `moon_masked == 0` on the frame,
   fatal-bit mask pixels — per the NEOWISE Explanatory Supplement
   §III.1.c. Metadata is already in the coarse snapshots; no re-query.
2. Two alternates run on the development set only: _strict_ (also drop
   `qual_scan < 10`, Moon separation < 30°) and _loose_ (v1's
   `qual_frame ≠ 0` only). Report the change in per-cell R, candidate
   count, and m90 between the three; the primary is the reported result
   and the spread is quoted as a systematic.
3. Null ensemble construction (§3) applies the same mask to controls —
   automatically true for spatial/time-scramble nulls since they share
   frames; assert it in the trajectory-randomisation null.

## 8. Reproducibility and bookkeeping

(Review §9 — moderate.)

1. `scripts/manifest.py`: content hash (sha256 of bytes + schema
   version) of every cutout, mask, tensor, PRF file, config, registry
   and snapshot used by a run → `runs/wise/v2/manifest_<run>.json`; the
   AnalysisRun's `observation_set_hash` / `intersection_set_hash` become
   real hashes of those manifests.
2. Tensor files carry the hash of their cutout inputs in their header;
   `--only-missing` refuses a tensor whose recorded input hash differs
   from the current cutouts (stale-product detection mandatory).
3. `tests/` (pytest, run by `run_v2.sh` before each batch): §4.8
   invariants, plus record-count reconciliation, supersession-link
   resolution, and hypothesis-hash pinning.
4. Report tables are generated from the ledger by `scripts/report_tables.py`;
   the build fails on count mismatch.

### 8.5 Changes to shared code in `sglsurvey/`

`sglsurvey/vetting.py` is consumed by the PS1, ZTF, and joint scripts,
so every change below is additive with defaults that leave the other
surveys' behaviour unchanged; v2 of those surveys (project plan §10)
adopts the new paths deliberately, not by accident.

| Module                          | Change                                                                                                                                                                                                                                   | Kind     |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- |
| `records.py` `Constraint`       | add `completeness_kind` (threshold / final_candidate), `ci_68`, `n_injections`, `prf_model`, `spectrum_model`, `motion_bound_norm`; keep them out of the `constraint_id` hash so v1 IDs are unchanged                                    | additive |
| `records.py` `Candidate`        | add `annotations`, `global_p_value`, `rank_statement`, `rejection_test`                                                                                                                                                                  | additive |
| `photometry.py`                 | seam between cutout read and `matched_filter` for image-level injection (`inject=` callable or a `read_cutout` / `flux_map_from_arrays` split); `matched_filter` already takes an arbitrary `kernel`, so the empirical PRF needs nothing | additive |
| `vetting.py`                    | `static_source_test`: add flux-consistency (catalogue mag → matched-filter S at the test epochs) and return the phase S0/S1 signature; add `parallax_phase_test()`; `track_position` gains an optional cross-track offset for §2.3 cells | additive |
| `geometry.py` `GeometryContext` | optional per-epoch observer (`observer_for(mjd)`) for the topocentric WISE state; reuse whatever ZTF did for Palomar rather than adding a second mechanism                                                                               | extend   |
| **new** `inject.py`             | PRF loading, DN placement with subpixel phase, spectrum → DN with colour corrections, stamp-response linearity trick                                                                                                                     | new      |
| **new** `nulls.py`              | promoted from `surveys/ztf/scripts/look_elsewhere.py`: leave-one-out R, pooled null ensembles, survey-wide max-R FWER / BH; ZTF's script becomes a thin caller                                                                           | promoted |
| **new** `manifest.py`           | content-hash manifests and stale-product detection (§8.1–8.2), usable by every adapter                                                                                                                                                   | new      |

Stays WISE-specific under `surveys/wise-v2/`: PRF files, colour
corrections, quality masks, `adjudicate_v2.py`, `null_ensemble.py`
configuration, the corridor table.

**Cross-survey dependency to resolve before v1 deletion:**
`surveys/{ztf,panstarrs,spherex}/scripts/*_corridors.py` and the
overlay builders import `CORRIDOR_OF` / `MEMBERS` from
`surveys/wise/scripts/wise_corridors.py`. Promote the corridor table
to `sglsurvey/corridors.py` (or `targets/`) during step B so nothing
outside `surveys/wise/` depends on it.

## 9. Reporting

1. **Now (before v2 computation):** `report/wise_survey_v3.md` erratum
   v3.1 — withdraws the four claims listed in §0, fixes the 17 mas
   observer arithmetic, relabels W3/W4 tables, replaces "FAR < 1/8"
   with the rank statement. Links `scientific_review.md`.
2. **`report/wise_survey_v4.md`** after v2: separates, in this order,
   (i) the confirmatory-set result under the frozen rule with the global
   p-value; (ii) the development-set result (labelled exploratory);
   (iii) threshold vs final-candidate completeness with CIs per band and
   temporal model; (iv) W3/W4 sensitivity with the explicit no-exclusion
   label; (v) retained-ambiguous candidates, if any, with their
   annotations; (vi) the systematic spread from §7; (vii) the positive
   control; (viii) a scope statement: targeted coverage of a frozen
   88-endpoint portfolio — no population or Picky-Network inference
   (review §8). Population inference is out of scope for v2 and is
   deferred to a separately pre-registered model.

## 10. Order of work and gates

| Step | Work                                                                                     | Gate to pass                                                        |
| ---- | ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| A    | §9.1 erratum; §1 hypothesis v2.0 freeze incl. stratified endpoint split (seed recorded); `configs/v2_0_freeze.json` | committed before any v2 script runs                                 |
| B    | §8.1–8.3 manifests, stale detection, invariant tests                                     | tests green on v1 products                                          |
| C    | §2 covariance propagation + observer fix; asteroid positive control                      | asteroid recovered blind; envelope table reviewed                   |
| D    | §7 primary mask; §3 null ensemble on the **development set**                             | three null constructions exchangeable per cell (or flagged)         |
| E    | §4 image-level injections on the development set; §6 selection functions                 | CI widths acceptable (m90 ± 0.15 mag at 400/cell)                   |
| F    | Decision rule sanity on development set: candidates, annotations, rejections             | rule frozen; any change → back to A and the hold-out is re-declared |
| G    | Confirmatory set: null, injections, real data, **blind**                                 | run once                                                            |
| H    | §5.3 confirmation procedure on any exceedance; post-hoc held-out-epoch test with its null/injection-calibrated pass rates (§1.8) | —                                                                   |
| I    | §9.2 report v4 from the ledger                                                           | tables build from records                                           |

Rough cost: A–B a day; C two days (sglseti MC + control); D one
CPU-day; E ~1 CPU-day per half of the portfolio (the stamp-response
trick is what makes this affordable — without it, full re-convolution
of 280k × ~150 cutouts is ~10× more); F–I two to three days of
analysis and writing. Disk: v2 tensors with 48 offsets are ~6× v1's
19 GB → reuse the PS1 purge pattern (`purge_products.py`) and keep only
the pooled null statistics per cell after the ensemble runs.

## 11. Decisions

Resolved 2026-08-21 and folded into the sections cited:

1. **FWER α = 0.05 is the candidate rule; BH q-values are informational**
   (§1.4, §3.2).
2. **Hold-out split: random endpoint split stratified by confusion
   class**, seed recorded in the freeze — not by batch (§1.7).
3. **Epoch hold-out: calibrated post-hoc.** All epochs searched; the
   2022–24 prediction test is applied to retained candidates with its
   pass rates measured on the null ensemble and injections, and run as a
   true hold-out on the development set only to measure power (§1.8).

Still open `[DECISION]`:

4. Empirical PRF vs a per-band mean PRF. Full focal-plane-position PRFs
   are more faithful but add a dependency on the frame (x, y) of every
   injection; a mean PRF with a measured scatter may be adequate at the
   6″ scale. Decide after measuring the throughput spread on 1 corridor.

## 12. Learnings to carry into the other surveys' v2

Recorded here so the WISE v2 work product includes its own transfer
list; the sequencing lives in `notes/project_plan.md` §10. Everything
in §1, §3, §4, §6, and §8 applies to PS1, ZTF, SPHEREx, and the joint
stage with the substitutions below; §2 (covariance) and §5 (W3/W4) are
partly specific.

| Item                                                   | PS1 / ZTF                                                                                                                                                                        | SPHEREx                                                      | joint PS1+ZTF(+WISE)                                                  |
| ------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ | --------------------------------------------------------------------- |
| 1/9 threshold problem (§3)                             | identical 8-control construction → same null ensemble + FWER; `look_elsewhere.py` already exists for ZTF                                                                         | 16 controls → 1/17 per search, same issue                    | the joint stage multiplies the family; needs one FWER across archives |
| Image-level injection (§4)                             | PS1 warp / ZTF sci cutouts retained? if purged (`purge_products.py`), re-fetch per batch as `run_common_t0.sh` does; PSFs: PS1 per-skycell, ZTF per-quadrant from the sci header | PSF cube is already in the product (`spherex_kernel`)        | inject the _same_ physical source into all archives                   |
| Physical spectrum + colour terms                       | AB system; reflected-solar and flat-Fν families                                                                                                                                  | per-pixel wavelength → inject an SED, not a magnitude        | colour consistency 0.5–22 µm becomes a calibrated veto                |
| Temporal models (§1.3)                                 | ZTF nightly cadence supports all four; PS1 single-phase cadence makes visit-scale ≈ persistent                                                                                   | ~2–3 phases: visit-scale model dominates                     | hold-out epochs: last ZTF year                                        |
| Calibrated vs annotation vetoes (§6)                   | phase-split, season, role-coincidence rules → annotations; catalogued-star test (already in `vetting.py`) → flux-consistent calibrated veto                                      | template-coverage and bright-static-neighbour rules likewise | same                                                                  |
| Hold-out split (§1.7–1.8)                              | stratified random endpoint split as WISE; epoch hold-out can be pre-registered for ZTF (future epochs keep arriving) rather than post-hoc                                        | northern vs southern run is a natural split; QR3 epochs are a true hold-out | —                                                                     |
| Gap-free intervals, L∞ motion bound, log-uniform draws | identical                                                                                                                                                                        | identical                                                    | identical                                                             |
| Covariance propagation (§2)                            | same check; optical PSFs (1–2″) make the envelope test _stricter_ — more endpoints may need the cross-track dimension                                                            | 6″ pixels: as WISE                                           | —                                                                     |
| Observer                                               | Palomar / Haleakalā already topocentric                                                                                                                                          | spacecraft state from the product header                     | —                                                                     |
| Positive control                                       | `asteroid_control.py` exists (ZTF); add for PS1                                                                                                                                  | add (asteroid through the spectral-image path)               | —                                                                     |
| Manifests + invariant tests (§8)                       | identical                                                                                                                                                                        | identical                                                    | identical                                                             |
