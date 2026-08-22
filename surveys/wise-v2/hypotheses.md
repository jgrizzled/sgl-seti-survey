---
title: "WISE/NEOWISE survey — hypothesis freeze v2.0 (decision rule and geometry specification)"
status: "v2.1 — frozen 2026-08-22 after the development-set rule check (step F); v2.0 → v2.1 change recorded in §8; no v2 script has run on the confirmatory set; hash in configs/v2_0_freeze.json"
date: 2026-08-22
---

# Hypothesis freeze v2.0

A new version, not an edit: `surveys/wise/hypotheses.md` v1.0–v1.5
(the physics cell) stands unchanged and is incorporated by reference —
targets and registry v1.5 (88 endpoints / 77 corridors), Rx and Tx as
separate hypotheses, relay distance 550–10,000 AU with a log-uniform
prior, unresolved point source, duty cycle ≥ 0.5, station-kept relay
with bounded residual motion. v2.0 fixes everything v1 §9 left to "the
run config once the pipeline exists": the decision rule, the null
model, the completeness definition, the geometry specification and the
hold-out. It responds to `surveys/wise/scientific_review.md` through
`surveys/wise/v2_plan.md`, whose section numbers are cited below.

Principle: **one frozen decision rule, applied blind to one end-to-end
injection set, one exchangeable null ensemble, and the real data — in
that order.**

## 1. Geometry and observer (plan §1.1, §1.6, §2)

1. **Motion bound: component-wise L∞.** |µ*α| ≤ 1″/yr and |µ*δ| ≤ 1″/yr
   (the 5 × 5 grid actually searched; corner nodes reach √2″/yr in L2).
   Constraint records carry `motion_bound_norm = "linf"`.
2. **Distance prior vs grid.** The log-uniform prior on [550, 10,000] AU
   governs injection placement and the weighting of completeness
   statements; the 64-node grid uniform in 1/z is the numerical
   tabulation of the search and carries no prior meaning.
3. **Observer.** Topocentric WISE spacecraft position from the L1b
   frame headers (SUN2SCX/Y/Z plus the Sun's barycentric position from
   the pinned astropy ephemeris), tabulated per frame and registered as
   the sglseti programmatic observer `wise-l1b-spacecraft` whose
   identity is the table's content hash. The v1 Earth-centre
   approximation was bounded by 6,900 km / 550 AU ≈ 17 mas (v1 stated
   2.5 mas in error; the measured locus shift for a 6,900 km
   displacement is ≤ 0.03″).
4. **Covariance propagation.** Per endpoint × role, the sglseti seeded
   Monte Carlo (`propagate_locus_uncertainty`, N = 2,000, seed
   recorded) of the full registry covariance at z = 550 and 10,000 AU
   and three epochs (2010.0, 2017.0, 2024.0) gives the 99 % cross-track
   envelope half-width `sigma_xt_99`. Where max `sigma_xt_99` ≤ 0.5 ×
   FWHM(W1) = 3.05″ for all epochs, the nominal locus is the search
   position and the envelope is a verified approximation. Where it is
   not, the cell gains a cross-track dimension (offsets 0, ±σ_xt99 or
   0, ±½, ±1 × σ_xt99 for > 2 FWHM) searched and null-calibrated as an
   extra grid axis. An independent astropy-only propagation of the
   Gaia DR3 catalogue solution for 50 endpoints is reported as an
   empirical coverage check.
5. **Positive control.** A numbered main-belt asteroid near the W1
   single-exposure limit, JPL Horizons ephemeris (observer `@-163`,
   the WISE spacecraft) as the trajectory model, driven blind through
   cutout → flux map → tensor → stack → the decision rule of §3. It
   must be recovered as a candidate; its recovered flux versus the
   Horizons-predicted H, G brightness is the chain's first empirical
   throughput check.

## 2. Temporal models (plan §1.3)

Four injection families, each at duty 0.5 unless persistent:
`persistent`; `flicker` (independent Bernoulli p = 0.5 per exposure —
v1's model); `visit` (one Bernoulli draw per NEOWISE visit, visits
defined by gaps > 5 d in the cell's epochs); `block` (on for one
contiguous span covering 50 % of the mission, random start). A cell is
"covered at duty ≥ 0.5" only where the **worst** of the four reaches
the completeness target.

## 3. Detection statistic, null ensemble and decision rule (plan §1.4–1.5, §3)

1. **Cell.** endpoint × role × band; 176 pair × 4 bands = 704 cells
   (more if §1.4 adds cross-track cells).
2. **Stack statistic.** S(z, µ) = A/√B with inverse-variance weights
   capped at 20 × the per-node median (v1 v0.2.0 effective-epoch
   floor), usable samples = finite with PSF good-fraction ≥ 0.7, nodes
   with < 5 epochs undefined; no single-epoch clip. S_max = max over the
   grid.
3. **Local threshold.** T = max of the grid maxima of the **8 designated
   spatial controls**: offsets (±20″, 0), (±30″, 0), (±40″, 0),
   (0, ±30″) from the locus, same cutouts. R = S_max / T (leave-one-out
   for the designated controls themselves).
4. **Null ensemble per cell.** Three constructions are computed, all
   preserving local coverage and contamination:
   - *spatial ring*: 48 offsets, radii 20/30/40″ × 16 position angles
     (the designated 8 are members), same cutouts — **the exchangeable
     ensemble**: under the null the real track is one more ring member;
   - *phase-coherence scramble*: 200 draws; within each parallax
     phase the (z, µ) node labels of the real-trajectory stack are
     relabelled by a common random cyclic z-shift and µ permutation,
     independently per phase, before the phases are recombined. Every
     epoch keeps its values, variance, cadence, phase and static-sky
     contamination (the annual parallax returns a track to the same
     static sources at the same node every year); only the alignment
     between phases is randomised;
   - *trajectory randomisation*: 50 draws; the cutouts sampled on the
     (z, µ) sky-motion pattern of another endpoint (same role, other
     corridor, seeded draw) re-centred on this corridor at T0.
   Measured on the development set before this freeze (2026-08-22):
   the ring null is locally matched (median R ≈ 0.5–0.8, q95 ≈ 1.0–1.5
   in nearly every cell); the trajectory-randomised tracks sample
   *other* parts of the corridor and have corridor-dependent tails
   (q95 up to 4–18 where a donor track crosses a bright source) — they
   differ from the ring at KS α = 0.01 in two thirds of cells; the
   phase scramble is degenerate in single-phase-contaminated cells
   (its draws sit at R ≈ R_real). Therefore **the candidate rule uses
   the ring only**, and the other two constructions are carried as
   per-cell annotations: p_trajectory = P(donor R ≥ R_real) and
   p_phase = P(phase-scramble R ≥ R_real). A per-epoch scramble (node
   labels permuted independently per epoch) is a diagnostic noise-only
   floor, reported and never pooled (gj-625/rx/W1: scramble maxima ≈ 12
   vs spatial controls 26–67).
5. **Null properties that void a cell** (judged on the null alone,
   never on the real R): `null_unstable` if the inner (20″) and outer
   (40″) ring members differ at KS α = 0.01 (the local null depends on
   the offset radius — a strong gradient such as a bright-star halo);
   `null_heavy_tail` if the ring maximum exceeds 2.5 × the ring's 95th
   percentile (the local null is dominated by one offset position).
   Such cells are excluded from the family, cannot be candidates and
   carry no constraint (`not_constrainable`); their number is
   reported. KS tests between all constructions are reported.
6. **Survey-wide statistic and error target.** Per cell R̃ = R / q95,
   q95 the 95th percentile of the cell's ring null (so heavily
   contaminated cells do not set the family-wide scale). The family
   statistic is max R̃ over the cells of the set being tested; its null
   is 10,000 pseudo-experiments drawing one ring R̃ per cell (cells
   treated as independent — conservative for positively correlated
   cells). **Candidate rule: R̃ ≥ R̃_FWER, the 95th percentile of the
   pseudo-experiment maxima (family-wise α = 0.05).** R̃_FWER is a
   deterministic function of the family's own null ensemble, computed
   before any real R̃ of that family is looked at. BH q-values on the
   per-cell p = (k + 1)/(N + 1) are reported as information only. Every
   cell carries its rank statement (rank among the 48 exchangeable
   controls with Wilson 68/95 % intervals) and its global p-value
   P(max null R̃ ≥ R̃). On the partial development set the rule implies
   a required excess of R ≈ 2.5 × T in a typical cell (≈ 3.9 without
   the heavy-tail exclusion); the price of a controlled family-wise
   error over ~500 heavy-tailed local nulls, to be quantified by the
   injections.
7. **Calibrated vetoes** (plan §5–6) — the only grounds for rejecting a
   candidate, each with an injection-measured selection function:
   a. *flux-consistent catalogued static source / halo*: CatWISE2020
      sources (≥ 3 detections) within 30″ of the track, pushed through
      the band's PRF-vs-Gaussian radial response at every epoch's
      separation and stacked with the cell's weights, account for the
      measured S within a factor 2 overall and in the dominant phase;
   b. *(v2.1: annotation, not a veto — see §8)* held-out-epoch
      prediction test: refit (z, µ) on epochs ≤ MJD 59579 (2021-12-31),
      forced photometry at that node on the 2022–24 epochs; pass iff
      S_late ≥ 3 and f_late/f_early ≥ 0.3; its false-pass rate on null
      trajectories and pass rate on injections per temporal model are
      quoted with every candidate; applicable only with ≥ 5 late epochs;
   c. *W3/W4 confirmation procedure* (plan §5.3): W1/W2 forced
      photometry at the same (z, µ) on the cryo frames against the
      300 K flux-ratio prediction (W2/W3, W1/W2, 3σ), post-cryo W1/W2
      persistence along the track, catalogue counterpart at the peak;
      external imaging only if all three are inconclusive.
   **Annotations** (reported, never grounds for rejection): parallax-
   phase split and ratio, single-epoch dominance, W1:W2 significance
   ratio, cryo visit count, bare catalogue proximity, grid-edge fit.
   A candidate that survives is `retained-ambiguous`.
8. **Quality mask** (plan §7). Primary: `qual_frame > 0`, `qual_scan ≥ 5`
   where present, `saa_sep > 0`. Strict: additionally `qual_frame = 10`
   and `moon_sep ≥ 30°`. Loose: `qual_frame ≠ 0` (v1). Primary is the
   reported result; strict/loose run on the development set and their
   spread is quoted as a systematic.

## 4. Completeness (plan §4)

Injections are image-level: the band's empirical PRF (IRSA second-pass
9 × 9 focal-plane templates, element chosen by the frame position,
sub-pixel phase uniform) is added in DN to each `-int` cutout before
background estimation, masking and the matched filter, exactly via the
stamp-response linearity of the estimator (invariant-tested). Per cell
**400 injections**, 100 per temporal model; z continuous from the
log-uniform prior, (µ*α, µ*δ) uniform on the L∞ box, cross-track offset
from the §1.4 envelope, Vega magnitude uniform on [m90_v1 − 2, m90_v1
+ 2] where m90_v1 is the cell's v1 median depth (a fixed mid-scale
window where v1 found none). Spectrum: flat Fν for W1/W2, 300 K
blackbody for W3/W4, with the Explanatory Supplement colour
corrections; records carry F_ν at λ_iso and the Vega magnitude.

Two curves, always separate: **threshold completeness** (the injected
cell's R̃ ≥ R̃_FWER with the maximum within 2 z-nodes and 1 µ-node of the
injection) and **final-candidate completeness** (and it survives every
calibrated veto of §3.7 run blind on it). Each as a function of
magnitude, in 8 z-intervals whose edges are reciprocal-distance
midpoints tiling 550–10,000 AU with no gap; m90 and m50 per interval
and temporal model from a logistic fit with bootstrap 68 % intervals
and the per-interval injection count; worst temporal model defines the
cell's coverage at duty ≥ 0.5. W3/W4 constraints carry
`completeness_kind = threshold` and `exclusion_claim = false` until
§3.7c has been exercised on a real exceedance and on injections.

## 5. Hold-out (plan §1.7–1.8)

A random **corridor** split stratified by the v1 overlay confusion
class (low / mid / high / mid+bright-star / high+bright-star), seed
20260822, ≈ 30 % of corridors (every class represented) to the
**development set** — all rule tuning, exchangeability, veto selection
functions, mask sensitivity, and the true hold-out power test of §3.7b
(2022–24 excluded from the search) — and the remaining ≈ 70 % to the
**confirmatory set**, run once, blind, under this rule. Both roles of an
endpoint and all components of a system share a corridor and therefore
a side. The lists and their hash are in `configs/v2_0_freeze.json`.
If the rule changes after the confirmatory set is touched, the
confirmatory claim is void and the report must say so. All epochs are
searched (NEOWISE ended 2024-08; a pre-registered epoch hold-out would
cost ≈ 25 % of depth unconditionally); the 2022–24 prediction test is
the calibrated post-hoc test of §3.7b. Cross-archive confirmation (ZTF /
PS1 / SPHEREx at the predicted position) remains the primary
independent evidence for any retained candidate.

## 6. Scope statement

Targeted coverage of the frozen 88-endpoint portfolio under the v1
physics cell. No population or Picky-Network inference is made; that
requires a separately pre-registered generative model.

## 7. Freeze gate

`configs/v2_0_freeze.json` records this file's hash, the registry
hash, the split, the seed and every numerical parameter above; its own
hash is pinned in every v2 AnalysisRun. No script in plan §3–§8 runs
against the confirmatory set before that file is committed.

## 8. v2.0 → v2.1 (step F, development set only, 2026-08-22)

The development-set selection functions showed that the held-out-epoch
prediction test passes 59 % of null trajectories above their cell's
q95 (43 % of all null trajectories) while rejecting 74–76 % of in-scope
long-block sources (W1/W2; they are simply off during 2022–24) and
8–11 % of visit-scale sources. By §6 of the plan a veto may not reject
a present, detectable in-scope signal, so the test does not qualify as
a rejection test. **v2.1 demotes it to an annotation** carried on every
candidate with its measured rates; the calibrated rejections are the
flux-consistent static/halo test (3.7a) and the W3/W4 confirmation
procedure (3.7c). Final-candidate completeness is therefore threshold
completeness minus the static-test loss. No other parameter changed;
the hold-out split is unchanged and the confirmatory set has not been
touched. The freeze file is re-hashed (`configs/v2_0_freeze.json`,
version field `wise-hypotheses-v2.1`).
