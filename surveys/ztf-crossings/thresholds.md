# ZTF crossings threshold freeze v1.0

Frozen 2026-08-23 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json`
(`sha256:45e51b46…76fa8de`), bound to content hashes of
`hypotheses.md`, `results/coverage_v1_events.ecsv`, and
`results/saturation_cut_v1.ecsv`. No pixel data touched. Conventions
inherited from the v2 engine: 8 designated controls, exceedance ratio
R = S/T with T = max control statistic (per-search crossing probability
1/9), WEIGHT_CAP 20× effective-epoch floor, split seed 20260822,
dev fraction 0.30, fwer α = 0.05, KS α = 0.01.

## Statistics

**Channel A** (blended star, searchable bands only): PSF forced
photometry on ZTF *difference* images at the star position — the star
lives in the reference, so difference flux is window excess by
construction. Per-window weighted excess S_w; primary statistic
S_A = stack of S_w over covered windows, one-sided positive, requiring
≥ 2 covered windows. Single-window units are constraint-only. The
≤ 10″ station annulus uses the channel-B construction minus the track.

**Channel B** (antipode track): PSF forced photometry on difference
images along the per-(event, z) apparent relay track; the z grid is a
*nested family* sharing epochs, so the unit statistic is the max over z
and controls take the same max — one trial per (event, radius, band).
Units with a single in-window epoch form the reportable "single-epoch
exceedance" class; promotion to candidate requires recurrence at a
second covered window.

## Controls and rule

- A: 8 temporal pseudo-window sets at frozen offsets ±23/47/71/97 d
  (avoiding lunation multiples and the semiannual crossing spacing;
  overlaps with any real window re-draw to ±113/127 d).
- B: 8 spatial ring trajectories, radii 20/30/40″, identical designated
  construction as ztf-v2, offsets following the track.
- Exceedance: R = S/T > 1. Expected control crossings across all 77
  units: **8.6**. Every exceedance is individually adjudicated by the
  frozen veto ladder (census, chromatic, rate, recurrence) — no silent
  drops.
- Quality: primary mask for the search; strict-mask re-run reported for
  any exceedance; any ZTF reference containing an in-window epoch is
  flagged and its node masked.

## Search units (from the frozen coverage × saturation intersection)

| channel | units | constraint-only |
|---|---|---|
| A (target × band × radius, ≥2 windows) | 29 | 9 single-window |
| B (event × radius × band, ≥2 epochs) | 48 | 23 single-epoch |

## Dev / confirmatory split (unit = target, seed 20260822)

- **A dev (6):** ez-aqr, gj-1087, gj-1111, gj-11068, gj-518, gj-9193 —
  stratified by band-searchability class; remaining 21 confirmatory.
- **B wide-rung dev (1):** ross-154. gj-1276 and teegarden are
  *excluded from dev eligibility*: their grazing-rung windows are
  temporal subsets of their wide-rung windows, and the grazing rung
  (5 covered events, the survey's headline hypothesis) is fully
  confirmatory with procedures locked from A + B-wide dev experience.

## Completeness

Injections per search unit (100 per event-window cell) on the
difference-image substrate, chord temporal profile from the unit's
(b, v⊥, radius), line spectra 532 nm (g) / 650 nm (r), v2 recovery
window [2,1].

## Amendment v1.1 (2026-08-23, channel A only)

Adopted after the dev findings in `results/dev_search_v1.md`
(`configs/threshold_freeze_v1_1.json`). Channel B unchanged.

1. **Systematics template (finding A2).** Per (target, band), the
   per-epoch forced flux is modeled as
   c0 + c1·Δt + cx·P_α(t) + cy·P_δ(t) + cs·(seeing − median), where the
   parallax factors P are fully fixed by registry astrometry (no
   free-phase annual term); only linear couplings are fitted, on a
   deterministic ≤150-epoch uniform-in-time off-window subsample with
   3×3σ clipping. The window statistic is formed on residuals. Quality
   gate: ≥30 surviving fit epochs and |median control statistic| ≤ 1,
   else constraint-only.
2. **Pseudo-window exclusion (finding A1)** relaxed to same-rung-only;
   wide-beam leakage into narrow-rung controls is conservative (inflates
   T). Units still lacking 8 valid offsets are constraint-only. Rung
   verdicts per target are declared correlated.
3. **Exceedance rule:** S > max(T, 0), margin S − T reported; the R
   ratio is retired (unstable near T → 0⁺, sign-invalid for T ≤ 0).

## Amendment v1.2 (2026-08-23, channel A only)

`configs/threshold_freeze_v1_2.json`. Adds empirical variance rescaling
to v1.1: per (target, band), per-epoch variances are scaled by
k = median(r²/v)/0.4549 over the clipped off-window fit epochs
(robust χ²/dof; floored at 1), standardizing S by the empirical
residual scatter — the v1 "confusion, not noise, sets the floor"
precedent. The exceedance comparison is scale-invariant; the change
makes the |median control| ≤ 1 gate meaningful. Everything else in
v1.1/v1.0 unchanged.
