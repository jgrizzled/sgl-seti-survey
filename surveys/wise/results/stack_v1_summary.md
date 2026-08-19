---
title: "stack_v1 — trajectory-aware forced-photometry stack (plan §4.6)"
date: 2026-08-18
run_dir: "runs/wise/stack_v1 (gitignored; regenerate with scripts/fetch_cutouts.py + forced_stack.py)"
analysis_run: "run-2f0f4c6fd24a"
---

# stack_v1 results

Stage 2 of the search: matched-filter forced photometry on every usable
frame, stacked along predicted SGL trajectories over a joint
(relay distance × residual motion) grid — the plan §4.6 "fit relay
distance and bounded residual motion jointly", implemented as
shift-and-stack over 64 log-z × 5×5 µ (|µ| ≤ 1 "/yr frozen bound)
per endpoint × role × band.

## Inputs and method

- 4,280 int+unc cutout pairs (720 MB, 0 fetch errors) sized to each
  frame's corridor arc; masks reused from precise_v1 (600 MB).
- Per cutout: Gaussian-PSF matched-filter flux + variance maps
  (FFT cross-correlation, fatal-bit masked pixels excluded, sigma-
  clipped-median background), fluxes rescaled to a common zero point
  (magzp → 20.0). Declared pilot approximations: Gaussian PSF,
  constant background — they cost sensitivity, not validity.
- Stacks: inverse-variance-weighted sums over all usable epochs, with
  (a) parallax-phase-split accumulators (the static-star veto validated
  in screen_v1) and (b) a +35" offset-trajectory negative control
  through identical machinery.
- AnalysisRun `run-2f0f4c6fd24a` (config + input-set hashes pinned);
  outputs `stacks.npz`, `epoch_fluxes.npz` (14 MB).

## Outcome across all 56 endpoint × role × band stacks

| Verdict | Count |
| --- | --- |
| max S ≤ offset-control max (field-star contamination dominates) | 36 |
| exceeds control but single-phase data only (cryo-era W3/W4 — one visit, no recurrence power) | 12 |
| exceeds control, two-phase data, fails the phase test (single-phase ⇒ static field star) | 6 |
| marginal cells carried to stage-7 adjudication | 2 |

The W1/W2 stacks are decisively confusion-dominated: control maxima run
up to S≈830 (bright field stars anywhere along a trial trajectory), and
every strong real-trajectory peak is single-phase — the static-star
signature — or below its control. **No detection.**

The two flagged marginal cells (both Lalande 21185 corridor):

- `rx/W1`: S=30.1 at z≈2630 AU, µ=(+0.5,−1.0), vs control 27.6 (+9%);
  phase split (16.5, 26.5). W2 at the same cell = 8.0 — a W1:W2
  significance ratio ~3.8:1, typical of an ordinary faint star.
- `tx/W2`: S=15.2 at z≈1660 AU vs control 12.9 (+18%); only ~7
  effective epochs carry the weight.

Both are consistent with placement variance against a **single** control
draw in non-Gaussian statistics with ~90,000 correlated cells searched;
neither meets any predeclared threshold (none exists yet — that is
step 7's job). They are recorded here for adjudication by the empirical
false-alarm machinery, not as candidates (plan §3.2: a Candidate
requires competing-model tests against a declared threshold).

## Interpretation caveats

- Per-pixel `unc` excludes confusion noise, so S values are inflated
  relative to true significance; only the control/empirical comparisons
  are meaningful. Injection calibration (§4.7) will measure the real
  throughput and false-alarm scale.
- One offset control per stack is a 1-sample null; §4.7 adds matched
  control corridors (sky-rotated / time-scrambled) for real false-alarm
  rates.
- Cryo W3/W4 stacks are effectively single-visit forced photometry —
  no recurrence or phase power; their excesses over one control draw
  carry little meaning. Waste-heat constraints from W3/W4 will be
  sensitivity statements, not recurrence tests.
- No flux limits are quoted yet: converting stack depth into
  completeness-qualified constraints requires the injection stage.

## Next (plan §4.7)

Injection calibration: blind synthetic point sources injected into the
real cutouts across (flux × z × µ × duty cycle × background), recovered
through this exact pipeline, plus a battery of matched control
corridors → empirical false-alarm rates, recovery curves, and the first
ledger-ready Constraint records — including adjudication of the two
flagged Lalande cells.
