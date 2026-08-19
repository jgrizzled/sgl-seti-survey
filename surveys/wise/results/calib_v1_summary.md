---
title: "calib_v1 — injection calibration, thresholds, and first constraints (plan §4.7)"
date: 2026-08-18
run_dir: "runs/wise/calib_v1 (gitignored; regenerate with scripts/sample_tensor.py + injection_calibrate.py)"
analysis_run: "run-88887a3ac8d9"
---

# calib_v1 results

Injection calibration over the full pilot: per-epoch sample tensors for
the real trajectory plus **8 offset-control trajectories** per
endpoint × role (14 tensors, 713–916 epochs each, 1.3 GB), searched on
an improved relay-distance grid, with analytic injection-recovery and
the pilot's first ledger-ready Constraint records.

## Grid fix (supersedes stack_v1 geometry)

The stack_v1 log-spaced z-grid left ~17" gaps between adjacent nodes at
z≈550 AU (locus position ∝ 1/z) — nearly 3 PSF widths, so low-z sources
between nodes could be missed. calib_v1 uses a grid **uniform in 1/z**:
constant ~5.6" node spacing (just under the W1 PSF FWHM) across the
whole 550–10,000 AU range with the same 64 nodes.

## Predeclared thresholds and negative controls

Per endpoint × role × band: T = max of the 8 control-trajectory grid
maxima (empirical FAR < 1/8 per grid search by construction). Controls
in crowded corridors reach T ≈ 75–230+ — the search is field-star
contamination limited, not noise limited, and the thresholds inherit
that honestly.

**Adjudication of the two stack_v1 marginal cells (Candidate records,
both `vetoed`):**

- lalande-21185/rx W1: S=30.4 vs T=76.6
- lalande-21185/tx W2: S=17.0 vs T=74.1

Both sit far inside the 8-control null on the finer grid. **The pilot
search ends with no surviving candidate.**

## Injection recovery → constraints

Injection model: Gaussian-PSF point sources at grid nodes, exploiting
matched-filter linearity (injected flux adds f·g_e per epoch), with
per-band grid- and µ-mismatch factors applied; duty cycle realizations
randomized (32 repeats, seed 20260818); recovery criterion S ≥ T at the
true cell; 90%-recovery flux = 90th percentile of per-repeat minimum
recoverable fluxes. Depths in WISE Vega mag via per-frame magzp.

**448 Constraint records** (14 endpoint-roles × 4 bands × 8
z-intervals; every one `recovery_curve`, all 64 nodes valid), each:
duty ∈ [0.5, 1] (bounded by the duty-0.5 curve), |µ_resid| ≤ 1 "/yr,
90% recovery, band-specific physical interpretation, traceable to
AnalysisRun `run-88887a3ac8d9` → tensors → usable precise-pass pixels.

Median 90%-recovery depths (duty ≥ 0.5), by corridor character:

| Endpoint/role (example) | W1 | W2 | W3 | W4 |
| --- | --- | --- | --- | --- |
| lalande-21185 rx (sparse field) | 16.9 | ~15 | ~12.6 | ~10.4 |
| barnard-star rx | 15.6 | ~14.5 | ~12.5 | ~10.4 |
| ross-154 rx | 15.5 | 14.4 | 12.8 | 10.5 |
| alpha-cen-a rx (Galactic plane) | 15.0 | ~14 | ~11 | ~10 |
| sirius-a rx (bright-star field) | 12.9 | 13.5 | 10.6 | 10.2 |

(Full 64-node curves per pair × band × duty in `m90_curves.npz`;
per-interval values in the Constraint records.)

Reading: stacked depths land near or somewhat beyond single-exposure
depth rather than the ideal √N ≈ 3.5-mag gain, because the max-over-grid
thresholds are set by field-star contamination of control trajectories,
not by pixel noise. Sirius A's shallow W1 reflects its bright-star-
artifact corridor. W3/W4 depths rest on the single 2010 cryo visit.

## Plan §4 success criteria status

- Reproducible joins from frozen metadata and manifests — **yes**
  (snapshots + content-addressed records at every stage).
- Injected sources achieve a declared recovery probability over
  nontrivial flux×distance regions — **yes** (90% over all 448 cells).
- Negative controls within predeclared thresholds — **yes** (threshold
  IS the declared control statistic; both marginal cells vetoed).
  Positive controls: synthetic injections only — recovery of known
  moving objects (e.g. via the snapshotted SSO association lists) is an
  open follow-up, noted as a gap.
- Constraints traceable to usable pixels and an analysis run — **yes**
  (record chain, never a bare footprint).
- At least one ledger-ready constraint — **yes, 448 across all 7
  endpoints** (or an explicit not-constrainable statement — none
  needed).

## Caveats (carried into the report)

- Gaussian-PSF injection measures pipeline throughput, not PSF-shape
  mismatch against the true WISE PSF (wings/spikes); an empirical-PSF
  upgrade would tighten depths modestly.
- FAR resolution is 1/8 per grid search (8 controls); a production run
  wants O(100) sky-rotated/time-scrambled corridors.
- Depths assume the frozen morphology/persistence cell (point source,
  duty ≥ 0.5, |µ| ≤ 1 "/yr); other grid cells are separate hypotheses.
- Background estimation is cutout-global; crowding (α Cen) biases it
  slightly — shared by injections, so recovery curves absorb it to
  first order.

## Next (plan §4.8)

Report: assemble observation/intersection/analysis/constraint/candidate
records, coverage by epoch and z-interval, completeness curves, and the
physical interpretation of the flux limits into the per-target,
per-archive shakedown report.
