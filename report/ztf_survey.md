---
title: "ZTF survey (frozen decision rule, blind confirmatory run)"
date: 2026-08-23
status: "complete — confirmatory 45 endpoints / 228 cells: 0 candidates at family-wise α = 0.05; development 24 / 102: 0 candidates; 5,312 injection-calibrated constraints; supersedes the v1 report (removed in the v1 retirement, 2026-08-24; git history)"
---

# ZTF survey

**Pins:** registry v1.5 · hypotheses `surveys/ztf/hypotheses.md` v2.0
(freeze `configs/v2_freeze.json` `sha256:9daed5f2…`) · physics cell as
v1.0 · engine `sglsurvey/v2/` · AnalysisRuns `run-ab66811d2d93`
(confirmatory) / `run-d8f0905483fb` (development) · all numbers from
`surveys/ztf/results/report_tables.md` (ledger-generated).

## Summary

The v1 ZTF survey (62 corridors, "FAR < 1/8", 2,640 threshold
constraints) was re-analysed under the WISE v2.1 design: a 48-offset
exchangeable null per cell with a family-wise error rate, image-level
Moffat injections into the difference/science hybrid search image with
four temporal models, the flux-consistent catalogued-static veto (with
the difference-image fraction taken into account), covariance
propagation with a cross-track dimension where needed, and a blind
confirmatory hold-out. **No cell of the confirmatory set (228 cells)
reaches the family-wise threshold** R̃_FWER = 1.658; 19 cells have
R > 1 against a ring expectation of 26 — the v1 exceedances were the
threshold's rank statistics. The development set (102 cells, R̃_FWER
1.872) likewise has 0 candidates. Median 90 %-completeness for a
persistent source is zg 23.2 / zr 23.1 / zi 20.9 AB, and for the worst
of the four temporal models (duty ≥ 0.5 coverage) zg 22.3 / zr 22.1 /
zi 19.4 — ZTF's ~1,000-epoch cadence makes the family-wise price small
(R̃_FWER 1.66 vs WISE's 2.44, with only 4 void cells).

## Rule (frozen)

As WISE v2.1 (`surveys/wise/hypotheses.md` §3–§5, §8) with ZTF
bindings: 192 × 5 × 5 grid (T0 = 59800), single-epoch clip |S_e| ≤ 5,
per-frame weight cap, ring 20/30/40″, R̃ = R / q95(ring), heavy-tail and
inner/outer-ring void flags, R̃_FWER from 10,000 pseudo-experiments over
the family; the only veto is the flux-consistent catalogued static
source (ZTF DR objects, ≥ 3 good observations) with the predicted flux
scaled by the science-image fraction at each epoch (static sources are
absent from the difference-image regime — without that scaling the
veto "explained" injected sources, a defect found and fixed on the
development set before the confirmatory completeness was computed;
the confirmatory null result was unaffected). The held-out-epoch test
(split 2024-09-01) is an annotation: it passes 12 % of null
trajectories (26 % of those above q95) and 98 / 97 / 92 / 27 % of
persistent / flicker / visit / block injections.

## Results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | candidates |
| --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 45 / 228 | 4 | 1.658 | 19 (26.3) | **0** |
| development | 24 / 102 | 1 | 1.872 | 11 (12.0) | 0 |

Mask sensitivity (development, 98 common cells): R̃_FWER 1.87 / 1.70 /
1.92 (primary / strict / loose), 0 candidates under each, ΔR median
−0.01.

## Completeness (confirmatory set, AB, median over cells × intervals)

| band | persistent | flicker | visit | block | worst |
| --- | --- | --- | --- | --- | --- |
| zg | 23.21 | 22.44 | 22.35 | 22.43 | 22.4 |
| zr | 23.08 | 22.29 | 22.14 | 22.21 | 22.1 |
| zi | 20.96 | 20.34 | 20.26 | 20.32 | 20.3 |

Final-candidate curves are within 0.05 mag of threshold curves (static
veto fired on 1.7 % of threshold-recovered injections). The
Gaussian-kernel throughput on the Moffat PSF is 0.675 (0.43 mag),
consistent with the v1 asteroid control's 0.41–0.59 mag. Constraints:
5,312 (confirmatory 3,680, of which 3,356 are recovery curves; the
rest `not_constrainable` — void null, or a fit whose m90 would sit
brighter than the cell's single-epoch clip limit, i.e. in the
catalogue layer's regime). Sources brighter than
the cell's single-epoch 5σ limit are clipped from the stack by the
layered-search rule and belong to the catalogue layer; each constraint
carries that bright limit.

## Geometry and control

Cross-track dimension needed for 16 of 138 cells (threshold 0.95″,
max σ_xt,99 17.5″ for the ε Ind Ba/Bb photocentre), searched and
null-calibrated with 3–5 offsets. Positive control: asteroid (60000)
re-scored under the v2 rule — R̃ = 3.9 (zg) / 3.4 (zr), rank 1/48, at
the Horizons position; because it is individually detected in single
frames (S/N ≈ 9) the frozen clip removes those frames and the stack
recovers only its sub-clip share (1.6 mag faint); without the clip the
recovered magnitude is 0.35 / 0.53 mag fainter than predicted (zg / zr),
the Moffat-vs-Gaussian throughput. A fainter mover in the stack regime
is the right control and is left for the next cycle.

## Scope

Targeted coverage of the 69-endpoint ZTF-visible subset under the
v1.0 physics cell; no population inference.
