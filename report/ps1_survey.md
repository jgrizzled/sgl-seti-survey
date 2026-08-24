---
title: "Pan-STARRS1 survey (frozen decision rule, blind confirmatory run)"
date: 2026-08-23
status: "complete — confirmatory 45 endpoints / 450 cells: 0 candidates at family-wise α = 0.05; development 24 / 240: 0 candidates; 11,040 injection-calibrated constraints; supersedes the v1 report (removed in the v1 retirement, 2026-08-24; git history)"
---

# Pan-STARRS1 survey

**Pins:** registry v1.5 · hypotheses `surveys/panstarrs/hypotheses.md`
v2.0 (freeze `sha256:be6f3a1c…`) · common T0 = 59800 · engine
`sglsurvey/` · AnalysisRuns in `runs/panstarrs/v2/records/` · numbers
from `surveys/panstarrs/results/report_tables.md`.

## Summary

The v1 PS1 survey was re-analysed under the WISE v2.1 design with the
purged warp cutouts re-fetched per batch, the v1 per-warp star zero
points, skycell duplicates collapsed, Moffat injections at each warp's
seeing, and the PS1 DR2 catalogued-static veto. **No cell of the
confirmatory set (450 cells) reaches the family-wise threshold**
R̃_FWER = 1.470; 49 cells have R > 1 against a ring expectation of 55.
Development set: 240 cells, R̃_FWER 1.504, 0 candidates. PS1's local
nulls are the best behaved of the four archives (no heavy-tail cell;
8 void cells from the inner/outer-ring test), Where the stack can constrain at all
(see §Completeness), persistent-source m90 is g 21.7 / r 21.6 / i 21.0 /
z 20.5 / y 19.4 AB — but the calibrated analysis shows that in most
cells PS1's ~20-warp stack has no 90 %-complete regime between the
single-epoch clip and the contamination threshold (v1's uniform
"~21 AB" claim did not account for either).

## Rule (frozen)

As WISE v2.1 with PS1 bindings: 360 × 5 × 5 grid (1″), clip |S_e| ≤ 5,
per-frame cap, ring 20/30/40″, R̃ = R/q95, void flags, R̃_FWER over the
family; veto = flux-consistent PS1 DR2 mean source (≥ 3 detections;
sentinel magnitudes filtered). No epoch hold-out (mission over); the
3π single-phase cadence leaves the phase split an uninformative
annotation — the joint stage supplies the second phase.

## Results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | candidates |
| --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 45 / 450 | 8 | 1.470 | 49 (54.9) | **0** |
| development | 24 / 240 | 4 | 1.504 | 27 (29.4) | 0 |

Mask sensitivity (development): R̃_FWER 1.50 / 1.54 / 1.50, 0 candidates
under each.

## Completeness (confirmatory, AB)

| band | persistent | flicker | visit | block |
| --- | --- | --- | --- | --- |
| g | 21.71 (107) | 21.53 (9) | 21.73 (36) | 21.49 (24) |
| r | 21.56 (64) | 21.14 (15) | 21.50 (25) | 21.41 (4) |
| i | 20.99 (71) | 21.01 (8) | 21.01 (17) | 20.97 (30) |
| z | 20.46 (60) | 20.38 (32) | 20.52 (44) | 20.03 (10) |
| y | 19.38 (104) | 19.58 (3) | 19.29 (42) | 19.35 (12) |

(numbers in parentheses: cell × interval count with a valid fit)

**The structural result of the calibrated re-analysis:** with ~20
warps per band, a source bright enough for a single-warp 5σ detection
is clipped from the stack by the layered-search rule, so the stack's
completeness lives only in the narrow window between the single-epoch
limit and the contamination-driven family threshold. A fit whose 90 %
point would land brighter than the clip is an extrapolation into the
catalogue layer's regime and is voided: only **1,336 of 7,200**
confirmatory constraint records are recovery curves (the surviving,
cleaner cells quoted above); everywhere else the PS1 *stack* adds no
90 %-complete depth of its own and the constraint belongs to the
catalogue-screening layer and the joint stage. The v1 claim of
uniform ~21 AB stack depth is thereby localised to the cells that can
actually support it: with 100 injections per
model and the bright cut, the logistic fit often lacks support. Static
veto fired on 0.2 % of threshold-recovered injections; Moffat throughput
0.73.

## Geometry

Cross-track dimension needed for 18 of 138 cells (threshold 0.5″ at
the 1″ PSF; max σ_xt,99 11.5″), as the WISE retrospective predicted for
optical PSFs. No new positive control was run (the v1 (60000) control
at F51 stands; a v2 re-score through the stack regime is left for the
next cycle).

## Scope

Targeted coverage of the 69-endpoint PS1 subset; no population
inference.
