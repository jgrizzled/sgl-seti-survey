---
title: "Joint Pan-STARRS1 + ZTF stage — Report v2 (one family-wise rule across archives)"
date: 2026-08-23
status: "complete — confirmatory 45 endpoints / 230 joint cells: 0 candidates at family-wise α = 0.05; development 24 / 102: 0 candidates; 5,280 joint constraints; supersedes joint_ps1_ztf_v1.md (status note 2026-08-22)"
---

# Joint PS1 + ZTF stage — Report v2

**Pins:** PS1 v2.0 freeze `sha256:be6f3a1c…` + ZTF v2.0 freeze
`sha256:9daed5f2…` (identical corridor splits by construction) · joint
freeze `surveys/joint-v2/configs/v2_freeze.json` `sha256:9bb5791f…` ·
common T0 = 59800, AB zero point 25 · `surveys/joint-v2/joint.py` ·
AnalysisRuns in `runs/joint-v2/records/` · numbers from
`surveys/joint-v2/results/report_tables.md`.

## Summary

The joint stage combines the two archives' v2 accumulated stack sums
(per trajectory, per mask, capped and clipped per archive; the ZTF
192-node sums interpolated in 1/z onto the PS1 360-node grid) into
S_joint = (A_ztf + A_ps1)/√(B_ztf + B_ps1) per joint cell (endpoint ×
role × band pair g/r/i). Because the two archives share the identical
48-offset ring, the joint ring null is exact and **one family-wise
threshold is calibrated across archives**. The j-th injection of a cell
is the same physical source in both archives (per-injection seeded
draws on the union magnitude window), so joint completeness is the
recovery of one source through both chains. **No joint cell of the
confirmatory set (230 cells) reaches R̃_FWER = 1.778**; 14 cells have
R > 1 against a ring expectation of 25. Development set: 102 cells,
R̃_FWER 1.54, 0 candidates. Persistent-source m90 g 23.1 / r 22.9 /
i 21.1 AB (worst-of-four ≈ 22.2 / 22.1 / 21.0): the joint stage adds
little depth to ZTF (PS1 contributes ~20 epochs against ~1,000) — its
value is the second parallax phase for PS1's single-phase cadence and
the cross-archive consistency it now provides for free on every cell.

## Rule (frozen)

Joint R = S_max/T over the 8 designated controls, R̃ = R/q95(ring),
heavy-tail / inner-outer-ring void flags, R̃_FWER from 10,000
pseudo-experiments over the set's joint cells; the only veto is the
flux-consistent catalogued static source with the PS1 DR2 and ZTF DR
predictions summed against the joint peak (an earlier version of the
stage OR-ed the per-archive flags, which rejected injections whose
ZTF share alone was explained; fixed on the development set before the
confirmatory completeness was computed — the confirmatory null result
was unaffected). Phase split, p_phase and p_trajectory are annotations.

## Results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | candidates |
| --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 45 / 230 | 5 | 1.778 | 14 (25.2) | **0** |
| development | 24 / 102 | 1 | 1.540 | — | 0 |

## Completeness (confirmatory, AB)

| band pair | persistent | flicker | visit | block |
| --- | --- | --- | --- | --- |
| g (PS1 g + ZTF zg) | 23.06 / 22.83 | 22.36 / 22.22 | 22.16 / 22.09 | 22.22 / 22.10 |
| r | 22.93 / 22.86 | 22.19 / 22.15 | 22.05 / 21.98 | 22.09 / 22.06 |
| i | 21.07 / 21.09 | 21.17 | 21.02 | 20.97 |

(threshold / final-candidate; the static veto fired on 2.9 % of
threshold-recovered joint injections). Injections brighter than the
fainter of the two archives' single-epoch clip limits are excluded
from the fits (layered-search rule). Constraints: 5,280, with
supersession links to the v1 joint records.

## Provenance note

The confirmatory null-ensemble file was re-serialised once after the
blind run to fix a field format in the report generator; the
computation is deterministic (seeded) and every number (230 cells,
5 void, R̃_FWER 1.778, 0 candidates) is identical to the first run.
