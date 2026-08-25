---
title: "Joint Pipeline A stage v3 — PS1 + ZTF stack with the WISE colour axis"
date: 2026-08-25
status: "complete — confirmatory 45 endpoints / 230 joint cells (blind, once): 0 candidates at family-wise α = 0.05; optical statistic reproduces the v2 run bit-identically; colour axis: 230/230 cells annotated, 0 would-fire vetoes, measured false-veto rate 0 in the fitted population; supersedes the v2 report (2026-08-23, git history)"
---

# Joint Pipeline A stage (PS1 + ZTF + WISE)

**Pins:** PS1 v2.0 freeze `sha256:be6f3a1c…` + ZTF v2.0 freeze
`sha256:9daed5f2…` + WISE v3.0-joint-conventions freeze
`sha256:14f681a3…` (identical corridor splits by construction) · joint
v3.0 freeze `surveys/joint/configs/v3_freeze.json` `sha256:bb68869c…`
(supersedes v2 `sha256:9bb5791f…`) · common T0 = 59800, AB zero point
25 across all three archives · `surveys/joint/joint.py` · AnalysisRuns
in `runs/joint/v3/records/` · numbers from
`surveys/joint/results/report_tables.md`.

## Summary

v3 extends the two-archive joint stage with the WISE W1/W2 colour
axis (plan §4 open item; learnings §10 item 3). The candidate
statistic is **unchanged from v2**: joint cell = endpoint × role ×
band pair (g/r/i), S_joint = (A_ztf + A_ps1)/√(B_ztf + B_ps1) on the
PS1 grid, exact shared 48-offset ring null, one family-wise threshold
across archives. WISE enters only through (i) per-cell W1/W2 stack
annotations at every joint peak and (ii) a **calibrated
colour-consistency veto**: under the frozen flat-Fν spectrum on the
common AB ZP-25 flux scale, a candidate's WISE stack flux must not
fall ≥ 5σ below the joint optical flux estimate in every usable W
band, with σ_w floored by the empirical 48-ring scatter (the WISE
confusion floor). The WISE tensors were rebuilt for this under the
joint conventions (T0 59800; Vega→AB W1 +2.699 / W2 +3.339; 9-node µ
grid whose [::2] subgrid is the optical grid; `runs/wise/v3/`).

**The blind confirmatory set (45 endpoints / 230 cells) returned 0
candidates** at R̃_FWER = 1.778 — and every optical number (cells,
voids, threshold, 14 R > 1 vs 25.2 expected, all 3,648 constraints
and completeness medians) is **bit-identical to the v2 run**, as the
freeze declared it must be. The confirmatory content new in v3 is the
colour axis: all 230 cells have both W bands usable (n ≥ 5), the
flat-Fν deficit statistic at the real peaks spans median −0.33 to max
3.73 (no cell approaches the ν = 5 veto line), and the top optical
cells now carry two-band colour consistency for free.

## Colour veto calibration (measured selection function)

The same per-injection seeded draws were pushed through the WISE chain
(one run per optical band family with that family's exact PS1/ZTF
union magnitude window, `runs/wise/v3/injections_{g,r,i}` — verified
draw-identical to the frozen optical injection products). On
threshold-recovered injections the veto fired **5 / 46,302
(confirmatory)** and **1 / 23,117 (development)** — and every firing
was an injection at mag 16.8–18.8, brighter than its cell's optical
single-epoch bright limit (21.0–21.9), i.e. in the catalogue-layer
regime already excluded from the fitted constraints. **In the fitted
population the measured false-veto rate is 0/46,302 (confirmatory) and
0/23,116 (development)** — the veto's completeness cost on every
quoted constraint is zero. Mechanism of the bright firings: the
uncorrected Gaussian-on-PRF throughput (W1 0.80 / W2 0.75) plus, for
block/visit models, the different overlap of the WISE 2010–2024 epochs
with the injected on-window — real effects confined to the
above-clip regime.

## Results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | colour usable / would-fire | candidates |
| --- | --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 45 / 230 | 5 | 1.778 | 14 (25.2) | 230 / 0 | **0** |
| development | 24 / 102 | 1 | 1.540 | — | 102 / 0 | 0 |

## Completeness (confirmatory, AB; unchanged from v2)

| band pair | persistent | flicker | visit | block |
| --- | --- | --- | --- | --- |
| g (PS1 g + ZTF zg) | 23.06 / 22.83 | 22.36 / 22.22 | 22.16 / 22.09 | 22.22 / 22.10 |
| r | 22.93 / 22.86 | 22.19 / 22.15 | 22.05 / 21.98 | 22.09 / 22.06 |
| i | 21.07 / 21.09 | 21.17 | 21.02 | 20.97 |

(threshold / final-candidate; static veto on 2.9 % of
threshold-recovered joint injections, colour veto on none in the
fitted region; injections brighter than the fainter of the two
archives' single-epoch clip limits excluded from the fits.)
Constraints: 3,648 confirmatory + 1,632 development under the v3
hypothesis version, each carrying the colour-veto census in `extra`.

## Stack-regime positive control (learnings §10 item 1)

Asteroid (220000) on nights with predicted V 21.9–23.05 — below the
single-frame limits, so the frozen |S_e| ≤ 5 clip removes (almost)
nothing and only the ephemeris-weighted stack can recover it (the
Horizons track playing the SGL trajectory's role):

- **ZTF** (35 frames, all sub-clip, max single-frame S/N 4.9):
  zr S = 9.6, R̃ = 2.90; zg S = 5.3, R̃ = 1.95; recovered − predicted
  = +0.23 / +0.29 mag (Horizons V + solar colours).
- **PS1** (25 warps; its one 5.6σ frame removed by the frozen clip):
  i R̃ = 1.77 at −0.14 mag; z an honest non-detection (median
  single-frame S/N 0.34, 8 epochs); g/r/y too few epochs standalone.
- **Joint g/r/i families through the exact v3 statistic:
  R̃ = 2.21 / 3.04 / 1.77** — every family above the development
  threshold (1.54): the family rule promotes the control, as a
  positive control must.

Scripts `surveys/{ztf,panstarrs}/scripts/asteroid_stack_fetch.py`,
`surveys/joint/scripts/asteroid_control_stack.py`; summaries under
`runs/{ztf/v2,panstarrs/v2,joint/v3}/control/220000_stack/`.

## Scope and provenance notes

- The colour veto rejects only flat-Fν-inconsistent *deficits*; excess
  W flux (red blends) can never fire it. Sources with steeply red or
  line SEDs remain detectable optically but sit outside the
  colour-veto cell, as frozen.
- W3/W4 take no part (scientific-review rule; threshold statements
  stay with the standalone WISE survey).
- WISE depth limits the veto's power to candidates brighter than
  roughly the optical single-epoch clip — measured, not assumed; for
  the quoted constraint space the veto is free insurance on
  adjudication rather than a depth contribution.
- The v2 report (2026-08-23) is superseded by this document; its run
  remains fully archived (`runs/joint/v2/`, freeze
  `configs/v2_freeze.json`). No optical number changed.
