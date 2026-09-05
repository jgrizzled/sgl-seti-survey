---
title: "SPHEREx deferred controls on QR2: six-detector joint cell and template-absorption control"
date: 2026-09-04
status: "complete — joint cell (hypotheses v2.1): confirmatory 61 endpoints / 122 cells, 0 candidates at family-wise α = 0.05, dev 27 / 54, 0 candidates; 2,816 joint constraints; template-absorption control on all 1,056 per-detector cells: median 49 % of a slow source's flux absorbed by the static template (Δm 0.73), corrected limits published report-level; plan §5.15 item O3"
---

# SPHEREx deferred controls (QR2)

**Pins:** SPHEREx v2.0 freeze `sha256:96fbfadb…` unchanged · amendment
v2.1 `surveys/spherex/hypotheses.md` (dev-driven amendment (a)
recorded) · freeze `surveys/spherex/configs/joint6_freeze.json`
(`sha256:6af3c9cd…`, supersedes the pre-amendment `5ba1ff30…`) ·
code `surveys/spherex/joint6.py`,
`surveys/spherex/scripts/template_{control,absorb}.py` · products
`runs/spherex/joint6/`, `runs/spherex/v2/injections_joint/`,
`runs/spherex/template_control/` · tables
`surveys/spherex/results/joint6/report_tables.md`,
`surveys/spherex/results/template_absorption{,_corrected}.md`
(ledger-generated).

## Summary

The two controls deferred from the v2 SPHEREx survey
(`report/spherex_survey.md`) were run on QR2 without waiting for QR3.

1. **Six-detector joint cell — 0 candidates.** One cell per endpoint ×
   role, S_J = Σ_b A_b / √Σ_b B_b over D1–D6 from the stored v2
   accumulators (the declared flat-Fν SED), ring null and R̃_FWER as v2.
   Confirmatory (blind, once): 122 cells, 3 void, R̃_FWER 2.249, 11
   cells with R > 1 against a ring expectation of 12.7, **0
   candidates**. Development: 54 cells, 4 void, R̃_FWER 1.589, 0
   candidates under all three masks. Median final-candidate m90 for a
   persistent flat-Fν source **20.93 AB** (confirmatory; dev 21.39),
   ~0.7–0.9 mag deeper than the best per-detector cell (D4 20.36), the
   expected √6 gain less the family threshold; worst-of-four temporal
   models 19.15 (visit-scale, as in v2). The v1 headline 20.75 was a
   joint stack without error control; the calibrated joint figure is
   20.93 before, and **20.02 after**, the absorption correction below.
2. **Template-absorption control — the v2 depths are overstated for
   slow sources by a median 0.7 mag.** For a ladder of persistent
   sources on every real track (1,056 cells, 29,537 measurements) the
   static template was refitted with the source present: it absorbs a
   median **49 %** of the source's stacked flux (10th–90th percentile
   25–72 %), i.e. Δm = 0.73 mag (p90 1.3–1.9 by detector), with 7.6 %
   of cells above 75 % and 0.6 % above 90 %. The effect is
   cadence-driven, not distance-driven: a node's epochs come from only
   2–3 SPHEREx visits and a two-parameter (a, b·λ) fit across them
   soaks up a source present during one, so absorption falls with the
   number of visits (deep-field gj-687 5–10 %; single-visit corridors
   Ross 128 / Ross 154 / GJ 581 80–90 %; rank correlation with epoch
   count −0.65) and does not depend on residual motion (µ = ±1″/yr:
   0.42 vs 0.42) or on magnitude below the clip (0.44 / 0.47 / 0.47 at
   m90 −1.5 / 0 / +1). Validated end to end on van Maanen (source added
   to the images, templates regenerated, tensor rebuilt): module vs
   image-level absorbed fraction agree to < 0.01 in all six detectors.
   Visit- and block-scale sources are absorbed like persistent ones
   (what matters is the node's weight in the source's visit, not the
   duty cycle); only exposure-flicker is absorbed proportionally less.

**Corrected limits (report-level; frozen records unchanged, flagged).**
Median final-candidate m90 over cells × z intervals, confirmatory set:

| cell | persistent m90 | corrected | worst model m90 | corrected | median Δm | p90 Δm |
| --- | --- | --- | --- | --- | --- | --- |
| D1 | 20.08 | 19.30 | 18.41 | 17.60 | 0.70 | 1.50 |
| D2 | 20.07 | 19.20 | 18.29 | 17.54 | 0.65 | 1.16 |
| D3 | 20.24 | 19.43 | 18.51 | 17.71 | 0.66 | 1.38 |
| D4 | 20.36 | 19.44 | 18.63 | 17.78 | 0.76 | 1.62 |
| D5 | 19.58 | 18.75 | 17.81 | 17.05 | 0.73 | 1.41 |
| D6 | 19.05 | 18.24 | 17.22 | 16.42 | 0.70 | 1.58 |
| J6 (joint) | 20.93 | 20.02 | 19.15 | 18.35 | 0.69 | 1.31 |

Per cell × z interval values: `surveys/spherex/results/template_absorption_corrected.json`.

## Joint cell — rule (frozen, v2.1)

S_J on the common 96 × 3 × 3 grid from the per-detector accumulators
(per-frame cap, |S_e| ≤ 5 clip and template subtraction inherited per
detector); n_J = Σ n_b ≥ 5; a cell needs ≥ 2 detectors with ≥ 5 usable
epochs (all 176 cells had 6); cross-track variants combined per
trajectory by the larger S_max (3 confirmatory cells). Ring null (48
offsets), R = S_max/T (8 designated), R̃ = R/q95, heavy-tail /
inner–outer KS void; R̃_FWER at α = 0.05 over the joint family, declared
a second dependent family on the same data (survey-wide FWER across the
per-detector and joint families ≤ 0.10). No catalogue veto (v2.0 rule);
a survivor would be `retained-ambiguous` with per-detector
decomposition and flat-Fν χ² annotations — none occurred. Injections:
the v2 chain re-run once with a common window per pair (v1 joint m90 ±
2) so the j-th injection is one source in all six detectors, 400 per
pair; six window sums combined as the accumulators; bright limit = the
brightest of the six single-epoch clip limits (median 18.6).

**Dev-driven amendment (a).** The first dev pass gave R̃_FWER = 5.26:
proxima-cen/rx lies in an over-subtracted Galactic-plane field where
every trajectory's joint S_max is negative, so T = 0.057 and the ring
q95 = −12.3, and dividing by a negative normaliser turned the most
negative ring values into the family's largest R̃. Rule added before
the confirmatory run: q95 ≤ 0 → void (`null_degenerate`). The
per-detector v2.0 engine has no such guard (its degenerate cells
happened to be caught by the heavy-tail rule); recorded for QR3.

## Joint cell — results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | candidates | persistent m90 | worst model |
| --- | --- | --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 61 / 122 | 3 (heavy-tail) | 2.249 | 11 (12.7) | **0** | 20.93 | 19.15 |
| development | 27 / 54 | 4 (2 heavy, 2 degenerate) | 1.589 | 5 (5.6) | 0 | 21.39 | 19.46 |

Top confirmatory cells: sigma-dra/rx R̃ 1.34 (global p 0.56), 82-eri/rx
1.29 (0.66, flat-Fν χ² 15.2/5 — the only top cell with an inconsistent
detector decomposition, consistent with a static residual), gj-908/rx
1.12. Flat-Fν χ² over all cells: median 2.9 for 5 d.o.f.; 3.3 % of
confirmatory cells above the p < 0.05 point (annotation only). Mask
sensitivity (dev): strict R̃_FWER 1.502 vs 1.589, 0 candidates in all.
Constraints: 2,816 (1,888 confirmatory recovery curves of 1,952; 772
of 864 dev). Caveat: the visit / block on-patterns of the injections are
drawn per detector (D1/D4, D2/D5, D3/D6 share exposure times in
reality), an inherited property of the chain.

## Template-absorption control — method

`surveys/spherex/scripts/template_absorb.py`: per corridor × detector the
epoch cache of `static_template.py` (node samples f, w, λ) is rebuilt at
the nodes within 15″ of the ladder tracks; for each source the injected
per-epoch flux is rendered with the exposure PSF plane and the v2
matched-filter response window at those nodes, the template is refitted
(full two-pass 3σ clipped fit) with the injection added, and
f_abs = Σ_e w_e T_inj(track_e, λ_e) / Σ_e w_e δ_e(track_e) with the v2
stack's own capped, clipped weights. Ladder per cell: the 8 z-interval
centres × (m90 −1.5, m90, m90 +1) at µ = 0, plus µ = (±1, ±1)″/yr at the
two most distant z, persistent, cross-track 0. Products
`runs/spherex/template_control/absorb/<corridor>__<band>.json`; summary
`surveys/spherex/results/template_absorption.md`.

## Templates regenerated

The v3 static templates cited by the v2 report were lost in the v1
retirement (2026-08-24): `calib_v4/templates` was a symlink into the
deleted `calib_v2`. They were regenerated from the retained cutouts with
the unchanged algorithm (`template_control.py fit`, 462 groups) into a
real `runs/spherex/calib_v4/templates/` directory. Verification
(`verify`): a van Maanen tensor rebuilt with the regenerated templates
matches the stored v2 tensor at every valid node to < 0.05σ per epoch
and to 2.5 × 10⁻⁴ in S_max — equivalent for every statistic, not
bit-identical (the lost set was built from a slightly earlier
photometry state).

## Scope and forward fix

Targeted coverage of the frozen 88-endpoint portfolio in QR2; no
population inference. The absorption is a property of fitting a static
template with 2–3 visits per node and will shrink as SPHEREx visits
accumulate; for the QR3 re-run the plan is a source-excluded template
(per node, drop the epochs during which the hypothesised trajectory is
within ~2 FWHM), which removes it by construction, plus the
degenerate-null guard in the per-detector engine.
