---
title: "SPHEREx survey — Report v2 (frozen decision rule, blind confirmatory run)"
date: 2026-08-24
status: "complete for QR2 — confirmatory 61 endpoints / 732 detector cells: 0 candidates at family-wise α = 0.05; development 27 / 323: 0 candidates; 16,896 injection-calibrated constraints; supersedes spherex_survey_v1.md (status note 2026-08-22); re-run planned when the next quick release adds a parallax phase"
---

# SPHEREx survey — Report v2

**Pins:** registry v1.5 · hypotheses `surveys/spherex-v2/hypotheses.md`
v2.0 (freeze `sha256:96fbfadb…`) · QR2 only · engine `sglsurvey/v2/` ·
v3 static templates (`runs/spherex/calib_v4/templates`) · numbers from
`surveys/spherex-v2/results/report_tables.md` (ledger-generated).

## Summary

The v1 SPHEREx survey (77 corridors, 16-control thresholds, template-
subtracted stacks) was re-analysed under the WISE v2.1 design: cells
per detector D1–D6, the exposure's own PSF plane as the matched filter
and as the injected source, the v3 static-sky template subtracted at
sampling (with full SIP-WCS sampling — a local-linear approximation
mis-samples the 20′ cutouts at the pixel level and defeats the
template), a 48-offset ring null with a family-wise threshold, and a
blind confirmatory hold-out. **No cell of the confirmatory set (732
cells) reaches R̃_FWER = 2.314**; 74 cells have R > 1 against a ring
expectation of 81. Development set: 323 cells, R̃_FWER 2.353, 0
candidates. Median 90 %-completeness for a persistent flat-Fν source:
D1–D4 20.1–20.4, D5 19.6, D6 19.1 AB; worst-of-four temporal models
17.3–18.7 (the visit-scale model dominates: QR2 has only 2–3 observing
seasons). The injected-PSF throughput is 1.03 (the kernel is the PSF).

## Rule (frozen)

As WISE v2.1 with SPHEREx bindings: 96 × 3 × 3 grid (T0 = 61000, ZP
23.9 µJy AB), single-epoch clip |S_e| ≤ 5, per-frame cap, ring
20/30/40″, R̃ = R/q95(ring), heavy-tail / radius-dependent void flags
(18 confirmatory cells void), R̃_FWER over the family. **No catalogue
flux-consistency veto:** the static template is the calibrated
treatment of catalogued static sources — a catalogue test on
template-subtracted fluxes double-counts (measured on the van Maanen
test cell before the freeze). A surviving candidate would be
`retained-ambiguous` with annotations; none occurred. No epoch
hold-out (one quick release).

## Results

| set | endpoints / cells | void | R̃_FWER | R > 1 (expected) | candidates |
| --- | --- | --- | --- | --- | --- |
| confirmatory (blind, once) | 61 / 732 | 18 | 2.314 | 74 (81.3) | **0** |
| development | 27 / 323 | 9 | 2.353 | 28 (33.2) | 0 |

Mask sensitivity (development): primary = loose = all exposures and
strict (deep-field collection excluded) give identical R̃_FWER (2.353)
and 0 candidates.

## Completeness (confirmatory, AB, median over cells × intervals)

| detector | persistent | flicker | visit | block |
| --- | --- | --- | --- | --- |
| D1 | 20.08 | 19.26 | 18.44 | 18.79 |
| D2 | 20.07 | 19.14 | 18.38 | 18.73 |
| D3 | 20.24 | 19.38 | 18.61 | 18.71 |
| D4 | 20.36 | 19.51 | 18.64 | 18.97 |
| D5 | 19.58 | 18.69 | 17.87 | 18.19 |
| D6 | 19.05 | 18.17 | 17.34 | 17.34 |

16,276 of 16,896 constraint records are recovery curves; the rest are
`not_constrainable` (void null, too few epochs on, or a fit brighter
than the single-epoch clip limit). v1's headline 20.75 AB was the
six-detector joint stack without error control; the v2 per-detector
figures carry the family-wise threshold and the temporal-model floor.
The six-detector joint cell can be formed from the stored per-tensor
accumulators and is deferred to the QR3 re-run.

## Geometry and caveats

Cross-track dimension for 3 of 176 cells (ε Ind Ba/Bb rx/tx, EZ Aqr
tx), as WISE. Caveats: the static template was fitted without the
injected sources, so partial absorption of a slow (z ~ 10,000 AU) real
source into the template is not captured by the injections; no moving
positive control has been run through the spectral path (the v1 star
control validated photometry only) — both are queued for the QR3
cycle, when the added parallax phase also makes the phase-coherence
scramble informative.

## Scope

Targeted coverage of the frozen 88-endpoint portfolio in QR2; no
population inference.
