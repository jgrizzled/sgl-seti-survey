---
title: "Batch 3 — five mid-confusion corridors through the full pipeline"
date: 2026-08-19
analysis_run: "run-04d6098fb26e (supersedes run-2d18a690aa95)"
---

# Batch 3 results

The five mid-confusion clean-gate corridors from the WISE overlay queue
(ε Eri, Lacaille 9352, GJ 1061, GJ 12724, Wolf 1061), curated into
registry v1.3 (`sha256:a476a3e3…`, 26 targets / 23 corridors) and run
through the complete chain. Hypotheses v1.3 addendum: endpoint set
only; physics parameters unchanged from v1.0.

## Curation notes

ε Eri gets the τ Cet treatment (RUWE 2.72 → uncertainties ×3; debris
disk and RV planet noted, reflex negligible at locus scale); GJ 1061's
Gaia RV (+1.49) cross-checked against Medina 2022 (+1.0 ± 0.5),
consistent; GJ 12724's RV from the CNS5 compilation (−18.0 ± 2.0);
Lacaille 9352 brings 6.9 ″/yr corridor drift (best of the batch);
Wolf 1061's small planets noted as reflex-negligible.

## Pipeline outcomes

12,696 new frame-band rows discovered (GJ 1061's corridor the densest
query yet at 4,035); 21,862 new ScreenMatch records (120,911 total);
recurrence peaks all chance/static. Calibration rerun over all **52
endpoint-roles → 1,664 Constraint records** (all valid recovery
curves).

**Threshold exceedances: 13 of 208 pair-bands** — 8 carried over with
verdicts unchanged, 5 new, all vetoed:

| Cell | Excess | Verdict |
| --- | --- | --- |
| gj-1061/tx W1 (+W2, same cell) | 60.5 vs 57.4 | single-phase static star (−4.6/+99.6); star-like W1:W2 (18.3); the W1/W2 pairing is the star's cross-band consistency, not a relay's |
| lacaille-9352/tx W2 | 36.0 vs 34.0 | single-phase static star (43.1/−0.7); W1 same cell 32.3 |
| lacaille-9352/rx W4, tx W4 | +0–3% | single-visit cryo cells; placement variance |

**No surviving candidate across all 26 endpoints / 52 hypotheses.**
Candidate records now total 13, all vetoed with recorded reasons.

## Depths (median 90%-recovery W1, Rx)

Lacaille 9352 17.0, ε Eri 16.7, Wolf 1061 16.4, GJ 1061 15.6,
GJ 12724 13.7 (dense southern field, b = −25° but crowded) — the
mid-confusion tier lands between batch 2 (≈17.5) and the crowded pilot
corridors (≈13), as the overlay grades predicted.

## Status

22 of the universal list's 32 searchable systems now carry completed
WISE searches (23 corridors, 52 endpoint-role hypotheses, 1,664
constraints). Remaining queue: 10 systems — 8 high-confusion (Luhman 16,
61 Cyg, Procyon, Groombridge 34, Luyten's Star, LP 145-141, GJ 1221,
GJ 9193) + 2 bright-star corridors (GJ 1111, Kapteyn's Star); 5
deferred with named unblocking conditions.
