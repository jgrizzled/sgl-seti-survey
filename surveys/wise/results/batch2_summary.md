---
title: "Batch 2 — nine low-confusion corridors through the full pipeline"
date: 2026-08-19
analysis_run: "run-2d18a690aa95 (supersedes run-20040d16f249)"
---

# Batch 2 results

The first nine corridors of the WISE overlay work queue (all
low-confusion, clean-gate: Ross 128, ε Ind A, τ Cet, GJ 54, Teegarden's
Star, Lacaille 8760, van Maanen's Star, GJ 908, GJ 784) curated into
registry v1.2 (`sha256:f03b821f…`, 21 targets) and run through the
complete chain. Hypotheses v1.2 addendum: endpoint set only; physics
parameters unchanged from v1.0.

## Curation notes (documented per entry)

- τ Cet adopted with ×3-inflated uncertainties (Gaia RUWE 2.63,
  saturation-bright); GJ 54 adopted linear despite literature binarity
  hints (RUWE 1.00, noted); van Maanen's Star is the first white-dwarf
  wildcard searched — its catalog RV is flagged GR-contaminated with a
  ±50 km/s allocation (locus impact < 0.1″).
- **ε Ind Ba/Bb deferred** (RUWE 4.3: photocenter rides the ~11-yr
  T-dwarf binary orbit → ±1.5–3″ if treated as linear). Unblock by
  curating the published Ba/Bb orbit.

## Pipeline outcomes

18,000 new frame-band rows discovered; 31,675 new ScreenMatch records
(99,049 total); recurrence peaks in all nine corridors resolved as
chance/static patterns. Calibration rerun over all **42 endpoint-roles →
1,344 Constraint records** (every cell a valid recovery curve).

**Threshold exceedances: 8 of 168 pair-bands** — the 5 previously
adjudicated (verdicts unchanged) plus 3 new, all vetoed:

| Cell | Excess | Verdict |
| --- | --- | --- |
| lacaille-8760/rx/W1 | S=399 vs T=42 | **single-epoch dominance**: one frame (MJD 57855) carries S=397.8 of 399.3; single-phase; W2=4.8. Transient/artifact — fails persistence by construction |
| ross-128/tx/W1 | S=225 vs T=183 | single-phase static star (10.9/377.5); star-like W1:W2 |
| van-maanen/tx/W1 | S=95 vs T=88 | single-phase static star (158.4/−2.4); star-like W1:W2 |

**No surviving candidate across all 21 endpoints / 42 hypotheses.**
Candidate records now total 9, all vetoed with recorded reasons.

## Depths: the queue ordering worked

The low-confusion corridors delivered the survey's deepest limits —
median 90%-recovery W1 depths (Rx): **Lacaille 8760 17.5, GJ 54 17.4,
GJ 908 16.8, van Maanen 16.6, τ Cet 16.3**, GJ 784/Ross 128 15.6–15.7,
Teegarden 14.3 (bright target-adjacent field). Compare the pilot's
crowded corridors (Sirius ≈ 13). Full curves in
`runs/wise/calib_v1/m90_curves.npz`.

## Pipeline lesson

Single-epoch dominance (the Lacaille 8760 excess) is a new veto class:
an inverse-variance stack can be captured by one frame with tiny formal
variance. Adopted as a manual check this round; a production fix is an
effective-epoch floor (or per-epoch weight cap) inside the stack — filed
as an improvement for the next calibration version.

## Status

17 of the universal list's 32 searchable systems now carry completed
WISE searches (8+9 system nodes = 18 corridors incl. Proxima/α Cen
split). Remaining queue: 15 systems (mid/high confusion + bright-star
corridors); 5 deferred with named unblocking conditions.
