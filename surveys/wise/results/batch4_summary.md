---
title: "Batch 4 + unblocked deferred — portfolio completion"
date: 2026-08-20
analysis_run: "run-76251a7d4ad3 (calibration v0.2.0; supersedes run-04d6098fb26e)"
---

# Batch 4 + deferred results — the WISE portfolio is complete

21 endpoints across 16 corridors: the full high-confusion +
bright-star queue (61 Cyg A/B, Procyon A/B, Groombridge 34 A/B,
Luhman 16 A/B, Luyten's Star, LP 145-141, GJ 1221, GJ 9193, GJ 1111,
Kapteyn's Star) **plus every formerly deferred system** (Struve 2398
A/B and GJ 783 unblocked by clean per-component Gaia solutions with
negligible orbital curvature; ε Ind Ba/Bb and EZ Aqr as photocenter
tracks whose wobble sits inside the |µ| ≤ 1″/yr cell; WISE 0855 and
GJ 11068 on CNS5 compilation astrometry; Luhman 16 via the Garcia
2017 orbit; Procyon via the Bond 2015 orbit with hip2-as-barycenter
validated against 2MASS to 0.31″). Registry v1.4
(`sha256:09366624…`, 47 targets / 38 corridors); hypotheses v1.4
(physics unchanged from v1.0).

## Calibration v0.2.0 (effective-epoch floor)

Per-cell epoch weights capped at 20× the median positive weight;
N_eff reported at every peak. Validation: the batch-2 single-epoch
artifact (lacaille-8760/rx/W1, one frame carrying S=398 of 399) **no
longer exceeds threshold** — the capture mechanism is closed at the
estimator level. (One axis bug in the N_eff reporting crashed the
first calibration attempt and was fixed; tensors were unaffected.)

## Outcomes

- ~50k new frame-band products (98,982 total; near-ecliptic-pole
  corridors like GJ 1221 and Struve 2398 the richest, 5.5–11k rows);
  19,295 cutout pairs fetched with 8 frame failures; 72,409 new
  ScreenMatch records (193,320 total); recurrence peaks all
  chance/static.
- **3,008 Constraint records** over 94 endpoint-roles under
  `run-76251a7d4ad3`. Depth highlights: Kapteyn's Star W1 m₉₀ = 17.2
  (its bright-star flag proved conservative), EZ Aqr 17.1, WISE 0855
  16.8, GJ 783 16.8; the galactic-plane GJ 11068 corridor bottomed at
  13.6 as its confusion grade predicted.
- **27 of 376 pair-bands exceeded threshold; all vetoed** (19 new
  Candidate records; 34 cumulative, all vetoed). New verdict classes
  this batch: single-visit W3/W4 cells coincident with W1-identified
  static sources (Luhman 16, Luyten, Groombridge 34), and
  dense-corridor contamination where the parallax-phase test loses
  power because *distinct* static sources occupy both phase positions
  (gj-11068/rx/W1: balanced phases 754/904 in a b=2.3° field with a
  S≈1000 control null; the catalog-screening null is the backstop —
  no track-consistent source exists at that catalog-bright level).

**No surviving candidate in any of the 94 endpoint-role hypotheses.**

## Status

The universal target portfolio's WISE overlay queue and deferred list
are both empty: **32 systems, 38 corridors, 47 endpoints — done.**
Report: `report/wise_survey_v2.md` (updated to the final state).
Growing coverage further means raising the universal-list horizon or
adding baskets (config change), or moving to the second archive.
