# DECam southern survey — hypotheses (DRAFT v0.1 — NOT FROZEN)

Status: **draft**. This document accompanies the geometric stages
(coarse discovery, precise pass), which are hypothesis-light and
reusable. Before any screening threshold, stack, injection, or
candidate rule runs, this document must be frozen as v1.0 with the v2
discipline (dev/confirmatory split with declared seed, FWER α = 0.05,
null-ensemble thresholds, image-level injections) — no v1-style
exploratory analysis pass exists for this survey. The user owns the
freeze (scientific parameters).

## 1. Physics cell (inherited)

Same baseline local-artifact hypothesis as WISE v2.1 / ZTF v2 /
PS1 v2: compact station-kept relay on the anti-star focal ray,
log-uniform relay-distance prior 550–10,000 AU, Rx and Tx as separate
hypotheses, duty cycle ≥ 0.5, |µ_resid| ≤ 1"/yr (L∞ box). Bands grizY
carry the reflected-sunlight and self-luminous-optical
interpretations; u and narrowband frames are recorded but excluded
from the search cell (draft; freeze decides).

## 2. Geometry

- Model/ephemeris: Tusay2022Eq57V1 + astropy built-in (as all v2
  surveys); 99% confidence, +10" padding, 5" tolerance.
- Observer: CTIO Blanco 4-m terrestrial site (`ctio-blanco-decam`,
  lon −70.80655°, lat −30.16928°, h 2207 m; MPC W84), via
  `GeometryContext.decam_default()`. Validation target: topocentric
  shift vs Earth-centre ≤ 0.02" agreement with sglseti's site model,
  as the ZTF pilot did for Palomar (≤ 0.016").

## 3. Archive and products

NOIRLab Astro Archive, DECam `instcal` (community pipeline) image +
dqmask + wtmap triplets joined on EXPNUM; latest CP version per
exposure. Discovery via the advanced-search API (snapshotted). Exact
astrometry and usable pixels from the dqmask (per-CCD TPV WCS).
Single-CCD retrieval via `?hdus=` (the archive ignores HTTP Range).
Heterogeneous PI programmes are the norm: per-exposure quality cuts
(draft: obs_type = object, EXPTIME ≥ 30 s, grizY only, PHOTFLAG
recorded) carry more weight than in survey-uniform archives and are a
freeze-level decision.

Data-quality convention (draft): fatal dqmask codes = all except
INTERPOLATED (4); mask-sensitivity (primary/strict/loose) reported as
a systematic per the v2 common rules.

## 4. Flux scale

Per-frame star calibration through the identical matched filter
(PS1 lesson 1): primary catalog ATLAS RefCat2 (VizieR
J/ApJ/867/105), cross-checked against Gaia DR3 and the header
MAGZERO; NSC DR2 (Data Lab, `astro-datalab` client, anonymous sync)
where reachable. Verified to ≤ 0.2 mag on catalogued stars before any
depth is quoted (SPHEREx rule).

## 5. Screening and vetoes

Layer-1 screening against NSC DR2 `meas` (per-exposure detections;
object-cone-first then `objectid IN (...)` — direct meas cones time
out) with the caveat that NSC DR2 is time-partial (ends ~2017-2019 in
probed fields). Flux-consistent catalogued-static veto per
`sglsurvey/vetting.py` on Gaia DR3 + NSC objects (≥ 3 detections).
Parallax-phase test is expected to carry weight single-archive
(corridors see 4–12 calendar months). All other rules are annotations
unless injection-calibrated (v2 common rules).

## 6. Controls

Positive control: a known main-belt asteroid recovered through the
full chain with Horizons magnitudes (ZTF/PS1 pattern). Negative
controls and thresholds: v2 null ensemble (extended spatial offsets +
time scrambling + trajectory randomisation), survey-wide max-R null,
FWER ≤ 0.05.

## 7. Known survey-specific risks (to resolve at freeze)

- gj-1221 / gj-687 corridors sit near the LMC: crowding class +
  exposure-selection policy needed (12,652 / 1,158 exposures).
- Coverage is inhomogeneous (35–318 exposures elsewhere): per-cell
  epoch floor (ZTF lesson) and single-epoch clip conventions.
- No archive difference images: static-sky template route (SPHEREx
  `static_template.py` pattern) vs direct stack decided at freeze.
- Pre-2012 Mosaic II epochs (Data Lab SIA) are out of scope for this
  survey (separate instrument/adapter if ever wanted).
