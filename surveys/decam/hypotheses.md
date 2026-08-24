# DECam southern survey — hypotheses v1.0 (FROZEN 2026-08-24)

Status: **frozen** (user approved the §8 decisions 2026-08-24; the
hash of this document is pinned in `configs/v2_freeze.json`). The v2
discipline applies from the first analysis: stratified
dev/confirmatory split with a declared seed and the pilot corridors
forced into the development set, null-ensemble thresholds,
FWER α = 0.05 with BH q-values informational, image-level injections.
There is no v1-style exploratory pass for this survey. Version history:
v0.1 draft (2026-08-24) accompanied the geometric stages; v1.0 freezes
the §8 decisions and the positive-control selection.

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
  `GeometryContext.decam_default()`. **Validated 2026-08-24**
  (`results/observer_validation.json`): site effect 7–16 mas at
  z = 550 AU, residual vs an independent astropy computation
  ≤ 1.7 mas (typical 0.3 mas) — inside the ZTF pilot's 16 mas bar.

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
a systematic per the v2 common rules. Known archive-side losses (recon
addendum): ~1% of files need the adapter's recovery paths (hdus-500
fallback, per-file HDU order, content-verified md5 drift); 1 of 214
pilot exposures has an unrecoverably corrupt image and is excluded as
product-corrupt, not as a scientific veto.

## 4. Flux scale

Per-frame star calibration through the identical matched filter
(PS1 lesson 1): primary catalog NSC DR2 (same instrument, verified
reachable anonymously), with ATLAS RefCat2 / Gaia DR3 via VizieR as
cross-checks. Verified to ≤ 0.2 mag on catalogued stars before any
depth is quoted (SPHEREx rule). Measured evidence
(`results/flux_scale_check.json`, 24 frames): header MAGZERO is in
*counts* convention (not counts/s) and is per-frame unreliable
(outliers to +3.3 mag; u systematically +2.7) — it serves only as a
cross-check annotation. Per-frame star ZP scatter (aperture MAD)
0.014–0.10 mag with 35–6,300 NSC stars per CCD.

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
- Positive-control asteroid: the corridors sit at ecliptic latitudes
  −28° to −66°, so no main-belt object crosses them — the control uses
  separate near-ecliptic DECam exposures through the same adapter and
  estimator (target to be picked at freeze).

## 8. Frozen decisions (approved 2026-08-24)

1. **Quality cut**: EXPTIME ≥ 30 s, bands grizY (u and narrowband
   recorded, excluded from the search cell), obs_type = object;
   PHOTFLAG recorded only. Overlay v1 under this cut: 10 ok /
   5 crowded / 0 sparse corridors.
2. **Fatal dqmask codes**: all-but-INTERPOLATED (4) as primary;
   strict (all codes fatal) / loose (CR = 16 also allowed) reported as
   the mask-sensitivity systematic.
3. **Dev/confirmatory split**: stratified by the WISE v2 confusion
   class (shared across surveys), seed declared in the freeze, with
   the three pilot corridors (lalande, sigmadra, hd219134) **forced
   into the development set** — their data shaped the pipeline. The
   engine's `forced_dev` mechanism records them in the freeze.
4. **gj-1221 exposure selection**: deepest exposure per night per band
   (cap against the 2,033-exposure LMC pile-up); the cap is part of
   the frozen configuration, applied at tensor build.
5. **Estimator**: direct stack with ring-offset controls (PS1/WISE v2
   pattern); a static-sky template is deferred to a later version if
   the LMC corridors demand it.
6. **Positive control**: **(60000) Miminko** — the PS1 control, reused
   for cross-survey continuity. 12 instcal exposures over 12 nights
   (2013–2019), bands g/r/i/z, Horizons V 19.1–20.6, all triplets
   re-discovered through the survey adapter
   (`configs/asteroid_control_v1.json`; SSOIS + Horizons responses
   snapshotted under `runs/decam/asteroid_control_v1/`).
7. **Multi-CCD stamps**: mosaic in the flux-map builder — per-CCD
   matched filtering (injections enter each CCD image before the
   filter), then nearest-neighbour pasting onto a common TAN canvas
   (≲ 0.15" placement error vs the 1" locus tolerance).
