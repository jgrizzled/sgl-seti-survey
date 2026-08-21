---
title: "SPHEREx pilot — baseline hypothesis freeze"
status: "v1.0 — frozen 2026-08-20 (3-star pilot)"
date: 2026-08-20
---

# Baseline hypothesis freeze — SPHEREx v1.0

Third-adapter pilot (plan §6). Parameters are inherited from the WISE
freeze (`surveys/wise/hypotheses.md` v1.0) wherever the physics is
archive-independent, so that a SPHEREx null, a WISE null and a ZTF null
on the same endpoint close the *same* cell of the §3.1 grid in
different bands. Only items marked **SPHEREx-specific** are new. Run
configs record this file's version and content hash; any change
creates a new version.

## 1. Target endpoints and target-state models

Pilot corridors (3), chosen from the 15 universal-list systems whose
corridors lie at Dec < −28° (unreachable by ZTF), reusing registry
`registries/pilot_wise_2026.yaml` entries **unchanged**:

| endpoint | why it is in the pilot |
|---|---|
| `lalande-21185` | nearest star ZTF cannot reach (2.55 pc); corridor at RA 345.8°, Dec −35.9°, b −68° (clean, high-latitude field); 349 Level-2 exposures |
| `gj-687` | corridor at RA 84.0°, Dec −68.3°, 3° from the south ecliptic pole, inside the SPHEREx deep field: 7,308 exposures with continuous cadence — the parallax-phase and spectral-sampling stress case; engineering-backbone basket member |
| `sigma-dra` | engineering-backbone top pick (5.8 pc, quiet G9V); corridor at RA 113.0°, Dec −69.7°, 9° from the pole and outside the deep field; 643 exposures |

All three are linear-astrometry endpoints. Component endpoints of the
multiple systems (61 Cyg, Struve 2398, Groombridge 34, GJ 338) follow
in the scale-up, as in WISE and ZTF.

## 2. Role

Rx and Tx as separate hypotheses, both evaluated (unchanged).

## 3. Relay-distance prior and sampling

550–10,000 AU, log-uniform prior, sglseti adaptive-locus sampling with
exact covered-interval reporting (unchanged).

## 4. Geometry model and observer-state versions — SPHEREx-specific observer

- Geometry: sglseti `tusay2022_eq5_7_v1` (unchanged).
- Observer: **Earth centre**. SPHEREx is in a ~650 km sun-synchronous
  LEO; the spacecraft-vs-geocentre locus shift is ≤ 7,000 km / 550 AU
  = 0.017" at the nearest relay distance and is carried as an accuracy-
  budget term, not modelled (the Level-2 headers carry the geocentric
  spacecraft state `X_SC, Y_SC, Z_SC` should a later version need it).
- Epoch: ObsCore `t_min`/`t_max` bracket the 113 s integration; loci
  are evaluated at the midpoint. Corridor drift over one exposure is
  ≤ 0.009" (plan §3.3 midpoint model adequate).

## 5. Uncertainty confidence level and search padding

99% propagated locus confidence + **10" fixed padding** (unchanged).
At 6.15"/pix this is 1.6 pixels; the coarse stage adds the polygon-vs-
SIP footprint allowance (§9).

## 6. Source morphology and spectral model — SPHEREx-specific interpretation

- Morphology: unresolved point source at SPHEREx resolution (PSF FWHM
  ≈ 5.3" median, under-sampled on 6.15" pixels). Extended or trailed
  morphologies are separate cells.
- Bands: the six detectors D1–D6 (0.75–1.12, 1.10–1.64, 1.62–2.42,
  2.40–3.82, 3.80–4.42, 4.40–5.00 µm) are the reporting bands. Each
  detector is a linear-variable filter, so every exposure samples one
  wavelength (R ≈ 40–130) at a given sky position; **the per-epoch
  wavelength and bandwidth at the track position are carried on every
  sample** and every Constraint quotes the wavelength range it covers.
  A stack across epochs of one detector is therefore a broadband
  (detector-wide) stack; the spectral axis is used for candidate
  discrimination, not for the baseline detection statistic.
- Physical interpretation of a near-IR flux limit, in order of
  relevance by detector:
  1. **Reflected sunlight** (D1–D3 primarily; D4–D6 secondarily):
     H-D convention with a solar-coloured reflector,
     H_V = m_AB(λ) − (M_⊙,AB(λ) − 4.81), D_km = 1329 p^{−1/2} 10^{−H_V/5}.
     Limits quote p explicitly (p = 0.1 reference).
  2. **Hot components** (D4–D6, 3.8–5 µm, the W1/W2 regime): a
     blackbody-equivalent emitting diameter at T = 400, 700 and
     1000 K is quoted per constraint. A 550 AU solar-equilibrium body
     (12 K) is invisible at all SPHEREx wavelengths; these limits are
     never to be described as cold-thermal coverage.
  3. **Self-luminous non-thermal emission** — reported as µJy.

## 7. Persistence / duty-cycle assumption

Duty cycle ≥ 0.5; injections sample [0.5, 1] (unchanged).

## 8. Stationkeeping / residual-motion bounds

|µ_resid| ≤ 1"/yr, fit jointly with relay distance (unchanged). Over
the 14.5-month QR2 baseline this is ≤ 1.2" = 0.2 pixel: residual
motion is an *unresolved* nuisance parameter in this pilot (as in
WISE), retained in the grid for uniformity with the other surveys.

## 9. Detection pipeline and decision threshold — SPHEREx-specific inputs

Layered per plan §3.4. Frozen data-quality inputs:

- Frame selection: `spherex.obscore` rows with `data_rights='public'`,
  `dataproduct_type='image'`, `calib_level=2`, collections
  `spherex_qr2` and `spherex_qr2_deep` (QR1 rows are superseded by the
  QR2 reprocessing and are not served). No exposure-level quality cut
  exists in ObsCore; quality enters per pixel through FLAGS.
- Usable pixels: FLAGS fatal template **708343** = bits {0 TRANSIENT,
  1 OVERFLOW, 2 SUR_ERROR, 4 PHANTOM, 5 REFERENCE, 6 NONFUNC,
  7 DICHROIC, 9 MISSING_DATA, 10 HOT, 11 COLD, 14 PHANMISS,
  15 NONLINEAR, 17 PERSIST, 19 OUTLIER}; bits 12 (FULLSAMPLE,
  informational) and 21 (SOURCE — a catalogued source is present) are
  **not** fatal.
- Search image: IMAGE − ZODI (pipeline zodiacal-light model) with a
  robust constant residual background removed per cutout. No
  reference-subtracted product exists; the static field is handled by
  the parallax-phase test and the layer-1 screen.
- Matched filter: the exposure's own PSF-zone plane (11 × 11 detector
  zones, 10× oversampled) binned to detector sampling — not a
  Gaussian (lesson of the ZTF asteroid control).
- Variance: pipeline VARIANCE plane rescaled per cutout to the robust
  (MAD) scatter of the residual image (factor recorded per epoch; it
  is > 1 in sparse fields where confusion and zodi-model residuals are
  missing from the model and ≈ 0.3 in the deep fields where the model
  is conservative). Depths are invariant to a uniform rescale; the
  single-epoch clip is not, hence the empirical choice.
- Flux scale: matched-filter amplitude (MJy/sr summed over the PSF) ×
  per-exposure median pixel solid angle (`OMEGA_MEDIAN`) → µJy;
  AB = 23.9 − 2.5 log10(F/µJy). **Verified against catalogued
  2MASS/CatWISE stars in the same cutouts before any depth is quoted
  (§6 success criterion: ≤ 0.2 mag).**
- Screening: SPHEREx QR2 publishes no source catalogs. Layer 1 is a
  per-exposure peak search (S_e > 5) within 10" of the track; peaks
  recurring at a fixed sky position are static sources, cross-matched
  against CatWISE2020 (snapshotted). Catalog absence is neither
  detection nor null (plan §3.4).
- Parallax-phase test, offset-trajectory controls (8), effective-epoch
  weight cap (20× median), per-cell epoch floor (5) and single-epoch
  clip (|S_e| ≤ 5) carried over from ZTF calibration.
- Positive control: flux-scale recovery of catalogued stars (above).

Specific thresholds are frozen in the run config when each stage is
first executed.

## Out of scope for this freeze

Barycenter/planetary endpoints, swarms, off-axis infrastructure,
inactive relics, duty cycle < 0.5, trailed sources, spectral-line
(single-channel) searches along the track, Level-3 mosaics.

## Amendments (v1.0, same freeze)

- **2026-08-20, calibration stage — role-coincidence veto.** For
  linear-astrometry endpoints the Rx and Tx loci coincide to < 1" at
  every relay distance, so a threshold exceedance in one role's stack
  is vetoed when the other role's stack at the same (z, µ) cell is
  below 2σ. Added after the first calibration run showed one
  phase-passing exceedance (σ Dra rx D4) absent from the coincident Tx
  stack. Coarse/precise run configs carry the pre-amendment content
  hash of this file; the calibration run config carries the amended one.
- **2026-08-21, calibration stage — v2 estimator and two adjudication
  rules.** (i) Epoch variances are the per-cutout MAD variance scaled
  by the static template's per-node reduced χ² (the MAD includes
  static confusion the template removes; depths improve by a median
  0.44 mag). (ii) A joint six-detector stack (band `ALL`, flat-F_ν
  injection) is searched and constrained alongside the per-detector
  stacks, with the same controls and vetoes. (iii) An exceedance is
  vetoed when the template covered < 80 % of the real-trajectory
  epochs at the peak cell (the stack layer's cell is "static sky
  removed"), or when a layer-1 static cluster ≥ 5 mag brighter than
  the stack depth lies within 3 PSF FWHM of the track position (PSF-
  wing residuals; the GJ 13157 D2 case at z = 10,000 AU, 10–19″ from
  a W1 = 7.7 star in a b = −3° field). Calibration run `calib_v3`.
- **2026-08-21, calibration stage — 16 controls and the cross-detector
  phase test (amendment 3).** (i) Sixteen offset-control trajectories
  (eight added: the four ±50″ diagonals, ±75″ Dec, and two skew
  offsets) so the per-search empirical false-alarm rate of the max-of-
  controls rule is < 1/16. (ii) Where the exceeding stack's own
  parallax-phase split is undefined (fewer than 3 epochs in one
  phase — the deep-field case, in which detectors observe in
  different seasons), the lacking phase is stacked across *all six
  detectors* at the same cell; the exceedance is vetoed if that
  complement has ≥ 3 epochs and S < 2, and as single-season if it has
  < 3. Calibration run `calib_v4`.
