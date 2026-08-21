---
title: "Pan-STARRS1 pilot — baseline hypothesis freeze"
status: "v1.0 — frozen 2026-08-20 (3-corridor pilot)"
date: 2026-08-20
---

# Baseline hypothesis freeze — PS1 v1.0

Fourth-adapter pilot (plan §3.6 priority 4; first non-IRSA archive,
started while IRSA throughput is saturated by the ZTF scale-up and the
SPHEREx pilot). Parameters are inherited from the WISE freeze
(`surveys/wise/hypotheses.md` v1.0) and the ZTF freeze
(`surveys/ztf/hypotheses.md` v1.0) wherever the physics is
archive-independent, so that PS1, ZTF and WISE nulls on the same
endpoint close the *same* cell of the §3.1 grid. Only items marked
**PS1-specific** are new. Run configs record this file's version and
content hash; any change creates a new version.

## 1. Target endpoints and target-state models

The **same three corridors as the ZTF pilot**, registry
`registries/pilot_wise_2026.yaml` unchanged, so that PS1 (2009–2014)
and ZTF (2018–) give a direct 5–10-year baseline extension per corridor:

| endpoint | Dec of corridor | why |
|---|---|---|
| `ross-128` | −0.8° | ecliptic corridor, clean high-latitude field |
| `eps-ind-a`, `eps-ind-b` | +56.8° | near-circular parallax ellipse; southern star reached from the north |
| `proxima-cen` | +62.6° | Galactic-plane (b +1.9°) crowding stress case at 1.2" seeing |

All three are north of the PS1 3π limit (δ > −30°).

## 2. Role

Rx and Tx as separate hypotheses, both evaluated (unchanged).

## 3. Relay-distance prior and sampling

550–10,000 AU, log-uniform prior, sglseti adaptive-locus sampling with
exact covered-interval reporting (unchanged).

## 4. Geometry model and observer-state versions — PS1-specific observer

- Geometry: sglseti `tusay2022_eq5_7_v1` (unchanged).
- Observer: **PS1 Haleakalā site** from the warp headers
  (`FPA.LONGITUDE` 10.4171 h W = −156.2559°, `FPA.LATITUDE` +20.7071°,
  `FPA.ELEVATION` 3048 m) via `Observer.from_geodetic`. Topocentric
  shift ≤ 0.016" at 550 AU (as for Palomar), carried in the accuracy
  budget.
- Epoch: warp `MJD-OBS` is exposure start (fitscut labels it TAI; the
  ≤ 35 s TAI−UTC offset moves the locus < 0.003" and is carried in the
  budget). Loci are evaluated at mid-exposure using the nominal 3π
  exposure time per filter (g 43, r 40, i 45, z 30, y 30 s) at discovery
  and the header `EXPTIME` from the precise pass on. Drift over a 60 s
  exposure ≤ 0.005".

## 5. Uncertainty confidence level and search padding

99% propagated locus confidence + **10" fixed padding** (unchanged; 40
PS1 pixels).

## 6. Source morphology and spectral model — PS1-specific interpretation

- Morphology: unresolved point source at PS1 resolution (median 3π
  seeing ≈ 1.0–1.3"; warps are *not* PSF-homogenised, so the per-warp
  header `CHIP.SEEING` is used).
- Bands: g (0.48 µm), r (0.62), i (0.75), z (0.87), y (0.96); every
  public warp used, constraints reported per filter.
- Physical interpretation identical to ZTF §6: reflected sunlight is the
  primary cell (m_⊙,r = −26.93; limits quote albedo explicitly);
  self-luminous optical emission as a flux limit; never thermal coverage.

## 7. Persistence / duty-cycle assumption

Duty cycle ≥ 0.5; injections sample [0.5, 1] (unchanged).

## 8. Stationkeeping / residual-motion bounds

|µ_resid| ≤ 1"/yr fit jointly with z (unchanged). Over the ~4.5-year PS1
baseline this is ±9 pixels, a resolved nuisance parameter as in ZTF.

## 9. Detection pipeline and decision threshold — PS1-specific inputs

Layered per plan §3.4. Frozen data-quality inputs:

- Frame selection: all warps listed by `ps1filenames.py` (`badflag` is
  recorded; none were flagged at the pilot antipodes). No hard cuts on
  seeing or airmass; they enter the stack as per-epoch weights.
- Usable pixels: mask template **16255** = IPP `MASK.VALUE` 8575
  (DETECTOR|FLAT|DARK|BLANK|CTE|SAT|LOW|CR|CONV.BAD) | SPIKE 512 |
  GHOST 1024 | STREAK 2048 | STARCORE 4096. SUSPECT (128) and CONV.POOR
  (16384) are **not** fatal. Pixels the exposure did not touch are NaN
  in the image and BLANK/CONV.BAD in the mask.
- Search image: the sky-subtracted warp (PS1 publishes no single-epoch
  difference images). The static field is handled by the single-epoch
  clip (|S_e| ≤ 5, as ZTF) and the parallax-phase test. Variance from the
  `.wt` plane (verified to be variance, recon note), Gaussian kernel at
  the header seeing.
- **Flux scale (PS1-specific, motivated by the 2026-08-20 flux-scale
  erratum):** every warp's zero point is calibrated *empirically* by
  recovering DR2 `mean`-table stars (15 < mag < 19.5, nDetections ≥ 5,
  unsaturated) through the identical matched filter in the same cutout:
  ZP_star = median(m_DR2 + 2.5 log10 F_mf). This absorbs the
  Gaussian-vs-true-PSF throughput of the estimator and any header ZP
  convention; the header `FPA.ZP` is recorded alongside. Frames with < 5
  calibrators fall back to the median ZP_star − FPA.ZP offset of their
  filter.
- Screening catalogs: DR2 `detection` (per-epoch PSF photometry) matched
  to warps by `obsTime`, plus the DR2 `mean` table as recurrence
  context. Catalog absence is neither detection nor null (plan §3.4).
- Parallax-phase test, 8 offset-trajectory controls, effective-epoch
  weight cap (20× median), per-cell epoch floor N ≥ 5 and 1/z-uniform
  z-grid carried over from WISE calibration v0.2.0 / ZTF v1.0. z-grid:
  **360 nodes** (≈ 1" spacing, under the seeing).
- Duplicate epochs: the same exposure warped onto two overlapping
  skycells is one epoch; the sample tensor keeps the warp with the
  larger usable fraction along the locus.
- Positive control: a numbered main-belt asteroid recovered through the
  shift-and-stack code along its JPL ephemeris (site code F51).

Specific thresholds are frozen in the run config when each stage is
first executed.

## Out of scope for this freeze

Barycenter/planetary endpoints, swarms, off-axis infrastructure,
inactive relics, duty cycle < 0.5, trailed sources, PS1 stack images
(static-sky products), proprietary/other PS1 surveys (MD fields).
