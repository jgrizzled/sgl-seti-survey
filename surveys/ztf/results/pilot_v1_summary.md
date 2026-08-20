# ZTF pilot v1.0 — stage summaries (2026-08-20)

Hypotheses `ztf-hypotheses-v1.0`; registry `pilot_wise_2026.yaml` v1.5;
corridors `ross128`, `epsind`, `proxima`. All products regenerable from
the tracked scripts; bulky data under `runs/ztf/` (1.6 GB).

## Coarse discovery (`coarse_v1`)

IBE `ztf/products/sci`, public (`ipac_gid=1`), `imgtype=object`,
MJD 58150–61300. 7,930 quadrant-exposures discovered (g 3,888 /
r 3,924 / i 118; 313 flagged bad-quality via infobits bit 25),
epochs 2018-03-25 → 2026-06-20. 15,860 coarse evaluations, 1,896 hits.

| endpoint | frames | rx hits (g/r/i) | tx hits |
|---|---|---|---|
| ross-128 | 1,785 | 174 / 183 / 0 | 149 / 152 / 0 |
| eps-ind-a | 3,460 | 5 / 12 / 0 | 6 / 14 / 0 |
| proxima-cen | 2,685 | 287 / 322 / – | 281 / 311 / – |

**Grid-gap finding.** All three antipodes fall in or beside inter-CCD
gaps of the fixed primary ZTF field grid: ε Ind A at the 4-CCD
junction of field 786 (only secondary field 1821, 31 frames, covers
it); Proxima in the ~8′ gap between CCDs 9 and 10 of field 808 (only
the low-z half of the locus, z ≲ 1,900 AU, ever reaches a live
quadrant, and only at one parallax phase); Ross 128 87″ outside
quadrant (446, 9, 2). With ≈13 % of sky in primary-grid gaps and a
6′ locus, this will be common; the scale-up overlay must compute
per-corridor grid position.

## Precise pass (`precise_v1`)

600-pix mskimg cutouts (TPV WCS + ZSDS mask, fatal template 6141),
tolerance 1″, seed step 10″. 1,896 evaluations: 1,471 usable, 9
partial, 216 unusable, 200 unknown (product directory not served —
≈10 % of public metadata rows).

| endpoint / role | usable | disjoint z-cover | covered z [AU] | epochs |
|---|---|---|---|---|
| ross-128 rx / tx | 259 / 225 | 3 / 1 | 550–10,000 | 2018-06 → 2026-06 |
| eps-ind-a rx / tx | 10 / 14 | 0 | 550–10,000 | 2018-05 → 2020-10 |
| proxima-cen rx / tx | 496 / 476 | 15 / 20 | 550–1,434 / 550–1,913 | 2018-06 → 2026-02 |

## Catalog screening (`screen_v1`)

Per-exposure `psfcat` matched within 10″ of each locus; 7,786
ScreenMatch records (Proxima 6,660, Ross 128 1,104, ε Ind 22); DR24
objects snapshotted per corridor (Proxima 43,623 within 0.2°).
Bin-occupancy saturates (a 10″ strip in a 1″ survey always contains
stars). Fixed-z point filter (R = 2″): Ross 128 tx best z = 892 AU
resolved as a static background star (sky rms 0.2″ vs track rms 0.3″);
Proxima rx/tx best z ≈ 557–585 AU with 13–14 of 65 visits, support
spread along the track (rms ≈ 62–137″) — many different field stars,
not one follower; the stage-2 stack at those cells is below threshold.
In ZTF the difference image removes static stars at the pixel level,
so layer-1 screening is a veto input, not the discriminator.

## Forced photometry, stack, calibration (`calib_v1`, AnalysisRun `run-1f02c325fd66`)

- 797 sci+diff+msk cutout sets (787 with difference images). Search
  image is the pixelwise hybrid (difference image inside the reference
  footprint, sky-subtracted science image in the edge strip outside
  it); the matched-filter weight fraction from difference pixels is
  carried per cell (`diff_image_weight_frac`).
- Per-pixel variance model validated: S/N histogram N(0,1) to 2 % at
  1σ and 2.33σ in both regimes.
- 192-node 1/z-uniform grid × 5×5 µ grid (|µ| ≤ 1″/yr, T0 = MJD
  59800), 8 offset controls, weight cap 20× median, **per-cell epoch
  floor N ≥ 5** (new: gap geometry leaves cells with 1–2 epochs where
  a single residual sets S), **single-epoch clip |S_e| ≤ 5** (new:
  defines the stack cell as sub-threshold sources; brighter sources are
  the screening layer's cell; bright edge ≈ 20.1–20.4 AB recorded on
  every constraint).
- Thresholds T = 4.3–7.9; every real-track maximum ≤ T after the clip.
  Two earlier exceedances (Proxima rx g, tx g) were phase-vetoed; the
  final run has **0 candidates**. Proxima has single-phase coverage
  only (phase_n ≈ {38, 1}) — a Proxima candidate could not be
  adjudicated by the phase test at all; Ross 128 r has both phases
  (28/25, 16/23).
- **Positive control** (`control_v1`, asteroid (60000), V 19.1–20.2):
  22 frames; stack peak at exactly the Horizons position (0″, 0″),
  S = 28 (g) / 41 (r) vs T = 3.0 / 2.3; 7 of 10 (g) and 11 of 12 (r)
  single frames > 5σ. Recovered magnitude was 1.8 (r) / 2.8 (g) mag
  too faint on the first run — the shared matched filter returned the
  Gaussian *peak amplitude*, not total flux (2.5 log10 2πσ²_pix). Fixed
  in `sglsurvey/photometry.py` (unit-sum kernel); residual offset after
  solar colours +0.41 (g) / +0.59 (r) mag = Gaussian-vs-Moffat aperture
  loss, applied as a 0.5 mag throughput correction to every depth.
- **80 Constraint records** (59 recovery curves, 21 not constrainable):
  m90 (90 % recovery, duty ≥ 0.5, |µ| ≤ 1″/yr, control-corrected):

| endpoint / role | band | N | T | m90 median / best / worst (AB) |
|---|---|---|---|---|
| ross-128 rx | g / r | 130 / 129 | 7.0 / 7.5 | 21.39 / 21.70 / 20.54 — 21.30 / 21.60 / 20.67 |
| ross-128 tx | g / r | 114 / 111 | 6.9 / 7.9 | 21.25 / 21.71 / 19.93 — 21.14 / 21.78 / 20.40 |
| proxima-cen rx | g / r | 232 / 264 | 5.1 / 5.7 | 21.47 / 22.04 / 17.56 — 21.36 / 22.21 / 18.21 (z ≤ 1,050 AU only) |
| proxima-cen tx | g / r | 226 / 250 | 4.5 / 6.3 | 21.19 / 22.21 / 17.31 — 20.80 / 21.74 / 18.33 (z ≤ 1,350 AU only) |
| eps-ind-a rx / tx | r | 6 / 10 | 4.3 / 4.4 | 19.01 / 19.69 / 18.42 — 19.44 / 20.04 / 17.41 (z ≥ 1,350 AU only) |

Depth is set by T (faint static-star confusion in the science-image
regime and residuals in the difference regime), not by photon noise:
the Gaussian-noise stack limit would be ≈ 24 AB.

- Reflected-light conversion (albedo 0.1, asteroid H–D convention) at
  the mid-z of each interval: m90 ≈ 21.3 excludes only D ≳ 7×10⁴ km at
  550 AU, rising to ≳10⁷ km at 3,000 AU — optical limits constrain
  planet-scale reflectors, as anticipated in plan §5.

## Lessons for the scale-up

1. Compute primary-grid pixel position per corridor in the overlay;
   weight secondary-grid fields; expect single-phase coverage for
   gap-adjacent corridors and record it.
2. Reference-image strips: the hybrid search image works but the
   science-image regime is confusion-limited; a custom per-corridor
   reference (median of the 8-year cutout stack) would extend
   difference imaging to the edge strip.
3. The psfcat screening is cheap (3–4 exposures/s) but not
   discriminating at 1″; keep it as the bright-source layer.
4. The asteroid control must stay in the chain: it caught a factor-of-
   2πσ² flux-scale error in the shared estimator that three WISE
   batches had not.
