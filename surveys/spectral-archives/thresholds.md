# Spectral-archive family threshold freeze v1.0

**FROZEN 2026-09-05** together with `hypotheses.md` (D1–D9; user
approval of the recon recommendations). Machine form
`configs/threshold_freeze_v1.json`, sha256 in
`results/threshold_freeze_v1_sha.txt`. Written before any spectrum is
resampled or any statistic formed; the three recon downloads (NIRPS
wolf-359 2025-03-03T05:20, HARPS ross-128 2021-03-17T02:26, SPIRou
wolf-359 2021-03-03T11:51) were opened for headers only and are
declared here.

## 1. Units

Unit = target × instrument × event, over in-window public spectra of
the reduced products (hypotheses §2), restricted to star × instrument
combinations with ≥ 20 out-of-window public spectra. Population fixed
by `scripts/select_units.py` from `results/recon_rows_v0.json` →
`results/units_v1.json` before any download. Cells per unit: {532,
1064, 1550, generic} ∩ band. Rungs per unit: union over its spectra.

## 2. Statistics

Grid: barycentric vacuum log-λ, Δln λ = 1/(3R). Normalisation: running
median over 201 px. Template: leave-one-out median of the unit's
null ensemble (60 seeded same-star same-instrument out-of-window
spectra; all when fewer). Scale σ: max(MAD scatter of null residuals
per pixel, photon error). Matched filter: Gaussian, instrumental FWHM
(2.5 px; 3.0 px X-shooter).

1. **S_line** (primary): max over the unit's in-window spectra of the
   cell maximum of S. T_line = max over the null ensemble of the same
   per-spectrum quantity. Trial per unit × cell.
2. **S_coadd** (secondary, n ≥ 2 in-window spectra): cell maximum of S
   on the mean spectrum. T_coadd = max over 60 seeded draws of n nulls.
   Trial per unit × cell.

Exceedance ⇔ S > T. Expected exceedances under the null =
Σ_trials 1/(N_null + 1).

## 3. Gates

- Spectrum usable: file parses; ≥ 70 % of the cell's grid pixels
  unmasked; header time inside the window (mid-exposure); SNR keyword
  or median photon SNR ≥ 5 in the cell.
- Unit searchable: ≥ 1 usable in-window spectrum and ≥ 20 usable nulls.
- Masks: QUAL ≠ 0, non-finite, SPIRou `Recon` < 0.3 or `OHLine` > p99
  of the order; ±3 px around grid gaps between orders.

## 4. Adjudication

Frozen ladder (hypotheses §8): shape → stellar_flare (list in the
config, ±150 km/s) → telluric/sky (OH list, observer frame; σ-band
rule) → recurrence. Dispositions: `defect/cosmic`, `stellar_flare`,
`telluric`, `retained-ambiguous`, `candidate`.

## 5. Split (seed 20260905)

Dev: ross-154 (HARPS, CARMENES-VIS), gj-908 (HARPS) — analysed first;
amendments, if any, are frozen as v1.x before any confirmatory
spectrum is normalised. Confirmatory: all wolf-359, ross-128,
teegarden units, analysed once.

## 6. Completeness

Injections per unit × cell (hypotheses §9): 12 amplitude steps × 50
random positions per in-window spectrum; A90 by linear interpolation;
F_λ,cont from Gaia DR3 + 2MASS photometry (Data Lab, snapshotted),
NIRPS `FLUX_CAL` cross-check; P = 1.0645 A90 F_λ FWHM_λ π b_rung².

## 7. Order of work

Freeze (sha) → `select_units.py` → `fetch.py` (dev, then
confirmatory; snapshots + sha256 of every file) → machinery checks
(BERV sign, resampling) → `search.py dev` → amendments → `search.py
confirmatory` once → `adjudicate` → `completeness.py` → report
`report/spectral_archives.md` → plan / history / learnings.

## Amendment v1.1 (2026-09-05, dev-driven; frozen before any confirmatory spectrum was normalised)

Machine form `configs/threshold_freeze_v1_1.json` (sha256 in
`results/threshold_freeze_v1_1_sha.txt`). Found on the four dev units
(gj-908 HARPS + X-shooter-NIR, ross-154 HARPS + CARMENES-VIS; three dev
passes, `runs/spectral-archives/v1/search_dev_v1_1*.log`; details in
`notes/dev_machinery_log_2026-09-05.md`). None changes a hypothesis;
all are gates or the literal implementation of frozen wording:

- **A1 wavelength-solution gate**: stellar-region cross-correlation
  shift vs the ensemble template ≤ 3 px for every spectrum (one
  CARMENES ross-154 frame, 2017-09-11, carries a wrong header BERV).
- **A2 SNR gate at spectrum and pixel granularity**: spectra with band
  median flux/error < 5 unusable (a 118-s HARPS frame, SNR 0.95);
  pixels whose local continuum/error < 5 masked (zero-flux blue orders
  of the M dwarfs, normalised to ±10³); pixels with ensemble-median
  SNR < 5 excluded.
- **A3 PSF-consistent statistic**: cell maxima taken over local
  maxima of S whose Gaussian fit has FWHM within [0.6, 2.0]× the
  instrumental value (1–3 px cosmic hits and normalisation blow-ups set
  T at 10³–10⁵ otherwise); isolated single-pixel despike retained.
- **A4 loaders**: X-shooter WAVE in nm; SPIRou Recon/OHLine optional;
  HARPS s1d has no error vector → photon-shaped error scaled to the
  header SNR.
- **A5** D8 BERV check in stellar-region form (whole-band RMS is
  telluric-dominated and misleading).
- **A6 σ floor** = the spectrum's own photon error (frozen wording),
  not the ensemble median.
- **A7 bookkeeping**: expected S_line exceedances n/(N+n) per cell.
- **A8** X-shooter excluded (dev finding; no confirmatory unit).
- The frozen ±150 km/s stellar-line list is masked out of every cell
  (declared unconstrained), for nulls and in-window alike.

Dev outcome under v1.1: 3 searched units (X-shooter recorded, not
counted), 10 trials, 1 exceedance (ross-154 CARMENES S_coadd generic at
760.03 nm, the O₂ A band, σ-ratio 20.7 → `telluric`), expected 0.6.
