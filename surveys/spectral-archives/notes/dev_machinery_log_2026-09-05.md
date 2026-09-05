# Spectral-archive family — dev machinery log (2026-09-05)

Chronological record of the dev stage (thresholds §5/§7; hypotheses
D8/D9). Dev units: gj-908 HARPS 2011-09-21 (3 spectra), gj-908
X-shooter-NIR 2010-09-21 (36 products), ross-154 HARPS 2017-07-03 (6),
ross-154 CARMENES-VIS 2018-07-03 (3). Null ensembles 60/36/60/53 drawn
with seed 20260905 (`results/units_v1.json`). No confirmatory spectrum
was opened before the v1.1 amendment was frozen
(`results/threshold_freeze_v1_1_sha.txt`).

## Fetch

961 files wanted for both families (dev + confirmatory) — 706 HTTP
files (ESO phase 3, CADC SPIRou `t`) fetched with 0 failures, 4
CARMENES DR1 zips (ross-128 247 MB, ross-154 251 MB, teegarden 1.15 GB,
wolf-359 336 MB) with 255 members extracted; sha256 + size per file in
`runs/spectral-archives/v1/fetch_manifest.jsonl`. Two fetch defects
fixed before any analysis: the Karmn `+` must be percent-encoded
(the server returns an HTML page otherwise); zip members carry an
`_A` fibre suffix (matched by TIMEID stem). 12 GB on the runs SSD.

## D8 checks

- **BERV sign** (v1.1 form, stellar-region cross-correlation of the
  two nulls with the largest |ΔBERV|, three treatments):
  HARPS (already barycentric, `SPECSYS=BARYCENT`): none 60 px / applied
  0 / flipped 60 (gj-908, ΔBERV 59 km/s = 68 px) — consistent.
  CARMENES caracal (observer frame): applied brings the stellar
  features to 0 px for 51 of 52 nulls (per-spectrum table in the
  session record); `none` shows ΔBERV; **λ_bary = λ_obs (1 + BERV/c)
  confirmed**. The whole-band RMS form of the check (first
  implementation) is telluric-dominated and gave the opposite verdict —
  replaced (A5).
- **Wrong header BERV**: CARMENES `car-20170911T*` (ross-154) states
  +5.45 km/s inside a run at −27; its stellar lines sit 31 px off →
  the per-spectrum wavelength-solution gate (A1) drops it.
- **Resampling / σ**: σ/photon ratio 1.0 (HARPS, ESPRESSO, NIRPS),
  1.19 (CARMENES-VIS: telluric residuals in the red).
- **Cell mask fractions** (after A2): HARPS 532 cell 16 % (blue-order
  gaps + low SNR), generic 4–19 %; CARMENES generic 6 %.

## What set the null thresholds, pass by pass

| pass | T_line (HARPS gj-908 532 / generic; CARMENES ross-154 generic) | dominant cause |
|---|---|---|
| v1.0 | 64 / 2856; 102,607 | cosmic hits (1–3 px) in HARPS s1d, flare Hα/Ca lines in nulls, zero-flux blue orders normalised to ±10³ |
| + despike + stellar mask | 64 / 1588; 102,620 | 2-px cosmics survive the single-pixel despike; blue-edge blow-ups |
| + A3 PSF-consistent statistic | 24 / 1588; 10,060 | PSF-shaped blow-ups next to masked runs in zero-flux regions (HARPS frame with header SNR 0.95; CARMENES 520–540 nm at flux 0 ± 0.01) |
| + A2 spectrum/pixel SNR gate | 24 / 491; 208 | low-SNR spectra over-weighted: robust z scale 0.6–3.4, corr(scale, SNR) = −0.9 (σ floor used the ensemble-median error) |
| + A6 own-error floor (final) | 19 / 49; 70 | genuine residual floor: null S medians 6 (532) / 17 (generic) HARPS, 35 CARMENES generic |

## Final dev outcome (v1.1)

| unit | usable | S_line (cell = S / T) | S_coadd | exceedances |
|---|---|---|---|---|
| gj-908 HARPS 2011 | 3/3 | 532 4.9/19.3; generic 34.3/48.6 | 3.8/14.0; 20.5/45.2 | 0 |
| ross-154 HARPS 2017 | 6/6 | 532 14.5/15.6; generic 17.0/38.5 | 7.8/10.3; 10.7/24.9 | 0 |
| ross-154 CARMENES 2018 | 3/3 | generic 57.5/69.5 | 98.7/92.7 | 1 → `telluric` (760.03 nm, O₂ A band, σ-ratio 20.7, FWHM ratio 1.53) |
| gj-908 X-shooter-NIR 2010 | 22/36 | 1064 5.3/4.3; 1550 8.1/3.9; generic 16.8/32.0 | 6.6/6.2; 28.1/7.8; 60.9/33.2 | 5 (recorded, not counted — A8) |

Counted dev: 10 trials, 1 exceedance vs 0.6 expected (A7 formula).
X-shooter: 4 phase-3 product variants per exposure (36 products = 9
exposures), 14/36 fail the wavelength gate (no BERV keyword; frames not
all in one system), and the single-night null pool (2010-08) is blind to
the night-to-night sky residual that appears in all 2010-09-18 products
at 1538.57 nm (OH 3-1 band region; not in the compact sky list) →
excluded from the family (A8).

## Notes for the confirmatory run and the report

- The generic-cell T of CARMENES-VIS (~70) is systematics-dominated
  (telluric residuals 900–960 nm); the 532 cell is not in CARMENES's
  band ([515, 545] vs 520–960). CARMENES units are generic-only.
- A `retained-ambiguous` disposition in confirmatory will need the
  manual OH/telluric look the compact list cannot give in the NIR
  (SPIRou `OHLine` p99 mask + σ-band rule are the automatic guards).
- Bookkeeping: units with many in-window spectra (ross-128 SPIRou 40,
  HARPS 18; wolf-359 SPIRou 20) have null exceedance probabilities of
  0.4/0.23/0.25 per S_line cell by construction (A7); the coadd
  statistic is the sharper test for them.
