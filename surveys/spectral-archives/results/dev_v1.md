# Spectral-archive family — dev results (v1.1, run 2026-09-05T15:33 UTC; frozen config threshold_freeze_v1_1_sha.txt sha 020baba13ecd)

Units 4 (4 searched), trials 16, exceedances 6 vs 1.61 expected; dispositions {'telluric': 1, 'retained-ambiguous': 5}.

## Combos (null ensembles)

| combo | nulls used | T_line per cell | null S median / p90 | σ/photon | stellar mask | status |
|---|---|---|---|---|---|---|
| gj-908|HARPS | 60 | 532 19.3, generic 48.6 | 532 6.0/13.0, generic 17.5/30.1 | 532 1.00, generic 1.00 | 0.023 | ok |
| gj-908|XSHOOTER-NIR | 36 | 1064 4.3, 1550 3.9, generic 32.0 | 1064 2.8/3.4, 1550 2.8/3.6, generic 5.5/15.1 | 1064 1.00, 1550 1.00, generic 1.00 | 0.006 | ok |
| ross-154|CARMENES-VIS | 51 | generic 69.5 | generic 35.4/60.4 | generic 1.19 | 0.015 | ok |
| ross-154|HARPS | 58 | 532 15.6, generic 38.5 | 532 4.9/10.0, generic 17.6/27.4 | 532 1.00, generic 1.00 | 0.023 | ok |

## Units

| unit | b (R☉) | rungs | usable / n | S_line (S / T per cell) | S_coadd (S / T per cell) | exceed. | expected |
|---|---|---|---|---|---|---|---|
| gj-908|HARPS|2011-09-21 | 12.32 | 0.1AU | 3 / 3 | 532 4.9/19.3; generic 34.3/48.6 | 532 3.8/14.0; generic 20.5/45.2 | 0 | 0.13 |
| gj-908|XSHOOTER-NIR|2010-09-21 | 12.32 | 0.1AU | 22 / 36 | 1064 5.3/4.3**; 1550 8.1/3.9**; generic 16.8/32.0 | 1064 6.6/6.2**; 1550 28.1/7.8**; generic 60.9/33.2** | 5 | 1.19 |
| ross-154|CARMENES-VIS|2018-07-03 | 3.40 | 0.1AU | 3 / 3 | generic 57.5/69.5 | generic 98.7/92.7** | 1 | 0.07 |
| ross-154|HARPS|2017-07-03 | 3.40 | 0.1AU | 6 / 6 | 532 14.5/15.6; generic 17.0/38.5 | 532 7.8/10.3; generic 10.7/24.9 | 0 | 0.22 |

## Exceedance ledger (automatic ladder; manual notes in the report)

| unit | statistic | cell | S / T | λ_bary (nm) | spectrum (UTC) | FWHM ratio | ladder evidence | disposition |
|---|---|---|---|---|---|---|---|---|
| gj-908|XSHOOTER-NIR|2010-09-21 | S_line | 1064 | 5.3 / 4.3 | 1058.048 | ADP.2014-05-16T19:19:01.763.fits 2010-09-18T04:40 | 0.97 | recurs in 0 other spectra | retained-ambiguous |
| gj-908|XSHOOTER-NIR|2010-09-21 | S_line | 1550 | 8.1 / 3.9 | 1538.572 | ADP.2014-05-16T19:19:01.763.fits 2010-09-18T04:40 | 0.79 | recurs in 18 other spectra | retained-ambiguous |
| gj-908|XSHOOTER-NIR|2010-09-21 | S_coadd | 1064 | 6.6 / 6.2 | 1058.048 | ADP.2014-05-16T19:19:01.763.fits  | 0.97 | recurs in 0 other spectra | retained-ambiguous |
| gj-908|XSHOOTER-NIR|2010-09-21 | S_coadd | 1550 | 28.1 / 7.8 | 1538.572 | ADP.2014-05-16T19:19:01.763.fits  | 0.79 | recurs in 2 other spectra | retained-ambiguous |
| gj-908|XSHOOTER-NIR|2010-09-21 | S_coadd | generic | 60.9 / 33.2 | 1985.468 | ADP.2014-05-16T19:19:03.190.fits  | 0.95 | recurs in 0 other spectra | retained-ambiguous |
| ross-154|CARMENES-VIS|2018-07-03 | S_coadd | generic | 98.7 / 92.7 | 760.026 | car-20180701T01h21m33s-sci-gtoc-vi  | 1.53 | σ-band ×20.7 | telluric |

## Completeness and power limits (A90 = 90 % recovery amplitude in units of the local pseudo-continuum; P through the rung cone at the cell reference wavelength; range over the cell in brackets)

| unit | cell | T_line | A90 best / median | λ_ref (nm) | F_λ (erg s⁻¹ cm⁻² Å⁻¹) | P per rung (W) |
|---|---|---|---|---|---|---|
| gj-908|HARPS|2011-09-21 | 532 | 19.3 | 0.268 / 0.326 | 530 | 9.19e-13 | 0.1AU 8.5e+03 [7.7e+03–9.4e+03] |
| gj-908|HARPS|2011-09-21 | generic | 48.6 | 2.000 / 2.000 | 511 | 8.41e-13 | 0.1AU 5.59e+04 [4.1e+04–1.4e+05] |
| gj-908|XSHOOTER-NIR|2010-09-21 | 1064 | 4.3 | not recovered ≤ 2.0 | 1060 | 4.23e-12 | — |
| gj-908|XSHOOTER-NIR|2010-09-21 | 1550 | 3.9 | not recovered ≤ 2.0 | 1550 | 1.48e-12 | — |
| gj-908|XSHOOTER-NIR|2010-09-21 | generic | 32.0 | not recovered ≤ 2.0 | 1567 | 1.94e-12 | — |
| ross-154|CARMENES-VIS|2018-07-03 | generic | 69.5 | not recovered ≤ 2.0 | 707 | 7.24e-13 | — |
| ross-154|HARPS|2017-07-03 | 532 | 15.6 | 0.513 / 0.546 | 530 | 2.50e-13 | 0.1AU 4.41e+03 [3.7e+03–5.2e+03] |
| ross-154|HARPS|2017-07-03 | generic | 38.5 | 1.865 / 2.000 | 511 | 2.08e-13 | 0.1AU 1.29e+04 [9.5e+03–5.8e+04] |
