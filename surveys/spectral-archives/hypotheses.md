# Spectral-archive family hypothesis freeze v1.0 (Pipeline B, channel A)

**FROZEN 2026-09-05**: the user approved the six recon recommendations
("Looks good, proceed") — this document encodes them as decisions
D1–D9 before any spectrum is analysed. **Amendment v1.1** (dev-driven
gates A1–A8, plus A9 found in confirmatory adjudication) is recorded in
`thresholds.md` and `configs/threshold_freeze_v1_1.json`; no hypothesis
changed. Survey complete the same day: `report/spectral_archives.md`. Drafted after the reachability
recon (`notes/spectral_archives_recon_2026-09-05.md`; scan
`results/recon_scan_v0.json`, rows `results/recon_rows_v0.json`) and
three verification downloads (one NIRPS, one HARPS, one SPIRou
in-window spectrum, headers read, no statistic formed). First survey of
the spectral-archive family (plan §3.6 row 10, §5.15 item O4); the
crossings discipline (frozen hypotheses → dev → blind confirmatory →
adjudication → injections → report) is inherited from the ZTF/TESS/
ATLAS/LASCO/HI-1 surveys; substrate-specific substitutions are marked
**[SPEC]**.

Input: `crossings/universal_v1` (`xng-a09e2db7681d`, Earth-center
observer — valid: every spectrograph here is on the ground), channel A
(`link_direction = inbound`, `axis_distance_au > 0`), the seven targets
with b_min ≤ 0.1 AU. Windows: flat chords from b_min and v⊥ at the
three rungs (1.2 R☉, 2.5 R☉, 0.1 AU), as in every crossings survey.

## 1. Physical hypothesis **[SPEC]**

The target system transmits a **narrow optical/near-IR line** (laser;
width ≤ the instrumental resolution element) toward the Sun's focal
line. Earth intercepts the beam while it sits within b_rung of the
Sun–star axis on the star side (star at opposition). In a spectrum of
the star taken inside the window the line appears as an unresolved
emission feature **at rest in the star's barycentric frame** (the
systemic RV is known to ≪ 1 km/s; a ±50 km/s tolerance covers a
transmitter in orbit), superposed on the photospheric spectrum. Cells:

| cell | vacuum wavelength interval | rationale |
|---|---|---|
| 532 | 515–545 nm | frequency-doubled Nd/Yb solid-state and fibre lasers |
| 1064 | 1030–1090 nm | Nd:YAG/YVO₄/YLF, Yb fibre (1030–1090) |
| 1550 | 1530–1570 nm | telecom C-band, Er fibre |
| generic | the instrument's full usable band | any line |

Cells are evaluated where the instrument's band contains them.
Continuous (d = 1) emission through the exposure is the primary
hypothesis; pulses are integrated by the 90–2700 s exposures and are
carried as a duty-cycle × peak-power reinterpretation (declared). The
photospheric flux in one resolution element is the noise floor, so
thresholds scale with the star's brightness (bright stars give weaker
limits — the opposite of the imaging surveys).

## 2. Substrate and archives (D1) **[SPEC]**

Reduced 1D (or order-by-order) spectra with sub-minute time stamps,
fetched anonymously:

- **ESO phase 3** (`ivoa.ObsCore`, `dataportal.eso.org/dataportal_new/file/<ADP>`):
  HARPS (378–691 nm, R 115k), ESPRESSO (378–789 nm, R 140k), NIRPS
  (966–1923 nm, R 70k; `FLUX_TELL_CAL_SKYSUB` telluric-corrected,
  sky-subtracted, absolutely flux-calibrated), X-shooter (per-arm
  UVB/VIS/NIR, R 5–10k). Wavelengths barycentric (`SPECSYS='BARYCENT'`),
  Å; HARPS/ESPRESSO in air → converted to vacuum (Morton 2000).
- **CADC/CFHT SPIRou APERO `t` product** (`ws.cadc-ccda…/data/pub/CFHT/<odometer>t.fits`):
  49 orders × 4088 px, FluxAB/WaveAB/BlazeAB + `Recon` (telluric) and
  `OHLine` (sky) models; wavelengths vacuum nm in the observer frame,
  `BERV` (km/s, barycorrpy) in the extension header; 955–2515 nm,
  R 73k. Barycentric shift λ_bary = λ_obs (1 + BERV/c), sign verified
  at dev machinery (§10).
- **CARMENES GTO DR1 VIS** (`carmenes.cab.inta-csic.es/gto/getDR1DataPublic.action?id=<Karmn>_VIS.zip`):
  caracal per-order spectra, vacuum wavelengths in the observer frame,
  BERV in header; 520–960 nm, R 94.6k; the DR1 CSV is the exposure
  list (TIMEID = UT start).
- **SOPHIE** (OHP, s1d if served anonymously): the wolf-359 2022-03-03
  unit is included iff the archive serves the s1d and its header UT
  puts the exposure inside the 2.5 R☉ window — decided at fetch, before
  any statistic (recorded in the machinery log).

Excluded from v1 (hand-off list, `report/spectral_archives.md`): raw
frames that need échelle extraction — KOA HIRES (wolf-359 2010 1.2 R☉),
Gemini MAROON-X (teegarden 2021, ross-128 2024 2.5 R☉), CRIRES 2009
(K band), BL APF (teegarden 2016); the proprietary SPIRou teegarden
2025 set (pre-registered re-run when public, ~2026-11); FEROS (< 20
nulls); EFOSC/LRIS/LAMOST low-resolution frames; every 0.1 AU spectrum
of instruments without a reduced-product route.

## 3. Unit population (D2)

**Unit = target × instrument × event**, over in-window public spectra
of the D1 products, restricted to star × instrument combinations with
**≥ 20 out-of-window public spectra** of the same product (the null
pool). Membership is fixed from `results/recon_rows_v0.json` by
`scripts/select_units.py` → `results/units_v1.json` (seed 20260905),
before any download. Expected from the recon (per-minute
deduplicated; exact counts in `units_v1.json`):

| target | instrument | in-window spectra | events | grazing (≤ 2.5 R☉) | null pool |
|---|---|---|---|---|---|
| wolf-359 | SPIRou | 24 | 3 | 4 (1.2 R☉, 2021) | 675 |
| wolf-359 | NIRPS | 8 | 1 | 4 (1.2 R☉, 2025) | 219 |
| wolf-359 | CARMENES-VIS | 9 | 3 | 0 | 70 |
| wolf-359 | HARPS | 1 | 1 | 0 | 107 |
| wolf-359 | X-shooter | 3 | 1 | 0 | 25 |
| ross-128 | SPIRou | 40 | 1 | 4 (2.5 R☉, 2022) | 216 |
| ross-128 | HARPS | 28 | 6 | 4 (2.5 R☉, 2021) | 279 |
| ross-128 | ESPRESSO | 8 | 4 | 0 | 77 |
| ross-128 | CARMENES-VIS | 3 | 2 | 0 | 55 |
| teegarden | CARMENES-VIS | 12 | 4 | 2 (2.5 R☉, 2017) | 248 |
| teegarden | ESPRESSO | 4 | 1 | 0 | 32 |
| ross-154 | HARPS | 6 | 1 | — | 120 |
| ross-154 | CARMENES-VIS | 3 | 1 | — | 53 |
| gj-908 | HARPS | 3 | 1 | — | 84 |
| gj-908 | X-shooter | 46 files | 1 | — | 46 files |

Rungs per unit = the union of the rungs of its spectra; a unit at a
grazing rung is also a 0.1 AU unit. van-maanen and gj-1276 have no
unit (recorded as uncovered; van-maanen → future-observation
recommendation).

## 4. Null ensemble (D3) **[SPEC]**

For each star × instrument combination the **null ensemble** is a
seeded random draw (seed 20260905, stratified by year) of **N_null = 60**
out-of-window public spectra (all of them where fewer than 60 exist).
The same star, the same instrument, the same product: photospheric
lines, flares, telluric residuals, detector defects and reduction
artefacts are all represented. Null spectra double as the **template**
source (leave-one-out median, §6). This replaces the imaging surveys'
8 control patches: the empirical false-alarm rate per trial is
1/(N_null + 1) instead of 1/9 (§7).

## 5. Common grid and normalisation (D4)

Per instrument: log-λ grid, vacuum, barycentric, step Δln λ = 1/(3 R)
(≈ 1 pixel). Each spectrum is linearly resampled; pixels with QUAL ≠ 0,
non-finite flux, non-positive error, or (SPIRou) `Recon` < 0.3 or
`OHLine` above the 99th percentile of the order are masked. Each
spectrum is divided by its **running median over 201 grid pixels**
(pseudo-continuum: removes throughput, blaze residuals and airmass
curvature, keeps features narrower than ~70 resolution elements). Line
amplitudes are therefore in units of the local pseudo-continuum.

## 6. Detection construction (D5) **[SPEC]**

For a spectrum s (normalised): template T = pixel-wise median of the
null ensemble (excluding s if s is a null — leave-one-out); residual
r = s − T; per-pixel scale σ = max(MAD-based scatter of the null
ensemble residuals at that pixel, propagated photon error); z = r/σ.
**Matched filter**: Gaussian kernel of the instrumental FWHM
(2.5 grid px for R-sampled grids: HARPS/ESPRESSO/NIRPS/SPIRou/CARMENES;
X-shooter 3 px), S(λ) = Σ z·k / √Σ k² over the kernel footprint.

Statistics per unit and cell:

1. **S_line** (primary) = max over the unit's in-window spectra of
   max_λ∈cell S(λ). Null value per null spectrum: the same quantity.
   T_line = max over the null ensemble.
2. **S_coadd** (secondary; units with ≥ 2 in-window spectra) =
   max_λ∈cell S(λ) evaluated on the **mean** of the unit's in-window
   normalised spectra (z formed with σ/√n). Null: 60 seeded random
   draws of n null spectra, same construction; T_coadd = max over the
   draws.

Exceedance: S > T. A unit with S ≤ T yields a constraint; S > T enters
adjudication (§8). Trials = Σ_units Σ_cells (1 + [n ≥ 2]).

## 7. Error-rate bookkeeping

Under the null each trial exceeds T with probability 1/(N_null + 1)
(= 1/61 for full ensembles), so the expected number of exceedances is
Σ_trials 1/(N_null,unit + 1). The dev family establishes machinery
behaviour; the confirmatory family is analysed **once**.

## 8. Adjudication ladder (D6; frozen order)

For every exceedance, at the peak wavelength λ_p (barycentric vacuum):

1. **Shape**: fit a Gaussian to r around λ_p; FWHM < 0.6× or > 2.0× the
   instrumental FWHM → `defect/cosmic` (single-pixel spikes, order
   edges, bad columns). Also `defect` if the peak pixel is masked in
   ≥ 20 % of null spectra.
2. **Stellar activity**: |λ_p − λ_line(1 + RV_star/c)| ≤ 150 km/s for
   any line in the frozen M-dwarf emission list (H Balmer α–δ, Ca II
   H/K and IRT, He I D3 587.6 and 1083.0 nm, Na I D, K I 766.5/769.9,
   Mg I b, Fe II 501.8/516.9, Paschen β/γ/δ 1282/1094/1005 nm, Brackett
   γ 2166 nm, Pfund lines, Ca I 422.7) → `stellar_flare` (M-dwarf
   flares fill these lines; the same-star null ensemble usually shows
   them too).
3. **Sky/telluric**: |λ_p − λ_OH| ≤ 1 resolution element in the
   **observer** frame for the Rousselot et al. (2000) OH list, or the
   pixel inside a telluric band where the null-ensemble σ exceeds 3×
   its median → `telluric`.
4. **Recurrence**: the same barycentric λ_p (± 1 resolution element)
   exceeds T in another in-window spectrum of the same unit and in
   none of the nulls → `retained-ambiguous` if 1–3 alone cannot explain
   it; a single-spectrum feature that passes 1–3 → `retained-ambiguous`
   with the held-out prediction (the next window) as the designated
   follow-up.
5. Anything surviving 1–4 with recurrence → **candidate**, subject to
   plan §3.4 (held-out epochs).

## 9. Injections and power conversion (D7) **[SPEC]**

Completeness: Gaussian lines of the instrumental FWHM are injected
into copies of each in-window spectrum at random λ within each cell
(50 per amplitude step; amplitudes A ∈ {0.02 … 2} × local
pseudo-continuum, log-spaced, 12 steps); recovered iff S_line > T_line
of the unit. A90 = the amplitude at 90 % recovery (linear
interpolation). Line flux F = 1.0645 × A90 × F_λ,cont(λ) × FWHM_λ, with
F_λ,cont from the star's photometry (Gaia G/BP/RP, 2MASS J/H/Ks zero
points, log-log interpolated to λ; NIRPS's absolute `FLUX_CAL` is the
cross-check for wolf-359). **Power through the rung cone**
P = F × π b_rung² (flat-top footprint of radius b_rung, the convention
of every crossings report: 1.2 R☉ 8.35e8 m, 2.5 R☉ 1.74e9 m,
0.1 AU 1.496e10 m), reported per unit × cell × rung. Declared: no
pulse cell; duty cycle d re-scales P by 1/d.

## 10. Machinery checks before the dev search (D8)

- BERV sign (SPIRou, CARMENES): the residual between two null
  spectra with |ΔBERV| > 20 km/s must be minimised by the applied
  shift, not its opposite (cross-correlation of the normalised spectra).
- Resampling: the template's residual scatter vs photon error per
  instrument (report the ratio; expected 1–3 at telluric bands).
- Wavelength cells: fraction of masked pixels per cell per unit.
Recorded in `notes/dev_machinery_log_2026-09-05.md`.

## 11. Dev / confirmatory split (D9; seed 20260905)

- **Dev**: ross-154 (HARPS 6, CARMENES 3 — flare star, the
  `stellar_flare` rule test) and gj-908 (HARPS 3, X-shooter 46 files —
  bright star, the multi-arm loader) → 4 units.
- **Confirmatory** (analysed once after dev, amendments frozen as
  v1.x before any confirmatory spectrum is normalised): every
  wolf-359, ross-128 and teegarden unit (SPIRou, NIRPS, HARPS, ESPRESSO,
  CARMENES-VIS, X-shooter; SOPHIE if §2 admits it).

## 12. Decisions (all approved 2026-09-05 with the recon recommendations)

| D | decision |
|---|---|
| D1 | substrate = ESO phase 3, CADC SPIRou `t`, CARMENES DR1 VIS (+ SOPHIE if served); raw-frame hits → hand-off |
| D2 | unit = target × instrument × event; combos need ≥ 20 out-of-window spectra |
| D3 | null ensemble = 60 seeded same-star same-instrument out-of-window spectra; leave-one-out template |
| D4 | log-λ grid at 1/(3R), running-median (201 px) pseudo-continuum |
| D5 | matched-filter S_line (primary) and S_coadd (secondary); T = max over nulls / 60 coadd draws |
| D6 | adjudication ladder shape → stellar line → sky/telluric → recurrence |
| D7 | injections at instrumental FWHM, A90, photometric F_λ, P = F π b² |
| D8 | BERV-sign and resampling checks before dev |
| D9 | dev = ross-154 + gj-908 units; confirmatory = wolf-359, ross-128, teegarden units, once |
