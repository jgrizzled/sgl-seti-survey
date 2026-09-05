# Spectral-archive family — laser-line search in the crossing windows (channel A)

**0 candidates.** First survey of the spectral-archive family (plan
§3.6 row 10, §5.15 item O4): archived high-resolution spectra of the
seven deep-family stars taken **inside** their channel-A crossing
windows (Earth on the Sun–star axis, star at opposition), searched for
an unresolved emission line at rest in the star's frame. Substrate:
ESO phase-3 HARPS/ESPRESSO/NIRPS, CFHT SPIRou (APERO, via CADC) and
CARMENES DR1 VIS spectra of wolf-359, ross-128 and teegarden
(confirmatory) and ross-154, gj-908 (dev), 2006→2025, with 60-spectrum
same-star same-instrument out-of-window null ensembles. Blind
confirmatory: **24 units, 77 trials, 18 exceedances vs 4.7 expected
under the null, all adjudicated as sky/telluric (12), artificial-lamp
contamination of one spectrum (4) or a two-pixel cosmic hit (2); no
retained-ambiguous feature, no candidate.** The excess over the
expectation is a known weakness of the construction (observer-frame
airglow smeared across a barycentric-grid ensemble; §6). Injection
completeness gives the programme's **first line-SED constraints, and
the first at 1064 and 1550 nm**: through the **1.2 R☉ cone a
continuous 1550 nm laser ≳ 12 W and a 1064 nm laser ≳ 22 W** (wolf-359,
NIRPS, 2025-03-03; ≳ 190–220 W with SPIRou 2021-03-03); through the
2.5 R☉ cone ≳ 37 W at 532 nm (ross-128, HARPS, 2021-03-17), ≳ 53–95 W
at 1550/1064 nm (wolf-359), ≳ 400 W at 1550 nm (ross-128, SPIRou 2022);
through the 0.1 AU cone 0.1–7 kW at 532/1064/1550 nm across 13 units.
CARMENES-VIS units and the SPIRou/CARMENES generic cells are
systematics-limited beyond the 200 % injection cap (unconstrained at
v1). Docs `surveys/spectral-archives/` (hypotheses D1–D9, thresholds
v1.0 + amendment v1.1, dev machinery log, results tables, adjudication,
completeness); recon `notes/spectral_archives_recon_2026-09-05.md`.

## 1. Hypotheses and geometry

- **Signal**: a narrow (≤ resolution element) emission line from the
  target system, at rest in the star's barycentric frame (RV known;
  ±50 km/s orbital tolerance), continuous through the exposure
  (pulses integrated; duty cycle re-scales power). Cells: 532
  (515–545 nm), 1064 (1030–1090), 1550 (1530–1570), generic (whole
  band), evaluated where in band. Power convention as in every
  crossings report: P = F_line × π b_rung² (flat-top footprint of
  radius b_rung: 1.2 R☉, 2.5 R☉, 0.1 AU).
- **Geometry**: `crossings/universal_v1`, channel A (inbound, star
  side), b ≤ 0.1 AU: seven targets, one event per sidereal year on a
  fixed date, flat-chord windows ±5.8 d (0.1 AU), ±11–16 h (2.5 R☉),
  ±2.5–7.5 h (1.2 R☉; van-maanen, wolf-359, gj-1276, teegarden only).
  No elongation gate — opposition is the easiest geometry for a
  ground spectrograph.
- **Unit** = target × instrument × event over in-window public
  reduced spectra; star × instrument combinations need ≥ 20
  out-of-window public spectra. 29 units fixed before any download
  (`results/units_v1.json`, seed 20260905); van-maanen and gj-1276
  have none.
- **Frozen decisions** D1–D9 (`hypotheses.md`): substrate; unit;
  60-spectrum seeded null ensemble with leave-one-out median template;
  log-λ grid at 1/(3R) with 201-px running-median normalisation;
  matched-filter S_line (max over the unit's spectra) and S_coadd
  (mean spectrum) with T = max over the ensemble / 60 complement-
  template draws; adjudication ladder shape → stellar line → sky/
  telluric → recurrence; injections at the instrumental FWHM; dev =
  ross-154 + gj-908, confirmatory = wolf-359 + ross-128 + teegarden,
  once.

## 2. Chain

Recon (2026-09-05: 11 archives, 7,644 on-star spectra, 633/30/14 in
0.1 AU / 2.5 R☉ / 1.2 R☉ windows) → freeze v1.0 (user-approved) →
`select_units.py` → `fetch.py` (961 files, 12 GB, 0 failures, sha256
manifest) → D8 machinery checks → dev (three passes) → **amendment
v1.1** frozen before any confirmatory spectrum was normalised
(`thresholds.md` §Amendment; `notes/dev_machinery_log_2026-09-05.md`):
A1 wavelength-solution gate (one CARMENES frame with a wrong header
BERV), A2 the frozen SNR ≥ 5 gate at spectrum and pixel granularity
(zero-flux blue orders of M dwarfs normalised to ±10³), A3
PSF-consistent statistic (1–3 px cosmic hits and normalisation
blow-ups set T at 10³–10⁵ otherwise), A4 loader fixes, A5 D8 check on
stellar regions, A6 σ floor = each spectrum's own photon error, A7
exceedance bookkeeping n/(N+n), A8 X-shooter excluded (4 product
variants per exposure, single-night nulls, no BERV) → confirmatory
once → adjudication → completeness → this report. One machinery
defect surfaced *during* adjudication: **A9** — ESPRESSO phase-3 `WAVE`
is vacuum (HARPS is air); the v1.0 flag double-converted the ESPRESSO
grid by +84 km/s. The two ESPRESSO combos were re-run with the
corrected frame (the statistic is translation-invariant; the
first-pass file is kept as `results/confirmatory_search_v1_espresso_airflag.json`;
T_line changed by up to 12 %, S values by < 5 %, no exceedance
appeared or vanished). The
SOPHIE wolf-359 2022 unit was not admitted: the public s1d header is
date-stripped (BJD rounded to the day), so the frozen "header UT inside
the window" condition cannot be met (hand-off).

Dev outcome (v1.1): gj-908 HARPS 2011, ross-154 HARPS 2017, ross-154
CARMENES 2018 — 10 trials, 1 exceedance (ross-154 CARMENES S_coadd,
760.03 nm, O₂ A band → `telluric`) vs 0.6 expected;
`results/dev_v1.md`.

## 3. Results (blind confirmatory, 2026-09-05)

`results/confirmatory_v1.md`, `results/confirmatory_search_v1.json`,
adjudication `results/adjudication_confirmatory_v1.json`, evidence
`results/confirmatory_exceedance_evidence_v1.json`.

**Combos and thresholds** (null S medians / T_line): HARPS ross-128
532 cell 5.5/16.0, generic 19/43; ESPRESSO ross-128 4.4/14.9,
39/127; ESPRESSO teegarden (27 usable nulls) 4.0/8.5, 94/171;
NIRPS wolf-359 1064 3.3/7.3, 1550 3.5/7.5, generic 42/93; SPIRou
wolf-359 1064 6/46, 1550 11/127, generic 68/127; SPIRou ross-128
1064 5.9/56, 1550 6.1/26, generic 51/106; CARMENES generic 35–45 /
67–163. wolf-359 HARPS: only 18 usable nulls (V 13.5, short
exposures) → not searchable.

| unit | b (R☉) | rungs | spectra | outcome |
|---|---|---|---|---|
| wolf-359 NIRPS 2025-03-03 | 0.76 | 1.2 R☉ | 8 | 1550 S_line 8.4/7.5, S_coadd 17.4/9.0 → **telluric** (CO₂ line, transmission 0.79); generic coadd 62/58 → telluric (1.8 µm band); 1064 cell null (5.5/7.3, 10.3/10.8) |
| wolf-359 SPIRou 2021-03-03 | 0.77 | 1.2 R☉ | 20 | 1064 coadd 37.7/28.6 → telluric (OH-subtraction plateau, 2021-03-02 night); generic coadd 281/71 → telluric (1.35 µm band); 1550 null |
| wolf-359 SPIRou 2019-03-03 | 0.78 | 0.1 AU | 4 | null (all cells) |
| wolf-359 CARMENES 2016/17/18 | 0.79 | 0.1 AU | 1/7/1 | 2017 coadd 97.9/90.3 → telluric (O₂ A band); others null |
| ross-128 HARPS 2021-03-17 | 1.81 | 2.5 R☉ | 18 | 532 null (14.7/16.0; coadd 4.3/9.1); generic coadd 69.5/35.9 → telluric ([O I] 5577 airglow, in all 18) |
| ross-128 SPIRou 2022-03-17 | 1.81 | 2.5 R☉ | 40 | 1064/1550 null; generic coadd 139/127 → telluric (OH 1456.4 nm, model 208) |
| ross-128 HARPS 2006/13/14/15/19 | 1.80 | 0.1 AU | 2/2/2/3/1 | 2015 generic coadd 37.6/28.8 → telluric ([O I] 5577); others null |
| ross-128 ESPRESSO 2019/20/21/22 | 1.81 | 0.1 AU | 2/2/3/1 | 2020: 4 exceedances (532 20.3/14.9, generic 465/127, both statistics) → **defect/contamination**: the 2020-03-16 spectrum carries 31 peaks above threshold incl. the Hg I pentad 404.8/436.0/546.2/577.1/579.2 nm (vacuum; the +84 km/s offset in the first pass exposed A9) — lamp light in the science fibre; the 2020-03-22 spectrum is clean (532 S 6.0/14.9, generic 44/127); others null |
| ross-128 CARMENES 2016/17 | 1.81 | 0.1 AU | 1/2 | null |
| teegarden CARMENES 2017-11-08 | 1.12 | 2.5 R☉ | 7 | S_line 83/67 and coadd 93/73 → **telluric**: OH 7-3 8827.10 Å (882.952 nm) and OH 9-4 7750.64 Å (775.277 nm) airglow; the same spectra carry the OH 8-3/9-4/6-2/7-3 forest at z 9–50; amplitude anticorrelates with SNR (additive sky) |
| teegarden CARMENES 2016/18/19 | 1.12 | 0.1 AU | 1/3/1 | 2019 S_line 82/67 → telluric (the same OH 8827 line, same observer wavelength); others null |
| teegarden ESPRESSO 2024-11-08 | 1.08 | 0.1 AU | 4 | 532 S_line 37.9/8.5, coadd 31.9/7.7 → **defect/cosmic** (2-px spike, z 23/29/8, one spectrum, absent in the other three); generic null |

**Tally**: 77 trials, 18 exceedances (12 sky/telluric, 4 lamp
contamination = one spectrum × two cells × two statistics, 2 cosmic =
one feature × two statistics) vs 4.7 expected; **0 retained-ambiguous,
0 candidates.** The only feature that briefly looked like the
hypothesis — a weak PSF-shaped excess at 1569.499 nm present in all
eight wolf-359 NIRPS spectra of 2025-03-02/03 and in none of 60 nulls
— sits exactly on a CO₂ telluric line (NIRPS `ATM_TRANSM` = 0.79; the
CO₂ series 1569.01/.25/.50/.75/1570.01 nm) and is absent at the same
barycentric wavelength in the 24 SPIRou in-window spectra of 2019 and
2021 and in the SPIRou nulls: a telluric-correction residual of those
two nights.

## 4. Completeness and power limits

Injections (12 log-spaced amplitudes 0.02–2.0 × local pseudo-continuum,
50 positions each, per in-window spectrum, recovered iff the
PSF-consistent S_line > T_line); A90 = best spectrum of the unit
(continuous emission is recovered if any spectrum recovers it);
FWHM_λ = λ/R; P = 1.0645 A90 F_λ FWHM_λ π b² (F_λ source below the
table).
`results/completeness_confirmatory_v1.json`; full table in
`results/confirmatory_v1.md`. Ranges in brackets = variation of F_λ
across the cell.

| unit | cell | A90 | P through 1.2 R☉ | 2.5 R☉ | 0.1 AU |
|---|---|---|---|---|---|
| wolf-359 NIRPS 2025 | 1550 | 0.066 | **12 W** | 53 W | 3.9 kW |
| wolf-359 NIRPS 2025 | 1064 | 0.10 | **22 W** | 95 W | 7.0 kW |
| wolf-359 NIRPS 2025 | generic (966–1923) | 1.62 | 203 W | 882 W | 65 kW |
| wolf-359 SPIRou 2021 | 1064 | 0.93 | 193 W | 839 W | 62 kW |
| wolf-359 SPIRou 2021 | 1550 | 1.24 | 221 W | 960 W | 71 kW |
| wolf-359 SPIRou 2019 | 1064 / 1550 | 1.25 / 1.26 | — | — | 83 / 72 kW |
| ross-128 HARPS 2021 | 532 | 0.56 | — | **37 W** | 2.8 kW |
| ross-128 SPIRou 2022 | 1550 | 0.24 | — | 404 W | 30 kW |
| ross-128 SPIRou 2022 | 1064 | 1.02 | — | 2.9 kW | 214 kW |
| ross-128 SPIRou 2022 | generic | 1.82 | — | 4.2 kW | 309 kW |
| ross-128 HARPS 2006–2019 (5 units) | 532 | 0.55–0.93 | — | — | 2.7–4.6 kW |
| ross-128 ESPRESSO 2019–2022 (4 units) | 532 | 0.21–0.35 | — | — | 0.87–1.4 kW |
| teegarden ESPRESSO 2024 | 532 | 0.98 | — | — | **129 W** |
| all CARMENES-VIS units; SPIRou/HARPS/ESPRESSO generic cells | — | > 2.0 | unconstrained at v1 | | |

F_λ: for the NIR cells of wolf-359 the **measured NIRPS absolute SED**
(telluric-corrected `FLUX_TELL_CAL`, 10-nm bins outside absorption
bands, checked against 2MASS J/H to 0.12–0.17 mag — a log-log
interpolation of Gaia+2MASS photometry under-estimates the 1.0–1.1 µm
peak of this M6 dwarf by 2.3×); for ross-128's NIR cells the
photometric interpolation times the wolf-359 ratio (2.7 at 1064 nm,
1.5 at 1550 nm), a conservative, limit-weakening choice; for the 532
cells the BP–G interpolation with a declared ±30 % systematic.

Reading: the star's photospheric flux in one resolution element is
the noise floor, so faint stars give deeper limits (teegarden ESPRESSO
107 W at 0.1 AU vs ross-128 1 kW), and telluric-corrected NIR products
with sky fibres (NIRPS) beat products without (CARMENES). These are
the programme's first constraints on a line SED: broadband photometry
cannot see a 10 W line against an M dwarf at any depth (10⁻⁹ of the
bolometric flux). Declared unconstrained: lines coincident with the
frozen stellar emission-line list (±150 km/s), pulses shorter than the
exposure (integrated), the CARMENES generic cell, every generic cell
beyond 200 % of the local continuum, and wavelengths inside the
per-pixel masks (order gaps, telluric bands with transmission < 0.3,
strongest 1 % OH pixels per SPIRou order).

Dev completeness (`results/completeness_dev_v1.json`): gj-908 HARPS
532 A90 0.27 → 8.5 kW (0.1 AU); ross-154 HARPS 532 0.51 → 4.4 kW
(both bright stars; photometric F_λ).

## 5. What this survey adds

- The **first laser-line (line-SED) constraints in the programme**,
  and the first at 1064 and 1550 nm — cells structurally open in every
  imaging survey (`report/joint_crossings.md` §3): 12 W (1550 nm) and
  22 W (1064 nm) through the 1.2 R☉ cone (wolf-359), 37–95 W through
  2.5 R☉ cones on two targets, 0.1–7 kW through 0.1 AU cones on 13
  units.
- The **archival on-star channel is deep and cheap**: reduced 1D
  products from four archives, no WCS/PSF/mask machinery, a
  same-star null ensemble for free. The 2010 HIRES, 2021/2024 MAROON-X,
  2009 CRIRES and 2016 APF grazing-window frames (raw) and the 2025
  SPIRou teegarden 1.2 R☉ set (proprietary to ~2026-11) are the
  next-cheapest additions.
- The Hg-lamp contamination of ESPRESSO ross-128 2020-03-16 and the
  airglow-dominated CARMENES statistics are archive facts worth
  recording for any laser-line search on these products.

## 6. Open items and hand-offs

- **v2 design notes** (`notes/learnings.md` §14): (i) an
  observer-frame sky-emission mask (OH Meinel bands, [O I], Na, Hg
  lamps) applied before the barycentric resampling — the barycentric
  ensemble smears observer-frame lines over ±30 km/s so the per-pixel
  σ does not see them, which is why 12 of 18 exceedances were sky;
  (ii) a spectrum-level contamination veto (count of PSF-consistent
  peaks above T across the band; > 10 → contaminated); (iii) for
  CARMENES a sky model or the NIR channel with its OH-corrected
  reduction; (iv) NIRPS/SPIRou telluric-transmission mask at
  T < 0.9 for the 1550 cell (CO₂ band); (v) 2-px despike.
- **Hand-offs** (raw frames / non-public): KOA HIRES wolf-359
  2010-03-03 (1.2 R☉, I₂ in; Tellis & Marcy 2017 may already carry a
  per-spectrum null — cross-match their table), Gemini MAROON-X
  teegarden 2021-11-08 and ross-128 2024-03-17 (2.5 R☉, +3.3 h /
  +1.4 h), CRIRES teegarden 2009-11-08 (K band, generic cell), BL APF
  teegarden 2016 (raw 2D, 0.1 AU), SOPHIE wolf-359 2022-03-04 night
  (2.5 R☉; hour requires the observatory), SPIRou teegarden 2025-11-08
  ×4 at +0.2 h (1.2 R☉; pre-registered re-run when public), HPF
  teegarden 2019/2020 (8 epochs at 0.1 AU, 810–1280 nm; data ask),
  CARMENES NIR twins of the 2017 grazing pair (Calar Alto form).
- **van-maanen** (deepest graze, b 0.32–0.46 R☉) has no échelle
  spectrum in any grazing window in any archive → a
  future-observation recommendation: one 20-min HARPS/ESPRESSO or
  NIRPS exposure at t_ca ± 7 h on Oct 6–7 of any year opens the
  deepest channel-A grazing cell of the programme. gj-1276 is
  uncovered.
- **Yearly refresh** (§5.7): new NIRPS/SPIRou/ESPRESSO public epochs
  (ESO 1 yr, CFHT 1 yr proprietary) land in the 2026 and 2027 windows;
  the ensembles and thresholds are cached
  (`runs/spectral-archives/v1/ensembles/`).
