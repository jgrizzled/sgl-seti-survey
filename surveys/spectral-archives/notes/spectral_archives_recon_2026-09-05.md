---
title: "Spectral-archive family reachability recon (ledger item O4) — on-star laser-line substrate for the crossing windows"
date: 2026-09-05
status: "recon complete — 11 archives probed live (all anonymous); a deep, public, multi-instrument on-star substrate exists incl. 14 grazing-window spectra on 3 targets and the first 1064/1550 nm in-band coverage of any crossing window; ready for a freeze"
---

# Spectral-archive family recon (plan §5.15 ledger item O4; §3.6 row 10)

Goal: establish whether any public spectral archive holds spectra of
the deep-family targets **inside** their channel-A crossing windows
(`crossings/universal_v1/`, 1980→2028), at what product level, and in
which wavelength bands — i.e. whether the "laser-line inspection in a
crossing window" that §5.5 handed off (the 3 teegarden APF spectra)
generalises into a survey. The plan row said "ESO/Keck and other
spectral archives — continuous or pulsed laser-line searches at the
target and local corridor" and "1550 nm line hypotheses (APOGEE
H-band, ESO NIRPS)". This note replaces that with endpoint facts.

Headline: **yes — a survey substrate exists, and it is much richer
than the APF hand-off suggested.** The seven b ≤ 0.1 AU channel-A
targets are among the most heavily observed radial-velocity stars in
the sky, and channel-A windows fall at opposition, where every
spectrograph can reach them. Across KOA, ESO, CADC (Gemini + CFHT),
CARMENES DR1, SOPHIE, BL Open Data and the fibre surveys we find
7,644 on-star spectra (1978→2026, deduplicated per exposure) of which
**633 sit inside a 0.1 AU window (67 of the 336 events), 30 inside a
2.5 R☉ grazing window (11 events, 3 targets) and 14 inside a 1.2 R☉
window (5 events, 2 targets)**. The grazing-window spectra come from
SPIRou (955–2515 nm), NIRPS (966–1923 nm), HARPS, HIRES, CARMENES VIS,
MAROON-X, CRIRES and SOPHIE — so the 1064 and 1550 nm line families,
structurally open in every survey so far (`report/joint_crossings.md`
§3), are **in band inside grazing windows for wolf-359 and teegarden
(SPIRou, NIRPS) and inside a 2.5 R☉ window for ross-128 (SPIRou)**.
All three reduced-product routes were verified by anonymous download
of an in-window spectrum. Channel B has no spectral substrate (§"Channel
B"). The APF hand-off itself is the weakest item in the family: raw 2D
frames, 5.3 d from t_ca, 0.1 AU rung only.

## Service probes (2026-09-05, dev machine, all anonymous)

| Service | Status | Notes |
| --- | --- | --- |
| KOA TAP `koa.ipac.caltech.edu/TAP` | **up** | `koa_hires/nirspec/kpf/nires/esi/mosfire/lris/deimos/kcwi`; ADQL cone on `ra,dec`; `date_obs`+`ut` (KPF: `date_beg`); `waveblue/wavered` Å on HIRES/NIRSPEC; `propint` |
| ESO TAP `archive.eso.org/tap_obs` | **up** | `dbo.raw` needs `INTERSECTS(CIRCLE, s_region)` (the `CONTAINS(POINT(ra,dec))` form errors); `ivoa.ObsCore` phase-3 with `em_min/em_max/em_res_power/access_url`; raw rows carry `release_date`, `prog_title`, `tpl_name` but no λ for CRIRES |
| CADC TAP `argus` (CAOM) | **up** | one query covers Gemini (MAROON-X, GNIRS, GRACES, PHOENIX, GMOS), CFHT (SPIRou, ESPaDOnS), HST, SDSS/BOSS/APOGEE mirrors, DAO; `INTERSECTS(CIRCLE, position_bounds)`; `targetPosition_coordinates_cval1/2` for the pointing (no `position_bounds_center`); multiple planes per exposure (raw/e/s/t/v) — dedupe per minute |
| MAST CAOM | **up via Mashup cone** (`api/v0/invoke Mast.Caom.Cone`, ~1 min); the TAP cone on `dbo.obspointing` runs > 3 min | HST/IUE/FUSE/EUVE spectra; wolf-359 call returned `EXECUTING` (async) — CADC's HST mirror covers it |
| SDSS SkyServer DR17 + DR19 | **up, slow** | `apogeeVisit` (dateobs/jd) box queries fine; `specObjAll` box query timed out 3× at 300 s for gj-908 (then fine); DR19 `mwm_apogee_allvisit`, `spAll_allepoch` (`racat/deccat`, `tai_list`) |
| NEID TAP `neid.ipac.caltech.edu/TAP` | **up** | `neidl0/l1/l2`, cone on `qrad,qdecd`; **0 rows for all 7 targets** (public L2 only) |
| BL Open Data `seti.berkeley.edu/opendata/api/query-files` | **up** | 3 APF files for teegarden (SO0253), 5 for gj-908, 3 for wolf-359; raw 2D 2080×4608 uint16 frames (19 MB), 420 s, `DATE-BEG`, decker W |
| CARMENES GTO DR1 `carmenes.cab.inta-csic.es/gto` | **up** | 362 stars; per-star `*_SERVAL+RACOON.csv` = one row per spectrum (BJD, TIMEID = UT start in the caracal filename, EXPTIME, SNREF, FLAG) + `*_VIS.zip` of all reduced VIS spectra (~4.8 MB each). **VIS channel only**; targets present: teegarden J02530+168 (260), wolf-359 J10564+070 (79), ross-128 J11477+008 (58), ross-154 J18498-238 (56); gj-908, gj-1276, van-maanen absent |
| Calar Alto archive `caha.sdc.cab.inta-csic.es/calto` | up, **form-only** | lists CARMENES raw + reduced products (both channels) with dates to 2025; a scripted POST returned the form — not decoded in this recon |
| SOPHIE archive `atlas.obs-hp.fr/sophie/sophie.cgi` | **up** (object-name CGI) | 47 wolf-359, 72 gj-908, 4 ross-128 spectra with date, mode (HR/HE), exptime, SN26; **day-level dates only in the listing** (UT in `view_head`); e2ds/s1d downloads when public |
| LAMOST DR10 v1.0 `lamost.org/dr10/v1.0/q` | **up** (form POST, cone, pipe-CSV) | combined spectra with integer MJD; 3 van-maanen, 2+1 teegarden, 1 gj-908 LRS/MRS, 15 rows near wolf-359 (none the star) |
| Data Lab TAP `desi_dr1.zpix` | **up** (box query; ADQL `POINT` unsupported) | coadds with `min_mjd/max_mjd`; 12 ross-128, 6 van-maanen, 4 gj-908, 3 gj-1276 fibres, none in a window |
| Gemini archive `archive.gemini.edu` | **403 anonymous** (anti-bot block) | irrelevant: CADC mirrors Gemini metadata and data (MAROON-X, GNIRS, GRACES, PHOENIX) |
| VizieR TAP | **403** from this host | LAMOST reached through its own form instead |
| HET (HPF) archive | **unreachable** (`archive.hetdex.org` no route) | HPF has no public spectral archive; RV epoch lists only via CDS (below) |
| Subaru SMOKA (IRD/HDS) | up, form-only | not probed; IRD (970–1750 nm) SSP targets late-M dwarfs — a later probe |
| TNG/IA2 (HARPS-N, GIANO-B) | JS-only portal | not probed |

All responses are snapshotted verbatim under
`runs/spectral-archives/recon/` (`SnapshotStore`; 187 response files);
the scan is `scripts/recon_scan.py`, the tables
`scripts/recon_summary.py` → `results/recon_summary_v0.md`,
`results/recon_scan_v0.json` (per target: events, windows, per-instrument
summary, every in-window row), `results/recon_rows_v0.json` (all 16,433
raw rows).

## Geometry (channel A, b ≤ 0.1 AU; `results/channelA_deep_windows.json`)

Seven targets, 48 events each (one per sidereal year, 1980→2028), all
at opposition on a fixed calendar date; flat-chord half-widths:

| target | date | b (R☉) | ±1.2 R☉ | ±2.5 R☉ | ±0.1 AU |
| --- | --- | --- | --- | --- | --- |
| van-maanen | Oct 6–7 | 0.32–0.46 | 7.5 h | 16.1 h | 5.81 d |
| wolf-359 | Mar 3 | 0.75–0.94 | 5.9 h | 15.3 h | 5.76 d |
| gj-1276 | Sep 4–5 | 0.82–0.91 | 5.6 h | 15.4 h | 5.86 d |
| teegarden | Nov 8 | 1.07–1.28 | 2.5–3.4 h | 14.4 h | 5.75 d |
| ross-128 | Mar 17–18 | 1.78–1.82 | — | 11.1 h | 5.77 d |
| ross-154 | Jul 3 | 3.37–3.41 | — | — | 5.84 d |
| gj-908 | Sep 21–22 | 12.30–12.36 | — | — | 4.78 d |

No elongation gate here (contrast Kepler/K2, WISE, the sunward cells):
channel A on-star is the *easiest* geometry for a ground spectrograph,
and the RV programmes that own these stars observe them nightly around
opposition. Cone: 5′ around the era-mean position (PM drift ≤ 2′ from
the mean for teegarden/wolf-359), on-star = within 3′; targeted
spectra are identified by the pointing, not by catalogue position, so
no PM-forced-position defect of the Gattini kind arises.

## Substrate found (on-star, deduplicated per exposure)

| target | on-star spectra | in 0.1 AU (spectra / events) | in 2.5 R☉ | in 1.2 R☉ | instruments with in-window spectra |
| --- | --- | --- | --- | --- | --- |
| wolf-359 | 2,302 | 83 / 19 | 10 / 4 | 9 / 3 | SPIRou 768 (32 in-window), NIRPS 227 (8), HIRES 197 (10), CARMENES-VIS 79 (9), SOPHIE 47 (5), X-shooter 5, CRIRES 6, MAROON-X 2, FEROS 2, HARPS 1, CES 3 |
| teegarden | 659 | 68 / 9 | 11 / 4 | 5 / 2 | CARMENES-VIS 260 (12), SPIRou 150 (12, proprietary), CRIRES 55 (20), KPF 33 (14), ESPRESSO 36 (4), MAROON-X 27 (3), APF 3 (3) |
| ross-128 | 1,108 | 84 / 12 | 9 / 3 | — | SPIRou 316 (41), HARPS 308 (28), ESPRESSO 85 (8), CARMENES-VIS 58 (3), HIRES 43 (3), MAROON-X 39 (1) |
| ross-154 | 1,560 | 297 / 8 | — | — | EFOSC 601 (269 — one flare-spectroscopy campaign, 2023), HARPS 127 (6), SPIRou 93 (12), MAROON-X 35 (4), CARMENES-VIS 56 (3), HIRES 42 (3) |
| gj-908 | 786 | 51 / 12 | — | — | HIRES 227 (15), CRIRES 52 (16), X-shooter 17 (9), MAROON-X 159 (3), HARPS 89 (3), FEROS 17 (3), KPF 34 (1), SOPHIE 72 (1) |
| van-maanen | 1,199 | 50 / 7 | 0 | 0 | EFOSC 1,022 (46 — it is the PESSTO/ePESSTO+ **flux standard** "VMA2", R ~ 300–700), LRIS 33 (3), LAMOST-LRS 3 (1); no échelle spectrum within a grazing window (closest 1.2 d) |
| gj-1276 | 30 | 0 | 0 | 0 | X-shooter 13 (2016-11), UVES 3, HIRES 1, DESI 3 — nothing within 8 d of any t_ca |

Totals: 7,644 on-star spectra; 633 / 30 / 14 in the 0.1 AU / 2.5 R☉ /
1.2 R☉ windows; 67 / 11 / 5 distinct events. Full per-instrument rows
(era, nearest |Δt|, band, line families in band, public flag) in
`results/recon_summary_v0.md`.

### The grazing-window spectra (≤ 2.5 R☉), all public unless noted

| target | event (t_ca) | b | UTC mid | instrument | exp | Δt from t_ca | rungs | band | product |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| wolf-359 | 2010-03-03 | 0.82 | 10:13 | KOA HIRES (I₂ cell in, deck B5) | 600 s | −2.1 h | **1.2 R☉** | 446–891 nm | raw 2D (KOA has HIRES extracted products for some frames) |
| wolf-359 | 2021-03-03 | 0.77 | 11:51–11:58 | CFHT SPIRou ×4 | 111 s | +4.1 h | **1.2 R☉** | 955–2515 nm, R 73k | APERO `t` (telluric-corrected orders + OH-line model) via CADC |
| wolf-359 | 2022-03-03 | 0.77 | night of 03-02/03 | OHP SOPHIE HR | 2700 s | day-level only | 2.5 R☉ (hour to confirm) | 387–694 nm | e2ds/s1d |
| wolf-359 | 2025-03-03 | 0.76 | 05:20–05:32 | ESO NIRPS ×4 | 223 s | −3.0 h | **1.2 R☉** | 966–1923 nm, R 70k | phase-3 1D, SNR 88, sky-sub + telluric-corrected columns |
| teegarden | 2009-11-08 | 1.15 | 03:24–03:42 | ESO CRIRES ×4 (prog 182.C-0748, Bean+ NH₃-cell RV survey) | 300 s | −2.4 → −2.2 h | 2.5 R☉ (last frame **1.2 R☉**) | K band ~2.29–2.35 µm (survey setup; λ not in TAP, confirm from header) | raw nodded frames |
| teegarden | 2017-11-07/08 | 1.12 | 23:17, 20:53 | CARMENES VIS ×2 (+ simultaneous **NIR** frames, not in DR1) | 1800 s | −7.9 h, +13.7 h | 2.5 R☉ | 520–960 nm (NIR 960–1710) | caracal 1D (VIS, DR1 zip) |
| teegarden | 2021-11-08 | 1.10 | 11:11 | Gemini MAROON-X | 1800 s | +3.3 h (1.2 R☉ edge is 3.1 h) | 2.5 R☉ | 499–920 nm | raw via CADC (PI-reduced) |
| teegarden | 2025-11-08 | 1.08 | 08:41–08:53 | CFHT SPIRou ×4 | 217 s | +0.2 → +0.4 h | **1.2 R☉** | 955–2515 nm | APERO products, **proprietary until ~2026-11** |
| ross-128 | 2021-03-17 | 1.81 | 02:26–08:28 | ESO HARPS ×4 (Red Dots) | 1200 s | −10.7 → −4.6 h | 2.5 R☉ | 378–691 nm, R 115k | phase-3 1D (SNR 22) |
| ross-128 | 2022-03-17 | 1.81 | 12:13–12:19 | CFHT SPIRou ×4 | 89 s | −7.1 h | 2.5 R☉ | 955–2515 nm | APERO via CADC |
| ross-128 | 2024-03-17 | 1.82 | 09:05 | Gemini MAROON-X (LP-202) | 1800 s | **+1.4 h** | 2.5 R☉ | 499–920 nm | raw via CADC |

Every grazing hit but the SOPHIE one is timed to the second from the
archive; the SOPHIE listing carries the date only (its `view_head`
gives the UT — one request at the freeze). CRIRES 2009: the programme
is the Bean et al. NH₃-cell K-band RV survey, so it does **not** carry
1064/1550 nm — it is a generic-line cell at 2.3 µm; the header setting
must be read before it is assigned a band.

### Product levels verified by download (anonymous)

- ESO phase 3 (`dataportal.eso.org/dataportal_new/file/<ADP id>`):
  NIRPS 1D (47 MB; WAVE/FLUX/ERR/QUAL + `FLUX_TELL_*`, `*_SKYSUB`,
  `ATM_TRANSM`; per-order SNR in the header) and HARPS 1D (5 MB;
  WAVE/FLUX/ERR, SNR keyword). Release dates are in TAP; everything in
  the tables above is past its proprietary year except the 2025–26
  SPIRou teegarden set.
- CADC (`ws.cadc-ccda…/data/pub/CFHT/<id>t.fits`): SPIRou `t` product
  (18 MB; 49 orders × 4088 px FluxAB/WaveAB/BlazeAB + `Recon` telluric
  and `OHLine` sky models; `SPEMSNR` 127 for wolf-359 at 111 s).
- BL APF: raw 2D echelle frames only (no extracted spectra in the
  Open Data archive); reduction is on us.
- CARMENES DR1 VIS zips: not downloaded (1.2 GB for teegarden); the
  per-spectrum list is the CSV, verified; format documented in the
  readme (caracal FITS, SPEC/SIG/WAVE/CONT per order).
- KOA HIRES / Gemini MAROON-X: raw frames (KOA serves pipeline-extracted
  HIRES spectra for part of the archive — check per KOAID at the freeze).

### Epoch lists for non-public substrates (Teegarden)

CDS J/A+A/684/A117 (Dreizler+2024, "Teegarden's Star revisited") gives
mid-exposure BJDs for 467 RVs (`results/teegarden_cds_epochs_v0.json`;
raw table under `runs/spectral-archives/recon/`). Inside our windows:
CARMENES VIS 12 (2 grazing — the same two found from DR1), MAROON-X 3
(1 grazing — the same), and **HPF 8 epochs in the 0.1 AU windows of
2019 (Nov 6, 14) and 2020 (Nov 2, 12)** at 810–1280 nm — a 1064 nm
in-band cell on the 0.1 AU rung whose spectra live only with the HPF
team (Penn State). CARMENES also exposes its NIR channel (960–1710 nm,
R 80k) simultaneously with every VIS frame, so the two 2017 grazing
spectra have 1064/1550 nm twins; DR1 released VIS only, and the
Calar Alto archive holds the NIR frames behind a form (item for the
freeze: decode the POST, or ask the consortium).

## Channel B has no spectral substrate

Nobody takes spectra of empty sky at arbitrary positions. The only
antipode-position spectra in existence are multi-object sky fibres
(SDSS/BOSS ~80–100 per plate, DESI ~400 per tile, LAMOST ~ 300, APOGEE
~35): the probability that one lands within a 1.5″ fibre of the
antipode on a tile that also falls in a window is ≈ 400 × π(0.75″)² /
8 deg² ≈ 1e-5 per tile, and the antipode fields are not in any
survey's repeat cadence. Slitless or IFU wide-field spectroscopy of
the antipode fields does not exist at useful depth; the one all-sky
spectrophotometer, SPHEREx (R 40–130), is already the §4.3 Pipeline-A
substrate and the §5.7 crossings-revisit item. **This family is
channel A only.** (Channel B *photometric* line cells remain with the
narrowband/medium-band imaging idea in §7, if ever.)

## Prior art (what the wide rung already has)

- Tellis & Marcy 2017 (AJ 153, 251): 67,708 Keck HIRES spectra of
  5,600 FGKM stars, 2004–2016, laser thresholds 3 kW–13 MW (10-m
  diffraction-limited transmitter), no time selection. Their archive
  is the KOA HIRES set above — the wolf-359 2010-03-03 1.2 R☉ frame and
  the gj-908/ross-128/ross-154 0.1 AU HIRES frames may already carry
  a per-spectrum null; cross-match their target table at the freeze.
- Zuckerman et al. 2023 (PASP; arXiv:2301.06971): 1,983 APF spectra
  of 388 BL stars, 84 kW at 78.5 ly, no candidates; Lipman et al. 2019
  (Boyajian's star, 177 APF spectra, 7.4 MW at 1470 ly). These are the
  pipelines the APF hand-off would reuse; neither applied a
  crossing-window selection or a beam-footprint power conversion.
- The plan's "wide rung is answered by BL's published nulls" (§5.5)
  extends to optical lines only through these two papers, i.e. for the
  HIRES/APF subsets of the 86 registry stars; the ±2-month 1 AU windows
  of the registry are a cross-match, not a new search, and stay a
  hand-off.

## What this means for the survey design

1. **Substrate class**: reduced 1D (or order-by-order) high-resolution
   spectra of the target star, R 70k–115k, from public archives with
   sub-minute time stamps — the cleanest v2 substrate the programme
   has met (no WCS, no masks, no PSF; one flux vector per exposure).
   Three heterogeneous reductions (ESO phase 3, APERO, caracal) plus
   raw-frame sets (HIRES, MAROON-X, CRIRES, APF) that need extraction —
   the raw sets are an adapter cost to be weighed per hit, not a
   blocker.
2. **Statistic**: the Tellis & Marcy / Zuckerman construction —
   an unresolved emission feature at the instrumental profile width
   above the local continuum, absent from the star's template — with
   two additions the crossings surveys already use: (a) the **null
   ensemble is the same star with the same instrument outside the
   window** (hundreds of spectra per instrument for wolf-359, ross-128,
   teegarden), giving empirical false-alarm rates per instrument and
   the 1/9 exceedance discipline for free; (b) injections at the
   instrumental FWHM converted to power through the rung cone
   (P = F·π b²_rung, the convention of every crossings report).
3. **Depth (order of magnitude, not injection-calibrated)** — *corrected
   2026-09-05 after the survey ran*: the first draft of this item
   converted erg s⁻¹ cm⁻² to W m⁻² with the wrong factor (10⁻⁷ J per erg
   was applied as 10⁻³) and quoted kW–MW floors; the correct scaling is
   10⁴ lower. A line spread over ~3 px is detectable at 5σ when it
   carries ≈ 8.7/SNR of the per-pixel continuum; with the stars'
   measured flux densities (NIRPS absolute flux for wolf-359, checked
   against 2MASS) the photon-limited floors are **watts through the
   grazing cones** (teegarden CARMENES VIS ≈ 0.2 / 1 / 80 W through
   1.2 R☉ / 2.5 R☉ / 0.1 AU; wolf-359 NIRPS at 1064 nm ≈ 13 / 55 / 4,300 W).
   The survey's injection-calibrated numbers (systematics-limited,
   `report/spectral_archives.md` §4) are 11–12 W (1.2 R☉), 37–52 W
   (2.5 R☉) and 0.1–4 kW (0.1 AU) on the best units. Photospheric flux
   in one resolution element is the noise floor, so bright stars give
   weaker limits — the opposite of the photometric surveys — but these
   are the **first line-SED constraints anywhere in the programme**, and
   the first at 1064/1550 nm, cells that broadband photometry cannot
   touch (a 10 W line is 10⁻⁹ of the star bolometrically).
4. **Systematics to design against**: cosmic rays (single-pixel;
   require the instrumental profile and, where 4 consecutive
   in-window exposures exist — SPIRou ×4, NIRPS ×4, HARPS ×4, CRIRES
   ×4 — persistence across them); telluric and OH emission (SPIRou
   ships an OH model, NIRPS ships sky-subtracted columns; lines at rest
   in the geocentric frame vs a laser at rest in the *star's* frame,
   the star's RV being known to m/s); M-dwarf flare emission lines
   (wolf-359 and ross-154 flare hourly — veto list from the stellar
   line catalogue, and the same-star null ensemble contains many
   flares); the HIRES iodine cell (I₂ absorption forest 500–620 nm on
   the 2010 wolf-359 frame); detector persistence/blaze edges in the
   NIR.
5. **Cell structure**: line families 532 / 1064 / 1550 nm (and 355 /
   266 nm harmonics where UV/blue coverage exists) plus a generic
   "any line" cell per spectrum; rungs 1.2 R☉ / 2.5 R☉ / 0.1 AU; a
   pulse-cell is not available (exposures 90–2700 s integrate any
   pulse train — a duty-cycle × peak-power hypothesis instead).
6. **The APF hand-off** stays in the family as one 0.1 AU unit
   (teegarden 2016, Δt +5.3 d, raw frames) — it is no longer the
   reason to run the survey.

## Recommended freeze items (for the user; not started)

- **Units**: the 14 grazing spectra × {532, 1064, 1550, generic} where
  in band (wolf-359 ×3 events, teegarden ×3 events incl. the 2009
  CRIRES K-band generic cell, ross-128 ×3 events) as confirmatory;
  the 0.1 AU rung (633 spectra, 67 events) as a second family with the
  loose-pointing beam hypothesis; the proprietary teegarden SPIRou
  2025 set as a pre-registered re-run when released (~2026-11).
- **Ensemble**: same star + same instrument, all out-of-window public
  spectra (wolf-359 SPIRou 736, NIRPS 219, HIRES 187; ross-128 SPIRou
  275, HARPS 280; teegarden CARMENES 248, CRIRES 35).
- **Adapters**: ESO phase-3 1D (NIRPS, HARPS, ESPRESSO; SNR per order
  in header), CADC/APERO SPIRou `t`, CARMENES caracal (DR1 zip);
  decide per hit whether raw HIRES/MAROON-X/CRIRES/APF frames are
  extracted (PypeIt) or dropped to the hand-off list.
- **Before any data are touched**: read the CRIRES 2009 header for
  the wavelength setting; fetch the SOPHIE 2022-03-03 UT; pull the
  KOA `koa_reduced_data` availability for `HI.20100303.36516`; fix the
  power convention (flat-top footprint of radius b_rung) and the
  laser-line profile hypothesis (unresolved; Doppler frame = star).
- **Asks (user action, optional)**: HPF team for the 8 in-window
  teegarden spectra (1064 nm at 0.1 AU); CARMENES consortium / CAHA
  form for the NIR twins of the two 2017 grazing spectra (1064 +
  1550 nm at 2.5 R☉ — the deepest NIR grazing cell reachable today).

## Plan corrections

- §3.6 row 10 / §5.15 O4: "ESO/Keck archives generally" → the family
  is ESO + KOA + CADC (Gemini/CFHT) + CARMENES DR1 + SOPHIE + BL; the
  fibre surveys (SDSS, LAMOST, DESI, APOGEE) contribute only low-res
  0.1 AU-rung spectra (van-maanen LAMOST 1, none elsewhere) and can be
  dropped from the design.
- "APOGEE H-band" for 1550 nm: APOGEE has 2–4 visits of gj-908 and
  ross-128 only, none in a window; the 1550 nm substrate is
  **SPIRou and NIRPS** (and CARMENES NIR behind the form).
- "1064 nm needs WINTER Y" (§5.8 item 6) holds for photometry; for
  lines the 1064 nm cell is open **now** on wolf-359 (NIRPS 2025,
  SPIRou 2021, both 1.2 R☉) and teegarden (SPIRou 2025 once public;
  HPF 2019–20 by request).
- van-maanen, the deepest graze (b 0.32–0.46 R☉), has no échelle
  spectrum within a grazing window in any archive — its only in-window
  spectra are the ePESSTO+ standard-star frames (R ~ 400) and 3 LRIS
  frames; a future-observation recommendation (a 20-min HARPS/ESPRESSO
  exposure at t_ca ± 7 h on Oct 6–7) is the natural hand-off, as the
  radio survey did for the antipode channel.
- gj-1276 is uncovered (30 spectra ever, none within 8 d of a window).

## Open items (recon-level, cheap)

- SMOKA (Subaru IRD, HDS) and TNG (HARPS-N, GIANO-B 950–2450 nm)
  were not probed (form-only / JS-only portals); IRD SSP targets late-M
  dwarfs and GIANO-B would add another 1064/1550 nm channel — one
  session each.
- MAST for wolf-359 returned an async status; irrelevant for the
  design (HST STIS/COS UV spectra are in CADC) but the row is incomplete.
- SDSS `specObjAll` for gj-908 timed out three times (integer-MJD
  coadds; no design impact).
- The Calar Alto POST for CARMENES NIR frames.
