---
title: "LASCO reachability recon (Pipeline B sunward channels)"
date: 2026-08-25
status: "recon complete — service probes + real-frame astrometry/star check; hypothesis freeze not started"
---

# LASCO recon (plan §5.8 item 9; follows the sunward geometry study)

Goal: DECam-style reachability recon for a LASCO crossings survey of
the sunward channels S1 (downlink post-lens, apparent source at the
solar-limb graze point) and S2 (uplink past the Sun, star near
conjunction), per `notes/sunward_geometry_study_2026-08-25.md`
(LASCO-era census: S1 grazing 128 events ≤ 1.2 R☉ / 160 ≤ 2.5 R☉;
0.1 AU rungs 224 events / 7 targets each).

## Service probes (2026-08-25, dev machine; all anonymous, no auth anywhere)

| Service | Status | Notes |
| --- | --- | --- |
| SDAC `umbra.nascom.nasa.gov/pub/lasco_level05/` | **works** | canonical level-0.5 tree, `YYMMDD/c2\|c3/*.fts`, plain HTTP dir listings; spans 1996-02 → **2025-02** (lags ~18 months) |
| NRL `lasco-www.nrl.navy.mil/lz/level_05/` | **works** (dated paths) | same layout; current through **2026-06** (2026-07 404; site says "up to 12 months behind realtime"). Bare top-level listing times out — always hit dated paths |
| NRL `lz/level_1/` | **works** | calibrated level-1 (`.fts.gz`, ~1 MB), same daily layout, spans 1996 → **2017-08-31** only; ~1:1 frame parity with level 0.5 (36 vs 37 on 2010-06-15 c3) |
| JPL Horizons `COMMAND='-21'` (SOHO) | **works** | SSB vectors returned at both era ends (1996-02, 2010) — `register_spacecraft_table_observer` route transfers from TESS unchanged |
| VSO SOAP (`sdac.virtualsolar.org` WSDL) | reachable | not needed — deterministic date paths cover discovery; keep as alternative |
| ESA SSA (`ssa.esac.esa.int/ssa/`) | reachable (landing) | not probed further; alternative only |
| Helioviewer API | works | JPEG2000 quick-look only (not photometric); confirmed **LASCO C3 still observing 2026-08-19** |

## Product facts (from real fetched frames, kept in `runs/heliospheric-crossings/recon/`)

- **Level 0.5** (`32227512.fts`, C3 2010-06-15; C2 `22335566.fts`):
  1024×1024 int16, ~2.1 MB plain FITS. C3 = Clear filter (~400–850 nm,
  contains the 532 nm doubled line), PLATESCL 56″/px; C2 synoptic =
  **Orange filter (~540–640 nm — 532 nm falls just outside)**,
  11.9″/px. Header WCS is **helioprojective only** (CTYPE SOLAR-X/Y
  arcsec, CRPIX = Sun center, CROTA 180 in the 2010 sample = SOHO roll
  state) — no celestial WCS, raw DN, bias in `OFFSET`. `MID_DATE`/
  `MID_TIME` = midpoint MJD + seconds-of-day. Some headers carry
  non-printable chars (astropy VerifyWarning, tolerated).
- **Level 1** (`35227512.fts.gz` = same exposure, filename 3rd digit
  2→5): calibrated to **BUNIT = MSB** (mean solar brightness;
  c3_calibrate + vignetting + ramp + distortion + mask per HISTORY),
  **derolled to solar north up** (RECTIFY=180) with residual
  CROTA −0.226°, and **time-corrected** DATE-OBS (−36 s vs level 0.5;
  original kept in HISTORY, residual ±15 s per TIME_DIFF.DAT) —
  negligible vs our 0.6 d+ windows.
- Frame rates (level 0.5, sampled): 1997: ~42 C2 + 26 C3 /day; 2002:
  157 + 87; 2010: 56 + 37 (low-activity year?); 2015: 121 + 111; 2024:
  124 + 112. Typical modern cadence ~12 min both cameras.
- **Failure mode seen live:** one C2 fetch silently truncated at 8 KB
  (server closed early; refetch of the same URL gave the full 2.1 MB).
  Every fetch must size-check + FITS-parse before acceptance (DECam
  lesson redux).

## Astrometry + star-visibility validation (the recon's go/no-go check)

`scripts/recon_star_check.py` → `results/recon_star_check_v1.json`, on
the level-1 C3 frame (2010-06-15 00:18 UT): built a synthetic celestial
WCS from geocentric Sun RA/Dec + solar P-angle (Meeus) + header CROTA,
extracted point sources (11-px median high-pass, local maxima), and
pattern-matched three bright Taurus stars:

- **Orientation resolved empirically** (the convention the adapter
  should hard-code): image x = −(East offset)/scale after rotating
  celestial offsets by +(P + CROTA); i.e. East to the left, solar
  north up. Wrong-sign alternatives fail outright (0–1 matches).
- **3/3 stars land at 0.0 / 0.4 / 0.7 px residuals** after one common
  translation: β Tau V 1.65 (det S/N 1124), ζ Tau V 3.01 (759),
  119 Tau V 4.32 (411) — in a *single 19 s exposure* with a crude
  high-pass. Crude single-frame depth estimate from the V 4.3 S/N:
  ~V 8–9 at S/N 8 (consistent with the known C3 stellar limit);
  window stacks (30–300 frames) should add ~2–3 mag effective.
- **The common translation is (−15.6, −2.6) px = 14.6′ ≈ 0.93 R☉
  projected** — the SOHO halo-orbit transverse offset measured
  directly. This is live confirmation of the geometry study's
  prerequisite: Earth-center impact parameters are invalid at R☉
  precision → the SOHO-observer crossing list (Horizons −21, probed
  working) is mandatory before any grazing-rung freeze; and per-frame
  star-fit WCS refinement (validated here at sub-px) absorbs the
  residual pointing error.

## What a survey needs (estimate)

1. **SOHO-observer crossing list** `crossings/soho_v1` (Horizons −21
   SSB vectors, 6 h sampling ≪ the ~178 d halo period; TESS pattern
   verbatim, raw response sha-pinned).
2. **Adapter** (new shape, but simple): discovery = dated directory
   listings (snapshot HTML), fetch = full frames (no cutout service;
   2 MB/frame is fine); synthetic celestial WCS (convention above) +
   per-frame star-fit refinement; per-frame star ZP (level-1 MSB and
   level-0.5 DN both go through star calibration — our standard
   anyway). Primary tree SDAC level 0.5 (+ NRL mirror for the
   post-2025-02 tail); level 1 (1996–2017-08) as the
   calibration-validation era.
3. **Rung accessibility check at pilot:** the S1 1.2 R☉ rung's
   apparent source likely sits **behind the C2 occulter** (usable C2
   inner edge ~2–2.2 R☉) — measure the real usable inner radius from
   the L1 mask files; 2.5 R☉ is in C2's sweet spot; both 0.1 AU rungs
   (ε ≤ 5.7°) live in C3. Wavelength cell: 532 nm doubled line only
   in C3 Clear; C2 Orange is broadband-leakage only.
4. **Data budget:** S1 grazing (160 windows × 0.6–1.3 d × ~100–230
   frames/day) ≈ 40k frames ≈ 80 GB transient — fine with per-batch
   purge. The 0.1 AU rungs (224 windows × 11.5 d) ≈ 280k frames ≈
   0.5 TB at full cadence → the freeze must set a cadence-subsampling
   or frame-selection policy (e.g. hourly picks ÷12).
5. **Era gaps to encode:** SOHO attitude loss ~1998-06 → 1998-10 (+
   gyro-loss recovery to 1999-02) and periodic 180° roll states
   (level 0.5 CROTA flips; post-2003 quarterly flips) — from general
   knowledge, verify against the trees during coverage.

## Open items before a hypothesis freeze

- Build `crossings/soho_v1`; re-run the sunward census on it (grazing
  b values will shift by up to ~0.9 R☉ — the family membership of
  individual events may change).
- Decide the sunward-channel statistics construction (recurrence
  stacks per the ATLAS freeze pattern; ring/temporal controls need a
  design pass — coronal background is structured and rotates).
- S1 post-lens annulus beam profile (does b_e < b_graze see flux?) —
  declared hypothesis, affects window/units definition.
- Verify C2/C3 usable occulter radii + vignetting-usable annulus from
  L1 masks; confirm the 1.2 R☉ rung's fate.
- Confirm the level-0.5-only era (2017-09+) star-ZP chain against the
  level-1 overlap era before trusting it.
- Volume policy for the 0.1 AU rungs (subsampling rule frozen, not
  improvised).

## Addendum 2026-08-25 (later) — `crossings/soho_v1` built + census re-run

The first pre-freeze open item is done. `sglsurvey/crossings.py` gained
`--observer soho` / `--fetch-observer soho` (Horizons −21, the TESS
pattern verbatim): table `crossings/observers/soho_sc_ephemeris.npz`
(44,925 rows, 6 h, 1996-01-01 → 2026-10-01 — **stop bounded by the
Horizons SPK end 2026-10-05**, not the universal 2028 window; extend at
the yearly refresh; sanity check vs Earth ephemeris: SOHO–Earth
0.0083–0.0113 AU ✓). Product `crossings/soho_v1` (`xng-298d55d0ce6b`):
10,714 events, 0 invalid, manifest parity.

Census + validation (`scripts/soho_census.py` →
`results/soho_census_v1.json`):

- **Grazing-family correction is TESS-scale, as predicted:** matched
  sunward grazing pairs (either frame ≤ 2.5 R☉) shift by |Δb| median
  0.124 / max 0.211 R☉ (both combos). Large-b events shift up to
  2.7 R☉ — there the 0.01 AU L1 radial offset projects transversely —
  so **any wide-rung unit definition must also use this list**, not
  just the grazing rungs.
- **Rung membership survives:** S1 ≤ 1.2 R☉ = 117 events / same 4
  targets (Earth-center same-era: 122), ≤ 2.5 R☉ = 153 / same 5
  targets (Earth: 153). S1 0.1 AU: 216 events / 7 targets.
- **The van-maanen annual October S1 family deepens** in the SOHO
  frame: the 2009–2015 events drop from b ≈ 0.27 R☉ (Earth) to
  0.13–0.23 R☉ — sub-quarter-photosphere grazes, annually recurring,
  during years with full LASCO level-1 coverage. This strengthens the
  headline cell (the family PS1 masks and TESS sector gaps lost).
- Era sub-censuses (level-1 era to 2017-08, WISPR era 2018+) are in
  the JSON for the freeze to draw on.

Remaining pre-freeze items are unchanged (C2 occulter fate of the
1.2 R☉ rung, control design, volume policy, level-0.5 star-ZP chain
validation, S1 annulus beam profile).
