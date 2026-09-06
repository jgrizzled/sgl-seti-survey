---
title: "PSP/WISPR reachability recon (Pipeline B sunward channels, plan §5.15 O14)"
date: 2026-09-05
status: "recon complete; hypotheses + threshold freeze drafted, awaiting approval"
---

# PSP/WISPR recon (plan §5.15 item O14; follows the O5 geometry pass, §5.18)

Goal: reachability recon for a WISPR crossings survey of the sunward
channels S1 (downlink post-lens) and S2 (uplink past the Sun) on the
**0.1 AU rungs**, the third heliospheric substrate after LASCO (§5.11)
and STEREO-A HI-1 (§5.16). Prerequisite done: the PSP observer list
`crossings/psp_v1` and the per-encounter visibility census
(`notes/psp_wispr_geometry_pass_2026-09-05.md`, `results/psp_census_v1.json`).

## 1. Service probes (2026-09-05, server; all anonymous)

| Service | Status | Notes |
| --- | --- | --- |
| **NRL WISPR tree** `wispr.nrl.navy.mil/data/rel/fits/{L1,L2,L2b,L3,LW}` | **works** | plain-HTTP Apache indexes; L2 by day (`L2/YYYYMMDD/`, **1,019 day dirs 2018-11-01 → 2026-03-17**), L3 by orbit and day (`L3/orbitNN/YYYYMMDD/`, orbits 01–27); files `psp_LX_wispr_YYYYMMDDTHHMMSS_VN_WXYZ.fits`, 3,954,240 bytes; full download ~1–2 s, HTTP Range honoured (header-only reads) |
| Users Guide v5 + Hess+2021 calibration paper | fetched | `runs/wispr-crossings/recon/` |
| JPL Horizons `-96` | works | 2018-08-13 → 2026-12-01 (the O5 table) |

Snapshot (`scripts/nrl_inventory.py` → `runs/wispr-crossings/coverage/
l2_inventory.json`, `l3_inventory.json`): **L2 106,752 files / 1,019
days; L3 84,732 files / 554 orbit-days.** Codes W X Y Z = telescope
(1 = I, 2 = O), readout microcode (XY), region (Z; 1–2 synoptic):
synoptic WISPR-I **55,918** (`1221` 31,974 — the E1–E5 and later
half-cadence series; `1211` 22,863; `1141` 851; `1301` 230), synoptic
WISPR-O 50,793 (`2222`, `2302`, `2212`). Versions V1–V3. The `2302`/
`2222` pairs share timestamps in E10–E16 (two products of one exposure
set — a dev item; WISPR-O is secondary here).

## 2. Product facts (real frames, `runs/wispr-crossings/recon/`)

Frames fetched: WISPR-I L2+L3 2024-12-21 23:30 (r 0.157 AU), 2024-12-24
00:00 (0.0585 AU, 12 h before the E22 perihelion), 2024-12-27 00:30
(0.157 AU); WISPR-O L2+L3 2024-12-24 00:02; plus headers of `1221`
(E1), `2302` (E3, E10), `2222` (E10), `1211` (E14), `1141` (E24).

- **Full-field synoptic frames**: 960 × 1024, 2 × 2 binned,
  **WISPR-I 0.04231°/px = 152″**, WISPR-O 0.05648°/px = 203″;
  `NSUMEXP` 5 (8 in E1) on-board-summed exposures, `XPOSURE` total
  9–47 s in WISPR-I (**15 s at r 0.058 AU, 47 s at 0.157 AU**; 9 s at
  the E24 perihelion), 24–41 s WISPR-O (300 s for `2302` at E3, the
  long-exposure product); `GAINMODE` LOW at perihelion, HIGH outside.
  Units **MSB** (mean solar brightness, exposure-normalised); 0.4 %
  NaN (edge rows/columns `DSTART/DSTOP`).
- **Two WCSs**: primary `HPLN/HPLT-ZPN` (helioprojective; centre
  CRVAL (31.9°, −7.4°) for I, (76.7°, −13.2°) for O; a rotation PC of
  ~11.5°) and the alternate **`RA---ZPN/DEC--ZPN` ('A')** celestial
  system (`get_wispr_pointing.pro` — star-fit pointing per the
  guide's lineage). Spacecraft HCI position/velocity and DSUN_OBS in
  the header.
- **Ephemeris caveat**: header `DSUN_OBS` agrees with the Horizons
  table to < 0.01 % on the E1/E3/E10/E14 frames but is **+0.65 % on the
  E22 frame** (a 12.5-min-equivalent offset; predicted-SPK processing)
  and −0.10 % at E24. Every geometric quantity of the chain (r, b_e,
  ε, encounter phase) comes from the Horizons table; only the
  star-verified celestial WCS is taken from headers.
- **L3 = L2 − background model, then × (r/0.2 AU)^2.3** (HISTORY:
  "Grand Minimum Model (0.10 perc.) created from Orbits 020-021-022";
  "Image brightness corrected by (S/C distance)^2.3, normalized at
  0.200 AU"). Star fluxes in L3 are the L2 fluxes times that scale
  (ratio 0.0594 measured vs 0.0592 computed at 0.0585 AU); undoing the
  scale gives L2 units with the F-corona removed: **L2/L3u aperture
  flux 0.95–0.98** on 100–415 calibrators. The 10th-percentile
  three-orbit model is built in the spacecraft frame, so it removes
  the F-corona and stray light and leaves the (moving) stars — the
  search substrate is **rescaled L3** with L2 as the parity check.

## 3. Field model verified from the header (the O5 recon item)

`scripts/recon_star_check.py` → `results/recon_star_check_v1.json`,
`field_model_check`: the frame's own celestial WCS mapped into the
`psp_geometry` ram frame (Sun direction, prograde tangential, orbit
normal from the Horizons velocity):

| | centre (ε, ram-lon, orbit-lat) | left edge | right edge | top / bottom orbit-lat |
| --- | --- | --- | --- | --- |
| WISPR-I | 32.8°, **+32.5°**, −5.1° | ram-lon **13.7°** (ε 13.8°) | 51.6° (ε 52.0°) | +14.9° / −25.2° (corners +16 / −26) |
| WISPR-O | 77.1°, **+76.9°**, −9.5° | 48.9° | 105.0° | +20.1° / −39.2° |

Both fields are on the **ram side** (positive ram longitude), as the
literature said; the inner edge of WISPR-I is at 13.5°–13.8° from Sun
centre; the fields are **offset south of the orbital plane by ~5°
(I) / ~10° (O)** and the vertical extents are asymmetric (−25°/+15°
rather than ±20°). Header speed 167.4 km/s vs table 168.0. The census
model (±20°, ±26.5° symmetric) is therefore slightly optimistic on
the north side and pessimistic on the south; the coverage intersect
below uses the census model — the survey's footprint test uses **each
frame's own WCS** (per-frame pixel position of the fixed direction),
so the model only affected which rows were counted.

## 4. Astrometry + photometry validation (the go/no-go check)

`scripts/recon_star_check.py` (first pass, 2-px aperture — superseded
for photometry) and `scripts/recon_depth.py` →
`results/recon_depth_v1.json` (3-px aperture, annulus 5–9 px, 15-px
median high-pass on rescaled L3, Hipparcos V ≤ 9.5, calibrators V
3–7.5 at S/N > 10, 3σ-clipped):

- **Header celestial WCS is ~1 px**: WISPR-I 191 Hipparcos matches
  within 3 px, rms (0.93, 1.03) px, mean offset (−0.14, +0.65) px
  (2.5′); WISPR-O 727 matches, rms 0.7 px, offset (0.38, 0.10). Good
  enough for a 3-px aperture; a per-frame translation refinement is
  cheap and frozen in (HI-1 pattern).
- **PSF**: radial profile 1.0 / 0.55–0.62 / 0.16–0.25 / 0.03–0.06 at
  0 / 0.5 / 1 / 1.5 px → FWHM ≈ 1.6 px = 4′; encircled fraction
  0.50–0.67 (r 1 px), 0.86–0.96 (2 px), **0.95–1.00 (3 px)**. Aperture
  r = 3 px frozen.
- **Photometry**: ZP (V = ZP − 2.5 log₁₀ ΣMSB) = **−21.36, −21.35,
  −21.32** on the three WISPR-I frames (r 0.157 / 0.058 / 0.157 AU:
  the exposure normalisation holds to 0.04 mag across a 3× change in
  exposure), colour coefficient **+0.23 to +0.31 mag/(B−V)**, scatter
  0.14–0.18 mag MAD (258 / 100 / 415 calibrators), per-V-bin residuals
  ≤ 0.1 mag from V 3 to 8, residual vs elongation −0.09 → +0.10 from
  ε 10° to 55° (a mild vignetting/flat residual at the outer edge; the
  2-D ZP handles it).
- **Single-frame 5σ depth (solar colour, 3-px aperture, star-masked
  high-pass of rescaled L3)** — the decisive number:

| r (AU) | XPOSURE | ε 15° | 22° | 30° | 38° | 46° | 50° |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.157 (E22 −2.2 d) | 47 s | V 7.2 | 7.3 | 7.6 | 7.8 | 7.75 | 7.6 |
| 0.0585 (E22 −0.5 d) | 15 s | **5.4** | 5.8 | 6.2 | 6.3 | 6.3 | 6.1 |
| 0.157 (E22 +2.5 d) | 47 s | 7.3 | 7.4 | 7.7 | 7.9 | 7.9 | 7.7 |

  L2 (F-corona in) is 0.2–0.4 mag shallower than rescaled L3. The
  noise is structure-limited (K-corona streamers, residual F-corona
  gradient), not photon-limited. So **per frame V 5.5–8 depending on
  r**, 2–3 mag shallower than HI-1's V 10–10.5 and ~1 mag shallower
  than LASCO C3's V ~8 — offset by 5–15 min cadence (90–380 frames
  per arc) and 20–27 encounters per unit.

## 5. Cadence and coverage (inventory × census, `scripts/coverage_intersect.py` → `results/coverage_v1.json`)

Per-encounter synoptic WISPR-I cadence: E1 34 min, E2 9, E3 30, E4
24, E5 18, E6–E13 15–16 (E7 30), **E14–E16 7.5, E17–E23 5.0, E24 2.6,
E25 5.0, E26 7.5, E27 5.0**; arcs span −5 → +5 d about perihelion
(−2.5 → +2.6 d for the high-cadence `1211` series, the `1221`
half-cadence series filling the outer days). Gaps > 3 h: 0–2 per
encounter. **E28–E29 not yet released** (E27 = 2026-03 is the last).
19,702 WISPR-I frames lie outside the r < 0.25 AU arcs (extended
campaigns).

Per-event intersect (0.1 AU rung rows with a WISPR-I arc, frames
within the geometric in-beam-in-field arc): **41 units with ≥ 10
covered events (988 covered events, 171,927 arc-frame measurements
over ~36k distinct in-encounter WISPR-I frames)**; per unit 10–27
covered encounters, median 86–384 arc frames per event, minimum 13.
The 7-system family: van-maanen S1 (22 events, median 302 frames) /
S2 (26, 208), wolf-359 S1 (26, 128), teegarden S1 (25, 167) / S2 (24,
179), gj-1276 S2 (26, 126), ross-128 S1 (26, 162) / S2 (20, 260),
ross-154 S2 (27, 89), gj-908 S1 (20, 377) / S2 (26, 164). In E22 the
perihelion-side arcs run t_ca + 2.6 h → + 13–17 h at b_e 2.6 → 14 R☉,
r 0.049–0.08 AU, ε 14° → 53.5°; the straddling arcs (gj-908 S1,
ross-128 S2, van-maanen S1) run 34–160 h after t_ca at b_e 8–17 R☉,
r 0.07–0.23 AU.

**No off-beam star-fixed baseline at perihelion crossings.** The
source enters WISPR-I at ε 13.5° already in-beam and leaves at 53.5°
with b_e = 1.35 r < 0.1 AU whenever r < 0.074 AU — the entire WISPR-I
transit is in-beam; the source only leaves the beam in WISPR-O at
ε ≈ 63°. Off-beam WISPR-I baselines exist only for the larger-r
crossings (10 units, e.g. teegarden S2 176 frames, gj-581 S1 210,
wolf-1061 S1 204). The HI-1 star-fixed differential is therefore not
available as the primary static-content remover; the LASCO **stellar
template** construction is (hypotheses §6), with the star-fixed
differential as a secondary check where a baseline exists.

## 6. Pre-freeze data contact

The four recon frames lie inside in-beam arcs of E22 (the busiest
encounter): WISPR-I 2024-12-21 23:30 in the arcs of gj-1002 S1,
gj-54 S1, gj-908 S1, ross-128 S2, wolf-437 S2; 2024-12-24 00:00 in
teegarden S1, gj-581 S2, wolf-1061 S2; 2024-12-27 00:30 in gj-581 S1,
wolf-1061 S1 (and the teegarden S2 baseline); WISPR-O 2024-12-24 00:02
in 12 arcs (list in `results/coverage_v1.json` `recon_contact`).
Nothing was measured at any of those positions (full-field Hipparcos
match and column-strip noise). Declared remedy (the HI-1 pattern):
the four frame epochs are **excluded from the series of every event
they touch** (1 of 86–553 frames each; recorded in the config) and
the recon frames are never search frames.

## 7. What a survey needs — drafted in `hypotheses.md`

1. Substrate: NRL L3 synoptic WISPR-I (rescaled to L2 units by the
   header DSUN factor of the L3 processing — the factor is a pure
   function of the header value and is inverted exactly), L2 as
   parity. Per-frame chain: celestial WCS refined by a Hipparcos
   match → 15-px median high-pass → top-hat photometry (r 3 px) at
   ICRS-fixed patches → 2-D colour-corrected ZP.
2. Constructions: stellar-template static-content removal (Tycho-2 +
   Hipparcos through the measured colour system), 8 parallel-track
   control patches at orbit-latitude offsets ±1.5°…±6°, same-frame
   differential, S_stack / S_event / S_pulse, T = max over controls;
   star-fixed differential reported where a baseline exists.
3. Masks: synoptic-region gate, XPOSURE/GAIN recorded, NaN validity,
   planet/comet proximity (Mercury, Venus, Earth, Jupiter, Saturn are
   all regular WISPR transients), Tycho VT ≤ 8 patch mask with the S2
   exemption, per-encounter product code recorded.
4. Volume: ~36k WISPR-I frames × 4 MB ≈ 145 GB transient (the whole
   in-encounter synoptic archive); one measurement pass per frame
   serves every patch mapped to it.
5. Positive controls: the 61-vir S2 patch (V 4.7 in the aperture —
   the template null at high S/N), gj-908 S2 (V 9.0), and a moving
   body through the arc region chosen at dev from Horizons `@-96`
   (Vesta/Ceres/Pallas or a bright comet; Mercury/Venus saturate).
