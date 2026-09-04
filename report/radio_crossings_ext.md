# Radio crossings coverage — extension v2: ASKAP (RACS/VAST/EMU) + LoTSS

Plan §5.8 queue item 7. Executed 2026-09-04 as a metadata-only
coverage intersection, repeating the VLASS construction of
`report/radio_crossings.md` on two archives that revisit fields and
publish per-observation dates. **The radio decision (§5.5, 2026-08-24)
is unchanged: geometry-only, no radio signal search in this
repository.** No pixel was touched and no statistic was formed. The
product is a coverage ledger and a hand-off list.

## 1. What was intersected

- **CASDA ObsCore** (`casda.csiro.au/casda_vo_tools/tap`, anonymous,
  live 2026-09-04): every Stokes-I continuum restored image (`cube`,
  `cont.restored.t0`, `image.i.*`) of *any* ASKAP collection whose
  `s_region` contains the target-channel position — RACS-low/mid/high
  and the low re-epochs (2019-04 →), the VAST pilot (2019-08 → 2021-11)
  and full survey (2022-11 →), EMU, FLASH, WALLABY, DINGO, commissioning
  and guest observations. Each row carries its scheduling block (SBID)
  with start/stop MJD; the VAST-pilot cubes have no `t_min` in ObsCore
  and were dated from `casda.observation` by SBID. Several cubes per
  SBID exist (pipeline v1/v2 re-runs, lowres/raw variants) — one
  observation per (SBID, field centre), best `quality_level` kept
  (GOOD > UNCERTAIN > NOT_VALIDATED > REJECTED > BAD; REJECTED = failed
  validation and re-observed). 91 cones (84 A + 7 B), 0 failures;
  archive span 2019-04-21 → 2026-08-28.
- **LoTSS DR3** (ASTRON VO TAP `lotss_dr3.pointings`): 2,551 HBA
  pointings, 5,458 observations 2014-05-23 → 2024-08-18, the mid-MJD of
  every 8-h run that fed each pointing (`dateallobs`). DR3 covers 88 %
  of the northern sky and supersedes DR2; only positions δ > −10° were
  tested.
- Windows from `crossings/universal_v1`, frozen ladder (A 0.1 / 1.0 AU;
  B 1.2 / 2.5 R☉ / 0.1 AU), events 2014–2028; an observation covers a
  window when its actual on-sky interval overlaps the flat-chord window
  (no scheduling tolerance, unlike VLASS).
- Declared geometry: ASKAP `in_image` = the archive's own polygon
  containment (RACS images extend ~4–5° from the field centre);
  `in_footprint` = separation from the field centre ≤ 2.25° + HPBW/2
  (closepack36 at 0.9° pitch; 3.1° at 888 MHz, 2.8° at 1.37 GHz, 2.7°
  at 1.66 GHz). LoTSS `in_image` ≤ 2.61° (the 30 %-of-beam mosaic
  trim), `in_beam` ≤ 1.98° (half power of the 3.96° HBA beam).
- Quality caveat: ASKAP `quality_level` is the survey team's
  validation flag on the image, not a per-position sensitivity;
  UNCERTAIN images are released, REJECTED ones are not.

## 2. Results

Window census — events whose whole window lies inside the archive's
date span; "covered" = ≥ 1 observation with the position inside the
footprint/beam; "validated" additionally excludes REJECTED/BAD images.

| arm | rung | events (targets) in span | covered | validated |
|---|---|---|---|---|
| ASKAP | **B 1.2 R☉** (grazing) | 29 (4) | **0** | 0 |
| ASKAP | **B 2.5 R☉** | 36 (5) | **0** | 0 |
| ASKAP | **B 0.1 AU** (antipode) | 50 (7) | **7** | **6** |
| ASKAP | A 0.1 AU | 50 (7) | 3 | 3 |
| ASKAP | A 1.0 AU | 606 (86) | 211 | 206 |
| LoTSS | B 1.2 / 2.5 R☉ | 10 (1) / 20 (2) | 0 / 0 | — |
| LoTSS | B 0.1 AU | 30 (3) | 0 | — |
| LoTSS | A 0.1 AU | 50 (5) | 1 | — |
| LoTSS | A 1.0 AU | 300 (30) | 13 | — |

The 2.25 R☉ grazing family is time-starved as expected: ±0.3–0.7 d
windows against 12–15-min visits at cadences of weeks (VAST) to years
(RACS) — 65 in-span ASKAP events, none touched. Every hit below is in
a ±5–6 d 0.1 AU window.

### 2.1 Antipode channel (B 0.1 AU) — first archival radio coverage

The 2026-08-24 report found the antipode channel "archivally virgin":
zero Breakthrough Listen pointings within 0.5° in any window, one
VLASS tolerance-edge tile. Wide-field ASKAP imaging changes that at
the 0.1 AU rung. In-footprint observations, one line per (event,
SBID):

| target | b | t_ca | collection / SBID / field | ν | obs start (UT) | offset of ±half-window | sep / footprint | quality |
|---|---|---|---|---|---|---|---|---|
| **wolf-359** | 0.71 R☉ | 2023-09-05 | VAST 52549 VAST_2257-06 | 888 MHz | 2023-09-03 16:24 | −2.1 / 5.9 d | 0.80° / 3.13° | UNCERTAIN |
| **wolf-359** | 0.71 R☉ | 2024-09-05 | VAST 65727 VAST_2257-06 | 888 MHz | 2024-09-10 15:24 | +5.6 / 5.9 d | 0.80° / 3.13° | GOOD |
| **ross-154** | 3.30 R☉ | 2022-01-02 | RACS-high 34957 RACS_0651+23 | 1656 MHz | 2021-12-29 17:23 | −3.4 / 5.6 d | 0.63° / 2.72° | GOOD |
| ross-128 | 1.85 R☉ | 2021-09-19 | VAST pilot 32330 VAST_2338+00 | 1368 MHz | 2021-09-21 15:13 | +1.7 / 5.8 d | 2.45° / 2.82° | GOOD |
| van-maanen | 0.25 R☉ | 2022-04-03 | RACS-low 38682 RACS_1237-06 | 888 MHz | 2022-03-29 16:54 | −4.8 / 5.8 d | 3.12° / 3.13° (edge) | GOOD |
| teegarden | 0.98 R☉ | 2026-05-06 | FLASH 84179 FLASH_389 | 856 MHz | 2026-05-03 11:44 | −3.2 / 5.9 d | 2.18° / 3.16° | UNCERTAIN |
| van-maanen | 0.24 R☉ | 2026-04-03 | FLASH 83234 FLASH_497 | 856 MHz | 2026-03-29 11:58 | −5.0 / 5.8 d | 3.15° / 3.16° (edge) | REJECTED |

Plus 3 in-image-only rows (ross-128 2021 second pilot field at 3.16°;
van-maanen 2022 second RACS field at 3.32°; gj-908 2022 RACS-low at
3.92°). The two wolf-359 rows are the clean cases: the antipode sits
0.8° from the VAST_2257-06 field centre, inside the primary-beam core,
in two consecutive years — the September wolf-359 B family that the
optical surveys (ZTF, PGIR visit list) also cover. Antipode visit
lists over all dates: gj-908 89, van-maanen 79, ross-128 57, gj-1276
45, wolf-359 44, teegarden 21, ross-154 13 dated ASKAP observations
each; LoTSS reaches only gj-1276 B (18), ross-128 B (4) and
ross-154 B (1), none in-window.

### 2.2 On-star narrow rung (A 0.1 AU)

| archive | target | b | t_ca | field | obs | offset | sep |
|---|---|---|---|---|---|---|---|
| VAST | gj-1276 | 0.83 R☉ | 2023-09-05 | VAST_2257-06 SB52549 (888 MHz) | 2023-09-03 | −1.6 / 5.9 d | 1.10° (UNCERTAIN) |
| RACS-mid | teegarden | 1.08 R☉ | 2024-11-08 | RACS_0248+18 SB67840 (1368 MHz) | 2024-11-11 | +3.6 / 5.8 d | 2.09° |
| VAST | van-maanen | 0.33 R☉ | 2024-10-06 | VAST_0037+06 SB66449 (888 MHz) | 2024-10-07 | +1.2 / 5.8 d | 3.11° (edge) |
| **LoTSS** | teegarden | 1.09 R☉ | 2023-11-08 | P043+19 (144 MHz) | 2023-11-04 22:09 | −3.9 / 5.8 d | 1.89° (in beam) |

Out-of-footprint (in-image only): gj-908 2021-09 in two VAST-pilot
fields at 3.6°; teegarden 2024-11 in two further RACS-mid fields at
3.7–4.2°. The 2026-08-24 report's BL narrow-rung on-star cell (zero
GBT/Parkes/MeerKAT pointings) therefore now has four wide-field
continuum epochs beside it — none of them a targeted or
high-time-resolution observation.

### 2.3 Wide rung (A 1.0 AU)

Abundant, as with BL: ASKAP 211 of 606 in-span events covered in
footprint (64 targets; over all in-era events VAST covers 124,
RACS 96), LoTSS 13 of 300 (10 targets). The rung's windows tile ⅓–⅔ of the calendar; the
constraint is the surveys' own published transient/variable
searches (VAST, RACS, LoTSS), cited as for BL — nothing to add here.

## 3. What changes and what does not

1. **Decision unchanged.** These are 12–15-min (ASKAP) and 8-h (LoTSS)
   wide-field continuum images at ~0.1–0.3 mJy/beam; continuum
   imaging dilutes a narrowband line by the ~10⁹ channel-to-band ratio
   (§5.5), so the cell a targeted narrow-rung analysis would open is
   still empty, and a broadband transmitter is what a cutout would
   test. That is a follow-up class the repo already hands off (VLASS
   six-cutout quick-look), not a search this repo performs.
2. **The antipode channel is no longer archivally virgin at the
   0.1 AU rung.** Six validated wide-field epochs (five targets) sit
   inside antipode windows, two of them well inside the primary beam
   (wolf-359, 2023 and 2024). The 2026-08-24 statement stands for
   *targeted* radio (BL) and for the grazing rungs (still zero
   everywhere); it should be quoted with that qualifier from now on.
3. **Grazing rungs stay open in radio.** 65 in-span ASKAP grazing
   events, 30 LoTSS, zero observations. The only route is scheduled
   observation during a predicted window (the geometry products give
   the times to the hour).
4. **VAST is the archive that will keep answering.** Its cadence
   (weeks) over the wolf-359/gj-1276 antipode field makes the 0.1 AU
   antipode windows a recurring test; the repeat intersection is a
   one-command, ~12-min TAP run.

## 4. Hand-offs (not executed here)

- (a) **Antipode cutout quick-look**: the 6 validated B 0.1 AU epochs
  (wolf-359 ×2, ross-154, ross-128, van-maanen 2022, teegarden 2026)
  — broadband transient check at the relay position, RACS/VAST
  catalogue cone-search first (CASDA `AS110` / `AS207` tables), image
  cutout second; same class as the VLASS quick-look, to be done
  together.
- (b) **On-star narrow-rung epochs**: gj-1276 VAST 2023, teegarden
  RACS-mid 2024 + LoTSS 2023, van-maanen VAST 2024 — nearby M dwarfs
  and a white dwarf; a catalogue detection at the star would itself be
  a flare, worth one look alongside (a).
- (c) **VLASS refresh** at the yearly crossing-list refresh (§5.7):
  epochs 3.2 / 4.x have accrued since 2026-08-24; the v1 arm re-runs
  unchanged.
- (d) BL-MeerKAT re-check stays a §5.7 item (unchanged).

## 5. Products

`surveys/radio-crossings/scripts/coverage_intersect_v2.py` (the
intersection; `census_v2.py` the tables above),
`results/{askap_inwindow_v2.ecsv, askap_targets_v2.json,
lotss_inwindow_v2.ecsv, lotss_targets_v2.json,
coverage_v2_summary.json}`; raw TAP responses (all 91 CASDA cones, the
SBID date look-ups, the LoTSS pointing table) sha-pinned under
`runs/radio-crossings/v2/` with `run.log` and `census.txt`; the CASDA
collection census and ASTRON table list under
`runs/radio-crossings/recon/`.
