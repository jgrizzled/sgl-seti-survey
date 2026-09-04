---
title: "Kepler prime + K2 campaign footprint intersect (plan §5.8 item 8)"
date: 2026-09-04
status: "complete — 0 channel positions on Kepler silicon during any crossing at any rung; structural (elongation gate), not a coverage accident; item 8 CLOSED, no freeze"
---

# Kepler / K2 footprint intersect (plan §5.8 item 8)

Goal: the plan row asked for a one-afternoon intersect of the K2
ecliptic campaign fields (2014–2018) and the Kepler prime field
(2009–2013) against the universal crossing list, proceeding to a
crossings survey only on a hit (30-min continuous cadence in the
pre-ZTF era would open the pulse cell on 2009–2015 windows).

Headline: **no hit, and none was possible.** The K2 pointing
constraint held every campaign boresight 61°–158° from the Sun as
seen from the spacecraft, i.e. never closer than 22° to the anti-sun
point (37° with the true campaign-14 start), while every channel
position at a searched rung sits within 3.5° of the anti-sun point
(b ≤ 0.1 AU). The Kepler prime field is at ecliptic latitude +65°
and stays 66°+ from anti-sun. The nearest any active field's
boresight came to a ≤ 0.1 AU channel position was 17.3° (2.3 FOV
half-diagonals). This is the elongation gate of
`notes/learnings.md` §8 — a quadrature surveyor cannot see channels
B or A-0.1 — confirmed a fourth time, now for a spacecraft.

Two corrections to the plan row fall out:

1. **The universal (Earth-center) list is the wrong observer.** Kepler
   flew an Earth-trailing heliocentric orbit (period 372.5 d) and was
   0.04 AU (2009) → 1.14 AU (2018) from Earth. Matched events shift
   by |Δt_ca| median 33 d / p90 61 d / max 117 d and |Δb| up to
   1 AU (median 0.006 AU = 1.4 R☉) — invalid at every rung, not just
   grazing ones as for TESS/SOHO. A spacecraft product
   `crossings/kepler_v1` (Horizons −227, `xng-07012bdf124d`, 3,352
   events, 0 invalid) was built and is what the intersect used.
2. **"Intersect needs no adapter" was right for the wrong reason** —
   the answer is set by the mission's pointing geometry, so any
   observer at ~1 AU gives the same null. The Earth-center list would
   have produced the same "no hit" with the wrong times.

## Inputs and method (`scripts/footprint_intersect.py`)

| Piece | Source | Notes |
| --- | --- | --- |
| Observer | JPL Horizons `-227` SSB vectors, 6 h, 2009-05-01 → 2018-11-01 (`crossings/observers/kepler_*`, raw + sha-pinned) | `--observer kepler` added to `sglsurvey/crossings.py`; ~11 min run |
| Crossings | `crossings/kepler_v1/events.ecsv` | 1,676 non-sunward events; era scope B 1.2 R☉ 18 ev / 2 targets (gj-1276, wolf-359), B 2.5 R☉ 28 / 4 (+ van-maanen, ross-128 ×1), B 0.1 AU 64 / 7, A 0.1 AU 64 / 7 (gj-1276, gj-908, ross-128, ross-154, teegarden, van-maanen, wolf-359 — the usual family, 9–10 annual events each) |
| Footprints | `K2fov` 8.0.1 (Kepler/K2 GO office; added to `pyproject.toml`) | silicon-level channel polygons, modules 3/7 dead for K2, module 4 dead after C10; 12-px edge padding (the package default, generous); prime field = campaign 1000 at roll 20° + 90° × season, envelope over the four seasons, no dead modules |
| Campaign dates | MAST CAOM TAP `dbo.obspointing` grouped by `sequence_number` (async UWS jobs; the sync endpoint 504s past 60 s) | *actual* timeseries data ranges, ±2 d slack, clipped to the K2fov planned dates ±15 d. One CAOM defect: a C12-era product is tagged sequence 14, so C14's raw range starts 168 d early (57737.9 vs planned 57904) — the clip leaves a conservative 17-d early start; C20 never flew (retired 2018-10-30). Kepler prime: 54953.0–56423.5 (2009-05-02 → 2013-05-11), 212,954 timeseries products |
| Test | per event: fields active within the 0.1 AU flat-chord window of t_ca (±0.5 d otherwise) → `isOnSilicon` of the channel position; diagnostics: boresight separation, anti-sun offset from the spacecraft, sky-only on-silicon for any field at any time | no signal statistic; nothing here is a search |
| Archive check | 14 MAST box counts (±8° = FOV half-diagonal) of Kepler/K2 products around each ≤ 0.1 AU unit position, by campaign | independent of K2fov; every unit box holds 1–4 campaigns' products — the fields pass through these ecliptic positions, at the wrong times |

Snapshots: `runs/kepler-crossings/recon/` (26 MAST responses).
Results: `results/footprint_intersect_v1.json` (per-event rows),
`results/onstar_k2_v1.json`.

## Campaign geometry as seen from Kepler

Sun elongation of the boresight over the actual data range, and the
minimum boresight-to-anti-sun angle. Backward-facing campaigns sweep
~142° → 61°, forward-facing (C9, C16, C17, C19) 61° → 143°; the
mid-campaign boresight is at quadrature (94°–118°).

| Campaign | RA, Dec (deg) | Data MJD | Elongation start → end | Min to anti-sun |
| --- | --- | --- | --- | --- |
| C0 | 98.3, +21.6 | 56728–56805 | 139 → 65 | 41.1° |
| C1 | 173.9, +1.4 | 56808–56890 | 138 → 63 | 42.1° |
| C2 | 246.1, −22.4 | 56893–56972 | 135 → 61 | 45.1° |
| C3 | 336.7, −11.1 | 56977–57046 | 142 → 72 | 37.6° |
| C4 | 59.1, +18.7 | 57061–57132 | 142 → 71 | 37.6° |
| C5 | 130.2, +16.8 | 57139–57214 | 131 → 61 | 49.0° |
| C6 | 204.9, −11.3 | 57217–57296 | 138 → 66 | 42.4° |
| C7 | 287.8, −23.4 | 57300–57382 | 142 → 61 | 38.4° |
| C8 | 16.3, +5.3 | 57392–57470 | 142 → 61 | 37.6° |
| C9 | 270.4, −21.8 | 57501–57572 | 76 → 143 | 36.7° |
| C10 | 186.8, −4.0 | 57582–57651 | 125 → 63 | 55.0° |
| C11 | 260.4, −24.0 | 57656–57730 | 132 → 61 | 48.2° |
| C12 | 351.7, −5.1 | 57738–57817 | 142 → 61 | 37.7° |
| C13 | 72.8, +20.8 | 57820–57901 | 141 → 61 | 38.7° |
| C14 | 160.7, +6.9 | 57889*–57985 | 158* → 70 | 22.2°* (≈ 38° from the true 2017-06-01 start) |
| C15 | 233.6, −20.1 | 57989–58077 | 142 → 61 | 37.6° |
| C16 | 133.7, +18.5 | 58095–58175 | 62 → 143 | 37.1° |
| C17 | 202.5, −7.7 | 58179–58246 | 75 → 143 | 37.1° |
| C18 | 130.2, +16.8 | 58251–58302 | 137 → 88 | 43.4° |
| C19 | 347.3, −4.2 | 58361–58387 | 106 → 130 | 49.9° |
| Prime | 290.7, +44.5 | 54953–56424 | 85–92 throughout | 66.0° |

\* C14 start is the conservative clip of the mis-tagged CAOM range.

## Result

- **On silicon during an active field: 0 events** (of 1,435 with a
  field active at t_ca, any b up to the 1.03 AU list maximum).
- Anti-sun offset of the channel position, max by rung: B 1.2 R☉
  0.26°, B 2.5 R☉ 0.66°, B 0.1 AU 3.45°, A 0.1 AU 3.46°. Kepler's FOV
  half-diagonal is ~7.6°, so a hit needs a boresight within ~11° of
  anti-sun; the closest campaign boresight ever was 22° (C14
  conservative) / 37°.
- Nearest active boresight to a ≤ 0.1 AU channel position: 17.3°
  (ross-128 A, C14, 2017-05-10 — inside the conservative C14 pad;
  with the true 2017-06-01 start the nearest is 29.6°, ross-154 B,
  C16 2018-03-04, then 39.1° van-maanen A C3). Per unit the nearest
  active boresight is 17–72°.
- 9 of the 128 ≤ 0.1 AU events fall in inter-campaign gaps (no field
  active at all).
- **Sky-only (time ignored): K2 did image 5 of the 14 ≤ 0.1 AU unit
  positions and 6 registry stars.** Channel positions on K2 silicon in
  *other* campaigns: gj-1276 B (C14), ross-128 A (C1), ross-154 A
  (C7), van-maanen A (C8), wolf-359 A (C14) — the fields pass through
  the ecliptic positions of these units at quadrature, 1–9 months from
  the crossings. Plus two on-star, wide-b-only: gj-876 A and ez-aqr A
  (C3). These are the K2 light curves of Wolf 359, Ross 128, Ross 154,
  van Maanen's star, GJ 876 and EZ Aqr (see below); EZ Aqr was on C3 silicon but not targeted.

## On-star K2 light curves of registry targets (channel A, wide rung)

`scripts/onstar_k2_check.py` confirms each against MAST (K2 products
within 1′ of the campaign-epoch star position) and computes the
impact parameter during the campaign, b = r sin(elongation):

| Target | Campaign | Data MJD | Elongation at mid | b at mid (AU) | K2 products |
| --- | --- | --- | --- | --- | --- |
| wolf-359 | C14 | 57905–57985 (2017-06-01 → 08-19) | 116° | 0.94 | `ktwo201885041` LC + SC (EPIC position is catalog-epoch, 80″ from the 2017 star position — a 1′ box misses it; PM 4.7″/yr) |
| ross-128 | C1 | 56808–56890 | 103° | 1.02 | `ktwo201518346` (+ `polar` HLSP) |
| ross-154 | C7 | 57300–57382 | 97° | 1.00 | `ktwo215632069`, `ktwo215632123` |
| van-maanen | C8 | 57392–57470 | 98° | 0.97 | `ktwo220434748` LC + SC, `ktwo229228408/520` |
| gj-876 | C3 | 56977–57046 | 112° | 0.91 | `ktwo206019387` LC + SC, `ktwo206019392` (+ `polar`) |
| ez-aqr | C3 | 56977–57046 | 109° | 0.93 | on silicon, **not targeted** (no product within 6′) |

These are 80-day 30-min (some 1-min) continuous light curves of five
registry stars at b ≈ 0.6–1.0 AU — the wide "A 1.0 AU" rung where the
plan (§5.5, radio; BL) already concludes that the archive's own
transient/flare searches are the constraint. Nothing in the crossing
design applies at that b (the beam-axis geometry is 0.6 AU off), and
they are recorded here as a fact, not a queue item.

## What this means

1. **Item 8 closes with no survey.** There is no Kepler/K2 pixel on
   any channel position during any crossing window at any rung;
   the null is geometric (pointing constraint), holds for any
   observer at ~1 AU, and cannot change with a data release.
2. **The 2009–2015 pulse cell stays open** on optical substrates. The
   fast-cadence archives that *could* reach it are those that pointed
   at opposition: the ground-based transient surveys already searched
   (PTF §5.10 is the era-parallel one; cadence minutes-to-days, not
   30 min continuous) and, prospectively, TESS (2018–, done §5.6).
   No further K2-class mission (continuous, 1 AU, ecliptic) observed
   near anti-sun: CoRoT (2007–2012) pointed at the Galactic
   centre/anti-centre eyes at elongation ~90° for the same
   solar-array reason and is expected to fail the same gate; not
   probed, noted for §5.7 only if a reason appears.
3. **Spacecraft crossing lists are cheap** (11 min, one Horizons
   fetch); the observer table pattern (`--fetch-observer`,
   `--observer`) now covers TESS, SOHO, Kepler. Any future spacecraft
   substrate (STEREO HI, PSP/WISPR, Gaia epoch photometry) should
   start from its own list, not Earth-center.
4. A CAOM data-quality note for anyone re-using MAST campaign ranges:
   group-by `sequence_number` ranges must be clipped to the planned
   campaign dates — at least one K2 product is mis-tagged (C12 data in
   sequence 14).
