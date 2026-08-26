# Sunward-channel geometry study (plan §5.8 item 9)

2026-08-25. Answers the question frozen at the ZTF crossings freeze
(2026-08-23): the two sunward (link direction × side-of-axis)
combinations were declared out of scope — are any of them **only**
visible at small solar elongation? If yes, heliospheric imagers are
the unique archival substrate and a recon is justified; if no, the
question closes for the record.

**Answer: yes.** Every frozen beam-radius rung of both sunward
combinations is confined to apparent-source elongation ≲ 6° — outside
every night-sky archive in the programme — except the outer skirt of
the widest (1 AU) uplink rung. Coronagraphs and heliospheric imagers
are the unique archival substrate. Recon recommended (LASCO first).

Inputs: `crossings/universal_v1` (`xng-a09e2db7681d`), Earth-center,
1980→2028, all 16,586 minima. Script
`surveys/heliospheric-crossings/scripts/sunward_geometry.py` →
`results/sunward_geometry_v1.json`.

## 1. Taxonomy and derived geometry

`side == target` ⇔ `axis_distance_au > 0` (Earth on the star side of
the Sun); the four combinations partition the list:

| combo | (direction, side) | n | apparent source | elongation at t_ca |
| --- | --- | --- | --- | --- |
| A (searched) | inbound, target | 4,169 | the star, near opposition | ≥ 92.2°, median 150° |
| B (searched) | outbound, anti_target | 4,124 | relay = antipode, near anti-Sun | ≥ 91.7°, median 152° |
| **S1** | outbound, target | 4,169 | **solar-limb graze point** | ε ≈ b_graze/1 AU |
| **S2** | inbound, anti_target | 4,124 | **the star, near conjunction** | ε ≈ arcsin(b/1 AU) |

- **S1 — downlink post-lens.** Earth on the star side sits in the
  relay's beam *after* it grazed the Sun at b_graze (1.2–2.5 R☉ frozen
  rungs) and was bent toward the target star. Post-lens the beam is a
  nearly collimated annulus of radius ~b_graze (a grazing ray's
  deflection displaces it only ~1,300 km per AU of travel), so the
  crossing condition reuses the same grazing-rung events, and the
  apparent source is the graze point at **1.2–2.5 R☉ from Sun center
  (19′–40′) regardless of Earth's in-beam offset**. The 0.1 AU
  pointing-margin rung passes the Sun effectively unlensed and arrives
  from ≤ ~6°. S1 can never appear in a night-sky archive.
- **S2 — uplink past the Sun.** Earth on the antipode side intercepts
  the star→relay beam after it passes the Sun; the apparent source is
  the star itself at ε ≈ arcsin(b_e/1 AU): ≤ 5.7° for the 0.1 AU rung.
  Occultation core: b < 1 R☉ puts the star behind the photosphere —
  144 events / 3 targets at t_ca. The 1 AU rung ("1 m-class
  transmitter, effectively any minimum") spans ε up to ~88°; its outer
  skirt (ε ≳ 30–40°) is ordinary evening/morning sky — the one
  non-sunward-exclusive corner (see §4).

Validation against the frozen list: computed source elongations match
the analytic forms with median |residual| 0.30° over all 16.6k events;
restricted to the sunward b ≤ 0.1 AU population (673 events) the
match is exact to ≤ 0.032°, all at ε ≤ 3.29° at t_ca. The A/B sources
never drop below 91.7° — confirming the two searched channels and the
two sunward ones cleanly partition the sky into night-side and
Sun-side populations.

## 2. Sunward event census (Earth-center, indicative)

Events / targets per frozen rung; eras are archive start → list end
(imager parameters from general knowledge, *unprobed* — recon items).

| rung | full 1980–2028 | LASCO era (1996–) | STEREO-HI era (2007–) | WISPR era (2018–) |
| --- | --- | --- | --- | --- |
| S1 ≤ 1.2 R☉ | 192 / 4 | 128 / 4 | 83 / 4 | 37 / 4 |
| S1 ≤ 2.5 R☉ | 240 / 5 | 160 / 5 | 103 / 5 | 46 / 5 |
| S1 ≤ 0.1 AU | 336 / 7 | 224 / 7 | 145 / 7 | 64 / 7 |
| S2 ≤ 0.1 AU | 337 / 7 | 224 / 7 | 145 / 7 | 64 / 7 |
| S2 ≤ 1 AU | 4,124 / 85 | 2,743 / 85 | 1,783 / 85 | 803 / 85 |

The grazing family is the familiar one — van-maanen, wolf-359,
teegarden, gj-1276 (+ ross-128 at 2.5 R☉). Note the scale: **128
LASCO-era S1 grazing windows at ≤ 1.2 R☉** vs channel B's 35 in the
whole ZTF era — 30 years of continuous minutes-cadence coverage over
the same family whose observable-side events were lost to PS1 chip-gap
masks and TESS sector gaps. The recurrence-stacked construction frozen
for ATLAS (`surveys/atlas-asassn-crossings/hypotheses.md`) transfers
directly.

## 3. Substrate map and sensitivity reality check

| imager | elongation window | era / cadence | point-source depth |
| --- | --- | --- | --- |
| SOHO/LASCO C2 | 1.5–6 R☉ (0.4°–1.6°) | 1996– / ~12–20 min | bright stars only (V ≲ 6–7) |
| SOHO/LASCO C3 | 3.7–30 R☉ (1°–8°) | 1996– / ~12–20 min | V ~ 8 |
| MLSO Mk3/Mk4/K-Cor | 1.05–3 R☉ | 1980– (ground, daytime) | coronal SB-limited |
| STEREO-A HI-1 | 4°–24° | 2007– / 40 min | stellar photometry to V ~ 13 |
| STEREO-A HI-2 | ~19°–89° | 2007– / 2 h | V ~ 12 |
| PSP/WISPR | 13.5°–108° from PSP | 2018– / encounters | V ~ 13–15 |

Mapping: the S1 graze-point source (1.2–2.5 R☉) falls in C2 / K-Cor
territory (the 1.2 R☉ rung sits at C2's occulter edge); the 0.1 AU
rungs of both combos (ε ≤ 5.7°) fall in C3 and the inner edge of HI-1.

Rough power scalings from the ZTF-crossings anchors (order-of-magnitude
only, **not** thresholds): channel B's median m90 ≈ 21.8 ↔ ≳ 130 W
through the 2.5 R☉ cone scales to ~40 MW at a C3-like V ≈ 8 limit and
~0.4 MW at an HI-1-like V ≈ 13; the channel-A uplink anchors
(56–400 kW, 10 m-class) scale to ~10–100 MW (C3) and ~0.1–1 MW (HI-1).
Framing per the plan: these archives **open the sunward cell at
MW-class power** with uniquely dense window coverage — not deep
exclusion.

## 4. Caveats recorded

1. **Observers.** Earth-center is exact for the geometry answer but
   not for a grazing-rung freeze: SOHO's halo orbit has ~R☉-scale
   cross-track amplitude (TESS lesson: |Δb| up to 0.33 R☉ mattered) —
   an L1/SOHO-observer crossing list via
   `register_spacecraft_table_observer` precedes any freeze. STEREO-A
   drifts along the 1 AU orbit (its crossing *dates* differ entirely;
   similar annual statistics); WISPR is an inner-heliosphere observer
   (0.05–0.7 AU) whose graze-point elongation ≈ b_graze/r_helio only
   enters its 13.5° inner FOV edge near closest perihelia — a
   dedicated PSP-observer geometry pass is prerequisite to adopting it.
2. **S1 beam profile.** Whether Earth *inside* the post-lens annulus
   (b_e < b_graze) sees flux depends on divergence filling the annulus
   interior — a hypothesis-freeze question (affects window/statistic
   definitions, not the sunward conclusion).
3. **S2 1 AU outer skirt.** The ε ≳ 30–40° portion of the widest
   uplink rung is night-sky visible in principle (twilight sky, star
   near conjunction) and is *not* exclusive to this substrate; recorded
   as a possible future extension of a channel-A-style blended search
   (e.g. ATLAS), low priority — widest rung, contrast-limited, worst
   airmass.
4. Imager FOVs, cadences, depths, and archive access above are from
   general knowledge; every number is a recon-stage verification item.

## 5. Recommendation

Positive answer → per plan §5.8 item 9, archive recon is justified.
Suggested order: **LASCO** (longest era — 30 yr × ~2 windows/yr/target
on the grazing family; near-Earth observer; NRL/SDAC archives well
served) → **STEREO HI-1** (depth V ~ 13 → sub-MW downlink cell;
needs its own observer list) → **WISPR deferred** pending the PSP
observer-geometry pass. Queue position relative to remaining §5.8
items is a planning decision; nothing here blocks on the server
migration.
