# Solar Orbiter / SoloHI observer-geometry pass (plan §5.23)

2026-09-06. The hand-off named by the PSP/WISPR geometry pass (§5.18)
and survey (§5.19): build the Solar Orbiter observer crossing list on
the generic spacecraft machinery and ask, per orbit, whether the
sunward apparent sources — the star (S2) or its antipode (S1) at
elongation ε = atan(b_e / r_along) — enter SoloHI's field, and what
the public archive covers. It decides whether SoloHI can be adopted as
the fourth heliospheric substrate after LASCO (§5.11), STEREO-A HI-1
(§5.16) and WISPR (§5.19), and on which cells.

**Answer.** SoloHI is adoptable, on the same cell as WISPR and no
other. The grazing cones never come near the field: from q ≥ 0.281 AU
(60 R☉) the 1.2 / 2.5 R☉ graze sources sit at ε ≤ 1.1° / 2.4° against
a 5.2° inner edge — `not_constrainable` by a factor 2–5 in angle, not
WISPR's knife-edge. The **0.1 AU rungs of both sunward channels are
open every orbit**: from 0.28–0.9 AU the rung subtends 6°–20°, and the
in-beam arc of a 0.1 AU window runs through the two inner detector
tiles for a median 42 h per event (the pre-t_ca half — the field is on
the anti-ram side, the opposite of WISPR), on 24 S1 and 17 S2 targets,
27–29 event-rungs per orbit. The public L2 archive (SOAR/NRL,
2021-12 → 2026-04, near-continuous, 12–24-min cadence in the inner
tiles) covers 227 of those arcs over nine orbits P04–P12 with a
median 114 frames per arc: **32 units (19 S1 + 13 S2) with ≥ 3 covered
events of ≥ 10 frames, 217 events, ~29,000 frames.** Per-frame depth
in the arc band (ε 7°–14.5°) is V 7–9 near perihelion and V 8–10
farther out — 1.5–2 mag deeper than WISPR per frame, LASCO C3 class,
with 4–5× fewer frames per event and 9 orbits against 27 encounters.
Two structural facts new to this substrate: the arcs see only the
**outer part of the cone** (b_e 7.5–21.5 R☉; the source clears the
5.2° edge only when b_e > 0.09 r_along), and the seam between the two
inner tiles lies on the orbital plane, so the deepest in-plane
crossings (ross-128, b 0.13 R☉, P05–P09) fall in a 0.49° gap.

Inputs: `crossings/solo_v1` (`xng-dea7100a725c`, Solar Orbiter
observer, 4,624 events, 2020-05-01 → 2028-01-01; 98 degraded),
`crossings/universal_v1`, the SOAR L2 inventory
(`runs/solohi-crossings/coverage/soar_l2_inventory.json`, 193,900
items, snapshot 2026-09-06), 48 + 125 L2 frames
(`runs/solohi-crossings/recon/`). Scripts `scripts/solo_geometry.py`,
`scripts/solo_census.py` → `results/solo_census_v1.json`;
`scripts/soar_inventory.py`, `scripts/sample_headers.py`,
`scripts/recon_star_check.py` → `results/recon_star_check_v1.json`,
`results/recon_depth_vs_eps_v1.json`.

## 1. The observer list

`python -m sglsurvey.crossings --fetch-observer solo` (Horizons −144,
SSB ICRF vectors; the ESA SPK `solo_ANC_soc-orbit-stp_20200210-20301120`
runs to 2030-11-20, tracking fit to 2026-08-30, prediction after) and
`--observer solo --coarse-step-days 0.5`, with the PSP settings
(10-min sampling in yearly chunks, 403,345 rows; 0.5-d coarse
bracketing). Era 2020-05-01 (SoloHI first light, cruise-phase
checkouts) → 2028-01-01 (the universal window end). At the 0.28 AU
perihelion the spacecraft moves 70 km/s and sweeps ~2° per 6 h; the
chord sagitta over a 10-min step is ~5e-8 AU. Build time 39 min.

The list is a different population from Earth's: 2,313 sunward events
vs 2,692 Earth-center in the same era; nearest same-target /
direction / side pairs within 120 d have median |Δt_ca| 48 d and
median |Δb| 0.20 AU (p95 0.62 AU), and 910 have no pair at all. The
dynamics are the PSP ones at larger radius: Solar Orbiter crosses
every Sun–star axis **twice per orbit** at whatever heliocentric
radius the orbit has at the star's longitude, and because the
Venus-resonant orbit keeps its line of apsides (the flybys raise the
inclination, they do not rotate the orbit), each target's crossing
geometry repeats orbit after orbit — until a flyby changes the
inclination, when the impact parameters step (ross-128: b 0.13 R☉ in
P05–P09, 11 R☉ from P10 after the 2025-02 Venus flyby; gj-1111 S2:
b 5.3 → 3.8 → 2.6 R☉ across P05–P09 / P10–P13 / P14–P16).

Perihelia from the trajectory (orbit = aphelion to aphelion):

| orbits | perihelion dates | q (AU) | Q (AU) | period (d) |
| --- | --- | --- | --- | --- |
| P01–P03 (cruise) | 2020-06-15, 2021-02-10, 2021-09-12 | 0.516, 0.495, 0.587 | 0.74–1.02 | 165–233 |
| P04–P05 | 2022-03-26, 2022-10-12 | 0.323, 0.293 | 1.02, 0.95 | 200, 190 |
| P06–P09 | 2023-04-10 → 2024-09-30 | 0.293 | 0.954 | 179.8 |
| P10–P13 | 2025-03-31 → 2026-08-19 | 0.294 | 0.90 | 176 → 168.5 |
| P14–P16 | 2027-02-06, 2027-07-06, 2027-12-03 | 0.281 | 0.82 → 0.58 (era end) | 162 → 150 |

## 2. Geometry and the field model (header-measured)

Same construction as `psp_geometry`: the apparent source is a fixed
ICRS direction and its elongation from Sun centre is
ε = atan(b_e / r_along). From 0.3 AU the 0.1 AU rung gives ε up to
18.4°; from 0.9 AU, 6.3°.

Frame: u = Sun direction, t = prograde tangential (Horizons
velocity), n = t × u (orbit normal); ram longitude λ = atan2(s·t, s·u),
orbit latitude β = asin(s·n). **The field model is measured, not taken
from the literature**: 48 L2 frames at 12 epochs (2021-12 → 2026-04,
r 0.29–1.01 AU, spacecraft roll −8° → +13°) mapped through the
celestial ('A', RA/DEC-ZPN) WCS into this frame. The four tile edges
are constant-λ / constant-β lines to ~0.1° and **stable in the
orbit-plane frame at every epoch** (in the solar-north /
helioprojective frame the same edges wander by ±5° with the roll: the
spacecraft rolls to keep the mosaic on its orbital plane). The field
is on the **anti-ram side** — λ < 0, helioprojective longitude −5° to
−45°, solar east — the opposite side from WISPR:

| tile | λ (deg) | β (deg) | descriptors | cadence (median gap) |
| --- | --- | --- | --- | --- |
| 1 (inner, south) | −25.6 → −5.2 | −19.3 → −0.2 | `solohi-1ft` | 24 min (p10 12, p90 48) |
| 2 (inner, north) | −24.3 → −5.1 | +0.2 → +20.5 | `solohi-2ft` | 24 min |
| 3 (outer, north) | −45.0 → −24.7 | +0.3 → +19.3 | `solohi-3ft` / `3fg` | 24–48 min |
| 4 (outer, south) | −45.1 → −26.0 | −20.5 → −0.1 | `solohi-4ft` / `4fg` | 24–48 min |

Seams between adjacent tiles are **0.47°–0.49° (23–24 binned px) with
no overlap** (tile-edge pixels mapped through the WCS: nearest
approach 0.487° between tile 1's top row and tile 2's bottom row, 0 of
240 edge points of one tile inside the other). The inner seam lies on
the orbital plane at β ∈ (−0.2°, +0.2°). Inner edge 5.2° from Sun
centre on the plane (rising to ~19° at the tile corners), outer edge
45°; 960 × 1024 2×2-binned pixels of 0.0206° = 74″; 5-exposure sums.

Pointing is not always nominal: of 125 `1ft` frames sampled every
~10 d, 104 have the tile centre within 1.5° of (λ −15.1°, β −9.9°);
~10 (2025-11 → 2026-03) are offset by 2–4°; 4 are roll campaigns
(SC_ROLL −244°, −113°, −67°, +88°) with the field elsewhere — in one
case on the ram side. The census uses the nominal tiles; **the survey
footprint test must use each frame's own WCS** (the WISPR rule).

A structural consequence checked on every visible row (§`results`,
`anti_ram_fraction`): as the observer moves prograde toward an axis
crossing, the source's λ runs from −ε to 0, so the **pre-t_ca half of
every window is on the anti-ram side (fraction 1.00) and the post-t_ca
half never is (0.00)**. SoloHI sees the first half of every sunward
window, WISPR the second — the two substrates are complementary in
phase.

## 3. Grazing rungs vs the 5.2° edge

Mission floor: q_min = 0.281 AU = 60.4 R☉ → ε(1.2 R☉) ≤ 1.14°,
ε(2.5 R☉) ≤ 2.37°; the edge needs r_along ≤ 13.2 R☉ (1.2) / 27.5 R☉
(2.5), a factor 2–5 inside any Solar Orbiter perihelion. Over the
era: 16 S1 events ≤ 1.2 R☉ (5 targets: ez-aqr, gj-667-c, ross-128,
ross-154, wolf-359; min r_along 0.295 AU), 25 ≤ 2.5 R☉ (7 targets,
adding gj-1276, gj-876); the best edge bound is gj-667-c
2025-03-30 (b 0.88 R☉, r_along 0.295 AU): ε_max 1.1° / 2.3°. **0** with
any in-field minute. Verdict: **grazing cones `not_constrainable`**;
LASCO C2 remains the unique substrate, as after WISPR.

## 4. The 0.1 AU rungs — the cell SoloHI opens

| channel, 0.1 AU rung | events (era) | targets | with a SoloHI arc | targets with any | arc hours / event (median, p10–p90) | covered by public L2 | units (≥ 3 events of ≥ 10 frames) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| S1 downlink post-lens | 268 | 24 | 242 | 24 | 42 (19–54) | 138 | 19 |
| S2 uplink past the Sun | 182 | 20 | 160 | 17 (never: gj-581, gj-682, wolf-1061) | 42 | 89 | 13 |

Whole list: 450 event-rungs, 402 with an arc (all in the **inner
tiles** — tile 1 97 / tile 2 130 of the covered rows; the source's
ε ≤ 20.9° never reaches tiles 3/4), 227 covered, 8,700 covered hours.
The arc is 28 % of the in-beam window (median window 150 h). The
source sits at ε = atan(b_e / r_along) with b_e shrinking toward t_ca,
so it **enters the field at the rung edge (b_e 21.5 R☉, ε 7°–20°) and
leaves through the 5.2° inner edge (b_e 7.5–19 R☉) on its way to
t_ca**: median arc −75 h → −30 h before t_ca, ε 14.5° → 7.2° (extremes
20.9° → 5.1°), b_e 21.5 → 13.4 R☉ (p10 7.5 R☉) at r_helio 0.30–0.76 AU
(median 0.40). The near-axis part of the cone (b_e < 0.09 r_along,
i.e. 6–8 R☉ near perihelion) is inside the inner edge from every
radius — SoloHI searches the outer half of the 0.1 AU cone, where
WISPR's perihelion-side arcs reached b_e 2.5 R☉.

Per orbit: 27–29 event-rungs with an arc and ~1,100 arc-hours in
P04–P16 (12–21 in the cruise orbits P01–P03); public L2 covers
**P04–P12** (2022-03 → 2026-03; 17–28 covered arcs and 670–1,100
covered hours per orbit; P13 = 2026-08 not yet released; the 2021-04
commissioning week and 2021-12 cruise data precede P04). Frames per
covered arc (inner-tile `Nft` files inside the geometric arc): median
114, p10–p90 26–254, 29,535 in all.

**Off-beam baselines exist for every arc** (the HI-1 star-fixed
pattern that WISPR lacked). Because the source enters the field from
larger elongation, it is in-field and off-beam (b_e > 0.1 AU) for the
days before the arc: in the 10 d preceding each covered arc, a median
155 h in-field (63 h in the same inner tiles at ε 15°–25°, the rest
in the outer tiles out to ε 45°), **135 h covered by L2 frames**
(p10 80 h; 224 of 227 arcs have ≥ 24 h). The baseline sits at larger
ε than the arc (a different F-corona level, the same detector), so
the star-fixed differential and the stellar template are both
available and check each other.

The seven-system family of the LASCO / HI-1 / WISPR surveys:

| target | S1 (outbound, star side) | S2 (inbound, antipode side) |
| --- | --- | --- |
| van-maanen | r 0.50–0.72 AU, b 3–9 R☉: 13/13 arcs, 9 covered, median 47 h, 8–307 frames | r 0.38–0.86, b 2–14 R☉: 11/12 arcs, 9 covered, 39 h, 17–137 frames |
| wolf-359 | r 0.51–0.97, b 0.8–8.8: 14/15, 8 covered, 44 h, 28–218 | r 0.34–0.52, b 0.6–6.6: 16/16, 9 covered, 52 h, 61–277 |
| teegarden | r 0.60–0.88, b 4–6: 9/9, 4 covered, 46 h, 39–114 | r 0.32–0.72, b 1.6–5.4: 12/13, 9 covered, 39 h, 68–231 |
| gj-1276 | r 0.34–0.52, b 1.3–5.9: 16/16, 9 covered, 53 h, 49–276 | r 0.51–0.97, b 1.7–7.9: 14/14, 8 covered, 46 h, 25–229 |
| ross-128 | r 0.44–0.93: 6/15 arcs (P05–P09 at b 0.13 R☉ in the seam), 4 covered, 38 h, 21–189 | r 0.38–0.67: 11/16 (P05–P09 in the seam), 4 covered, 49 h, 60–237 |
| ross-154 | r 0.28–0.59, b 1.1–4.2: 16/16, 8 covered, 43 h, 92–228 | r 0.81–0.93, b 3–10: 9/9, 4 covered, 8 h, 10–31 |
| gj-908 | r 0.38–0.67, b 5–12: 15/16, 9 covered, 55 h, 36–285 | r 0.43–0.92, b 4–10: 13/14, 7 covered, 48 h, 30–241 |

Every family member has covered arcs in both channels (14 of the 32
units); the other 18 units are gj-1002, gj-1111 (S2 8 events, S1
5 short 7-h arcs at r 0.77), gj-251 S2, gj-518 S1 (1), gj-581 S1,
gj-588 S1, gj-667-c S1, gj-674 S1, gj-682 S1, gj-783 S1, gj-876 S1+S2,
61-vir S1+S2, ez-aqr S1+S2, wolf-1061 S1 — 12 targets with their first
sunward 0.1 AU constraint from a heliospheric imager that WISPR did not
reach (gj-588, gj-674, gj-682, gj-784, gj-876 S1, ez-aqr S1, …), and
the WISPR-searched targets at a second phase of the window.

The S2 1 AU rung is out of scope as before (every observer position
inside 1 AU of the axis is in-beam): ledger of 36 stars with a covered
tile transit, 214 star-orbit pairs.

## 5. What the recon measured (go/no-go; `recon_star_check_v1.json`)

14 frames, every tile, r 0.30 / 0.385 / 0.62 / 1.01 AU:

- **Header celestial WCS ≈ 1 px**: 345–905 Hipparcos matches within
  3 px per frame, rms 0.8–1.2 px, mean offset (+1.2 … +1.5, ±0.8) px
  in the inner tiles (a constant ~110″ shift, refinable), < 1 px in
  the outer tiles. Header DSUN_OBS agrees with the Horizons table to
  1e-5; header HCI speed to 0.02 km/s.
- **Photometry**: ZP (V = ZP − 2.5 log₁₀ ΣMSB, 2-px aperture) −19.87
  to −20.02 across 16–145 s exposures and r 0.3–1.0 (exposure
  normalisation holds to 0.1 mag); colour coefficient +0.30 to +0.47
  mag/(B−V); scatter 0.09–0.20 mag MAD on 250–390 calibrators; per-V-bin
  residuals ≤ 0.1 mag from V 5 to 10 in the inner tiles. The outer
  tiles saturate bright stars (V < 6.5 residuals −0.6 to −1.4 mag at
  49–550 s exposures) — a bright-star class, irrelevant to the arcs.
- **PSF** FWHM 1.5–2.2 px (110–160″).
- **Single-frame 5σ depth (solar colour, 2-px aperture, 25-px median
  high-pass on L2)** binned by elongation over the whole inner tile
  (`recon_depth_vs_eps_v1.json`) — the decisive numbers, at the arc
  band ε 7°–14.5° and at ε 20°–26° for comparison:

| r (AU), mode, XPOSURE | ε 7° | 9° | 11° | 13° | 14.5° | 20°–26° |
| --- | --- | --- | --- | --- | --- | --- |
| 0.300, SYN_NEAR, 16 s | V 6.9 | 7.5 | 8.2 | 8.6 | 8.8 | 9.6–9.8 |
| 0.385, SYN_NEAR, 29 s | 7.4 | 8.1 | 8.8 | 9.2 | 9.4 | 10.2–10.4 |
| 0.620, SYN_FAR, 65 s | 8.0 | 8.9 | 9.6 | 10.0 | 10.2 | 11.0–11.2 |
| 1.007, SYN_FAR, 145 s | 10.2 | 9.7 | 10.0 | 10.3 | 10.3 | 10.8–11.4 |

  The noise is F-corona-gradient-limited inward of ε ~12° (background
  ×3–5 from ε 14° to 7°). Exposures scale with r (`HI_SYN_NEAR` 16–29 s
  inside 0.4 AU, `MID` 35–100 s, `FAR` 65–145 s outside 0.6 AU; all
  HIGH gain, 5-exposure sums; 125 sampled headers). So **per frame
  V 7–9 near perihelion, V 8–10 outside 0.5 AU** in the arc band —
  1.5–2 mag deeper than WISPR-I per frame (V 5.5–8), LASCO C3 class.
- **Products**: L2 only on SOAR (`v_sc_data_item`, anonymous TAP +
  data service, 4 MB per file) and on the NRL tree
  (`solohi.nrl.navy.mil/so_data/L2/YYYYMMDD/`, plain HTTP, Range
  honoured, 0.2 s header reads). The NRL `MOS/` "L3a" mosaics are
  8-bit display products; **no background-subtracted science product
  exists** — the WISPR L3 shortcut is unavailable and the static
  content has to be handled by the chain (LASCO pattern: high-pass +
  stellar template + same-frame differential). L2 is calibrated
  (bias, linearity, vignetting, MSB) with the exposure predicted
  (`shi_get_framexpdur`), a 0.1-mag-class systematic seen in the ZP.

## 6. Hand-offs

- `crossings/solo_v1` is the observer list for any Solar Orbiter
  instrument; refresh with the yearly cycle (SPK to 2030-11; the era
  can extend past 2028 with the universal window).
- The SoloHI survey (§5.23) is recon-complete: hypotheses and a
  threshold construction follow in `hypotheses.md` for the freeze.
- Complementarity: WISPR (ram side, post-t_ca, b_e 2.5–18 R☉, r < 0.25
  AU) and SoloHI (anti-ram, pre-t_ca, b_e 7.5–21.5 R☉, r 0.3–0.9 AU)
  sample the two halves of the same 0.1 AU windows from two vantage
  points; on the seven-system family both exist for every unit.
- The inner-tile seam on the orbital plane is a general SoloHI
  limitation for in-plane sources (any b_min ≲ 0.3 R☉ crossing).
