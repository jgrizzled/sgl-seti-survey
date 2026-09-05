# PSP/WISPR observer-geometry pass (plan §5.15 O5)

2026-09-05. The prerequisite study named by the sunward geometry study
(§4 caveat 1, 2026-08-25) and carried as ledger item O5: build the
Parker Solar Probe observer crossing list and ask, per encounter,
whether the sunward apparent sources ever enter WISPR's field. It
decides whether WISPR can be adopted as a third heliospheric substrate
after LASCO (§5.11) and STEREO-A HI-1 (§5.16), and on which cells.

**Answer.** WISPR is adoptable — but not for the cell the plan hoped
for. The grazing cones stay inside the 13.5° inner edge at every
encounter (1.2 R☉: ≤ 6.9°, never; 2.5 R☉: 13.3°–13.7° for two
targets in the final-orbit encounters E22–E29, i.e. the rung's outer
few hundredths of a solar radius for ≤ 3 minutes per window — a
knife-edge, `not_constrainable` in practice). The **0.1 AU rungs of
both sunward channels are, by contrast, wide open**: from an observer
at 0.05–0.25 AU the same transverse offsets subtend 14°–65°, so the
in-beam arc of a 0.1 AU window runs straight through WISPR-I (and, from E6 on, WISPR-O) for a median ~23 h per event, on 21 targets
per channel, every encounter, in public L2 data through E27. This is
the same cell STEREO HI-1 constrained at 12–61 MW, now at 5–50× shorter
range from the beam and at V ~ 13–15 depth — the recon item is whether
the inner-field F-corona and stray light at r < 0.1 AU leave that depth
on point sources.

Inputs: `crossings/psp_v1` (`xng-c339df63a608`, PSP observer, 7,448
events, 2018-08-15 → 2026-12-01), `crossings/universal_v1`, WISPR
Data Users Guide v5 (NRL, Sep 2025; `runs/wispr-crossings/recon/`).
Scripts `scripts/psp_geometry.py`, `scripts/psp_census.py` →
`results/psp_census_v1.json`.

## 1. The observer list

`python -m sglsurvey.crossings --fetch-observer psp` (Horizons −96,
SSB ICRF vectors) and `--observer psp --coarse-step-days 0.5`. Two
departures from the 1-AU-class recipe, both forced by the orbit:

- **10-min observer sampling in yearly chunks** (436,321 rows; Horizons
  caps a response at ~90k rows). At the 0.046 AU perihelion PSP moves
  191 km/s and sweeps ~34° of heliocentric longitude per 6 h; the
  6-h default would have been useless there. The chord sagitta over
  one 10-min step is ~2e-6 AU = 0.0003 R☉.
- **0.5-d coarse bracketing** (default 10 d). The observer's star-side
  and anti-star-side minima of one axis are separated by 180° of
  heliocentric longitude, which the final orbits sweep in ~2.2 d
  around perihelion; a 10-d scan would merge or miss them.

The list is a different population from Earth's, not a correction of
it: 3,724 sunward events vs 2,844 Earth-center in the same era;
nearest same-target/direction/side pairs within 120 d have median
|Δt_ca| 62 d and median |Δb| 0.35 AU, and 1,319 have no pair at all.
The dynamics are simple and worth stating because they fix everything
below: PSP crosses every Sun–star axis **twice per orbit** — once on
the star side, once on the antipode side — at whatever heliocentric
radius the orbit has at the star's longitude. Because the line of
apsides barely moves (the Venus flybys shrink the orbit, they do not
rotate it much), **each target's crossing geometry repeats orbit after
orbit**: a star whose direction lies near the perihelion longitude is
crossed near perihelion on one side and near aphelion on the other, at
nearly the same b every time. Impact parameters at the crossing are
r × (the star's offset from the orbital plane), so the deep family's
grazing membership is set by orbital latitude rather than ecliptic
latitude — the grazing family from PSP is wolf-359, ross-128, gj-1111,
gj-251, gj-518, wolf-437 (perihelion side) plus teegarden, van-maanen,
ross-154 at larger r.

Encounters (r < 0.25 AU arcs from the trajectory, the guide's
definition of the regular-observation periods): 29 through
2026-09-04, q = 35.7 R☉ (E1–E3) → 27.9 (E4–E5) → 20.4 (E6–E7) → 16.0
(E8–E9) → 13.3 (E10–E16) → 11.4 (E17–E21) → 9.85 R☉ (E22–E29); arcs
9.8–11.3 d. Dates reproduce the NRL release listing
(`runs/wispr-crossings/recon/nrl_release_listing_2026-09-05.txt`, 27
encounters public, E27 = 2026-03-06 → 03-17).

## 2. Geometry

Same construction as the LASCO tracks module and STEREO `hi_geometry`:
the apparent source is a fixed ICRS direction — the star (S2) or its
antipode (S1) — and its elongation from Sun center is
**ε = atan(b_e / r_along)**, where b_e is the observer's transverse
offset from the axis and r_along its distance from the Sun along it.
This is the whole story of the pass. From 1 AU the 0.1 AU rung gives
ε ≤ 5.7°; from 0.05 AU the same rung gives ε up to 63°.

WISPR frame (Sun direction u, prograde tangential t, orbit normal
n = t × u): ram longitude λ = atan2(s·t, s·u), orbit latitude
β = asin(s·n). Field model from the Users Guide v5 — inner telescope
sunward edge 13.5° from the pointing vector (Sun center whenever
r < 0.25 AU), 40° and 58° square fields, 3° overlap, outer edge
108.5°: WISPR-I λ ∈ [13.5, 53.5], WISPR-O [50.5, 108.5]; ram side and
vertical extent (±20°, ±26.5°) from Vourlidas et al. 2016 — **the side
and the half-heights are recon header items** (L2 headers carry the
S/C HCI velocity and a ZPN WCS). Visibility = in-beam ∧ in-field ∧
inside an encounter, at 1-min sampling.

A structural consequence checked directly (wolf-359 E27 window,
§`results`): the observer's transverse offset at t_ca is out-of-plane
(that is what sets b_min), and as PSP moves prograde the offset
acquires a +t component — so the **post-t_ca half of every window is
on the ram side and the pre-t_ca half is not** (ram-side fraction
0.99 post / 0.02 pre over all visible rows). WISPR sees the second
half of every sunward window, never the first.

## 3. Grazing rungs vs the 13.5° edge — per encounter

Mission floor: q_min = 9.85 R☉ → ε(1.2 R☉) ≤ 6.9°, ε(2.5 R☉) ≤ 14.2°;
the edge needs r_along ≤ 5.0 R☉ (1.2) / 10.4 R☉ (2.5). Only the
2.5 R☉ cone can touch the field, only at the E22+ perihelia, and only
if the axis crossing falls within hours of perihelion.

| target (S1, final orbits E22–E29) | b_min (R☉) | r_along (AU) | t_ca − perihelion | max ram-longitude of the 2.5 R☉ source | window | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| gj-1111 | 0.90 | 0.0472 | −0.16 d | **13.7°** (ε 14.6°, β −5.1°) | 4.9 h | in-field for 2–3 min at b_e 2.42–2.49 R☉, 8 windows |
| wolf-359 | 0.56 | 0.0468 | +0.13 d | 13.3° (ε 14.7°, β +3.2°) | 5.1 h | 0.2° short of the edge |
| ross-128 | 0.71 | 0.0492 | +0.24 d | 12.1° (ε 14.3°, β +3.8°) | 5.2 h | short |
| gj-251 / wolf-437 / gj-518 | 1.7–2.1 | 0.051–0.056 | ±0.4 d | 10.2° / 6.5° / 5.1° | 3–4 h | short |
| teegarden | 0.75 | 0.100 | −1.3 d | 7.2° | 10.6 h | crossing 1.3 d before perihelion, r too large |
| van-maanen | 2.41 | 0.207 | −3.7 d | 1.2° | 6.3 h | crossing 3.7 d before perihelion |
| ross-154 | 0.48 | 0.238 | +4.6 d | 2.5° | 26 h | crossing 4.6 d after perihelion |

Earlier encounters are strictly worse (larger q). Over the mission:
111 S1 events ≤ 1.2 R☉ (5 targets, 94 inside encounters), max edge
bound 6.8°, **0** in-field; 190 events ≤ 2.5 R☉ (9 targets, 170 inside
encounters), 16 with the rung-edge bound ≥ 13.5°, **8** with any
in-field minutes (all gj-1111, 0.03–0.05 h each). The graze point at
ε ≈ 14° from a 10 R☉ observer is also, physically, a point on the
corona at 2.4 R☉ transverse distance — the WISPR-I sunward edge is a
baffle-limited, F-corona-dominated strip. Verdict: **grazing cones
`not_constrainable` by WISPR**, as for HI-1; the sub-MW-through-the-
grazing-cone hope has no heliospheric-imager substrate and stays with
LASCO C2 (§5.11).

## 4. The 0.1 AU rungs — the cell WISPR does open

| channel, 0.1 AU rung | events | targets | events with WISPR-I arc | targets with any | WISPR-I hours / event (median) | WISPR-O events |
| --- | --- | --- | --- | --- | --- | --- |
| S1 downlink post-lens | 1,513 | 58 | 533 | 21 | 25 (p10–p90 9–58) | 479 |
| S2 uplink past the Sun | 1,701 | 64 | 561 | 21 | 22.5 both channels (9–50) | 627 |

Whole-list: 1,094 event-rung windows with a WISPR-I arc, 28,000 h.
Visible arcs sit at ε 14°–54° (WISPR-I) and 50°–108° (WISPR-O), at
b_e median 5.5 → 16.7 R☉ across the arc (the rung is 21.5 R☉), at
heliocentric distances 0.046–0.25 AU. Per encounter: 28 visible
event-rungs in E1–E3 (q 35.7 R☉, WISPR-I only) rising to 42 + 53
(I + O) from E22 on; ~680–1,110 WISPR-I hours per encounter.

The seven-system family of the LASCO/HI-1 surveys, final-orbit
geometry (repeats every encounter from E22; earlier encounters differ
only in r):

| target | S1 (outbound, star side) | S2 (inbound, antipode side) |
| --- | --- | --- |
| van-maanen | r 0.207 AU, t_ca −3.7 d: I 36 h + O 16 h (window straddles the arc) | r 0.054 AU, +0.4 d: I 21 h + O 8 h |
| wolf-359 | r 0.047 AU, +0.1 d: I 11 h + O 17 h | r 0.54 AU: outside encounters, never |
| teegarden | r 0.100 AU, −1.3 d: I 14 h + O 11 h | r 0.076 AU, +0.9 d: I 32 h (ε 14°–40°, no O) |
| gj-1276 | r 0.54 AU: never | r 0.047 AU, +0.1 d: I 11 h + O 17 h |
| ross-128 | r 0.049 AU, +0.2 d: I 14 h + O 14 h | r 0.34 AU, −7.8 d: I 74 h + O 23 h (far end of the window inside the arc) |
| ross-154 | r 0.238 AU, +4.6 d: none (source stays inside 13.5°) | r 0.053 AU, −0.4 d: I 7 h + O 13 h |
| gj-908 | r 0.31 AU, −6.8 d: I 68 h + O 22 h (far end of window) | r 0.050 AU, +0.3 d: I 14 h + O 13 h |

Every family member except ross-154-S1 / wolf-359-S2 / gj-1276-S1 has
a WISPR arc on the 0.1 AU rung at every encounter; 21 targets per
channel overall. Note the two kinds of arc: the **perihelion-side
crossings** (r ≈ 0.05 AU, t_ca within ±1 d of perihelion) give
10–30 h arcs at b_e 2.5–18 R☉ — the deep, near-axis part of the cone,
reached within hours of t_ca; the **straddling windows** (t_ca 4–8 d
before perihelion at r 0.2–0.4 AU) give 40–80 h arcs at b_e
12–20 R☉ — the outer part of the cone during the inbound leg.

## 5. Recon items (the survey stage, if adopted)

1. **Depth on point sources in-field vs r**: WISPR-I's F-corona and
   stray-light floor at r < 0.1 AU vs the guide's L2/L2b/L3 products
   (L3 = background-subtracted K-corona; LW an alternative
   background); the claimed V ~ 13–15 applies to the outer field.
   The relevant arcs sit at ε 14°–54°, r 0.05–0.25 AU.
2. **Header verification** of the field model: ram side, vertical
   half-heights, ZPN WCS, per-frame pointing (`HCIi_OBS`, `HCIi_VOB`
   in the headers; the pointing vector is Sun center within 0.25 AU).
3. **Cadence and campaign structure**: ~half the observations are
   high-cadence narrow-field programs (no L3), separate from the
   full-field synoptic series; coverage fraction of the ~23-h arcs.
4. **Star-fixed baseline**: the source is a fixed ICRS direction but
   PSP's perspective changes fast; the HI-1 19-d star-fixed baseline
   has no analogue — an encounter-local baseline design is needed.
5. **Positive control**: planets and comets transit WISPR every
   encounter (the guide lists them as standard sources); Mercury /
   Venus in-field passages replace Pallas.
6. Data: anonymous NRL tree (`wispr.nrl.navy.mil/wisprdata`, per-
   encounter FITS by level), public through E27 (2026-03).

## 6. Hand-offs

- `crossings/psp_v1` is the observer list for any PSP instrument, not
  only WISPR; refresh with the yearly list refresh (Horizons SPK
  currently ends 2026-12-01, the list stop).
- The ε = atan(b_e / r) result generalises: any inner-heliosphere
  imager (Solar Orbiter SoloHI, 5°–45° at 0.28–0.9 AU) opens the same
  0.1 AU cell; SoloHI is the natural fourth substrate and needs only
  `--observer solo` (Horizons −144) on this machinery.
- The grazing-cone question is closed for heliospheric imagers:
  LASCO C2 remains the unique substrate.
