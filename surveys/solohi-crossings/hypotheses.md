# Solar Orbiter / SoloHI sunward-channels hypothesis freeze v1.0 (Pipeline B)

**FROZEN 2026-09-06**: the user approved every §12 decision as recommended (D0 run, D1–D9).

Drafted 2026-09-06 after the Solar Orbiter observer crossing list
(`crossings/solo_v1`), the geometry pass + reachability recon
(`notes/solo_solohi_geometry_pass_2026-09-06.md`, one document — the
recon facts are its §5) and the census-level coverage
(`results/solo_census_v1.json`), before any frame is fetched at a
survey position beyond the recon frames declared in §11. Fourth
heliospheric substrate: the LASCO constructions
(`surveys/heliospheric-crossings/hypotheses.md`, D1–D9 + v1.1/v1.2),
the HI-1 substitutions (`surveys/stereo-hi-crossings/hypotheses.md`)
and the WISPR chain with its amendments v1.1–v1.5
(`surveys/wispr-crossings/`) are inherited; SoloHI substitutions
marked **[SoloHI]**. Freeze decisions D0–D9 in §12 for approval.

Input: **`crossings/solo_v1`** (`xng-dea7100a725c`, Solar Orbiter
observer, Horizons −144 at 10 min, 2020-05-01 → 2028-01-01, 4,624
events, 0 invalid, 98 degraded) — mandatory: Solar Orbiter is at
0.28–1.0 AU and crosses each axis twice per 150–230 d orbit; the
Earth-center list is a different population. This survey executes
plan §5.23: the two sunward combinations **S1 = downlink post-lens**
(outbound, target side) and **S2 = uplink past the Sun** (inbound,
anti-target side) on the **0.1 AU rung**, from an observer for which
that rung subtends 6°–20° (ε = atan(b_e / r)).

## 1. Archive and roles **[SoloHI]**

- **SoloHI level 2** (`solo_L2_solohi-{1,2}ft_*.fits`: the two inner
  tiles, MSB, calibrated — bias, linearity, vignetting, exposure
  normalisation; F-corona **in**) — the **search substrate**. Two
  identical public copies: the Solar Orbiter Archive (SOAR TAP
  `v_sc_data_item` + data service, the **inventory of record**,
  193,900 L2 items, snapshot 2026-09-06) and the NRL tree
  (`solohi.nrl.navy.mil/so_data/L2/YYYYMMDD/`, plain HTTP, Range
  honoured — the **fetch route**, 3,960,000 bytes per frame,
  size-checked and FITS-parsed on every fetch). **Decision D1.**
- **No background-subtracted science product exists** (the NRL `MOS/`
  "L3a" mosaics are 8-bit display products): the static content is
  handled by the chain (§6), with the **off-beam star-fixed
  baseline** that exists for every arc here.
- **Outer tiles** (`3ft/3fg`, `4ft/4fg`, ε 25°–45°) — *baseline only*:
  the 0.1 AU sources never reach them in-beam (ε ≤ 20.9°); the
  star-fixed baseline of each arc runs through them.
- L1, the 2022 `11t`–`42t` sub-field programs, `21s` — unused.

## 2. Eras and gaps

| item | value | note |
|---|---|---|
| survey era | P04 2022-01 → P12 2026-04-10 (release currency) | event list to 2028-01-01; P13 (perihelion 2026-08-19) unreleased — refresh item; the 2021-04 commissioning week and 2021-12 cruise data precede P04 |
| orbits | 9 covered (P04–P12); q 0.323 (P04), 0.293–0.294 (P05–P13); period 200 → 168 d; SoloHI observes at every r | geometry repeats orbit to orbit until a Venus flyby steps the inclination (2022-09, 2025-02) |
| inner-tile cadence | 24 min median (p10 12, p90 48) both tiles; gaps > 6 h: 68 per tile over the era | measured from the SOAR inventory; annual June–July gap (conjunction) |
| observing modes | `HI_SYN_NEAR` 16–29 s inside 0.4 AU, `MID` 35–100 s, `FAR` 65–145 s outside 0.6 AU; 5-exposure sums, HIGH gain | 125 sampled headers; the ZP absorbs exposure |
| pointing | 83 % of frames at the nominal orbit-plane field; 2–4° offsets 2025-11 → 2026-03; roll campaigns (4 of 125) elsewhere | **per-frame WCS footprint** decides membership (§7) |
| ephemeris in headers | DSUN_OBS agrees with Horizons to 1e-5; HCI speed to 0.02 km/s | all geometry from the Horizons table; only the celestial WCS from headers |

## 3. Channels, rungs, in-era scope (Solar Orbiter frame)

| channel / rung | events / targets (2020-05 → 2028-01) | SoloHI | note |
|---|---|---|---|
| S1 ≤ 1.2 R☉ | 16 / 5 | **not_constrainable** | graze source ≤ 1.1° from Sun centre vs the 5.2° edge; ledger rows |
| S1 ≤ 2.5 R☉ | 25 / 7 | **not_constrainable** | ≤ 2.4°; ledger rows |
| **S1 0.1 AU** | 268 / 24 | **searched**: 242 event-rungs with an arc on 24 targets, 138 covered; **19 units** (≥ 3 covered events of ≥ 10 frames) | arcs −75 → −30 h before t_ca, ε 14.5° → 7.2°, b_e 21.5 → 13 R☉ (p10 7.5) |
| **S2 0.1 AU** | 182 / 20 | **searched**: 160 / 17 targets, 89 covered; **13 units** | same kind |
| S2 1.0 AU | 1,186 / 85 | out of scope (D7) | 36 stars transit the tiles in 214 covered star-orbit pairs; the whole transit is in-beam |

**Apparent sources** (fixed ICRS directions, `scripts/solo_geometry.py`,
validated against the frame WCS): S2 = the target star
(PM-propagated), S1 = the star's antipode. In-beam ⇔ b_e(t) ≤ 0.1 AU
(filled cone, D2 inherited). **Anti-ram rule**: the pre-t_ca half of
every window is on the anti-ram side (SoloHI's side) and the post-t_ca
half never is — every arc is pre-t_ca, the complement of WISPR's.
**Seam rule**: sources within ±0.25° of the orbital plane fall in the
0.49° gap between tiles 1 and 2 (ross-128 P05–P09, both channels:
10 event-rungs, `seam` rows in the ledger).

**Unit population** (`results/solo_census_v1.json`; census-level —
the coverage intersect after the freeze re-derives it from per-frame
WCS footprints and may drop units; covered = ≥ 3 events with ≥ 10
inner-tile frames inside the geometric arc): **32 units** — S1:
61-vir, ez-aqr, gj-1002, gj-1111, gj-1276, gj-518, gj-581, gj-588,
gj-667-c, gj-674, gj-682, gj-783, gj-876, gj-908, ross-128, ross-154,
teegarden, van-maanen, wolf-1061, wolf-359 (20 with ≥ 3 covered
events; gj-518 has 1 — 19 units); S2: 61-vir, ez-aqr, gj-1002,
gj-1111, gj-1276, gj-251, gj-876, gj-908, ross-128, ross-154,
teegarden, van-maanen, wolf-359 (13). 217 events, ~29,200 inner-tile
frames (median 114 per arc, p10–p90 26–254). The LASCO/HI-1/WISPR
seven-system family is inside it (14 units).

## 4. Camera and wavelength **[SoloHI]**

Inner tiles: 960 × 1024 (2 × 2 binned), 0.0206°/px = 74″, PSF FWHM
1.5–2.2 px (110–160″), 5 summed exposures (XPOSURE 16–145 s scaling
with r), passband **~475–740 nm** (WAVELNTH 540 nm; the WISPR-class
broadband visible filter — the exact SoloHI curve from Howard et al.
2020 is recorded at dev). The 532 nm doubled line is **in band**;
1064/1550 nm out of band: hypothesis flux = monochromatic line via the
band effective width (W_eff ≈ 250 nm), broadband leakage declared for
the out-of-band lines. Colour system measured: +0.30 to +0.47
mag/(B−V) (recon, 14 frames), re-measured on the ensemble at dev.

## 5. Duty cycle and temporal models

- **d = 1 chord — primary.** 26–254 frames per arc at 12–48 min,
  **3–9 orbits per unit**: resolved light curves × recurrence. The
  per-orbit geometry repeats between flybys, so stacks are
  like-with-like within an inclination group (P04–P05 / P06–P09 /
  P10–P12).
- **Pulse — secondary, every unit**: per-frame max in-arc, periods
  ≥ 12–24 min (coarser than WISPR's 5 min; LASCO 12, HI-1 40).
- **Recurrence stack** over the unit's covered orbits — primary.
- **Declared unconstrained:** sub-exposure pulses (< 16–145 s), the
  post-t_ca (ram) half of every window, the near-axis part of the
  cone (b_e < 0.09 r_along, 6–8 R☉ near perihelion), the grazing
  cones, the seam rows, P13+ windows, the June–July gaps.

## 6. Detection constructions **[SoloHI]**

**Substrate chain (recon-validated, WISPR v1.1/v1.5 form):** per
frame — header celestial WCS (`RA/DEC-ZPN` 'A') **refined by a
Hipparcos match** (translation + 2 × 3 affine on ≥ 30 stars within
3 px; rms ≤ 1.5 px; frames failing → `no_astrometry`) → 25-px median
high-pass on L2 → top-hat aperture (r 3 px, annulus 5–9 px) at
ICRS-fixed positions → per-frame **2-D star ZP** (calibrators S/N ≥ 15,
≥ 30, linear surface below 60 / quadratic above; colour-corrected with
the ensemble coefficient; ZP-scatter gate 0.35 mag; ZP-uncertainty gate
per WISPR v1.5) → flux in ZP-normalised units. **Decision D3.**

**Patches.** The **source patch** (S2 star / S1 antipode) and **8
control patches at orbit-latitude offsets ±1.5°, ±3°, ±4.5°, ±6°**
from the source direction (parallel tracks: same ram longitude at the
same time, same frames, 70–290 px away; a one-sided ladder when the
source is within 6° of a tile's latitude edge, or of the seam —
recorded). Control interpolation **quadratic in orbit latitude** (the
WISPR v2 lesson: the median biases near-plane sources on the F-corona
ridge). **Decision D4.**

**Static-content removal — two routes, both available [SoloHI].**
(a) **Star-fixed differential (HI-1 pattern), primary**: every arc is
preceded by an off-beam in-field baseline (median 135 covered hours
in the 10 d before the arc, ≥ 24 h for 224 of 227 arcs, at ε 15°–45°
in the inner then outer tiles); E' = F_arc − ⟨F_base⟩ with the
baseline restricted to the same tile where ≥ 6 h exist there (median
63 h), else the outer-tile baseline with the tile-to-tile ZP offset
measured on the ensemble. (b) **Stellar template (LASCO/WISPR
route), secondary and check**: the predicted aperture flux of every
Gaia DR3 star (G ≤ 15, proper-motion propagated to the epoch, through
the measured colour system and encircled-fraction curve; the WISPR
v2 lesson in place of Tycho-2/Hipparcos; the S2 target star included
from the Simbad table) — E = F − (a + b·T) with the template response
(a, b) fitted on the controls (WISPR v1.4). The two excesses are
reported side by side; **S is formed on E' (route a)**; a unit whose
baseline fails the 6-h rule on > half its events falls back to (b),
declared per unit before any window-locked quantity.

**Same-frame differential (LASCO v1.1 rule 2):** D_src = E'_src −
Q(E'_k) over the 8 controls on the same frame, Q the quadratic-in-
latitude interpolation at the source latitude (controls: own E' minus
the interpolation of the other seven); per-event z = mean of in-arc D
over epochs / its empirical standard error (robust in-arc scatter;
≥ 10 usable epochs).

**Statistics per unit** (unit = target × channel; 3 trials/unit):
1. **S_stack** — Σ z_ev / √n over the unit's covered events (primary).
2. **S_event** — max z_ev (secondary).
3. **S_pulse** — max single-frame D_src in-arc, standardised by the
   in-arc robust epoch scatter; two-frame persistence (WISPR v1.5) at
   adjudication (tertiary; every unit).

Decision rule S > max(T, 0), **T = max over the 8 control patches** of
the identical statistic; expected control crossings 1/9 per trial.

**Veto ladder:** (1) **known-object census** — Mercury, Venus, Earth,
Mars, Jupiter, Saturn, bright asteroids and comets (Horizons `@-144`),
proximity veto to any of the event's nine patches (Mercury/Venus 10°,
Earth/Jupiter 5°, WISPR v1.3), annotated per exceedance epoch; (2)
**multi-frame persistence** — ≥ 2 consecutive frames on the fixed
patch, or sibling-orbit recurrence; (3) **chord shape** — TESS
split-half ramp + background-anticorrelation tests; (4) **stellar-
flare class (S2)** — an in-arc pulse against the star's own
baseline-subtracted distribution over its other orbits; (5) **coronal
transient annotation** — CME/streamer passages (SoloHI event lists,
HELCATS) at the epoch and position; (6) **recurrence** — the stack
demands it; any S_event exceedance is tested on the unit's other
orbits; (7) **static-content check** — a Gaia-predicted stellar flux
matching the excess to ~20 % is a systematic (WISPR adjudication
rule, now pre-registered).

## 7. Quality masks **[SoloHI]**

Frame level: inner-tile full-frame synoptic products (`1ft`, `2ft`;
960 × 1024, NBIN 4, NSUMEXP 5), OBS_MODE / XPOSURE / GAINMODE / SC_ROLL
recorded (the ZP absorbs them), astrometry refined (§6), **membership
by the frame's own WCS** (the patch inside DSTART/DSTOP with a 12-px
margin — the nominal field model is not used after coverage). Epoch
level: aperture ≥ 90 % finite pixels; proximity veto (§6); the recon
epochs dropped for the events they touch (§11). Patch level: Gaia
G ≤ 6 within 4 px → patch dropped for that event (bright-star wings),
S2 source patch exempt (the baseline/template models the star); a
dropped source patch removes the event; controls need ≥ 5 survivors.
S2 units whose target has 3 ≤ V < 8 use the self-calibrated template
route as the check (WISPR v1.3; 61-vir); V < 3 none in the
population. Fractions measured at coverage/dev; masks are not tuned
after any window-locked quantity exists.

## 8. Flux scale, injections, completeness

Per-frame 2-D star ZP with the ensemble colour system (recon: −19.87
to −20.02 across 16–145 s exposures, r 0.3–1.0 AU; the predicted
exposure `shi_get_framexpdur` is a ≤ 0.1 mag systematic the per-frame
ZP absorbs). Era stability measured orbit by orbit. Completeness =
**stamp-response injections into the real frames**
(`sglsurvey.inject.stamp_response`, empirical PSF), spanning flux ×
arc phase × orbit latitude × duty cycle, normalised per the C1 rule.
**Positive controls (decision D5):** (a) **61-vir S2** — V 4.7 in the
source aperture, the baseline/template null across 6 orbits (dev);
(b) **gj-908 S2** (V 9.0, dev) as in LASCO/HI-1/WISPR; (c) a **moving
body through the arc band** chosen at dev from Horizons `@-144` — the
brightest of Vesta/Ceres/Pallas or a catalogued comet with V ≤ 8 in
an inner-tile arc band, through the identical patch chain along the
predicted track, recovery ≤ 0.2 mag; Mercury/Venus saturate and are
veto bodies.

## 9. Expected sensitivity, declared before search

Recon single-frame 5σ depth in the arc band (solar colour, L2,
ε 7°–14.5°): **V 6.9–8.8 at r 0.30 AU (16 s), 7.4–9.4 at 0.385 AU
(29 s), 8.0–10.2 at 0.62 AU (65 s), ~10 at 1.0 AU (145 s)** —
F-corona-gradient-limited inward of ε ~12°; 1.5–2 mag deeper than
WISPR-I per frame, LASCO C3 class, 1–2 mag shallower than HI-1.
Declared expectations (systematics-bounded; gains measured): per-arc
stacks of 26–254 frames ~V 9–11; recurrence stacks over 3–9 orbits
~V 10–12 optimistic bound. Power scale via the HI-1 anchors
(recurrence m90 V_eq 11.9–12.9 ↔ downlink 12–29 MW, uplink 3–61 MW)
and the WISPR result (m90 V_eq 7.8–12.0 ↔ 0.08–3.9 GW downlink):
**~V 10–12 → downlink ~30–300 MW, uplink ~10–300 MW through the
0.1 AU cone**, LASCO-class. Honest framing: SoloHI does **not** deepen
the 0.1 AU sunward cells; proximity to the beam buys visibility, not
flux. What it adds: (i) the **pre-t_ca half of every window** (WISPR
saw only the post-t_ca half; LASCO/HI-1 the Earth-side windows at a
different phase entirely), (ii) the outer cone b_e 7.5–21.5 R☉ from
0.3–0.9 AU, (iii) an independent fourth instrument on the same cells,
with a star-fixed baseline WISPR lacked, (iv) 32 units — first
sunward 0.1 AU constraints for gj-588, gj-674, gj-682, gj-876 (S1),
ez-aqr (S1), gj-251 (S2) and the first from a heliospheric imager for
several more, (v) a second in-band 532 nm sunward substrate.

## 10. Volume policy (decision D6)

Full cadence: every inner-tile L2 frame inside the covered arcs
(~29k frames ≈ 115 GB transient, per-batch purge) plus the same-tile
baseline frames (~10 d × 60 frames per arc, capped at the nearest 200
frames ≈ 45k frames) and the outer-tile baseline where the same-tile
rule fails; one measurement pass per frame serves every patch mapped
to it. NRL Range/HTTP ~1–2 s per frame → ~20–40 h serial, ~8 h at
4 streams; measurement ~1 s per frame. The pulse cell is open on
every unit.

## 11. Pre-freeze data-contact declaration

Recon fetched 48 full frames (4 tiles × 12 epochs: 2021-12-27,
2022-03-26, 2022-10-10, 2023-04-18, 2023-10-06, 2024-04-02,
2024-10-01, 2025-03-20, 2025-03-28, 2025-09-15, 2026-02-01,
2026-04-05), 125 `1ft` frames at ~10-d spacing for header sampling
(`runs/solohi-crossings/recon/sample/`), and one L3a mosaic. Checked
against the census rows at coverage: any recon frame inside a covered
arc or baseline is **excluded from that event's series** (the config
lists them) and never re-used as a search frame; nothing was measured
at any survey position (full-field Hipparcos match, elongation-binned
noise, PSF from V < 6.5 stars only). Every survey frame is fetched
fresh under snapshot discipline.

## 12. Freeze decisions (recommendations marked; freeze on approval)

- **D0 — go/no-go.** Run the survey as the fourth substrate at
  LASCO-class depth, or table it as "geometry established; the same
  cell WISPR searched, at the other phase". **Recommended: run** —
  the pre-t_ca half of the windows, the first star-fixed baseline for
  a sunward inner-heliosphere substrate, and 6+ new targets are new
  cells; the cost is one chain on public data.
- **D1** — substrate: SoloHI L2 inner tiles, SOAR inventory of
  record, NRL fetch route; outer tiles baseline only. (Alternative:
  build our own F-corona model from the ensemble — rejected: the
  star-fixed baseline + high-pass make it unnecessary at this depth.)
  **Recommended.**
- **D2** — S1 beam model: filled-cone windows (inherited).
  **Recommended.**
- **D3** — chain: affine-refined header WCS + 25-px high-pass + 2-D
  colour-corrected star ZP (WISPR v1.1/v1.5 gates) + star-fixed
  differential primary, Gaia template secondary + same-frame
  quadratic-in-latitude differential. **Recommended.**
- **D4** — controls: 8 parallel-track patches at orbit-latitude
  offsets ±1.5 … ±6°, quadratic interpolation, one-sided ladders at
  edges and the seam. **Recommended.**
- **D5** — positive controls: 61-vir S2 and gj-908 S2 nulls (dev),
  plus a Horizons-selected moving body through an arc band (≤ 0.2
  mag); stamp-response injections. **Recommended.**
- **D6** — volume: full cadence over every covered arc plus capped
  baselines (~75k frames, ~300 GB transient); pulse cell on every
  unit. **Recommended.**
- **D7** — S2 1.0 AU out of scope (ledger: 214 star-orbit transits);
  S1 grazing rungs `not_constrainable`; seam rows `seam`.
  **Recommended.**
- **D8** — dev/confirmatory split (seed at threshold freeze): **dev =
  gj-908 S1+S2, 61-vir S1+S2, ross-154 S1** (5 units / 15 trials —
  the two template/baseline tests, the flare star at its best S1
  geometry, a near-perihelion unit); **confirmatory = the other 27
  units** (81 trials, 9.0 expected control crossings), with the
  seven-system family (11 remaining units) reported as the
  cross-substrate subset. (Alternative: confirmatory = family only —
  rejected as before.) **Recommended.**
- **D9** — era policy of §2: P04–P12 (release currency 2026-04-10);
  per-frame footprint from the header WCS; geometry from the Horizons
  table only; the nominal field model retired after coverage.
  **Recommended.**
