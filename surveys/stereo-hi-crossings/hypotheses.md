# STEREO-A HI-1 sunward-channels hypothesis freeze v1.0 (Pipeline B)

**FROZEN 2026-09-04**: the user approved every §12 decision as recommended (D1–D9). Drafted the same day after the STEREO-A observer
crossing list (`crossings/stereoa_v1`), the HI-1 reachability recon
(`notes/hi1_recon_2026-09-04.md`) and the census + coverage intersect
(`results/stereoa_census_v1.json`, `results/coverage_v1.json`), before
any frame is fetched at a survey position (three recon frames declared
in §11). Seventh crossings survey of the sunward family; the LASCO
constructions (`surveys/heliospheric-crossings/hypotheses.md`, D1–D9 +
amendments v1.1/v1.2) are inherited, HI-1 substitutions marked
**[HI-1]**. Freeze decisions are listed in §12 for approval.

Input: **`crossings/stereoa_v1`** (`xng-3dd5b764776e`, STEREO-A
spacecraft observer, Horizons −234, 2007-01-01 → 2026-12-01, 7,516
events, 0 invalid) — mandatory: the spacecraft leads Earth by up to
180° on a 0.96 AU orbit, so the Earth-center list mis-states every
crossing epoch by weeks to months (matched events: |Δt_ca| median
51 d, up to the 120-d match limit; 1,241 events have no Earth-center
counterpart within 120 d). This survey executes plan §5.15 item O2: the
second heliospheric substrate for the two **sunward** (direction ×
side) combinations — **S1 = downlink post-lens** (outbound, target
side) and **S2 = uplink past the Sun** (inbound, anti-target side) —
first searched with LASCO (§5.11, 0 candidates, MW–GW-class limits).

## 1. Archive and roles **[HI-1]**

- **RAL/UKSSDC level 1, `14h1A` product** (DN/s/CCDPIX, point-source
  calibration, desmeared, flat-fielded, saturated columns NaN'd; CGI
  tree `lz/L1/a/img/hi_1/YYYYMMDD/`, 2006-12-13 → 2026-08-31, anonymous
  POST fetch, ~4.2 MB float frames at ~2 s each) — the **search
  substrate** for the whole era. Every fetch size-checked (4,219,200
  bytes) and FITS-parsed. **Decision D1.**
- **NASA SSC level 2, `24h1A_br01`** (the same frames minus RAL's
  1-day running lower-quartile background; plain-HTTP daily
  directories, 2006-12-01 → currency) — *fetch fallback* if the RAL
  CGI throttles, *coverage inventory* (file stems identical, listings
  snapshotted, one header per day by HTTP Range), and *chain-validation
  substrate*: recon parity L2/L1 aperture flux = 1.005 ± 0.010 on 668
  stars. Not a second search substrate.
- Level 0 (NRL), L2 11-day, MSB/S10 products, HI-2 — unused. HI-2
  (ε 19°–89°) never contains a 0.1 AU-rung source (ε ≤ 6°).

## 2. Eras and gaps

| item | value | note |
|---|---|---|
| survey era | 2007-01-01 → RAL currency (2026-08-31 at recon) | event list capped 2026-12-01 (Horizons predicted trajectory end 2026-12-18); both extend at the yearly refresh |
| safe-mode / superior-conjunction gap | **2014-08-19 → 2015-11-16** (455 d, no synoptic frames; measured from the day listings) | 19 of the 295 event rows fall in it |
| other gaps ≥ 3 d | 2007-03-02→07, 2007-03-14→27, 2014-07-07→10, 2018-05-24→06-05, 2023-08-07→15 | measured; per-event frame counts carry them |
| pointing eras (per-day headers, 6,675 days) | **east** (CCD centre HPLN −14°) 2006-12 → 2014-08-18 and 2023-08-16 →; **west** (+14°, spacecraft rolled 180°) 2015-11-17 → 2023-08-15; HPLT centre drifts seasonally ±1.9° (excursions to +7° in early 2007, −4.6° in 2023-08) | 6,615 of 6,675 days nominal (|HPLN| within 0.5° of 14°, |HPLT| < 2.5°); the footprint test uses each frame's own centre; **in the west era the arc falls after t_ca** (the direction recedes) — 113 post-t_ca vs 160 pre-t_ca arcs |
| frame format | 234,689 frames listed; ~36/day; N_IMAGES 30, EXPTIME 1199.9 s on every sampled day header | synoptic-sum gate (§7) |

## 3. Channels, rungs, in-era scope (STEREO-A frame)

Beam-radius ladder inherited from the ZTF freeze; flat-chord windows
per event and rung from b_min, v⊥. In-era events
(`results/stereoa_census_v1.json`; covered counts in `results/coverage_v1.json`):

| channel / rung | events / targets | HI-1 | note |
|---|---|---|---|
| S1 ≤ 1.2 R☉ | 84 / 4 | **not_constrainable** | apparent source at ε ≤ 0.7°, inside HI-1's 4° inner edge; ledger rows only |
| S1 ≤ 2.5 R☉ | 105 / 5 | **not_constrainable** | same |
| S1 0.1 AU | 147 / 7 | **searched** | van-maanen (min b 0.09 R☉), gj-1276 (0.36), teegarden (1.06), wolf-359 (1.06), ross-128 (1.50), ross-154 (3.6), gj-908 (11.4); 21 events/target, 19–20 covered |
| S2 0.1 AU | 148 / 7 | **searched** | the same seven systems, star side of the link; 18–20 covered |
| S2 1.0 AU | 1,905 / 88 | **out of scope** (D7) | 16 targets transit HI-1 (289 transits, ~19 d each, ε 24° ↔ 4°): the whole transit is in-beam, so no window-locked construction exists; recorded as a ledger for a future persistent-excess design |

**Apparent sources (both fixed ICRS directions; derived for a general
observer in `scripts/hi_geometry.py`).** S2: the target star, PM-
propagated (parallax ≪ the 72″ pixel). S1 at the 0.1 AU rung: the
relay direction = the star's antipode (the beam passes the Sun
unlensed at b ≫ R☉; the LASCO graze-point construction reduces to it
exactly). In-beam ⇔ b_e(t) ≤ 0.1 AU with b_e the spacecraft's
transverse offset from the Sun–star axis (filled-cone, D2 inherited).

**Visibility geometry (census with the measured pointing eras).**
HI-1A looks 20° × 20° to one side of the Sun (|HPLN| 3.9° → 24.15°):
**east** in 2007–2014 and from 2023-08 (a fixed sky direction
approaches conjunction: the arc is **before t_ca**, the source enters
the beam at ε = 6.0° at t_ca − 5.7 d and leaves the CCD at the 4.05°
inner margin at t_ca − 3.85 d), **west** in the rolled 2015-11 →
2023-08 era (the direction recedes: the arc is t_ca + 3.85 d → + 5.7 d).
Either way **44 h ≈ 66 frames per window** (gj-908, b = 11.4 R☉:
24 h). 274 of 295 event rows are visible (the 21 others sit in the
gap); the coverage intersect against the real frame lists and per-day
headers finds **272 covered** (≥ 10 arc frames; median 65–66 arc
frames, 5th percentile 30). The same sky direction transits the
field for ~19 d on the same side of the arc (ε 24° ↔ 6°), which gives
every source a **star-fixed off-window baseline** (§6; median ~205
frames at |HPLN| 6°–12°) — the capability LASCO lacked.

## 4. Camera and wavelength **[HI-1]**

HI-1: 1024 × 1024 summed frames (30 × 40 s exposures, on-board
cosmic-ray scrubbed, 2 × 2 binned), 0.01998°/px = 72″, PSF FWHM 1.7 px
(recon), 40-min cadence, passband **630–730 nm**. The 532 nm doubled
line is **out of band**; 1064/1550 nm out of band. The survey
therefore tests **broadband leakage only** (declared); hypothesis flux
= monochromatic line converted via the band effective width
(W_eff ≈ 100 nm) as always. The band is red: HI-1 DN/s carry a
+0.60 mag/(B−V) colour coefficient (recon, 668 Hipparcos stars) — the
colour system is measured, not assumed; M-dwarf S2 sources are ~1
mag brighter in HI-1 than their V suggests.

## 5. Duty cycle and temporal models

- **d = 1 chord — primary.** ~66 frames per 44-h arc at 40-min
  cadence, over **21 annual windows per unit** (STEREO-A's synodic
  year is 346 d): resolved light curves × recurrence, as LASCO.
- **Pulse — secondary, every unit** (all windows are full cadence,
  D6): per-frame max statistic in-window, periods ≥ 40 min (a frame is
  a 20-min-effective sum of 30 × 40 s exposures; shorter pulses scale
  flat ×d, declared).
- **Recurrence stack** over the unit's 21 windows — primary statistic.
- **Declared unconstrained:** sub-exposure pulses (< 40 s), schedules
  avoiding the pre-conjunction arc, single-window transmissions below
  per-window depth, the grazing cones (not visible here).

## 6. Detection constructions **[HI-1]**

**Substrate chain (recon-validated):** per frame — header celestial
WCS (RA/DEC-AZP 'A' system, RAL star-fit pointing: 694/698 Hipparcos
stars within 3 px, 0.34 px rms, zero mean offset on the recon frame)
verified per frame by a Hipparcos match (≥ 50 stars, rms ≤ 0.6 px;
translation refinement applied; frames failing are `no_astrometry`)
→ per-frame large-scale background removal (31-px median filter,
frozen) → Gaussian-weighted aperture photometry (r 2.5 px, annulus
5–9 px) at ICRS-fixed positions → per-frame **2-D star ZP** (quadratic
in x, y; Hipparcos V ≤ 9.5 calibrators, colour-corrected with the
frozen 0.60 coefficient re-measured on the ensemble at dev; gate
≥ 100 calibrators, MAD ≤ 0.12 mag) → flux in ZP-normalized units.
**Decision D3.**

**Patches.** Every measured position is a fixed ICRS direction ("sky
patch"): the aperture content (field stars) is *static* through a
window — no star transits, unlike the corona-frame positions of the
LASCO chain. Per event, the **source patch** (S2 star / S1 antipode)
and **8 control patches** at HPLT offsets ±1°, ±2°, ±3°, ±4° from the
source's track (parallel ecliptic-ward tracks, ~50–200 px away,
same elongation at the same time, same frames). **Decision D4.**

**Per-patch series and star-fixed differential.** For each patch:
in-arc epochs (source in-beam, HPLN ∈ [−6.0°, −4.05°], inside the
frame's footprint) and **baseline epochs** = the same patch at
HPLN ∈ [−12°, −6°] (the 5.7 d before the arc, ~200 frames). The
per-patch excess is E = ⟨F_arc⟩ − g·⟨F_base⟩, where g(HPLN, HPLT) is
the **differential flat** — the ensemble ratio of a static star's
aperture flux in the arc region to its flux in the baseline region,
measured from thousands of Hipparcos/Tycho stars in the same frames
(machinery item; expected within a few %; frozen as a measured map,
identical for source and controls). E removes the static field-star
and target-star flux by construction — no stellar template is needed
(LASCO L3 retired here).

**Same-frame differential (LASCO v1.1 rule 2 transferred):** per
epoch, D_src = (F_src − g·B_src) − median_k (F_k − g·B_k) over the
control patches on the same frame (controls: D_k = own term minus the
median of the other seven); this cancels frame-common transients
(CMEs, exposure/ZP jitter, F-corona ridges). Per-event z = mean of
in-arc D over epochs divided by its empirical standard error (robust
scatter of the in-arc D epochs; ≥ 10 usable epochs required).

**Statistics per unit** (unit = target × channel; 3 trials/unit):
1. **S_stack** — Σ z_ev / √n over the unit's events (primary).
2. **S_event** — max z_ev (secondary).
3. **S_pulse** — max single-frame D_src in-arc, standardized by the
   in-arc robust epoch scatter (tertiary; every unit — full cadence).

Decision rule S > max(T, 0) with **T = max over the 8 control
patches** of the identical statistic; expected control crossings
1/9 per trial. S/MAD(controls) reported alongside.

**Veto ladder:** (1) **known-object census** — planets **including
Earth and the Moon** (both enter HI-1A when the spacecraft–Earth
separation is < 24°: 2007–08 and 2022–24), bright asteroids, comets
(Horizons observer `@-234`; SkyBoT with the MPC code C49 if it
supports it, else Horizons for the bright-asteroid list), annotated
per exceedance epoch; (2) **multi-frame persistence** — a real source
persists over consecutive 40-min frames on the fixed patch; single-
frame excesses are never candidates (the S_pulse cell requires ≥ 2
consecutive frames or a sibling-window recurrence, frozen); (3)
**chord shape** — TESS split-half ramp + background-anticorrelation
tests; (4) **stellar-flare class (S2)** — an in-arc pulse on an S2
unit is compared with the star's own off-window pulse distribution
over its 21 baseline transits (~4,000 epochs of the same star): a
candidate only if it exceeds the max off-window pulse; wolf-359 and
ross-154 are catalogued flare stars; (5) **coronal transient
annotation** — CDAW/HELCATS HI-1 CME catalogues at the epoch; (6)
**recurrence** — the stack demands it; any S_event exceedance is
tested on the unit's other 20 windows.

## 7. Quality masks **[HI-1]**

Primary (frame level): N_IMAGES = 30 and |EXPTIME − 1200 s| ≤ 60 s
(synoptic sums only), NMISSING = 0, 1024 × 1024, astrometry verified
(§6). Epoch level: aperture ≥ 90 % finite pixels (NaN'd saturation
columns, telemetry blocks), no planet/Moon/Earth within 10 px of any
patch of the event, and the **bright-body in-FOV veto** (Earth, Moon,
Venus, Jupiter, Mercury inside the footprint → epoch dropped for all
patches; the stray-light halo is frame-wide). Patch level: Tycho-2
VT ≤ 9 within 3 px → the patch is dropped for that event (source
dropped ⇒ event not searchable; controls require ≥ 5 survivors);
the S2 source patch is exempt (the star *is* the source — the
star-fixed differential models it). Fractions measured at the
coverage/dev stages; masks are not tuned after any window-locked
quantity exists.

## 8. Flux scale, injections, completeness

Per-frame 2-D star ZP with the ensemble colour system (§6) sets the
scale; the SSC level-2 product validates the chain on the overlap
(gate: L2/L1 ensemble flux agreement ≤ 0.02 mag median, recon 0.005).
Era stability measured from the ensemble year by year. Completeness =
**stamp-response injections into the real frames** (the established
`sglsurvey.inject.stamp_response` route with the empirical PSF),
spanning flux × window phase × HPLT × duty cycle, normalized per the
C1 rule. **Positive control (decision D5): a bright asteroid's passage
through the in-beam arc** — (4) Vesta / (1) Ceres / (2) Pallas as seen
from STEREO-A (Horizons `@-234` ephemeris + V), through the identical
patch chain along the predicted moving track; recovery ≤ 0.2 mag after
the colour term (asteroids B−V ≈ 0.7–0.8). Secondary in-situ
controls: gj-908 (V 9.0, dev) — its star-fixed differential must be
null across 21 transits; and the Hipparcos ensemble itself.

## 9. Expected sensitivity, declared before search

Recon single-frame 5σ depth (solar colour) V 12.2 at the outer edge,
**V 10.0–10.5 in the in-beam arc** (the F-corona is 5–10× brighter at
ε 4–6°). Declared expectations (systematics-bounded; the stack gains
will be measured): per-window 66-frame stacks ~V 12–12.5; recurrence
stacks (21 windows) ~V 13.5–14.5 optimistic bound. Power scale via the
LASCO anchors (S1 0.1 AU recurrence m90 V 6.2–8.4 ↔ 3–26 GW; S2 0.1 AU
V 6.2–7.6 ↔ 0.7–1.1 GW, 10-m class): **~5–7 mag deeper → downlink
~10–100 MW through the 0.1 AU cone, uplink ~5–20 MW**. Honest framing:
the plan's "sub-MW downlink cell" referred to the grazing cones, which
HI-1's 4° inner edge cannot see — this survey deepens the **0.1 AU
sunward cells by ~2 orders of magnitude** and opens the sunward pulse
cell at 40-min resolution across 21 recurrences; it does not reach the
grazing cones.

## 10. Volume policy (frozen, decision D6)

Full cadence everywhere: per event the arc (44 h, ~66 frames) + the
star-fixed baseline transit (5.7 d, ~205 frames) ≈ 270 frames; the
control patches share the same frames. 272 covered events × ~7.6 d ≈
2,070 frame-days (614 unique arc days) ≈ **70k frames ≈ 300 GB
transient** (per-batch purge, PS1 pattern; ~2 s/frame serial fetch →
~1.5 days at 1 stream, hours at 4). One measurement pass per frame
serves every patch mapped to it. The pulse cell is open on every unit.

## 11. Pre-freeze data-contact declaration

The recon fetched three frames of 2010-06-15 (L1 00:09 and 00:49 UT,
L2 00:09 UT) for service probing and the full-field star check.
Checked against `stereoa_v1` at census: **both L1 epochs lie inside
two in-beam arcs** — gj-1276 S1 (`evt-6a6141430884`, t_ca 2010-06-19,
source at pixel (981, 510), b_e 0.079 AU) and wolf-359 S2
(`evt-d39c468c8a27`, t_ca 2010-06-19, source at (957, 531), b_e
0.086 AU) — both confirmatory units under D8. No quantity was formed
at either position (the star check matched Hipparcos stars and
measured column-strip noise over 64 × 1024-px bands); no window-
locked series exists. **Declared remedy (the ATLAS wolf-359 pattern,
scaled to the contact):** the two frame epochs 2010-06-15 00:09 and
00:49 UT are **excluded from the in-arc series of those two events**
in every construction (2 of ~66 frames each; recorded in the config),
and the recon frames are never re-used as search frames. The frames
stay in `runs/stereo-hi-crossings/recon/` as recon artifacts; every
coverage/search frame is fetched fresh under snapshot discipline.

## 12. Freeze decisions (recommendations marked; freeze on approval)

- **D1** — substrate: RAL level 1 `14h1A` uniform for search; SSC
  level 2 `br01` as fallback mirror + validation. **Amended v1.1 at
  dev (`thresholds.md`): the search runs on level 2 — level 1 carries a
  ridge/edge background bias the controls cannot cancel.** (Alternative: search
  L2 — rejected: an extra processing layer whose running-background
  window we do not control.) **Recommended as stated.**
- **D2** — S1 beam model: filled-cone windows (inherited). **Recommended.**
- **D3** — chain: header WCS verified per frame + 31-px median high-
  pass + 2-D colour-corrected star ZP + star-fixed differential with
  the measured differential flat. **Recommended.**
- **D4** — controls: 8 parallel-track sky patches at HPLT ±1..±4°,
  same frames, same-frame differential. (Alternative: time-shifted
  patches along the source track — rejected: not simultaneous, triples
  the volume.) **Recommended.**
- **D5** — positive control: bright-asteroid arc passage (Vesta/Ceres/
  Pallas, Horizons `@-234`) ≤ 0.2 mag; gj-908 null differential and
  the Hipparcos ensemble as in-situ checks. **Recommended.**
- **D6** — volume: full cadence, arc + 5.7-d baseline per event
  (≈ 340 GB transient); pulse cell on every unit. **Recommended.**
- **D7** — S2 1.0 AU rung out of scope (no window-locked construction;
  ledger kept); S1 grazing rungs `not_constrainable` (ledger rows).
  **Recommended.**
- **D8** — dev/confirmatory split (seed at threshold freeze):
  dev = gj-908 (S1+S2; V 9.0 — the bright-star differential-flat
  test) + ross-154 (S1+S2; V 10.4, flare star, low galactic latitude);
  confirmatory = van-maanen, wolf-359, teegarden, gj-1276, ross-128 ×
  S1/S2 (10 units, the grazing family). **Recommended.**
- **D9** — era policy of §2 (RAL currency 2026-08-31; the 2014–15 gap
  and the east/west pointing eras measured from listings/headers;
  per-frame footprint from the header, so west-era post-t_ca arcs are
  searched identically to east-era pre-t_ca arcs). **Recommended.**
