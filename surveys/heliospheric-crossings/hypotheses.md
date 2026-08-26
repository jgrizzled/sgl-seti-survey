# LASCO sunward-channels hypothesis freeze v1.0 (Pipeline B)

**FROZEN 2026-08-25**: the user approved every §12 decision as
recommended (D1–D9). Drafted the same day after the sunward geometry
study
(`notes/sunward_geometry_study_2026-08-25.md`), the LASCO reachability
recon (`notes/lasco_recon_2026-08-25.md` incl. addendum), and the
SOHO-observer crossing list, before any archive data is touched at a
survey position (three recon frames declared in §11 — none in any
sunward window, verified). Sixth crossings survey; ZTF/PS1/TESS/ATLAS
constructions inherited, LASCO substitutions marked **[LASCO]**.

Input: **`crossings/soho_v1`** (`xng-298d55d0ce6b`, SOHO spacecraft
observer, 1996-01-01 → 2026-10-01) — mandatory; the Earth-center list
mis-states grazing b by up to 0.21 R☉ and wide-rung b by up to 2.7 R☉
(`results/soho_census_v1.json`). This survey executes plan §5.8
item 9: it searches the two **sunward** (direction × side)
combinations that every previous survey declared out of scope —
**S1 = downlink post-lens** (outbound, target side; apparent source =
the solar-limb graze point) and **S2 = uplink past the Sun** (inbound,
anti-target side; apparent source = the target star near conjunction).
The universal out-of-scope declaration is superseded *here only*.

## 1. Archive and roles **[LASCO]**

- **SDAC level 0.5** (`umbra.nascom.nasa.gov/pub/lasco_level05/`,
  1996-02 → 2025-02) + **NRL mirror tail** (`lz/level_05/`,
  → ~2026-06) — the uniform search substrate for the whole era.
  Discovery = dated `YYMMDD/c2|c3/` directory listings (HTML
  snapshotted); fetch = full 2 MB frames, every fetch size-checked +
  FITS-parsed (silent truncation observed live at recon).
- **NRL level 1** (1996 → 2017-08-31) — *validation substrate only*:
  the star-ZP + WCS chain is confirmed against it on the overlap era,
  and the C2/C3 occulter/vignetting masks are derived from its mask
  files. Not a second search substrate (uniformity of the recurrence
  stacks; our per-frame star calibration supplies the physical scale
  on level 0.5 everywhere). **Decision D1.**
- VSO / ESA SSA — alternatives only, unused.

## 2. Eras and gaps

| item | value | note |
|---|---|---|
| survey era | 1996-01-01 → NRL currency (≈ 2026-06; measured at coverage start) | event list capped at 2026-10-01 (Horizons SPK); both extend at the yearly refresh |
| SOHO attitude loss | ~1998-06-25 → 1998-10 (+ recovery to ~1999-02) | encoded as a coverage gap, verified against the trees |
| roll states | 180° flips (post-2003 quarterly + early-mission states) | handled per frame: the star-fit WCS refinement tests both roll hypotheses; header CROTA is a prior, never trusted alone |

## 3. Channels, rungs, in-era scope (SOHO frame)

Beam-radius ladder inherited from the ZTF freeze; flat-chord windows
per event and rung from b_min, v⊥ as always. In-era events
(`results/soho_census_v1.json`):

| channel / rung | events / targets | camera | note |
|---|---|---|---|
| S1 ≤ 1.2 R☉ | 117 / 4 | C2 | van-maanen (min b **0.127 R☉**, annual October family), wolf-359 (0.58), gj-1276 (0.61), teegarden (0.86) |
| S1 ≤ 2.5 R☉ | 153 / 5 | C2 | + ross-128 (1.61) |
| S1 0.1 AU | 216 / 7 | C3 | + ross-154 (3.18 R☉), gj-908 (12.0 R☉) — ~31 annual events/target |
| S2 0.1 AU | 217 / 7 | C3 | same seven systems, star side of the link; per-target min b mirrors S1 |
| S2 1.0 AU | — | — | **out of scope, frozen at draft**: the LASCO-visible portion (ε ≤ 8°) is already the 0.1 AU rung; the ε ≳ 30° twilight skirt is non-exclusive (geometry study §4.3) and not a heliospheric cell |

**S1 beam model (the annulus question, decision D2).** The post-lens
grazing beam is an annulus at the graze radius; whether Earth inside
it (b_e < b_graze) sees flux depends on divergence. Frozen primary =
**filled-cone**: in-beam ⇔ b_e(t) ≤ rung radius (conservative
superset in time; one contiguous window). The annulus alternative is
kept as a *shape discriminator*, not a mask: an annulus-wall beam
predicts a double-passage light curve (ingress arc / egress arc),
used at adjudication.

**Apparent-source positions.**
- S1: the limb graze point at radius = the ray's graze radius, at the
  position angle of the observer's offset from the axis (azimuth is
  preserved through the lens; derived in the geometry study). Both
  are computed per epoch from `soho_v1` geometry + the SOHO table.
  **C2 annulus gate:** the C2 usable inner radius (measured from L1
  masks at coverage stage, expected ~2.0–2.2 R☉) truncates the
  visible graze-radius range; rungs whose windows leave no epoch with
  a visible arc are reported `not_constrainable`, not silently
  dropped.
- S2: the star's PM-propagated position (parallax ≪ the 56″ pixel).
  **Elongation gate:** the star is unobservable when its elongation
  puts it behind the C3 occulter (3.7 R☉ ≈ 59′) — for grazing-depth
  events the search runs on the window wings where
  ε(t) ∈ [occulter, 5.7°]. Every target is V ≥ 9 (gj-908 V 9.0 the
  brightest): **fainter than the C3 single-frame limit → S2 is an
  empty-field excess search** (the wise-0855 shape) — no saturation
  cut exists in this survey.

## 4. Cameras and wavelength **[LASCO]**

C3 Clear (~400–850 nm, 56″/px) carries both 0.1 AU rungs and
**contains the frequency-doubled 532 nm line** — the line cell.
C2 Orange (540–640 nm, 11.9″/px) carries the S1 grazing rungs;
532 nm falls just outside → grazing rungs are **broadband-leakage
cells only** (declared). 1064/1550 nm remain outside every band in
the programme. Hypothesis flux = monochromatic line converted via
band effective width, as always.

## 5. Duty cycle and temporal models

- **d = 1 chord — primary.** At ~12-min cadence a 0.6–1.3 d grazing
  window holds ~70–300 frames: fully resolved crossing light curves
  (the TESS capability) over **~31 annual windows per target** (the
  ATLAS recurrence capability) — this survey is the first with both.
- **Pulse — secondary**: per-frame max statistic in-window, periods
  ~12 min → window length; **full-cadence (grazing C2) units only**
  (D6 subsampling removes the cell on 0.1 AU units — declared).
- **Recurrence stack** across all in-era windows of a unit — primary
  statistic (§6), the survey's headline cell.
- **Declared unconstrained:** sub-exposure pulses (< ~20 s) beyond
  flat ×d scaling; schedules avoiding Earth-crossing windows;
  single-window transmissions below per-window depth.

## 6. Detection constructions **[LASCO]**

**Substrate chain (recon-validated):** per frame — synthetic celestial
WCS (Sun position from the SOHO table + Meeus P-angle + CROTA,
orientation convention as measured: East left, solar north up) →
star-fit refinement (sub-px on 3/3 recon stars; also resolves the
roll-state hypothesis) → per-frame star ZP (field stars; gate ≤ 0.2
mag scatter, the DECam/TESS standard) → forced photometry at the
predicted position through an empirical field-star PSF → per-position
series.

**Detrending (frozen layer, TESS lesson):** per-frame radial-profile
high-pass (the structured F/K-corona is radial to first order) +
per-position temporal running median over off-window epochs (~30 d),
then the empirical variance rescale k = median(r²/v)/0.4549 (floor 1)
per position — the established recipe. **Decision D3.**

**Statistics per unit** (unit = target × channel × rung × camera;
3 trials/searched unit, ATLAS D4 shape):

1. **S_stack** — recurrence-stacked chord-weighted window-locked
   excess over all in-era events (primary).
2. **S_event** — max single-event chord amplitude (secondary).
3. **S_pulse** — max per-frame S/N in-window (tertiary; full-cadence
   units only).

**Controls (decision D4).** Primary = **8 same-radius PA-ring
controls**: the identical statistic at positions rotated about Sun
center at the same solar radius (radially matched to the coronal
background, computed from the *same frames* — no extra volume).
Secondary = temporal pseudo-windows ±23/47/71/97 d with same-rung
exclusion, run on the 1/day off-window series (geometrically
available: 0.6–11.5 d windows vs annual recurrence). Budget: expected
control crossings 1/9 per searched trial, tallied at threshold
freeze.

**Veto ladder:** (1) **known-object census** — near-Sun fields are
crossed by planets, bright asteroids at conjunction, and the Kreutz
sungrazer stream: JPL/Horizons planet ephemerides + SkyBoT where
usable + the Sungrazer comet catalog, annotated per exceedance epoch;
(2) **multi-frame persistence** — a real source persists over
consecutive ~12-min frames along the (slow) predicted track; cosmic
rays — dense in every LASCO frame — do not; single-frame excesses are
never candidates (frozen); (3) **chord/annulus shape tests** — the
TESS split-half ramp test and background-anticorrelation test,
plus the D2 double-passage discriminator; (4) **CME/streamer
annotation** — exceedances within a catalogued CME's position angle
and time (CDAW list) are re-adjudicated against the transient
corona; (5) **recurrence** — the stack demands it; any S_event
exceedance is tested on the unit's other ~30 windows; (6) rate test
vs the z-track prediction (S2; ~0.3–6.5″/day retrograde — sub-pixel
per window at 56″/px, so this discriminates only across windows;
declared weak here).

## 7. Quality masks **[LASCO]**

Primary (frame level): EXPTIME within 2× the camera synoptic norm;
NMISSING = 0 (level-0.5 telemetry-gap blocks excluded via the
level-1-derived block map on the overlap era, header DATAZER
heuristic after 2017-08); the frame's star-fit and ZP gates passed
(≥ 5 matched stars, scatter ≤ 0.2 mag); position outside the
occulter/pylon/vignetting mask with ≥ 0.7 relative response.
Strict (exceedance re-runs): ZP scatter ≤ 0.1, no CME at the
position angle, both roll hypotheses agree. Fractions measured at
coverage stage; masks must not be tuned after any window-locked
quantity exists.

## 8. Flux scale, injections, completeness

Per-frame star ZP (§6) sets the scale — level-0.5 DN and level-1 MSB
both pass through it, and the level-1 overlap era (1996–2017-08) is
the chain's validation set (**gate:** ZP chain agreement ≤ 0.1 mag
median across the overlap sample before any 2017-09+ depth is
quoted). Completeness = **stamp-response injections into the real
frames** (the established `sglsurvey.inject.stamp_response` route —
pixels are available here, unlike ATLAS) with the empirical
per-region PSF, spanning flux × window phase × position angle ×
duty cycle, normalized per the C1 rule (responses normalized at the
calibration reference through the same kernel chain). **Positive
control (decision D5): a Uranus conjunction passage** (V ≈ 5.9,
unsaturated, weeks inside the C3 FOV at a predicted ephemeris rate)
through the identical chain — recovery to ≤ 0.2 mag against its
predicted magnitude; Neptune (V 7.9) as the faint secondary control.

## 9. Expected sensitivity, declared before search

Recon-measured single-frame depth ~V 8.6 (S/N 8, C3 19 s). Declared
expectations (corona-systematics-bounded — the stack gains will be
*measured*, not assumed): per-window stacks (~70–300 frames) ~V
10.5–11.5; recurrence stacks (~31 windows) ~V 12.5–13.5 optimistic
bound. Power scale via the ZTF anchors: **downlink (S1) ~25 MW
single-frame → few-MW per window → 0.3–1 MW recurrence-stacked
through the grazing cones — including, for the first time in any
archive, the van-maanen 0.13 R☉ near-photosphere cone at ~12-min
resolved cadence**; uplink (S2, 10-m-class) ~GW single-frame →
~3–20 MW stacked. This opens the sunward cell at MW-class power; it
is not deep exclusion, and no depth is promised past what measured
coverage and the completeness gates deliver.

## 10. Volume policy (frozen, decision D6)

Grazing rungs (C2): full cadence — ~153 windows × ~160 frames ≈ 25k
frames. 0.1 AU rungs (C3): **deterministic 3-hourly subsampling**
(first frame per UTC 3-h bin per camera; ~92 points across an
11.5 d window — chord-sufficient, pulse cell declared closed on
these units). Off-window baseline: 1 frame/day over ±110 d per
event, shared across a target's events. Total ≈ 100k frames ≈
200 GB transient, per-batch purge (PS1 pattern); no server needed —
plain-HTTP fetch, hours not days.

## 11. Pre-freeze data-contact declaration

The recon fetched three frames of 2010-06-15 (C3 + C2 level 0.5, C3
level 1) chosen for service probing, and ran the star check on a
Taurus field. Verified against `soho_v1`: **2010-06-15 lies inside
zero sunward windows** (0 of the ≤ 0.1 AU event windows contain that
date; computed at census). No survey position was measured; no
window-locked quantity exists. No remedy required; the frames remain
in `runs/heliospheric-crossings/recon/` as recon artifacts and every
coverage-stage frame is fetched fresh under snapshot discipline.

## 12. Freeze decisions (recommendations marked; freeze on approval)

- **D1** — substrate: level 0.5 uniform for search; level 1 =
  validation + masks only. (Alternative: search L1 where it exists —
  rejected for stack uniformity.) **Recommended as stated.**
- **D2** — S1 beam model: filled-cone windows primary; annulus
  double-passage shape as an adjudication discriminator.
  **Recommended as stated.**
- **D3** — detrending: radial high-pass + 30 d temporal median +
  variance rescale k under all statistics. **Recommended.**
- **D4** — controls: 8 same-radius PA-ring controls primary (same
  frames, radius-matched); temporal pseudo-windows secondary on the
  1/day off-window series. **Recommended.**
- **D5** — positive control: Uranus conjunction passage (+ Neptune
  faint secondary) through the full chain, ≤ 0.2 mag gate.
  **Recommended.**
- **D6** — volume: full cadence grazing / 3-hourly 0.1 AU rungs +
  1/day baselines (≈ 200 GB transient); pulse cell restricted to
  full-cadence units. **Recommended.**
- **D7** — S2 1.0 AU rung out of scope (absorbed/non-exclusive, §3);
  S1 1.2 R☉ rung retained but gated on the measured C2 inner radius
  (`not_constrainable` reporting if occulted). **Recommended.**
- **D8** — dev/confirmatory split (seed at threshold freeze):
  dev = ross-154 + gj-908 (wide-b, rung-insensitive machinery
  validation) + ross-128 (shallowest graze, 2.5 R☉ rung only);
  confirmatory = van-maanen, wolf-359, teegarden, gj-1276 — all
  grazing-family units (TESS pattern). **Recommended.**
- **D9** — era policy of §2 (NRL currency measured at coverage
  start; 1998 gap encoded; per-frame roll resolution). **Recommended.**
