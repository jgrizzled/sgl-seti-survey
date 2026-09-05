# STEREO-A HI-1 sunward-channels crossings survey

**0 candidates.** Second heliospheric substrate for the sunward
(direction × side) beam combinations first searched with SOHO/LASCO,
using STEREO-A Heliospheric Imager HI-1 frames 2007–2026 on a
spacecraft-observer crossing list. Blind confirmatory: **27 searched
trials over 9 units, 3 exceedances vs 3.0 expected control
crossings, all adjudicated non-point-source systematics (two
coronal-transient fronts, one background-offset plateau); every
recurrence stack null.** Recurrence-stack depth V_eq 11.9–12.9 →
**downlink 12–29 MW through the 0.1 AU cone, 10-m-class uplink
3–61 MW** — two to three orders of magnitude below the LASCO limits
on the same cells. Survey docs `surveys/stereo-hi-crossings/`;
execution record `notes/project_history.md` §19; plan §5.16.

## 1. Hypotheses and geometry

The two sunward combinations excluded by declaration from every
night-sky survey (plan §5.11): **S1 — downlink post-lens** (outbound,
target side; the relay's beam after passing the Sun, apparent source
= the relay direction = the star's antipode, a fixed ICRS point) and
**S2 — uplink past the Sun** (inbound, anti-target side; apparent
source = the target star near conjunction). Observer:
`crossings/stereoa_v1` (Horizons −234, 7,516 events; the spacecraft
leads Earth by up to 180° on a 0.96 AU orbit, so the Earth-center
list mis-states every crossing epoch by weeks to months — 1,241
events have no Earth counterpart within 120 d). Era 2007-01 →
2026-08 (RAL currency). No z-grid.

HI-1A images a 20° × 20° band on one side of the Sun (|HPLN|
3.9°–24°, 72″/px, 40-min summed frames, 630–730 nm). The **grazing
rungs (ε ≤ 0.7°) are never visible** (`not_constrainable`, 84 / 105
ledger events); the **0.1 AU rungs** of both channels are visible for
an arc of ε 4.4°–6.0° — 44 h ≈ 66 frames per window, 21 windows per
unit at the 346-d synodic year — **before** t_ca when the camera looks
east (2007–2014-08, 2023-08–) and **after** t_ca in the rolled
west-looking era (2015-11 → 2023-08; found from the per-day header
inventory, 6,675 headers). The same sky direction transits the field
for ~19 d on the arc's side: every source has a **star-fixed
off-window baseline** the LASCO corona-frame construction could not
have. The 532 nm doubled line is out of band — a broadband-leakage
cell only. S2 1 AU (whole HI-1 transit in-beam; 16 targets, 289
transits) is out of scope by decision D7. Coverage: 272 of 295 event
rows (the 2014-08 → 2015-11 safe-mode/conjunction gap takes 19).

## 2. Chain

Level-2 `24h1A_br01` frames from the NASA SSC mirror (RAL's level 1
minus the per-pixel 1-day running lowest-quartile background — the
static F-corona exactly, while sources drift 54 px/day; amendment
v1.1), per frame: header celestial WCS verified against Hipparcos
(≥ 50 stars, ≤ 0.6 px rms; 0.44 px on ~700 stars typical) → 31-px
median high-pass → 2-D colour-corrected star zero point (+0.595
mag/(B−V); ~600 calibrators; MAD 0.054 mag) → top-hat aperture
photometry (r 2.5 px, annulus 5–9 px) at ICRS-fixed **sky patches**:
the source and 8 controls on parallel tracks at HPLT ±1°…±4°, same
frames. Per patch the static field is removed by the patch's own
baseline transit (|HPLN| 6°–12°, ~205 frames; ensemble differential
flat 0.999 from 693k calibrator records), then the same-frame
differential D_src = E_src − median_k E_k. Per event z = mean(D)/SE;
S_stack = Σz/√n (primary, 18–20 recurrences), S_event = max z,
S_pulse = max single-frame excess (every unit full cadence). Rule
S > max(T, 0), T = the max over the 8 control patches of the
identical statistic. Vetoes: planet/Earth/Moon 10-px proximity,
bright-body 2° proximity (v1.2), Tycho-2 VT ≤ 9 patch mask (S2
source exempt), synoptic-sum gate (30 × 40 s), the recon-frame
exclusions. 38,830 usable confirmatory frames.

**Amendments (dev, before any confirmatory pixel;
`surveys/stereo-hi-crossings/thresholds.md`).** v1.1: on level 1 the
median high-pass residual ramps to +16 DN/s across the last 25 CCD
columns and the annulus background is curvature-biased on the
F-corona ridge, where every source patch sits while its controls sit
on the flanks (145 vs 75 DN/s/px) — a ridge-dependent bias the
controls cannot cancel (first events z = −2 … −10); level 2 removes
it (same patch, same frame: −5.8 → +0.9 units); 24-px edge margin.
v1.2: the frame-wide bright-body veto emptied half of gj-908's
events; measured scatter with a body > 2° away equals clean-epoch
scatter (0.44 vs 0.41), ≤ 2° doubles it → 2° proximity veto. Recorded,
not amended: the ZP-scatter gate (MAD ≤ 0.12) rejects ross-154's
crowded west-era fields (0.12–0.17 mag blending scatter) — its dev
units keep 9–10 of 20 events.

## 3. Results (blind confirmatory, 2026-09-04)

**0 candidates** — 27 trials over 9 units, 3 exceedances vs 3.0
expected; all S_stack null (max 2.88 vs T ≥ 3.72;
`results/confirmatory_v1.md` has the unit table and adjudications):

| unit | n_ev | S_stack S/T | S_event S/T | S_pulse S/T |
|---|---|---|---|---|
| gj-1276 S1 | 20 | 0.11 / 9.66 | 4.06 / 7.75 | 5.73 / 14.65 |
| gj-1276 S2 | 18 | −0.21 / 17.64 | **8.22 / 7.75** | **15.71 / 11.16** |
| ross-128 S1 | — | constraint-only (antipode Tycho-blended, VT 7.2 at 1.7 px) | | |
| ross-128 S2 | 19 | −9.16 / 12.41 | 3.43 / 37.20 | 8.43 / 73.08 |
| teegarden S1 | 17 | −6.66 / 3.72 | 3.58 / 5.78 | **19.57 / 7.71** |
| teegarden S2 | 18 | 1.54 / 10.94 | 6.45 / 10.91 | 10.59 / 29.22 |
| van-maanen S1 | 18 | 1.51 / 21.08 | 4.81 / 15.73 | 9.68 / 106.35 |
| van-maanen S2 | 18 | −1.44 / 9.83 | 3.75 / 12.54 | 3.94 / 18.28 |
| wolf-359 S1 | 18 | −3.27 / 6.85 | 5.62 / 6.51 | 5.53 / 13.31 |
| wolf-359 S2 | 20 | 2.88 / 13.94 | 6.82 / 14.04 | 5.68 / 14.70 |

Adjudications (no planet or H < 7 asteroid within 0.5°, Horizons
`@-234`): **gj-1276 S2 S_pulse** — 2009-01-14, four consecutive
frames (+11.7, +17.6, +8.7, +3.0): the pixels show a ~10-px-wide
diagonal band sweeping across the stamp at 2–3 px/frame over 2.5 h,
an extended moving front (CME/streamer class), not a point source.
**teegarden S1 S_pulse** — 2010-08-14, the last two arc frames (+6.4,
+10.0, back to −2.6): the whole stamp lifts uniformly by ~0.8 DN/s
and drops back, an extended transient front. **gj-1276 S2 S_event** —
2012-10-22/24, a +0.39-unit plateau over 54 epochs (V_eq 12.1) with
no point source in the pixels at any epoch (gj-1276 is V 16): the
patch's own baseline sits at −0.27 units, a regional level-2
residual-background offset between the baseline and arc regions
that the flank controls do not share; the same patch is null in its
other 17 windows. Retained, non-promotable.

Dev stage (4 units, 12 trials, 2 exceedances vs 1.33): gj-908 S1
2022-10 extended-background plateau (same class); ross-154 S2
2010-10-06 single-frame +6.4-unit flare-class event on V1216 Sgr,
vetoed by the persistence rule. gj-908 (V 9.0, 13 units) repeats to
0.2 % across 18 transits — the in-situ bright-star control.

## 4. Completeness and controls

**Positive control (D5):** (2) Pallas through the in-beam arc as seen
from STEREO-A (Horizons `@-234`), two passages through the identical
moving-patch chain: 2008-02-22/26 (90 in-arc frames, V 9.7)
**+0.01 ± 0.04 mag** and 2012-11-16/20 (102 frames, HPLT −9.2°)
**−0.02 ± 0.07 mag** against the predicted V with the frozen colour
term — the ≤ 0.2 mag gate passed by an order of magnitude. (1) Ceres
2007-03 fell to the synoptic gate (early-mission 25 × 40 s sums).
Stamp response 0.997 ± 0.002 (5,382 stamps); level-2/level-1 star-flux
parity 1.005 ± 0.010; astrometry 0.44 px rms on ~700 stars per frame.

**Injections** into the real confirmatory null series against the
fixed thresholds (100 draws per flux, seed 20260904; recovery requires
S above max(T, observed null); `results/completeness_v1.json`,
`results/power_limits_v1.json`). Declared conversion: V_eq of a
B−V 0.65 source, line-equivalent through W_eff 100 nm at 680 nm,
±0.3 mag band systematic; S1 P_cone = F π (0.1 AU)²; S2 P_tx for a
10-m diffraction-limited transmitter at the registry distance.

| cell | m90 (V_eq) | power |
|---|---|---|
| S1 0.1 AU recurrence stacks (4 units) | 11.9–12.9 | **12–29 MW downlink through the cone** |
| S2 0.1 AU recurrence stacks (5 units) | 11.9–12.5 | **10-m uplink 2.9–61 MW** |
| S1 / S2 single-window (S_event) | 11.5–13.2 / 10.7–11.5 | 9–41 MW / 7–29 MW |
| pulse (≥ 1 frame, 40 min) S1 / S2 | 5.6–8.0 / 6.0–7.8 | 1–10 GW / 0.2–2.6 GW |

The recurrence stacks reach the declared expectation (per-window
~V 12–12.5, 18–20 recurrences); the pulse cell is systematics-limited
by the control patches' own single-frame excursions (bright static
stars, moving fronts) — T_pulse 7–106.

## 5. What this survey adds

The sunward 0.1 AU cells — LASCO's 3–26 GW (S1) and 0.7–1.1 GW (S2)
recurrence limits — are now constrained at **12–29 MW downlink and
3–61 MW uplink**, a 100–1000× deepening over 18–20 annual recurrences
per unit, with a validated moving-observer heliospheric-imager chain
(header-WCS astrometry, 0.05-mag ensemble photometry, asteroid
positive control at 0.01 mag). The grazing cones remain out of reach
of any heliospheric imager with a ≥ 4° inner edge: the plan's
"sub-MW downlink" hope for this substrate is closed by geometry, not
by data. The sunward pulse cell is open at 40-min resolution but
shallow (GW-class).

## 6. Open items and hand-offs

- The level-2 residual background carries regional ±0.3-unit offsets
  between the arc and baseline regions (the gj-908/gj-1276 plateau
  class); a v2 design would measure the baseline at the arc's own
  elongation band via the year's other transits, or fit the offset
  map from the calibrator dump.
- Pulse-cell thresholds are set by control patches containing bright
  static stars; a v2 could match control patches to the source's
  static content or mask VT ≤ 11 stars in control patches (the source
  patch is fixed).
- ross-128 S1 stays permanently blended (LASCO and HI-1); ross-154's
  crowded fields lose the west era to the ZP gate.
- PSP/WISPR (plan O5) remains the only imager with a smaller inner
  elongation at perihelion; the grazing cones need it or a
  coronagraph with a ≤ 2 R☉ occulter and stellar-depth photometry.
- Yearly refresh: extend `stereoa_v1` past the 2026-12 Horizons
  trajectory end; the archive currency (2026-08-31) adds one arc per
  unit per year.
