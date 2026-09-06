# PSP/WISPR sunward-channels hypothesis freeze v1.0 (Pipeline B)

**FROZEN 2026-09-05**: the user approved every §12 decision as recommended (D0 run, D1–D10). Drafted 2026-09-05 after the PSP observer crossing list
(`crossings/psp_v1`), the O5 geometry pass
(`notes/psp_wispr_geometry_pass_2026-09-05.md`), the WISPR reachability
recon (`notes/wispr_recon_2026-09-05.md`) and the coverage intersect
(`results/coverage_v1.json`), before any frame is fetched at a survey
position beyond the four recon frames declared in §11. Third
heliospheric substrate: the LASCO constructions
(`surveys/heliospheric-crossings/hypotheses.md`, D1–D9 + v1.1/v1.2) and
the HI-1 substitutions (`surveys/stereo-hi-crossings/hypotheses.md`) are
inherited; WISPR substitutions marked **[WISPR]**. Freeze decisions
D1–D10 in §12 for approval.

Input: **`crossings/psp_v1`** (`xng-c339df63a608`, PSP observer, Horizons
−96 at 10 min, 2018-08-15 → 2026-12-01, 7,448 events, 0 invalid) —
mandatory: PSP is at 0.046–0.73 AU and crosses each axis twice per
88–150 d orbit; the Earth-center list is a different population. This
survey executes plan §5.15 item O14: the two sunward combinations
**S1 = downlink post-lens** (outbound, target side) and **S2 = uplink
past the Sun** (inbound, anti-target side) on the **0.1 AU rung**, from
an observer for which that rung subtends 14°–65° (ε = atan(b_e / r)).

## 1. Archive and roles **[WISPR]**

- **NRL level 3 synoptic WISPR-I** (`psp_L3_wispr_*_1XY{1,2}.fits`, MSB,
  F-corona/stray-light "grand minimum" model subtracted, scaled by
  (DSUN/0.2 AU)^2.3 — the chain **divides that factor out exactly** from
  the header DSUN_OBS the pipeline used, restoring L2 units; recon
  parity L2/L3u 0.95–0.98) — the **search substrate**. Anonymous
  plain-HTTP tree `L3/orbitNN/YYYYMMDD/`, 3,954,240 bytes per frame,
  size-checked and FITS-parsed on every fetch. **Decision D1.**
- **NRL level 2** (same stems under `L2/YYYYMMDD/`) — parity check
  (gate: L2/L3u ensemble flux 0.90–1.05 per frame sample), fallback if
  an orbit's L3 is missing (L3 exists for orbits 01–27), inventory.
- **WISPR-O synoptic** (`2XYZ`) — *secondary*: the 0.1 AU arcs
  continue into WISPR-O at ε 50°–63° (b_e 13–21.5 R☉, the outer cone
  already covered by LASCO C3 / HI-1) and the source leaves the beam
  inside the O field, so an O off-beam baseline exists. Frozen as a
  **ledger + baseline substrate**, not a search substrate (its arcs
  add trials on the shallowest part of the cone; the `2302`/`2222`
  duplicate-stamp products need a dev look). Re-openable by amendment.
- L1, L2b, LW, high-cadence narrow-field regions (Z ∉ {1, 2}), PNG/
  movie products — unused.

## 2. Eras and gaps

| item | value | note |
|---|---|---|
| survey era | E1 2018-10-31 → E27 2026-03-17 (release currency) | event list to 2026-12-01 (Horizons SPK end); E28/E29 (2026-06, 2026-09) unreleased — refresh item |
| encounters | 27 released, r < 0.25 AU arcs 9.8–11.3 d; q 35.7 R☉ (E1–3), 27.9 (E4–5), 20.4 (E6–7), 16.0 (E8–9), 13.3 (E10–16), 11.4 (E17–21), 9.85 (E22–27) | geometry repeats within a q-group (line of apsides fixed) |
| WISPR-I synoptic cadence | 34 / 9 / 30 / 24 / 18 min (E1–E5), 15–16 (E6–E13; E7 30), 7.5 (E14–E16), **5.0 (E17–E23, E25, E27)**, 2.6 (E24), 7.5 (E26) | measured from the L2 inventory (`l2_inventory.json`, 106,752 files); two synoptic series `1211` (high cadence, ±2.6 d) and `1221` (half cadence, ±5 d) |
| gaps > 3 h inside arcs | 0–2 per encounter | per-event frame counts carry them |
| frames outside r < 0.25 AU | 19,702 WISPR-I | extended campaigns; used where an arc straddles the encounter start |
| ephemeris in headers | DSUN_OBS off by +0.65 % on E22, −0.10 % on E24 frames (predicted-SPK processing) | all geometry from the Horizons table; only the celestial WCS from headers |

## 3. Channels, rungs, in-era scope (PSP frame)

| channel / rung | events / targets (2018-08 → 2026-12) | WISPR-I | note |
|---|---|---|---|
| S1 ≤ 1.2 R☉ | 111 / 5 | **not_constrainable** | graze source ≤ 6.9° from Sun centre at any encounter; ledger rows |
| S1 ≤ 2.5 R☉ | 190 / 9 | **not_constrainable** | 13.3°–13.7° for gj-1111 / wolf-359 at E22+, minutes at the rung edge on the baffle strip; ledger rows |
| **S1 0.1 AU** | 1,513 / 58 | **searched**: 533 event-rungs with an arc on 21 targets; **17 units covered ≥ 10 events** | perihelion-side crossings (r ≈ 0.05 AU, t_ca + 2.6 → 17 h, b_e 2.6 → 14 R☉) and straddling windows (t_ca + 34 → 160 h, b_e 8–17 R☉) |
| **S2 0.1 AU** | 1,701 / 64 | **searched**: 561 / 21 targets; **24 units covered ≥ 10** | same two kinds |
| S2 1.0 AU | 1,964 / 86 | out of scope (D7) | 22 targets transit WISPR-I in 606 star-encounter pairs; the whole transit is in-beam |

**Apparent sources** (fixed ICRS directions, `scripts/psp_geometry.py`,
validated against the frame WCS): S2 = the target star (PM-propagated;
parallax from 0.05 AU ≪ the 152″ pixel), S1 = the star's antipode.
In-beam ⇔ b_e(t) ≤ 0.1 AU (filled cone, D2 inherited). **Ram-side
rule**: the post-t_ca half of every window is on the ram side (WISPR's
side) and the pre-t_ca half never is — every arc is post-t_ca.

**Unit population** (`results/coverage_v1.json`; covered = ≥ 10 arc
frames on ≥ 10 events): **41 units** — S1: 61-vir, gj-1002, gj-1087,
gj-1111, gj-251, gj-514, gj-518, gj-54, gj-581, gj-667-c, gj-908,
luyten-star, procyon-a, procyon-b, ross-128, teegarden, van-maanen,
wolf-1061, wolf-359, wolf-437 (20); S2: 61-vir, ez-aqr, fomalhaut,
gj-1002, gj-1087, gj-1276, gj-514, gj-518, gj-54, gj-581, gj-667-c,
gj-783, gj-876, gj-908, lacaille-8760, ross-128, ross-154, teegarden,
van-maanen, wolf-1061, wolf-437 (21). 986 covered events, 171,716
arc-frame measurements over ~36k distinct frames. The LASCO/HI-1
seven-system family is inside it (11 units).

## 4. Camera and wavelength **[WISPR]**

WISPR-I: 960 × 1024 (2 × 2 binned), 0.04231°/px = 152″, PSF FWHM
1.6 px (4′; encircled 0.95–1.00 at r 3 px), 5 summed exposures of
2–9 s (XPOSURE 9–47 s scaling with r), passband **490–740 nm**. The
532 nm doubled line is **in band** (the first heliospheric substrate
for which it is); 1064/1550 nm out of band. Hypothesis flux =
monochromatic line converted via the band effective width (W_eff ≈
250 nm) as always; broadband leakage declared for the out-of-band
lines. Colour system measured: +0.23 to +0.31 mag/(B−V) (recon, three
frames), re-measured on the ensemble at dev.

## 5. Duty cycle and temporal models

- **d = 1 chord — primary.** 90–380 frames per arc at 2.6–15 min
  (E1–E5: 34–9 min), **20–27 encounters per unit**: resolved light
  curves × recurrence. The per-encounter geometry repeats within a
  q-group, so stacks are like-with-like.
- **Pulse — secondary, every unit**: per-frame max in-arc, periods
  ≥ 5 min (E17+; 15 min E6–E13). The finest sunward pulse cell yet
  (LASCO 12 min, HI-1 40 min).
- **Recurrence stack** over the unit's covered encounters — primary.
- **Declared unconstrained:** sub-exposure pulses (< 2–9 s), the
  pre-t_ca (anti-ram) half of every window, the grazing cones, the
  E28/E29 windows, schedules avoiding perihelion ± 5 d.

## 6. Detection constructions **[WISPR]**

**Substrate chain (recon-validated):** per frame — header celestial
WCS (`RA/DEC-ZPN` 'A') **refined by a Hipparcos match** (≥ 30 stars
within 3 px, translation; rms ≤ 1.5 px; frames failing →
`no_astrometry`) → rescale L3 → L2 units → 15-px median high-pass →
top-hat aperture (r 3 px, annulus 5–9 px) at ICRS-fixed positions →
per-frame **2-D star ZP** (quadratic in x, y; Hipparcos V ≤ 8
calibrators, colour-corrected with the ensemble coefficient; gate
≥ 60 calibrators, MAD ≤ 0.25 mag) → flux in ZP-normalised units.
**Decision D3.**

**Patches.** The **source patch** (S2 star / S1 antipode) and **8
control patches at orbit-latitude offsets ±1.5°, ±3°, ±4.5°, ±6°**
from the source direction (parallel tracks: same ram longitude at the
same time, same frames, 35–140 px away; ±6° keeps the northern
controls inside the field's +15° edge for sources at orbit latitude
≤ +9°; for sources higher north the ladder flips to −1.5 … −12°,
recorded). **Decision D4.**

**Static-content removal — stellar template (LASCO route, generalised)
[WISPR].** No off-beam star-fixed baseline exists at the perihelion
crossings (recon §5), so the static content of every patch is
modelled: the predicted aperture flux of every Tycho-2 (VT ≤ 11) +
Hipparcos star within r + 1 px of the patch centre, through the
frame's ZP and the colour system (BT−VT → B−V), PSF-weighted by the
measured encircled fraction — F_template(frame, patch). The S2 target
star is part of the template (V/B−V table in the config, V ≤ 15).
E = F − F_template per epoch; the template's error (~0.1–0.15 mag
per star, plus the unmodelled VT > 11 background) is common to source
and controls and is calibrated by T. Where an off-beam WISPR-I
baseline exists (10 units), the **star-fixed differential** E' =
F_arc − ⟨F_base⟩ is reported as a secondary check, never as S.

**Same-frame differential (LASCO v1.1 rule 2):** D_src = E_src −
median_k E_k over the 8 controls on the same frame (controls: own E
minus the median of the other seven); per-event z = mean of in-arc D
over epochs / its empirical standard error (robust in-arc scatter;
≥ 10 usable epochs).

**Statistics per unit** (unit = target × channel; 3 trials/unit):
1. **S_stack** — Σ z_ev / √n over the unit's covered events (primary).
2. **S_event** — max z_ev (secondary).
3. **S_pulse** — max single-frame D_src in-arc, standardised by the
   in-arc robust epoch scatter (tertiary; every unit).

Decision rule S > max(T, 0), **T = max over the 8 control patches** of
the identical statistic; expected control crossings 1/9 per trial.

**Veto ladder:** (1) **known-object census** — Mercury, Venus, Earth,
Mars, Jupiter, Saturn, bright asteroids and comets (Horizons `@-96`;
WISPR's regular transients), annotated per exceedance epoch; (2)
**multi-frame persistence** — ≥ 2 consecutive frames on the fixed
patch, or sibling-encounter recurrence (single-frame excesses are
never candidates; cosmic rays survive the 5-exposure on-board sum
only partially); (3) **chord shape** — TESS split-half ramp +
background-anticorrelation tests; (4) **stellar-flare class (S2)** —
an in-arc pulse compared with the star's own template-subtracted
distribution over its other encounters; (5) **coronal transient
annotation** — CME/streamer-blob passages (WISPR encounter summaries,
HELCATS) at the epoch and position; (6) **recurrence** — the stack
demands it; any S_event exceedance is tested on the unit's other
encounters (same geometry within the q-group).

## 7. Quality masks **[WISPR]**

Frame level: synoptic full-field region (Z ∈ {1, 2}, `OBJECT
InnerFFV`), 960 × 1024, NSUMEXP ≥ 5 (8 in E1), XPOSURE and GAINMODE
recorded (the ZP absorbs them), astrometry refined (§6). Epoch level:
aperture ≥ 90 % finite pixels; no planet/comet within 12 px of any of
the event's patches; **bright-body in-FOV veto** (Venus, Mercury,
Earth, Jupiter inside the WISPR-I footprint → epoch dropped, the halo
is frame-wide — LASCO v1.2 lesson); the four recon epochs dropped for
the events they touch (§11). Patch level: Tycho-2 **VT ≤ 7** within
4 px → patch dropped for that event (bright-star wings exceed the
template's accuracy; ~1 % of patches), S2 source patch exempt
(template models the star); a dropped source patch removes the event;
controls need ≥ 5 survivors. Fractions measured at coverage/dev; masks
are not tuned after any window-locked quantity exists.

## 8. Flux scale, injections, completeness

Per-frame 2-D star ZP with the ensemble colour system; L2 parity on a
per-encounter sample (gate 0.90–1.05). Era stability measured from
the ensemble encounter by encounter (the exposure normalisation held
to 0.04 mag across a 3× exposure change at recon). Completeness =
**stamp-response injections into the real frames** (the established
`sglsurvey.inject.stamp_response` route, empirical PSF), spanning
flux × arc phase × orbit latitude × duty cycle, normalised per the C1
rule. **Positive controls (decision D5):** (a) **61-vir S2** — V 4.7
in the source aperture, the template null at S/N ~ 100 across 22
encounters (dev); (b) **gj-908 S2** (V 9.0, dev) as in LASCO/HI-1;
(c) a **moving body through the arc region** chosen at dev from
Horizons `@-96` — the brightest of Vesta/Ceres/Pallas or a catalogued
comet with V ≤ 7 in a WISPR-I arc band, through the identical patch
chain along the predicted track, recovery ≤ 0.2 mag; Mercury/Venus
saturate and are veto bodies, not controls.

## 9. Expected sensitivity, declared before search

Recon single-frame 5σ depth (solar colour, rescaled L3): **V 5.4–6.3
at r 0.058 AU (15-s exposures), V 7.1–7.9 at 0.157 AU (47 s)** —
2–3 mag shallower than HI-1 per frame, ~1 mag shallower than LASCO
C3, structure-limited (K-corona streamers). Declared expectations
(systematics-bounded; gains measured): per-arc stacks of 90–380 frames
~V 8–10.5; recurrence stacks over 20–27 encounters ~V 10–12
optimistic bound. Power scale via the HI-1 anchors (recurrence m90
V_eq 11.9–12.9 ↔ downlink 12–29 MW, uplink 3–61 MW): **~1–3 mag
shallower → downlink ~30–300 MW, uplink ~10–300 MW through the 0.1 AU
cone**, LASCO-class (40–530 MW). Honest framing: WISPR does **not**
deepen the 0.1 AU sunward cells; proximity to the beam buys
visibility, not flux (the beam footprint, not the observer's
distance, sets the received power). What it adds: (i) the **inner
cone** b_e 2.6–14 R☉ from a new vantage on the star-side perihelion
crossings, (ii) the **5-min sunward pulse cell** (12 / 40 min before),
(iii) an independent third instrument on the same cells and epochs,
(iv) 41 units instead of 14 — the first sunward constraints on 30
targets, (v) the first in-band 532 nm sunward substrate.

## 10. Volume policy (decision D6)

Full cadence: every synoptic WISPR-I L3 frame inside the covered arcs
(≈ 36k distinct frames ≈ 145 GB transient, per-batch purge) plus the
baseline transits where they exist; one measurement pass per frame
serves every patch mapped to it (41 units × 9 patches). ~1.5 s per
fetch → ~15 h serial, ~4 h at 4 streams; measurement ~1 s per frame.
The pulse cell is open on every unit.

## 11. Pre-freeze data-contact declaration

Recon fetched four full frames (WISPR-I L2+L3 2024-12-21 23:30,
2024-12-24 00:00, 2024-12-27 00:30; WISPR-O L2+L3 2024-12-24 00:02) and
five headers (E1, E3, E10 ×2, E14, E24). Checked against the coverage
rows: the frames lie inside the E22 arcs of gj-1002 S1, gj-54 S1,
gj-908 S1, ross-128 S2, wolf-437 S2 (12-21), teegarden S1, gj-581 S2,
wolf-1061 S2 (12-24 I), gj-581 S1, wolf-1061 S1 (12-27; + the
teegarden S2 baseline), and 12 WISPR-O arcs (12-24 O). No quantity
was formed at any of those positions (full-field Hipparcos match,
column-strip noise, PSF from V < 5 stars). **Declared remedy:** the
four frame epochs are excluded from the series of every event they
touch (1 of 86–553 frames each; listed in the config) and the recon
frames are never re-used as search frames; every survey frame is
fetched fresh under snapshot discipline.

## 12. Freeze decisions (recommendations marked; freeze on approval)

- **D0 — go/no-go.** Run the survey as a third substrate at
  LASCO-class depth with the five additions of §9, or table it as
  "geometry established, depth insufficient to improve on HI-1".
  **Recommended: run** — the inner-cone, 5-min pulse and 30-new-target
  cells are new, the cost is one chain on public data, and the
  programme's rule is that reachable substrates are searched and
  reported at their honest depth.
- **D1** — substrate: NRL L3 synoptic WISPR-I rescaled to L2 units;
  L2 as parity/fallback; WISPR-O ledger + baseline only.
  (Alternative: search L2 with our own background model — rejected:
  the NRL three-orbit grand-minimum model is the instrument team's
  and is reproducible from L2b.) **Recommended.**
- **D2** — S1 beam model: filled-cone windows (inherited).
  **Recommended.**
- **D3** — chain: Hipparcos-refined header WCS + rescale + 15-px
  median high-pass + 2-D colour-corrected star ZP + stellar-template
  static removal + same-frame differential; star-fixed differential
  secondary where a baseline exists. **Recommended.**
- **D4** — controls: 8 parallel-track patches at orbit-latitude
  offsets ±1.5 … ±6° (one-sided ladder for far-north sources), same
  frames. (Alternative: ram-longitude-shifted patches — rejected: not
  simultaneous in elongation, different coronal background.)
  **Recommended.**
- **D5** — positive controls: 61-vir S2 and gj-908 S2 template nulls
  (dev), plus a Horizons-selected moving body through an arc band
  (≤ 0.2 mag); stamp-response injections for completeness.
  **Recommended.**
- **D6** — volume: full cadence over every covered arc (~145 GB
  transient); pulse cell on every unit. **Recommended.**
- **D7** — S2 1.0 AU out of scope (ledger kept: 606 star-encounter
  transits); S1 grazing rungs `not_constrainable` (ledger rows).
  **Recommended.**
- **D8** — dev/confirmatory split (seed at threshold freeze): **dev =
  gj-908 S1+S2, ross-154 S2, 61-vir S1+S2** (5 units / 15 trials — the
  two bright-star template tests, the flare star, the low-latitude
  field); **confirmatory = the other 36 covered units** (108 trials,
  12.0 expected control crossings), with the LASCO/HI-1 family (8
  units) reported as the cross-substrate subset. (Alternative:
  confirmatory = family only, 8 units — rejected: it discards the
  first sunward constraints on 28 targets for no cost.) **Recommended.**
- **D9** — era policy of §2: E1–E27 (release currency 2026-03-17);
  per-frame footprint from the header WCS; geometry from the Horizons
  table only; `1211` + `1221` both synoptic. **Recommended.**
- **D10** — 532 nm in-band: the line hypotheses are evaluated in-band
  for 532 nm and as broadband leakage for 1064/1550 nm (declared per
  §4). **Recommended.**
