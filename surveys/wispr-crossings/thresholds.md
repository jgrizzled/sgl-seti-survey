# PSP/WISPR sunward-channels threshold freeze v1.0

**FROZEN 2026-09-05** (user approval of hypotheses D0–D10 and this construction in one decision). Drafted 2026-09-05 after the coverage intersect
(`results/coverage_v1.json`) and before any frame is fetched at a
survey patch or any window-locked quantity is formed. Machine form:
`configs/threshold_freeze_v1.json` (sha256 recorded in
`results/threshold_freeze_v1_sha.txt`). Constructions
inherit the established chain (ZTF v1.0–1.2 rules, TESS shape tests,
ATLAS recurrence stack, LASCO v1.1 same-frame differential + stellar
template, HI-1 parallel-track controls) with WISPR substitutions per
`hypotheses.md` (D0–D10).

## 1. Searchable-unit population (from coverage)

Unit = target × channel at the 0.1 AU rung, WISPR-I. **41 units**
(covered = ≥ 10 arc frames on ≥ 10 events; per unit in the config):

| units | covered events | median arc frames / event |
|---|---|---|
| S1 (20): 61-vir 25, gj-1002 20, gj-1087 27, gj-1111 26, gj-251 27, gj-514 26, gj-518 25, gj-54 20, gj-581 24, gj-667-c 11, gj-908 20, luyten-star 27, procyon-a 27, procyon-b 27, ross-128 26, teegarden 25, van-maanen 22, wolf-1061 22, wolf-359 26, wolf-437 26 | 499 | 86–384 |
| S2 (21): 61-vir 22, ez-aqr 26, fomalhaut 26, gj-1002 26, gj-1087 10, gj-1276 26, gj-514 22, gj-518 22, gj-54 26, gj-581 26, gj-667-c 27, gj-783 26, gj-876 26, gj-908 26, lacaille-8760 26, ross-128 20, ross-154 27, teegarden 24, van-maanen 26, wolf-1061 27, wolf-437 20 | 487 | 21–293 |

986 covered events, 171,716 arc-frame measurements over ~36k
distinct WISPR-I synoptic frames (E1–E27). Off-beam WISPR-I baselines
exist for 10 units (median 30–210 frames); the others rely on the
stellar template alone. S1 grazing rungs `not_constrainable`; S2 1 AU
out of scope (D7). No z-grid.

## 2. Statistics (per unit; S > max(T, 0) everywhere)

Series per ICRS-fixed patch (source + 8 controls), per frame: header
celestial WCS refined by a Hipparcos match (≥ 30 stars, ≤ 1.5 px rms)
→ L3 rescaled to L2 units → 15-px median high-pass → top-hat aperture
(r 3 px, annulus 5–9 px) → flux normalised by the frame's 2-D
colour-corrected ZP to ZP_ref −21.35 (V-equivalent of a B−V 0.65
source) → **E = F − F_template** (Tycho-2 VT ≤ 11 + Hipparcos stars
within r + 1 px, PSF-weighted, through the same ZP/colour system; S2
target star included).

Per event: arc epochs (in-beam, inside the frame footprint by the
header WCS with a 12-px margin, inside the encounter or a released
extended-campaign frame). Same-frame differential D_src = E_src −
median_k E_k (controls: E_k − median of the other seven); per-event
z = mean(D_arc) / [robust SD(D_arc)/√n_arc].

1. **S_stack** (primary): Σ z / √n over the unit's included events.
2. **S_event** (secondary): max z.
3. **S_pulse** (tertiary, every unit): max single-frame D_src in-arc
   divided by the robust epoch scatter; the frozen persistence rule
   (≥ 2 consecutive frames, or sibling-encounter recurrence) applies
   at adjudication, not to the statistic.

Trials: 3 per unit → **123 trials, expected control crossings 123/9 ≈
13.7** (dev 15 / 1.67, confirmatory 108 / 12.0).

## 3. Controls

- **Primary T: 8 parallel-track sky patches** at orbit-latitude
  offsets ±1.5°, ±3°, ±4.5°, ±6° from the source direction (same ram
  longitude at the same time, same frames, own templates; D4);
  sources north of +9° orbit latitude use the one-sided ladder −1.5°
  … −12° (the field's north edge is at +15°). T = max over the 8 of
  the identical statistic; ≥ 5 survivors per epoch.
- **Secondary (reported, never T):** the star-fixed differential
  (arc minus the patch's own off-beam WISPR-I baseline) for the 10
  units with ≥ 30 baseline frames; the WISPR-O continuation of the
  arc as a ledger.

## 4. Unit and epoch gates

- Frame usable: synoptic full-field region (Z ∈ {1, 2}, `OBJECT
  InnerFFV`, 960 × 1024), NSUMEXP ≥ 5, astrometry refined, ZP gate
  (≥ 60 Hipparcos V ≤ 8 calibrators, MAD ≤ 0.25 mag); XPOSURE and
  GAINMODE recorded, not gated.
- Epoch usable: aperture ≥ 90 % finite pixels; no planet, bright
  asteroid (H < 7) or catalogued comet within 12 px of any of the
  event's patches (Horizons `@-96`); bright-body in-FOV veto (Venus,
  Mercury, Earth, Jupiter inside the WISPR-I footprint → epoch dropped;
  to be re-examined at dev the way HI-1 v1.2 did if it starves units);
  **recon exclusion**: the four recon frame epochs (2024-12-21 23:30,
  2024-12-24 00:00 and 00:02, 2024-12-27 00:30) are dropped from the 17
  E22 events they touch (`recon_exclusions` in the config).
- Patch usable: no Tycho-2 VT ≤ 7 star within 4 px (S2 source patch
  exempt); a dropped source patch removes the event.
- Event included: ≥ 10 usable arc epochs on the source patch.
- S_stack searchable: ≥ 8 included events; else constraint-only.
- Exceedance adjudication: the frozen §6 veto ladder (known objects
  via Horizons `@-96` and the WISPR encounter summaries; persistence;
  TESS split-half ramp + background-anticorrelation; S2 stellar-flare
  rule against the star's own other-encounter distribution; coronal
  transient annotation; recurrence within the q-group).

## 5. Split (D8; seed 20260905 for all randomized draws)

- **Dev** (5 units, 15 trials): gj-908 S1 + S2 (V 9.0 — the LASCO/HI-1
  bright-star template test), ross-154 S2 (V 10.4, flare star,
  |b_gal| ≈ 13°), **61-vir S1 + S2** (V 4.7 — the template null at
  S/N ~ 100; its S1 unit is the matching antipode field).
- **Confirmatory** (36 units, 108 trials, analysed once after dev
  passes with any amendments frozen pre-pixel): the remaining covered
  units of §1; the LASCO/HI-1 cross-substrate subset (van-maanen S1 +
  S2, wolf-359 S1, teegarden S1 + S2, gj-1276 S2, ross-128 S1 + S2)
  reported as its own block.

## 6. Completeness and positive control (per D5/D6/§8)

Stamp-response injections into the real frames (empirical 1.6-px
PSF, C1-normalised) spanning flux × arc phase × orbit latitude ×
duty cycle against the fixed confirmatory thresholds; the moving-body
arc passage (Horizons `@-96`, ≤ 0.2 mag after the colour term); the
L2 product validates the rescaling on a per-encounter sample
(0.90–1.05); ensemble era stability encounter by encounter.

## 7. Order of work

Freeze (sha) → fetch + measure the dev frames (5 units: ~120 events,
~20k frames ≈ 80 GB transient) → dev search → amendments if dev demands
them (frozen before any confirmatory pixel, versioned v1.x) →
confirmatory fetch + measure (36 units; the remaining ~16k frames) →
search once → injections/completeness → adjudication → report
(`report/wispr_crossings.md`) → plan/history/learnings.

## Amendment v1.1 (2026-09-05, pre-dev machinery; frozen before any survey frame was measured)

Machine form `configs/threshold_freeze_v1_1.json` (sha256 in
`results/threshold_freeze_v1_1_sha.txt`). Found while validating the
chain on the four recon frames (no survey position measured):

1. **Affine WCS refinement.** Translation-only refinement left 1.8 px
   rms on the 2024-12-27 frame (a 0.3° residual rotation of the header
   solution: 4 px at the field edge); a 2 × 3 affine fitted to the
   ≥ 30 matched V ≤ 7.5 stars after the translation gives 0.62–0.75 px
   on all three WISPR-I frames. The rms gate (≤ 1.5 px) is unchanged.
2. **Calibrator S/N ≥ 15, ≥ 30 calibrators (linear ZP surface below
   60, quadratic above), ZP-scatter gate 0.35 mag** (were S/N 8, ≥ 60
   quadratic, 0.25). At perihelion exposures (15 s) only ~40 Hipparcos
   V ≤ 8 stars reach S/N 15, and admitting fainter ones biases the ZP
   by noise-boosted selection (−21.40 at S/N ≥ 15, −21.28 at ≥ 10,
   −21.05 at ≥ 8 on the 2024-12-25 04:50 frame, r 0.068 AU); the
   frozen gates would have dropped most perihelion frames. The MAD
   gate bounds the per-star scatter; the frame ZP itself is determined
   to MAD/√n ≈ 0.02–0.03 mag.
3. Everything else — controls, statistics, gates, split, seed, recon
   exclusions — unchanged.

## Amendments v1.2–v1.4 (2026-09-05, dev-driven on the partial dev reduce; frozen before any confirmatory pixel)

Machine forms `configs/threshold_freeze_v1_2.json` … `v1_4.json` (sha256
in `results/threshold_freeze_v1_{2,3,4}_sha.txt`); each item is
recorded with its evidence in the config `amendments` block and in
`notes/dev_machinery_log_2026-09-05.md`.

- **v1.2 — template completion and flux unit.** The Tycho-2 VT ≤ 11
  extract lacks the brightest and the high-proper-motion stars (61 Vir,
  V 4.7; gj-908, V 9.0, 1.4″/yr), so the dev S2 source patches carried
  the whole star as excess (61-vir S2 z ≈ +120 per event). Template =
  Tycho-2 ∪ Hipparcos (de-duplicated within 10″) + the S2 target star
  from a Simbad V/B−V table (`s2_stellar_template`, 21 targets).
  Flux unit set to 1e-12 of the MSB scale (ZP_REF 8.65) — the stored
  series had been rounded to zero at the 1e-13 scale.
- **v1.3 — proximity veto and the bright-star S2 class.** Control
  same-frame-differential MAD by class: clean 1.35; Mercury within
  10° 3.2–4.9; Venus 5–10° 1.42; Jupiter 2–10° 1.35–1.51; the
  frame-wide in-FOV veto had flagged 36 % of arc epochs and its
  large-separation excess is an encounter confound (E12–E14). Veto =
  proximity to any of the event's nine patches: Mercury/Venus 10°,
  Earth/Jupiter 5°. S2 units whose target has 3 ≤ V < 8 (61-vir,
  gj-783, lacaille-8760) use a self-calibrated template (the unit's,
  or the q-group's when ≥ 3 events, median ⟨F⟩/T — a signal at the
  same level in every encounter is absorbed, declared); V < 3
  (fomalhaut S2) is `bright_star_unsearchable`. gj-908 S2 (V 9.0)
  needs no self-calibration: ratios 0.65–1.04, z within ±2.6.
- **v1.4 — template response.** On the control-patch ensemble the
  measured arc flux is a + b·T with b ≈ 0.7–0.9 and a ≈ −0.07 to −0.13
  units (V 6–8 patches read 0.885 of the template; V 8–11 patches
  0.2–0.5 units low; blank patches slightly negative — the high-pass/
  annulus zero point in star fields). E = F − (a + b·T) with a, b
  fitted at reduce from the controls (source patches excluded) and
  applied to source and controls alike — the HI-1 differential-flat
  pattern. After it: gj-908 S2 z median +0.25, ross-154 S2 +0.8.
- Everything else — controls, statistics, gates, split, seed, recon
  exclusions — unchanged. Dev is reduced under v1.4.
