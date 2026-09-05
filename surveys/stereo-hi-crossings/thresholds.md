# STEREO-A HI-1 sunward-channels threshold freeze v1.0

**FROZEN 2026-09-04** (user approval of the hypotheses D1–D9 and this construction in one decision). Drafted after the coverage intersect
(`results/coverage_v1.json`) and before any frame is fetched at a
survey patch or any window-locked quantity is formed. Machine form:
`configs/threshold_freeze_v1.json` (sha256 in `results/threshold_freeze_v1_sha.txt`). Constructions
inherit the established chain (ZTF v1.0–1.2 rules, TESS shape tests,
ATLAS recurrence stack, LASCO v1.1 same-frame differential) with HI-1
substitutions per `hypotheses.md` (D1–D9).

## 1. Searchable-unit population (from coverage)

Unit = target × channel at the 0.1 AU rung. **14 units**: 7 targets ×
{S1, S2}. Coverage (arc ≥ 10 frames; per unit in the config):

| unit | covered / in-era events | median arc frames | median baseline frames |
|---|---|---|---|
| S1: gj-1276, van-maanen, teegarden, ross-154, gj-908 | 20 / 21 | 66, 66, 65, 59, 36 | 202–208 |
| S1: ross-128, wolf-359 | 19 / 21 | 65, 65 | 207–208 |
| S2: ross-128, wolf-359 | 20 / 21 | 66, 66 | 204–205 |
| S2: gj-1276, teegarden, van-maanen, ross-154 | 19 / 21–22 | 64–65, 58 | 201–208 |
| S2: gj-908 | 18 / 21 | 30 | 207 |

272 of 295 event rows covered; the 23 uncovered are the 2014-08-19 →
2015-11-16 safe-mode/conjunction gap (19), the 2007-03 gap (1), and
two events outside the archive era. Arc sides: 160 pre-t_ca (east-
pointing eras 2007–2014, 2023–), 113 post-t_ca (rolled west-pointing
era 2015-11 → 2023-08). S1 grazing rungs are `not_constrainable`
(ledger rows); S2 1 AU is out of scope (D7). No z-grid.

## 2. Statistics (per unit; S > max(T, 0) everywhere)

Series are per ICRS-fixed patch (source + 8 controls), per frame:
header WCS verified (≥ 50 Hipparcos stars, ≤ 0.6 px rms, translation
refined) → 31-px median high-pass → top-hat aperture (r 2.5 px,
annulus 5–9 px) → flux normalized by the frame's 2-D colour-corrected
ZP to ZP_ref 11.10 (V-equivalent of a B−V 0.65 source).

Per patch and event: arc epochs (in-beam, 4.05° ≤ |HPLN| ≤ 6.0°,
inside the frame footprint) and baseline epochs (same patch, 6° <
|HPLN| ≤ 12°, same side). Static-content removal
E_i = F_i(arc) − g·⟨F_i(base)⟩ with the measured differential flat g;
same-frame differential D_src = E_src − median_k E_k (controls:
E_k − median of the other seven); per-event z = mean(D_arc) /
[robust SD(D_arc)/√n_arc].

1. **S_stack** (primary): Σ z / √n over the unit's included events.
2. **S_event** (secondary): max z.
3. **S_pulse** (tertiary, every unit): max single-frame D_src in-arc
   divided by the robust epoch scatter; the frozen persistence rule
   (≥ 2 consecutive frames, or sibling-window recurrence) applies at
   adjudication, not to the statistic.

Trials: 3 per unit → **42 trials, expected control crossings 42/9 ≈
4.7** (dev 12 / 1.33, confirmatory 30 / 3.33).

## 3. Controls

- **Primary T: 8 parallel-track sky patches** at HPLT offsets ±1°,
  ±2°, ±3°, ±4° from the source track (same elongation, same frames,
  own star-fixed baselines; D4). T = max over the 8 of the identical
  statistic; controls need ≥ 5 survivors per epoch.
- **Secondary (reported, never T):** the source patch's own baseline
  transit cut into 1.9-d pseudo-arcs at |HPLN| 6°–12° (3 per event) —
  the temporal null of the same sky patch.

## 4. Unit and epoch gates

- Frame usable: N_IMAGES = 30, EXPTIME 1140–1260 s, NMISSING = 0,
  1024², celestial WCS present and verified (§2).
- Epoch usable: aperture ≥ 90 % finite pixels; no planet/Earth/Moon
  within 10 px of any of the event's patches; bright-body in-FOV veto
  (Earth, Moon, Venus, Jupiter inside the footprint → epoch dropped);
  **recon exclusion**: 2010-06-15 00:09 and 00:49 UT are dropped from
  the arcs of gj-1276 S1 `evt-6a6141430884` and wolf-359 S2
  `evt-d39c468c8a27` (hypotheses §11).
- Patch usable: no Tycho-2 VT ≤ 9 star within 3 px (S2 source patch
  exempt); a dropped source patch removes the event.
- Event included: ≥ 10 usable arc epochs and ≥ 30 usable baseline
  epochs on the source patch.
- S_stack searchable: ≥ 8 included events; else constraint-only.
- Exceedance adjudication: the frozen §6 veto ladder (known objects
  via Horizons `@-234`; persistence; TESS split-half ramp +
  background-anticorrelation; S2 stellar-flare rule against the
  star's own baseline-transit pulses; CME annotation; recurrence).

## 5. Split (D8; seed 20260904 for all randomized draws)

- **Dev** (4 units, 12 trials): gj-908 S1 + S2 (V 9.0 — the bright-
  star differential-flat test; shortest arcs, 24 h), ross-154 S1 + S2
  (V 10.4, flare star, |b_gal| ≈ 13°: `high_background` class).
- **Confirmatory** (10 units, 30 trials, analysed once after dev
  passes with any amendments frozen pre-pixel): van-maanen, wolf-359,
  teegarden, gj-1276, ross-128 × S1/S2.

## 6. Completeness and positive control (per D5/D6/§8)

Stamp-response injections into the real frames (empirical 1.7-px
PSF, C1-normalized) spanning flux × arc phase × HPLT × duty cycle
against the fixed confirmatory thresholds; the bright-asteroid arc
passage positive control (Horizons `@-234`, ≤ 0.2 mag after the
colour term); the SSC L2 product validates the chain (≤ 0.02 mag
ensemble agreement); ensemble era stability year by year.

## 7. Order of work

Freeze (sha) → fetch + measure the dev frames (4 units: ~76 events ×
~7.6 d ≈ 20k frames) → dev search → amendments if dev demands them
(frozen before any confirmatory pixel, versioned v1.x) → confirmatory
fetch + measure (10 units, ~190 events ≈ 50k frames) → search once →
injections/completeness → adjudication → report
(`report/stereo_hi_crossings.md`) → plan/history/learnings.

## Amendment v1.1 (2026-09-04, dev-driven; frozen before any confirmatory pixel)

Machine form `configs/threshold_freeze_v1_1.json` (sha256 in
`results/threshold_freeze_v1_1_sha.txt`). Found on the first five
complete dev unit-events measured on the level-1 product (≈ 2,500
frames; the level-1 dev measurements are discarded, kept as
`runs/.../measurements_dev_L1_discarded.jsonl`):

**Finding H1 — ridge/edge background bias on level 1.** Every source
patch sits on the F-corona ridge (the ecliptic; 145 DN/s/px at ε 4.9°)
while its HPLT-offset controls sit on the flanks (75 DN/s/px at +2°).
The 31-px median high-pass residual ramps from +1 to +16 DN/s across
the last 25 columns of the CCD (boundary handling on a 2 DN/s/px
gradient), and the local-annulus median is curvature-biased on the
ridge. The result is a ridge-dependent additive bias the same-frame
differential cannot cancel: first complete events z = −2 to −10, the
ross-154 S1 source arc flux falling to −5.8 units at |HPLN| 4.2°. On
the level-2 product (RAL's per-pixel 1-day running lowest-quartile
background removed — the static F-corona exactly, while sources drift
54 px/day) the same patch on the same frame reads +0.9 on a flat
residual (column medians 0.1–0.3 DN/s to the edge).

1. **Substrate (D1 amended):** the search runs on the SSC level-2
   `24h1A_br01` product — the same frames, the same header WCS, and the
   identical chain (astrometry verification, 2-D colour-corrected ZP,
   high-pass, top-hat aperture); level-2/level-1 star-flux parity
   1.005 ± 0.010 and chain ZPs 11.069/11.085 at recon. No level-1
   fallback: frames absent from the SSC tree (1.6 % in the parity
   sample) are recorded as uncovered.
2. **Footprint edge margin 24 px** (was 8) for every patch — the
   arc's inner limit becomes |HPLN| ≥ 4.37°.
3. Everything else — controls, statistics, gates, split, seed, recon
   exclusions — unchanged. Dev is re-run from scratch under v1.1.

## Amendment v1.2 (2026-09-04, dev reduce; frozen before any confirmatory pixel)

Machine form `configs/threshold_freeze_v1_2.json` (sha256 in
`results/threshold_freeze_v1_2_sha.txt`).

**Finding H2 — the frame-wide bright-body veto starves the units.** The
frozen in-FOV veto (Earth, Moon, Venus, Jupiter anywhere inside the
20° footprint → epoch dropped) removed every epoch of 11 of gj-908's
20 S1 events and 4 of 18 S2 events (2,244 Venus, 709 Jupiter and 1,682
Earth epochs in dev alone: the bodies sit somewhere in HI-1's field
for weeks to years). Measured on the dev series, epochs with a bright
body in the field but > 2° from the patches have the same same-frame-
differential scatter as clean epochs (gj-908 S1 arc MAD 0.44 vs 0.41
units; baseline 0.155 vs 0.153; ross-154 S1 arc 0.98 vs 0.92), while
≤ 2° roughly doubles it (ross-154 S1 arc with Venus: 1.85 vs 0.92).

4. **Bright-body veto = proximity**: an epoch is dropped when Earth,
   the Moon, Venus or Jupiter lies within **2.0°** of any of the
   event's nine patches. The 10-px all-body proximity mask is
   unchanged; the frame-wide veto is retired. Applied at reduce.

**Recorded, not amended.** The frozen ZP-scatter gate (MAD ≤ 0.12 mag)
rejects the west-era frames of ross-154's fields (both the star at
b = −13° and its antipode: blending-dominated scatter 0.12–0.17 mag,
vs 0.09–0.11 in the east era and 0.05 for gj-908's high-latitude
fields) — 4,129 of 18,298 dev frames; ross-154's units keep 9–10 of
19–20 events. The gate stands as frozen (the `high_background`
class); the confirmatory fields are all high-latitude.

**Dev result (v1.2; `results/dev_v1.md`).** 4 units, 12 trials, **2
exceedances vs 1.33 expected**, both adjudicated under the frozen
ladder:

- **gj-908 S1 S_event 7.80 vs T 7.46** (2022-10-13/14 arc, z = 7.8,
  20 epochs): a persistent +0.6-unit plateau (V_eq ≈ 11.7) with both
  split halves positive, no planet or bright asteroid within 0.5°
  (Horizons `@-234`, H < 7 list). The pixels show the patch's faint
  static source plus, at mid-arc, a broad +0.1 DN/s/px elevation over
  the whole stamp (an extended coronal background transient, not a
  point source), and the same-frame control medians jump by ±0.5
  between consecutive frames. The same sky patch is null in its other
  18 windows (z ≤ +1.7). **Adjudicated extended-background transient,
  non-recurrent; retained, non-promotable.**
- **ross-154 S2 S_pulse 9.41 vs T 5.13** (2010-10-06 17:43, one
  frame, D = +6.4 units with neighbours +1.5 / −1.3 — a 2.6×
  brightening of the star in a single 40-min sum). The flare rule's
  amplitude test is passed (the star's own off-window maximum is 1.6
  units), but the **persistence rule vetoes it**: single-frame
  excesses are never candidates. Class: single-frame stellar flare
  (ross-154 = V1216 Sgr) or particle hit. Non-promotable.
- gj-908's S2 differential is null across 18 transits (S_stack 6.6 vs
  T 9.4; the 13-unit star repeats to 0.2 %) — the in-situ bright-star
  control passes.

No further amendment; confirmatory proceeds once under v1.2.
