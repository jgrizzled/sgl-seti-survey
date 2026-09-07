# Solar Orbiter / SoloHI sunward-channels threshold freeze v1.0

**FROZEN 2026-09-06** (user approval of hypotheses D0–D9 and this
construction in one decision). Drafted 2026-09-06 after the coverage
intersect (`results/coverage_v1.json`) and before any frame is fetched
at a survey patch or any window-locked quantity is formed. Machine
form: `configs/threshold_freeze_v1.json` (sha256 in
`results/threshold_freeze_v1_sha.txt`). Constructions inherit the
established chain (ZTF v1.0–1.2 rules, TESS shape tests, ATLAS
recurrence stack, LASCO v1.1 same-frame differential, HI-1
parallel-track controls and star-fixed differential, WISPR v1.1–v1.5
gates) with SoloHI substitutions per `hypotheses.md` (D0–D9).

## 1. Searchable-unit population (from coverage)

Unit = target × channel at the 0.1 AU rung, inner tiles. **32 units**
(covered = ≥ 10 arc frames on ≥ 3 events; per unit in the config):

| units | covered events | median arc frames / event |
|---|---|---|
| S1 (19): 61-vir 9, ez-aqr 9, gj-1002 6, gj-1111 4, gj-1276 9, gj-581 6, gj-588 7, gj-667-c 7, gj-674 3, gj-682 7, gj-783 8, gj-876 9, gj-908 9, ross-128 4, ross-154 8, teegarden 4, van-maanen 8, wolf-1061 5, wolf-359 8 | 130 | 18–256 |
| S2 (13): 61-vir 6, ez-aqr 6, gj-1002 4, gj-1111 8, gj-1276 8, gj-251 8, gj-876 5, gj-908 7, ross-128 4, ross-154 4, teegarden 9, van-maanen 9, wolf-359 9 | 87 | 20–253 |

217 covered events, 29,243 arc frames (inner tiles, nominal-tile
index), same-tile baselines of median 74–422 frames per unit (capped
at the 100 nearest the arc at index). Not covered (ledger): gj-784 S1
(2 events), gj-251 S1, gj-518 S1, wolf-437 S1, gj-667-c S2 (1 each),
lacaille-8760 S1, gj-518/gj-588/wolf-437 S2 (0). S1 grazing rungs
`not_constrainable`; S2 1 AU out of scope (D7); seam rows `seam`.

## 2. Statistics (per unit; S > max(T, 0) everywhere)

Series per ICRS-fixed patch (source + 8 controls), per frame: header
celestial WCS refined by a Hipparcos match (translation + affine,
≥ 30 stars, ≤ 1.5 px rms) → 25-px median high-pass on L2 → top-hat
aperture (r 3 px, annulus 5–9 px) → flux normalised by the frame's 2-D
colour-corrected ZP to ZP_ref 10.05 (V-equivalent of a B−V 0.65
source; one unit = 1e-12 of the MSB scale) → **E' = F − median(F_base)**
(star-fixed differential: the same patch on the same-tile off-beam
frames preceding the arc; ≥ 10 baseline epochs on the source and ≥ 5
controls, else the event is `no_baseline` and drops from the stack).
The Gaia-template excess E = F − (a + b T) is formed alongside and
reported as `z_template` per event, never as S.

Per event: arc epochs (in-beam, inside the frame footprint by the
refined WCS with a 12-px margin inside DSTART/DSTOP). Same-frame
differential D_src = E'_src − Q(E'_k), Q the quadratic-in-orbit-
latitude interpolation of the controls at the source offset
(controls: own E' minus the interpolation of the other seven at their
own offset); per-event z = mean(D_arc) / [robust SD(D_arc)/√n_arc].

1. **S_stack** (primary): Σ z / √n over the unit's included events.
2. **S_event** (secondary): max z.
3. **S_pulse** (tertiary, every unit): max over consecutive frame
   pairs of min(x_i, x_i+1), x = (D_src − median)/robust scatter
   (two-frame persistence in the statistic, WISPR v1.5).

Trials: 3 per unit → **96 trials, expected control crossings 96/9 ≈
10.7** (dev 15 / 1.67, confirmatory 81 / 9.0).

## 3. Controls

- **Primary T: 8 parallel-track sky patches** at orbit-latitude
  offsets ±1.5°, ±3°, ±4.5°, ±6° from the source direction (same ram
  longitude at the same time, same frames, own baselines); one-sided
  ladder [1.5 … 12] toward the tile interior when the two-sided ladder
  would cross the seam or the outer latitude edge (recorded per
  unit-event in the index). T = max over the 8 patches of the
  identical statistic; ≥ 5 controls per epoch.
- **Positive controls** (D5): 61-vir S2 (V 4.7 in the source aperture;
  the star-fixed null at S/N ~ 100) and gj-908 S2 (V 9.0) in dev; a
  Horizons-selected moving body through an inner-tile arc band at
  dev (≤ 0.2 mag recovery); stamp-response injections on every 12th
  frame.

## 4. Unit and epoch gates

- Frame: `1ft`/`2ft`, 1024 × 960, NBIN 4, NSUMEXP ≥ 5, L2, celestial
  WCS present; astrometry refined; ZP uncertainty MAD/√n ≤ 0.10 mag
  with ≥ 30 calibrators (S/N ≥ 15, 3 ≤ V ≤ 9).
- Epoch: aperture ≥ 90 % finite; no planet within 12 px of any of the
  event's patches; Mercury/Venus ≥ 10°, Earth/Jupiter ≥ 5° from every
  patch; the recon frames (48 recon + 125 sampled headers: 96
  frame-event contacts, listed in the config) dropped for the events
  they touch.
- Patch: Gaia G ≤ 6 within 4 px → dropped (S2 source exempt); a dropped
  source patch removes the event; controls need ≥ 5 survivors.
- Event: ≥ 10 usable arc epochs, ≥ 10 source baseline epochs.
- Unit: ≥ 3 included events for the stack; otherwise
  `constraint_only`.

## 5. Split (D8; seed 20260906 for all randomized draws)

**Dev** (5 units, 15 trials): gj-908 S1, gj-908 S2, 61-vir S1, 61-vir
S2, ross-154 S1. **Confirmatory** (27 units, 81 trials): the rest, run
blind once after the dev machinery is fixed; the seven-system family
(11 remaining units) reported as the cross-substrate subset.

## 6. Completeness and positive control (per D5/D6/§8)

Stamp-response injections (`sglsurvey.inject.stamp_response` form: a
Gaussian PSF of the frozen FWHM added to the raw frame, re-high-passed,
re-measured; ratio recorded) on every 12th frame; the recurrence-stack
completeness from the measured per-event scatter of every included
unit (m90 V_eq per unit → power via the HI-1/WISPR anchors); the
moving-body recovery.

## 7. Order of work

index dev → run dev → reduce dev → dev report + any pre-confirmatory
amendment (frozen, sha-recorded) → index confirmatory → run → reduce
once → adjudicate exceedances under the §6 ladder → completeness →
report → plan/history/learnings.

## Amendment v1.1 (2026-09-06, dev machinery; frozen before any confirmatory frame)

Machine form `configs/threshold_freeze_v1_1.json` (sha256 in
`results/threshold_freeze_v1_1_sha.txt`). Found on the v1.0 dev reduce
(`results/dev_search_vtest.json`, retained as evidence):

1. **Saturation mask.** Along the sunward edge of tile 1 (x ≳ 935,
   ε ≲ 9°) in the ≥ 45-s exposure regime (r ≳ 0.5 AU) the F-corona
   reaches a flat plateau at 0.90–0.92 × the header `DSATVAL` and
   stars vanish (61-vir read 14–18 units in the arc against 101–106 in
   the baseline; a ±80-px search found no peak); the 16-s perihelion
   frames show no plateau. Rule: pixels ≥ 0.85 × DSATVAL are NaN
   before the high-pass; the frozen epoch gate (≥ 90 % finite
   aperture) then drops the saturated part of every arc, and
   saturated stellar cores leave the calibrator set (valid ≥ 0.99).
   The masked fraction is recorded per frame (`sat_masked`).
2. Everything else — controls, statistics, gates, split, seed, recon
   exclusions — unchanged. Dev is re-measured under v1.1 from scratch.

## Amendment v1.2 (2026-09-06, dev reduce; frozen before any confirmatory pixel)

Machine form `configs/threshold_freeze_v1_2.json` (sha256 in
`results/threshold_freeze_v1_2_sha.txt`). Found on the v1.1 dev reduce
(`results/dev_search_v1_1.json`, retained as evidence):

1. **Structure-noise epoch gate.** The same-frame differential series
   are flat (scatter 0.5 units, V_eq 9.1 per frame) over the outer
   three-quarters of every arc and swing by ±10–100 units in the last
   quarter, inside ~120 px (2.5°, ε ≲ 8°) of the sunward edge, where
   the unresolved F/K-corona structure dominates (local scatter 6–11
   units, p90 30; the recon's V 7 there). The per-epoch annulus noise
   `err` tracks it (0.4 → 1.0 → 2.6–6.0 units). Rule: a patch epoch
   with err > 1.5 units is dropped (source and controls alike; the
   ≥ 5-controls rule then applies). On dev this removes ~17 % of arc
   epochs. Control thresholds had been S_stack T 31–125 under v1.1.
2. **Colour calibrators V ≥ 4.5.** The ensemble colour fit excludes
   V < 4.5 stars (V 3–4 read 0.3 mag faint: bright-star nonlinearity);
   the per-frame ZP surface is unchanged (3σ clip).
3. Everything else unchanged. Reduce-stage only: the v1.1 measurements
   stand.

## Post-blind diagnostic (2026-09-06; not an amendment — the frozen statistics stand)

`SOLOHI_INTERP=median python series.py reduce confirmatory
v1_diag_median` → `results/confirmatory_search_v1_diag_median.json`:
the identical reduce with the robust median of the controls in place
of the frozen quadratic-in-latitude interpolation, run once after
the blind result to adjudicate the gj-1111 S2 / gj-251 S2 exceedances
(a single control excursion extrapolated into the source on a
one-sided ladder). 4 exceedances vs 6.0 expected (the two teegarden
S2 movers, gj-876 S1 ×2), 0 candidates, S_stack T median 6.8 vs 8.1.
Recorded as a v2 lesson; no completeness or limit is derived from it.
