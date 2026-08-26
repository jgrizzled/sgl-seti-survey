# GALEX crossings dev stage v1 — CLOSED (2026-08-26)

Machinery development per threshold freeze v1.0 §Dev, **pseudo-units
only** (D8): every visit touched is off-window by the
`assert_off_window` guard against the union of each field's in-era
0.1 AU flat-chord windows — the five confirmatory in-window visits
were unreachable by construction. Scripts `scripts/dev_stage.py`
(+ `galexlib.py`, `flare_template.py`); machine record `dev_v1.json`;
all queries snapshotted under `runs/galex-crossings/`. All five
deliverables closed; two frozen-rule assessments resolved; one
amendment (v1.2).

## (i) Flag-64 astrometry — GATE PASS

16 source-blocks (deduped to ~4 physical MCAT stars × 2 petal
blocks, NUV 18.8–19.2, 480–2,540 photons each) measured during
flag-64 aspect seconds only: centroid offsets 0.08–1.02″ (median
0.37″), per-photon RMS 5.1–9.0″ ≈ the NUV PSF. The amendment-v1.1
premise is formally verified: flag-64 aspect seconds carry
sub-arcsec astrometry.

## (ii) Live-time/aperture calibration — GATE PASS

15 isolated calibrators (MCAT deduped at 3″ — the catalog carries
near-duplicate per-visit rows; the recon's ×2.6 "discrepancy" was
this plus crude live-time) over 3 off-window blocks:
**ZP_eff = 19.592 ± 0.088** (robust scatter ≤ 0.2 ✓) for the 8″
aperture chain vs nominal 20.08 — a 0.49 mag aperture+dead-time
term, now measured, closing the recon item. Blank-sky background
**0.209 cts/s per 8″ aperture** (24 positions).

## (iii) Control census — assessed

- **B pseudo-position rule**: valid with zero +5° rotations on all
  3 blocks tested (locus radius ~23′; no MCAT collisions).
- **A segment supply** (frozen < 8 → constraint-only rule):
  wolf-359 NUV **51** ✓, ross-128 NUV **35** ✓, ross-128 FUV **25**
  ✓; **gj-1276 NUV 4 / FUV 3** (its only long visit, 2010-10-23, is
  rim-only > 33′) and **wolf-359 FUV 0** (long visits NUV-only) →
  **gj-1276 A (both bands) and wolf-359 A FUV are constraint-only**.
  Confirmatory family: **5 band-units / 15 trials / expected
  control crossings 1.67**.

## (iv) Flare census — 3 flares at wolf-359; templates extracted

gj-1276 star (NUV ≈ 0.23 cts/s ≈ sky) and ross-128 (0.41 cts/s)
quiet in every off-window visit. wolf-359 (quiescent ≈ 2.1–2.8
cts/s): burst structure in all three long visits — S_burst 105.7 /
62.0 / 14.4. FRED templates (`iv_flare_template`): 2009-03-25 peak
59.4 cts/s (contrast 21.6×, rise ≤ 10 s, decay e-fold 30 s);
2010-02-23a contrast 12× (≤ 20 s / 40 s); 2010-02-23b a slow-rise
3× event — **not** FRED-like, so the veto's second prong (two-band
correlation) is load-bearing for slow flares, as frozen.

## (v) End-to-end pseudo-units + injections

- **B pseudo-unit** (the contiguous 1,381 s flag-64 segment,
  2009-03-08, z-550 locus apertures, 280 photons): null-clean —
  S_rate 0.203 vs T 0.274; S_burst 3.96 vs 4.87; S_period
  (jittered, v1.2) 31.2 vs 37.8. Controls 287–378 photons,
  rate-matched.
- **A pseudo-unit** (wolf-359 97 s segment vs 8 sibling segments):
  the drawn segment contained the 2009-03-25 flare — **all three
  statistics exceeded their control thresholds through the full
  chain** (S_burst 19.5 vs 5.5; S_rate 7.8 vs 2.9; S_period 7,489
  vs 2,733 unjittered). The conditional positive control is
  satisfied: a real astrophysical transient is detected by the
  frozen machinery, and the flare veto then classifies it (FRED ✓).
- **Injections** (real background + synthetic arrivals, identical
  chain): persistent 90 % recovery between 0.05 and 0.1 cts/s →
  **m90 ≈ NUV 22.3 (1.4 ks class)** via ZP_eff; 0.5 s single-pulse
  fluence threshold ~5 photons (3 ph: 22 %; 5 ph: 100 %); trains at
  P = 0.1 s recovered at every amplitude tried, P = 1 s from
  ~0.15 ph/cycle; P = 10 s dead under quantization and **restored
  by v1.2 jitter** (6/6 at 0.4 ph/cycle;
  `jitter_validation_v0.json`). Orbital-drift injections showed no
  loss at these (coarse) amplitudes; the confirmatory completeness
  grids sample finer.

## Amendment v1.2 (user approved)

S_period on tick-jittered times (U(0, 5 ms), seed 20260826 +
per-series offset): the 5 ms tick rate (200 Hz) makes every grid
frequency with a commensurate harmonic accumulate coherent
quantization power (null H 2,700–3,800 vs ~32 expected). Validated
before adoption; hypotheses.md §Amendment v1.2.

## Dev exit — cleared for blind confirmatory

All gates pass; the frozen chain for the blind run is
v1.0 + v1.1 + v1.2 with the threshold-freeze config and the
dev-assessed 15-trial family. Confirmatory inputs fixed: ZP_eff
19.592 (depth conversion), FRED templates (veto), control
constructions verified, injection grids per unit. The D1a
forced_dev lane (S_rate, gj-1276 B 2010) is computed alongside but
outside the family.
