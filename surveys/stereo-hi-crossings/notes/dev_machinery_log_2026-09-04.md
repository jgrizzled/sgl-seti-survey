# Dev machinery log — STEREO-A HI-1 sunward survey (2026-09-04)

Chronological record of the dev stage (freeze v1.0 → amendments v1.1,
v1.2), the positive control, and the machinery facts the report
draws on. Scripts: `scripts/series.py` (index / run / reduce),
`hi_lib.py` (chain), `hi_geometry.py`, `asteroid_control.py`,
`adjudicate.py`, `completeness.py`, `summarize.py`.

## 1. Fork-safety bug (before any science number)

The first dev run reported zlib "decompressing data" errors on ~30 %
of frames and empty record lists on most others. Cause: the Hipparcos
and observer tables were `np.load`ed lazily from compressed npz files
at import, and the 20 forked workers shared one file handle. Fix:
eager `np.array` copies at import. The first ~600 measured frames
were discarded and the run restarted; no survey number was formed
from the corrupted pass.

## 2. Finding H1 — level-1 ridge/edge bias → amendment v1.1

First five complete unit-events on level 1 gave z = −2 … −10 with
source arc fluxes to −5.8 units near the inner edge. Pixel inspection
(2007-06-19/20, ross-154 S1 antipode): the source patch sits on the
F-corona ridge at 145 DN/s/px while its HPLT-offset controls sit on
the flanks at 75; the 31-px median high-pass residual ramps +1 → +16
DN/s across the last 25 columns; the annulus median is curvature-
biased on the ridge. Random-position tests cannot measure this
(avoiding detected peaks selects troughs, biasing every estimator low
by ~0.5σ), so the diagnosis rests on the pixels and on the level-2
comparison: the same patch on the same frame reads −5.8 (L1) vs +0.9
(L2) with flat L2 column medians (0.1–0.3 DN/s) to the edge. v1.1:
search on level 2, 24-px margin. Level-1 measurements
(`runs/.../measurements_dev_L1_discarded.jsonl`) discarded; dev
re-run from scratch (18,298 frames, ~7,000 frames/h from the SSC
mirror, 20 workers, ~4 s per frame dominated by the median filter).

Residual patch-to-patch scatter on level 2: the per-patch arc-minus-
baseline offsets E scatter by ±0.2–0.5 units for source and controls
alike, driven by few-percent static-content leakage (a 2–3-unit
static star in a patch changes by 4–20 % between the arc and
baseline regions), with no ridge profile; the gj-908 star itself
(13 units) repeats to 0.2 %. This is the generic scatter the eight
controls calibrate; no further detrend change.

## 3. Finding H2 — frame-wide bright-body veto → amendment v1.2

See `thresholds.md` v1.2: the in-FOV veto emptied ~half of gj-908's
events; measured scatter with a body > 2° away equals clean-epoch
scatter; ≤ 2° doubles it. Replaced by a 2° proximity veto (Earth,
Moon, Venus, Jupiter). Applied at reduce — no re-measurement.

## 4. Recorded, not amended

- ZP-scatter gate (MAD ≤ 0.12) vs crowded fields: ross-154's fields
  (b = −13° star, +10° antipode) run 0.09–0.11 mag in the east era and
  0.12–0.17 in the west era (same fields, other side of the CCD; the
  blending-dominated scatter of 700+ calibrators, the median ZP itself
  being good to ~0.006 mag). 4,129 dev frames rejected, all ross-154
  west-era; the units keep 9–10 of 19–20 events. gj-908's fields:
  0.05 everywhere.
- Early-mission frames (≤ 2007-03) are 25 × 40 s, 600-s sums:
  `n_images_15/25` gate (140 dev frames); this also removed the Ceres
  2007-03 control passage.
- Arc epoch counts: gj-908 20–27 (its arc is 24 h), ross-154 42–60
  (of 58–60 listed); baselines ~200.

## 5. Positive control (D5) — `results/asteroid_control_v1.json`

(2) Pallas through the arc band as seen from STEREO-A (Horizons
`@-234`, 10-min ephemeris interpolated to each DATE-AVG), the
identical chain on the moving position: 2008-02-22/26 (119 frames,
90 in-arc, V 9.7, HPLT −5.7°) **+0.01 ± 0.04 mag**; 2012-11-16/20
(128 frames, 102 in-arc, HPLT −9.2°) **−0.02 ± 0.07 mag** against
V_inst = V − 0.595 (B−V − 0.65), B−V 0.66. Same-frame control patches
around the moving position read 0.0–0.8 units vs the asteroid's 3.5–4.3.
Stamp response through the chain 0.997 ± 0.003 (1,940 stamps).

## 6. Dev result (v1.2)

`results/dev_search_v1.json`, `results/dev_v1.md`: 12 trials, 2
exceedances vs 1.33 expected, adjudicated in `thresholds.md` (gj-908
S1 2022-10 extended background transient; ross-154 S2 2010-10-06
single-frame flare-class). Machinery validated; confirmatory started
the same evening on the frozen v1.2 construction.
