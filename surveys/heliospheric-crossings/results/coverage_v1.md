# LASCO sunward-channels coverage intersect v1 (2026-08-25)

Stage 2 of the frozen chain (`hypotheses.md` v1.0, FROZEN 2026-08-25).
Inputs: `crossings/soho_v1` sunward events × the frozen rung ladder;
occulter annuli measured from level-1 frames
(`occulter_radii_v1.json`, adopted: C2 usable [2.2, 6.3] R☉ —
diffraction ring at 1.9–2.0 R☉; C3 usable [4.4, 29] R☉ — occulted
core is a constant fill plateau, not zeros, in the post-1997 level-1
convention). Frame counts from 4,320 unique (camera, day) archive
directory listings, snapshotted under
`runs/heliospheric-crossings/coverage/listings/` (SDAC 3,882 / NRL
247 / empty 191). Day-granularity counts; per-frame times are a
search-stage product. Full rows: `coverage_v1.json`.

## Results (703 event-rung rows)

| channel · rung | events | covered | covered + visible | median frames/window | visible fraction |
| --- | --- | --- | --- | --- | --- |
| S1 1.2 R☉ (C2) | 117 | 111 | **0** | ~124 | 0.00 |
| S1 2.5 R☉ (C2) | 153 | 148 | **148** | ~230–244 | 0.12–0.25 |
| S1 0.1 AU (C3) | 216 | 208 | **208** | ~650–1,300 | 0.80–1.00 |
| S2 0.1 AU (C3) | 217 | 212 | **212** | ~920–1,270 | 0.80–1.00 |

- **S1 1.2 R☉ is formally `not_constrainable`** — every window's
  apparent source stays inside the C2 usable inner radius (2.2 R☉),
  exactly as freeze D7 anticipated. The 111 covered windows enter the
  ledger as nominal-covered / occulter-unusable; the frames are the
  same days as the 2.5 R☉ rung (no extra cost).
- **S1 2.5 R☉ searches its [2.2, 2.5] R☉ wings**: visible fraction
  0.12–0.25 of each ~1.3 d window ≈ 28–60 C2 frames per window ×
  ~29 windows per target — the recurrence stack is well fed even in
  the wing regime. All five grazing-family targets (van-maanen,
  wolf-359, teegarden, gj-1276, ross-128) at 29–31 covered+visible
  windows each.
- **The 0.1 AU rungs are the workhorse**: ~96–98 % of windows covered
  and visible, ~1,000+ full-cadence C3 frames per window (→ ~90 at
  the frozen D6 3-hourly subsample), visible fraction 0.80 for the
  grazing-depth targets (the 4.4 R☉ occulter core removes the
  deepest-elongation 20 %) and 1.00 for gj-908 (b = 12 R☉ never
  enters the occulter).

## Losses (23 uncovered event-rung rows, all explained)

1. **The 1998 SOHO attitude-loss gap** (1998-07 → 1999-03 incl. the
   gyro recovery): one full annual generation across the family —
   van-maanen 1998-10 (b = 0.35 R☉), gj-1276 1998-09, wolf-359
   1998-09 + 1999-03, ross-128 1998-09, ross-154 1998-07/1999-01,
   gj-908 1998-09. Encoded era gap, as frozen (§2/D9).
2. **The 2026-07 → 2026-10 tail** past the NRL mirror currency
   (~2026-06): 2 gj-1276 + 1 each gj-908/ross-154/van-maanen/
   ross-128/wolf-359 rows — recoverable at the yearly refresh.
3. Two isolated day gaps: teegarden 2009-11-08 (1.2 R☉ row only —
   moot) and 2001-11-08 (2.5 R☉ row).

## Handoff to the threshold freeze

Searchable-unit population (target × channel × rung × camera) from
the covered+visible counts above: 5 S1-grazing (2.5 R☉ C2) units,
7 S1-0.1 AU + 7 S2-0.1 AU C3 units; S_pulse applies only to the
full-cadence C2 units (D6). Dev/confirmatory split per D8 (dev =
ross-154, gj-908, ross-128; confirmatory = van-maanen, wolf-359,
teegarden, gj-1276), seed to be declared at threshold freeze.
Expected control-crossing budget tallied there (1/9 per searched
trial, 3 statistics per unit — 2 on subsampled units).
