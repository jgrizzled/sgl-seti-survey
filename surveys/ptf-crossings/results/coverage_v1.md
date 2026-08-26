# PTF crossings coverage intersection v1 — results

Run 2026-08-26 against hypotheses.md freeze v1.0 (script
`scripts/coverage_intersect.py`). Inputs: `crossings/universal_v1`
(`xng-a09e2db7681d`); IBE level-1 discovery snapshots under
`runs/ptf-crossings/` (93 boxes — 86 channel-A star boxes, 7
channel-B antipode boxes; era MJD 54891–57051). An exposure covers an
event at ladder radius r if t_mid ∈ t_ca ± √(r²−b²)/v⊥ and its
nominal CCD footprint (linear WCS from metadata) contains the
predicted source position (channel B: apparent relay position
recomputed per epoch per z ∈ {550, 1000, 2500, 5500, 10000} AU).
Per freeze D5 there is no listing-level quality gate — every
`imgtype='object'`, `fid ≤ 2` row counts; dmask and calibration-gate
attrition happens at the search stage, after coverage freezes.
Coverage counting only; no pixels touched, no signal statistic
formed.

Note vs `era_scope_v0.json`: this run applies the full frozen channel
definition (link direction **and** side of axis), which the scoping
scan omitted — event denominators here are roughly half the scoping
values (e.g. B 1.2 R☉: 24 events, not 48). Every scoped covered
event survives unchanged; the frozen coverage record is this product.

## Channel B — downlink pre-lens (antipode track)

| radius | events | covered | detail |
|---|---|---|---|
| 1.2 R☉ | 24 | **0** | structurally uncovered, as frozen (§2) — ledger entry only |
| 2.5 R☉ | 30 | **3** | van-maanen 2011 b = 0.28 R☉ (2 g epochs, 1 same-night pair); ross-128 2010-09 and 2012-09 b = 1.84 R☉ (2 and 7 R epochs, 1 and 6 pairs) |
| 0.1 AU | 42 | **7** | van-maanen ×4 (2009/2010/2011/2013 — the full April deep-graze recurrence family in-era), ross-128 ×2, wolf-359 ×1 |

Coverage is insensitive to z across the grid (apparent-position
spread ≪ CCD). Epochs per covered wide-rung event: 1–65 (ross-128
2012-09 has 65 epochs / 39 same-night pairs — the deepest-sampled
crossing window in the programme's optical record).

**The headline row:** van-maanen `evt-5c45890fa2` (2011-04-03,
b = 0.28 R☉ — beam axis crossing the photospheric disc) holds 2
g-band epochs forming one same-night pair *inside the 1.34 d
2.5 R☉ window*. And `evt-4c4ea2c397` (2010-04-03) — the event PS1
covered and then lost entirely to correlated chip-gap masks — has 2
independent PTF R epochs in its 0.1 AU window. Between PS1 (2010,
masked), PTF (2009/2010/2011/2013) and the later archives, the
van-maanen April family now has its first multi-instance optical
record.

## Channel A — uplink interception (star position)

| radius | events | covered | detail |
|---|---|---|---|
| 0.1 AU | 42 | **4** | teegarden 2010 (24 R epochs), van-maanen 2014 (11 g), ross-128 2013 (4 R), gj-1276 2010 (3 R) |
| 1.0 AU | 513 | 76 | constraint-only by frozen declaration; 34 targets, 1,206 epochs, 518 same-night pairs |

The A 0.1 AU units face the frozen saturation rule (§6) at the cut
stage: van-maanen and ross-128 are expected excluded-class in R,
teegarden marginal in R — the van-maanen event is g-band, where the
white dwarf (g ≈ 12.2) is also bright; expect the A channel to
thin considerably.

## Calibration and cadence facts

- `photcalflag = 1` epochs are essentially absent from covered
  narrow-rung rows (only ross-128 A's 4 R epochs) — the field-star
  primary calibration (freeze D3) is confirmed as the only viable
  chain.
- Same-night pairs exist in 11 of 14 covered narrow-rung events — the
  same-night repeat veto (§5) is available for most, not all: the
  three pairless events (wolf-359 B — single epoch; van-maanen B
  2009 — 4 isolated-night g epochs; teegarden A — 24 isolated-night
  R epochs) rely on the SkyBoT + rate rungs, recorded per
  exceedance as frozen.
- 12+ of the 93 discovery boxes returned zero exposures (southern
  targets outside the P48 sky and unobserved fields).

## Consequences for the threshold freeze

1. Searched-unit population (before the saturation cut and dmask
   attrition): B 2.5 R☉ — van-maanen (g), ross-128 (R); B 0.1 AU —
   van-maanen (g+R), ross-128 (g+R), wolf-359 (R, single-epoch);
   A 0.1 AU — teegarden (R), van-maanen (g), ross-128 (R),
   gj-1276 (R). Per freeze D4, S_stack applies to van-maanen B
   0.1 AU (4 events) and ross-128 B 2.5 R☉ / 0.1 AU (2 events each);
   everything else is S_event-only.
2. The dev split (D8: wolf-359 B 0.1, gj-1276 A, ross-128 A) holds
   9 covered epochs total — enough to exercise the machinery
   (fetch, calibration, dmask, statistic) but the saturation rule
   will be exercised on ross-128 A exactly as intended.
3. Off-window baselines are healthy at every covered position
   (25–346 epochs per position over the era) — the off-window
   null/variance constructions transfer without a thin-baseline
   gate, except possibly wolf-359 B (55 epochs).

Products: `coverage_v1_events.ecsv` (1,035 rows: per event × radius
× z, per-band epoch/photcal/pair counts),
`coverage_v1_summary.json`; discovery snapshots + query records in
`runs/ptf-crossings/`.
