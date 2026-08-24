# PS1 crossings coverage intersection v1 — results

Run 2026-08-24 against hypotheses.md freeze v1.0. Inputs:
`crossings/universal_v1` (`xng-a09e2db7681d`); ps1filenames discovery
snapshots under `runs/ps1-crossings/` (93 cones — 86 channel-A star
cones, 7 channel-B antipode cones; era MJD 54900–57300, DR2 warps;
skycell WCS cache seeded from the Pipeline-A corridor cache, 122 new
skycells fetched). A warp covers an event at ladder radius r if
t_mid ∈ t_ca ± √(r²−b²)/v⊥ and its skycell footprint contains the
predicted source position (channel B: apparent relay position
recomputed per epoch per z ∈ {550, 1000, 2500, 5500, 10000} AU).
"Primary" = badflag 0 (uniformly true; seeing is unavailable at
listing time). Footprints are nominal skycells: Pipeline-A experience
is ~75 % of nominal hits usable under the exact warp mask, applied at
the search stage. Coverage counting only; no pixels touched, no signal
statistic formed.

## Channel A — uplink interception (star position)

| radius | events | covered (≥1 primary epoch) | median primary epochs (covered) | targets covered |
|---|---|---|---|---|
| 0.1 AU | 46 | 18 | 4 | 7/7 |
| 1.0 AU | 563 | 294 | 12 | 56/86 |

4,098 primary epochs total at 1.0 AU (g 925 / r 1,072 / i 1,604 /
z 208 / y 289), 2,650 same-night TTI pairs. 27 of the 30 uncovered
targets are δ < −30° endpoints outside the PS1 footprint; the other
three (gj-687, gj-1221, sigma-dra) are the ecliptic-pole
annual-minimum family — b ≈ 1 AU sits at the rung edge, so their
windows collapse toward zero duration and the cadence misses them.

## Channel B — downlink pre-lens (antipode track)

| radius | events | covered | epochs per covered event (best z) |
|---|---|---|---|
| 1.2 R☉ | 27 | 1 | 2 |
| 2.5 R☉ | 34 | 1 | 2 |
| 0.1 AU | 47 | 14 | 1–14 (median 3.5) |

All seven narrow-rung targets are covered at the wide rung (gj-1276,
gj-908, ross-128, ross-154, teegarden, van-maanen, wolf-359 — the PS1
3π season is centered on opposition, which is where channel-B windows
sit by geometry). Coverage at the grazing radii is insensitive to z
across the full grid (apparent-position spread ≪ skycell). 65 of the
70 covered wide-rung (event, z) rows contain at least one TTI pair.

**The covered grazing event is the deepest graze in the programme so
far:** van-maanen `evt-4c4ea2c397`, b = 0.28 R☉ (the beam axis crosses
the photospheric disc), t_ca MJD 55289.4 (2010-04-03), windows 0.63 d
(1.2 R☉) / 1.34 d (2.5 R☉), 2 primary i-band epochs forming one
same-night TTI pair. It is a new window: the ZTF era starts eight
years later.

## Consequences for the threshold freeze

1. Grazing rungs are epoch-starved as the cadence predicts (1/27 and
   1/34) — but the one covered event is a two-epoch TTI pair at
   b = 0.28 R☉, so the grazing constraint is a single-event statement
   of unusual geometric quality rather than ZTF's five-window set. The
   0.1 AU rung (14/47 events, all 7 targets, mostly 2–6 epochs)
   carries the channel's statistical weight, with real multi-epoch
   stacking (unlike the single-epoch-dominant expectation).
2. Channel A at 1.0 AU has depth (median 12 primary epochs per covered
   event) but stays constraint-only by frozen declaration; the 0.1 AU
   rung has 18 covered events on all 7 targets, pending the
   saturation cut.
3. TTI pairs are nearly universal in covered rows, so the pair-motion
   veto is available essentially everywhere; PS1 bands are not
   simultaneous, so the chromatic discriminator is cross-epoch only.

Products: `coverage_v1_events.ecsv` (per event × radius × z, per-band
epoch/pair counts), `coverage_v1_summary.json`.
