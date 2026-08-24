# ZTF crossings coverage intersection v1 — results

Run 2026-08-23 against hypotheses.md freeze v1.0. Inputs:
`crossings/universal_v1` (`xng-a09e2db7681d`); IBE discovery snapshots
under `runs/ztf-crossings/` (93 cones — 86 channel-A star cones,
7 channel-B antipode cones; era MJD 58178–61275, public sci quadrants).
An epoch covers an event at ladder radius r if t_mid ∈ t_ca ±
√(r²−b²)/v⊥ and the nominal footprint contains the predicted source
position (channel B: apparent relay position recomputed per epoch per
z ∈ {550, 1000, 2500, 5500, 10000} AU). "Primary" = v2 quality mask
(archival bad-quality false, seeing ≤ 4″). Coverage counting only; no
pixels touched, no signal statistic formed.

## Channel A — uplink interception (star position)

| radius | events | covered (≥1 primary epoch) | median primary epochs | targets covered |
|---|---|---|---|---|
| 0.1 AU | 59 | 44 | 3 | 6/7 |
| 1.0 AU | 720 | 402 | 12 | 57/86 |

20,796 primary epochs total at 1.0 AU (g 10,211 / r 10,197 / i 675).
All 29 uncovered targets are δ < −30° endpoints outside the ZTF
footprint (verified against the events table); every reachable target
has at least one covered event.

## Channel B — downlink pre-lens (antipode track)

| radius | events | covered | median primary epochs |
|---|---|---|---|
| 1.2 R☉ | 35 | 1 | 0 |
| 2.5 R☉ | 43 | 5 | 0 |
| 0.1 AU | 60 | 33 | 1.5 |

Coverage at the grazing radii is insensitive to z across the full
550–10 000 AU grid (apparent-position spread ≪ quadrant), so the z
family survives to the search stage undiminished. The covered grazing
events, all with same-night g+r pairs where noted:

- gj-1276: t_ca 58910.6 (r only), 59275.9 (g+r), 59641.1 (g+r; also the
  single 1.2 R☉-covered event)
- teegarden: t_ca 58244.7 (g+r), 61166.7 (g+r)

## Consequences for the threshold freeze

1. The solar-grazing downlink hypothesis is epoch-starved as predicted:
   ~0.6–1.3 d windows vs ~2–3 d cadence → 1/35 and 5/43 events covered,
   ~1–2 epochs each. The B-channel constraint at grazing radii is a
   per-event statement about 5 windows on 2 targets, not a survey-wide
   depth. The 0.1 AU wide-beam rung (33/60 events, 7 targets) carries
   the channel's statistical weight.
2. Channel A has real depth (median 12 primary epochs per event at
   1 AU) but is contrast-limited on bright stars; the per-target
   saturation cut must be applied before any threshold is frozen.
3. Both channels retain enough same-night g+r pairs for the chromatic
   discriminator.

Products: `coverage_v1_events.ecsv` (per event × radius × z),
`coverage_v1_summary.json`.
