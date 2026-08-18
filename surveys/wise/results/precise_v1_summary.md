---
title: "precise_v1 — precise intersection pass over all coarse hits"
date: 2026-08-18
run_dir: "runs/wise/precise_v1 (gitignored; regenerate with scripts/precise_pass.py)"
---

# precise_v1 results

Precise stage of Pipeline A (plan §4.4) over every coarse_v1 hit:
exact per-frame WCS (read from the -msk FITS headers, SIP-aware) plus
usable-pixel bitmask tests, evaluated with sglseti
`covered_z_intervals` at frame mid-exposure (7.7 s frames move the
locus < 1 mas over the exposure; recorded as the interval model).

## Pins

- Source: `coarse_v1` hit set (4,451 unique frame-bands; 12,490
  endpoint × role × frame evaluations)
- Registry v1.0 `sha256:82743c09…`; hypotheses `wise-hypotheses-v1.0`
- Geometry: `tusay2022_eq5_7_v1` v1.1.0, Earth center, z 550–10,000 AU;
  precise tolerance 2", seed step 30" along the locus (usable stretches
  narrower than the seed step can be missed — declared cap), padding 0
  (precise stage tests the nominal locus; the 10" search padding
  belongs to the detection stage)
- Usable-pixel definition: mask fatal bits {0–4, 9, 10–18} (Explanatory
  Supplement IV.4.a NaN set) ∪ {21, 27, 28} (transients)
- Products: 4,443 -msk files, md5-verified, 600 MB under
  `runs/wise/products/msk/` (regenerable via the recorded IBE URLs)

## Outcomes (12,490 evaluations)

| Outcome | Count |
| --- | --- |
| usable (full locus crossing on good pixels) | 10,775 |
| partial (crossing with masked segments) | 1,181 |
| unusable (no usable-pixel crossing; no qual_frame=0 frames occurred) | 534 |

Covered relay-distance interval structure among the 12,312 geometric
hits: 6,898 single-interval (4,168 of them the full 550–10,000 AU),
5,414 split into 2–6 disjoint intervals by masks and frame edges —
34% of hits cover the full range in one piece. Interval-valued
coverage is not an edge case; it is the norm.

## Per endpoint × role × band: evaluated / geometric hit / usable, mean usable fraction

| Endpoint | Role | W1 | W2 | W3 | W4 |
| --- | --- | --- | --- | --- | --- |
| barnard-star | rx | 412/406/396 (0.93) | 413/408/398 (0.92) | 32/32/32 (1.00) | 16/16/16 (0.99) |
| barnard-star | tx | 408/401/393 (0.93) | 407/403/395 (0.93) | 31/31/31 (0.98) | 16/15/15 (0.94) |
| ross-154 | rx | 352/348/342 (0.96) | 354/347/341 (0.91) | 15/15/15 (1.00) | 15/15/15 (1.00) |
| ross-154 | tx | 351/347/341 (0.94) | 352/348/342 (0.94) | 15/15/15 (1.00) | 15/15/15 (1.00) |
| lalande-21185 | rx | 413/408/404 (0.95) | 415/411/407 (0.91) | 15/15/15 (1.00) | 15/15/15 (1.00) |
| lalande-21185 | tx | 416/403/400 (0.92) | 420/409/406 (0.91) | 15/15/15 (1.00) | 15/15/15 (0.97) |
| alpha-cen-a | rx | 462/458/428 (0.95) | 462/457/427 (0.93) | 40/38/38 (0.93) | 20/19/19 (0.95) |
| alpha-cen-a | tx | 465/457/428 (0.95) | 465/458/430 (0.91) | 38/37/37 (0.97) | 19/19/19 (1.00) |
| alpha-cen-b | rx | 462/455/426 (0.95) | 463/457/427 (0.94) | 40/38/38 (0.95) | 20/19/19 (0.95) |
| alpha-cen-b | tx | 467/459/430 (0.95) | 467/459/430 (0.92) | 38/37/37 (0.97) | 19/19/19 (1.00) |
| sirius-a | rx | 433/428/422 (0.94) | 435/429/423 (0.92) | 18/18/18 (1.00) | 18/18/18 (1.00) |
| sirius-a | tx | 433/428/422 (0.95) | 434/430/424 (0.93) | 18/18/18 (1.00) | 18/18/18 (1.00) |
| sirius-b | rx | 432/428/422 (0.96) | 434/429/423 (0.93) | 18/18/18 (1.00) | 18/18/18 (1.00) |
| sirius-b | tx | 435/428/422 (0.93) | 435/429/423 (0.92) | 18/18/18 (1.00) | 18/18/18 (1.00) |

Usable frame-band epochs per endpoint × role (all bands): 713–916.

## Reading

- ~98–99% of coarse hits confirm at the exact-WCS level (coarse padding
  admitted a few boundary frames, as intended), and ~95% of confirmed
  crossings land on usable pixels — WISE masks cost only a few percent
  of epochs, but they restructure coverage into disjoint z-intervals in
  44% of frames.
- Every endpoint × role retains 713–916 usable frame-band epochs across
  ~14 years, including all or nearly all cryo W3/W4 frames. The pilot
  enters the catalog-screening and photometry stages (plan §4.5–4.6)
  with no coverage-starved cell.
- These are geometric/usability products only: still NOT scientific
  coverage — that requires analysis runs with measured detection
  efficiency (plan §3.2). Next: catalog screening (§4.5) and forced
  photometry along the tracks (§4.6), then injection calibration
  (§4.7).
