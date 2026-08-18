---
title: "coarse_v1 — coarse discovery pass, all 7 pilot endpoints"
date: 2026-08-18
run_dir: "runs/wise/coarse_v1 (gitignored; regenerate with scripts/coarse_discovery.py)"
---

# coarse_v1 results

First end-to-end run of Pipeline A's coarse stage: sglseti discovery
cones → snapshotted TAP queries against `wise.neowiser_merge_p1bm_frm`
(full mission span MJD 55150–60550, both roles per endpoint) → per-frame
adaptive-locus vs nominal-WCS evaluation at visit-clustered epochs.

## Pins

- Registry: `pilot_wise_2026.yaml` v1.0,
  `sha256:82743c0976c2f7d8f3a994e56f117380f0454691c903807d8c0ee9c0e610969d`
- Hypotheses: `wise-hypotheses-v1.0` (hash recorded per record)
- Geometry: `tusay2022_eq5_7_v1` v1.1.0, `astropy_builtin` ephemeris,
  Earth-center observer, z 550–10,000 AU, tolerance 5", padding 10"
  (+2" cluster-drift allowance), 99% locus uncertainty folded into
  padding at this stage (≤ 0.1" for all pilot endpoints)
- Full config per invocation: `run_config_*.json` in the run dir

## Records emitted

| Record | Count |
| --- | --- |
| QuerySnapshot (raw TAP responses, verbatim) | 7 |
| Observation (unique frame-band products) | 10,545 |
| IntersectionEvaluation (coarse; hits AND misses) | 30,180 |

Run dir: 52 MB. α Cen A/B and Sirius A/B corridors overlap almost
completely, so their observations dedupe (10,545 unique vs 15,090
endpoint-row evaluations per role).

## Coarse hits / evaluated, by endpoint × role × band

| Endpoint | Role | W1 | W2 | W3 | W4 |
| --- | --- | --- | --- | --- | --- |
| barnard-star | rx | 412/998 | 413/997 | 32/81 | 16/39 |
| barnard-star | tx | 408/998 | 407/997 | 31/81 | 16/39 |
| ross-154 | rx | 352/860 | 354/859 | 15/32 | 15/32 |
| ross-154 | tx | 351/860 | 352/859 | 15/32 | 15/32 |
| lalande-21185 | rx | 413/1007 | 415/1006 | 15/39 | 15/39 |
| lalande-21185 | tx | 416/1007 | 420/1006 | 15/39 | 15/39 |
| alpha-cen-a | rx | 462/1113 | 462/1112 | 40/95 | 20/46 |
| alpha-cen-a | tx | 465/1113 | 465/1112 | 38/95 | 19/46 |
| alpha-cen-b | rx | 462/1111 | 463/1112 | 40/95 | 20/46 |
| alpha-cen-b | tx | 467/1111 | 467/1112 | 38/95 | 19/46 |
| sirius-a | rx | 433/1040 | 435/1041 | 18/51 | 18/51 |
| sirius-a | tx | 433/1040 | 434/1041 | 18/51 | 18/51 |
| sirius-b | rx | 432/1043 | 434/1043 | 18/51 | 18/51 |
| sirius-b | tx | 435/1043 | 435/1043 | 18/51 | 18/51 |

## Reading

- Every endpoint × role × band cell has coarse hits: 350–470 W1/W2
  frames per cell across ~14 years, and 15–40 cryo W3/W4 frames — the
  waste-heat-sensitive bands are exercised for all five systems.
- Hit fraction ~40% of evaluated rows is expected: the TAP cone is
  inflated by the frame half-diagonal, so many returned frames cannot
  contain the corridor; misses are retained as audit records (plan
  §3.2).
- Rx and Tx hit sets differ by a few frames per cell (locus epoch
  offsets of 2z/c), confirming the roles are evaluated independently.
- A/B component pairs differ at the single-frame level in W1/W2 — the
  arcsecond-scale component separation only matters downstream
  (photometry), as predicted by `notes/component_vs_barycenter.md`.

## Caveats / next

- Coarse stage only: nominal WCS footprints, no masks or usable-pixel
  tests; a hit here is NOT coverage. Next: precise pass — interval-aware
  swept loci × exact WCS (`int` header incl. distortion) + `msk`
  bitmasks + quality cuts via `covered_z_intervals`, on the coarse-hit
  set (~2,800 unique frame-bands after dedup).
- Frame quality (`qual_frame`, `qa_status`, `moon_sep`, `saa_sep`) is
  recorded on every Observation but not yet cut on.
