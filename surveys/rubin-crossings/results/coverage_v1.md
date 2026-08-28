# Rubin DP2 crossings coverage v1 (2026-08-26)

`scripts/coverage_intersect.py` under hypotheses.md v1.0; TAP queries
snapshotted (`runs/rubin-crossings/coverage_v1/`, 21 snapshots). No
DiaSource row touched; no signal statistic formed.

Era re-measured (D2): **28,698 visits, MJD 60790.117 → 61047.354** —
identical to the recon value; the freeze era stands.

## Channel B (outbound, axis < 0, validity = valid)

| rung | in-era events | with candidate visits | searchable |
|---|---|---|---|
| 1.2 R☉ (0.0056 AU) | 2 | 0 | — ledger |
| 2.5 R☉ (0.0116 AU) | 3 | 0 | — ledger |
| 0.1 AU | 4 | 1 | **ross-128 r** |

**ross-128 B 0.1 AU** (`evt-5197f21cd67b`, t_ca MJD 60937.94,
b_min = 1.855 R☉, v⊥ 29.66 km/s, window ±5.84 d): one in-window
visit, **2025091300612 (r, MJD 60932.28, Δt = −5.657 d)**. All five
z-grid apparent positions (550→10,000 AU; the family spans ~36″ →
2″ offsets from the axis point) fall on **detector 102**, magLim
r = 23.366, seeing 0.86 px-units as served. The visit samples the
beam at **b(t_visit) = 20.92 R☉ = 0.0973 AU** — just inside the
0.1 AU cone (a rim sample of the wide rung, not the 1.86 R☉ core;
the constraint statement must say so).

## Channel A (inbound, axis > 0)

5 in-era 0.1 AU events; one with candidate visits: **ross-154**
(`evt-04c01641edb6`, t_ca 60859.41, b_min = 3.406 R☉) — 16
on-detector visits (i 5 / z 6 / r 1 / y 4, Δt −1.4 → +4.6 d,
magLim 20.6–24.0), including two at Δt = −0.31 d sampling
b ≈ 3.6 R☉. Disposition belongs to the §6 saturation rule at the
threshold freeze (Ross 154 G = 9.13 in the frozen universal-list
input → excluded expected).

## Ledger (geometry-only / structurally uncovered)

- Grazing rungs: 0 candidate visits on all in-era events
  (B 1.2 R☉: wolf-359 b = 0.70 R☉, teegarden b = 0.99 R☉;
  B 2.5 R☉: + ross-128 b = 1.86 R☉).
- A 1.0 AU: 85 in-era events / 84 targets, out of scope (§5.7
  deferred wide-beam rung) — geometry-only.
- In-era 0.1 AU events with no candidate visits: B ross-154,
  teegarden, wolf-359; A gj-1276, gj-908, teegarden, van-maanen
  (van-maanen b = 0.32 R☉ — the deep-graze family's 2025 window
  fell outside DP2's visited sky; per-event rows in
  `coverage_v1_events.ecsv` and `coverage_v1_summary.json`).
