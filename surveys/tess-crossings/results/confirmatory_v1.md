# TESS crossings confirmatory run v1 — results

Run 2026-08-24 under threshold freeze v1.1
(`confirmatory_search.py` → `results/confirmatory_v1.json`;
adjudication `adjudicate_v1.py` → `results/adjudication_v1.json`).
6 channel-B units (wolf-359 s42 and teegarden s91, three nested rungs
each), two frozen statistics per unit, nested z family made explicit
(unit statistic = max over the 550–10,000 AU grid; controls take the
same max — the programme convention). **No candidates.**

## Results

| unit | b | S_c | T_c | margin | S_p | T_p | k |
|---|---|---|---|---|---|---|---|
| wolf-359 1.2 R☉ | 0.55 | 12.5 | 15.3 | −2.8 | 2.5 | 4.5 | 4.9 |
| wolf-359 2.5 R☉ | 0.55 | 21.4 | 20.9 | **+0.5** | 2.9 | 5.5 | 3.9 |
| wolf-359 0.1 AU | 0.55 | −19.4 | 927.6 | −947 | 3.9 | 137.2 | 1.7 |
| teegarden 1.2 R☉ | 1.04 | 17.2 | 17.3 | −0.03 | 3.0 | 4.8 | 1.3 |
| teegarden 2.5 R☉ | 1.04 | 48.3 | 24.9 | **+23.4** | 5.2 | 7.8 | 1.3 |
| teegarden 0.1 AU | 1.04 | 79.7 | 52.1 | **+27.6** | 6.9 | 7.9 | 1.0 |

**The pulse statistic — the survey's genuinely new cell — is null in
all six units (0 exceedances):** no pulse of ≥ 1 cadence (200–600 s)
stands above the ring-control maxima anywhere in the two fully
resolved grazing crossings or the wide windows.

**The chord statistic is systematics-dominated everywhere** (S and T
in the tens; one control at 927 where a ring trajectory crosses the
bright field star that also sets the PSF fit — that unit, wolf-359
0.1 AU, is thereby deeply insensitive but null). Sector-scale
scattered-light drift is not removed by the k-rescale (variance-only)
and only partially shared by the ring controls; the 1/9
control-crossing budget is approximate at best for this statistic —
recorded as the survey's main design lesson (a v2 would need a
detrending/differential layer under the chord filter). 3 chord
exceedances observed vs 1.3 expected (P(≥3 | 1.3) ≈ 0.14).

## Adjudication (frozen ladder; census N/A — no pulse exceedance)

- **teegarden 2.5 R☉ (margin +23.4): VETOED by the chord-shape/timing
  test** — split-half amplitudes 15.1 / 48.2 (ratio 0.31): a
  monotonic ramp, not a flat-topped chord; background correlation
  −0.785 (strong anti-correlation = the background-subtraction
  residual signature); the window sits in the sector's final day.
  Strict-mask re-run persists (it is a drift, not a flagged-cadence
  artifact) — the shape test, not the mask, is the calibrated
  discriminator here.
- **wolf-359 2.5 R☉ (margin +0.46): adjudicated control-crossing
  statistic** — the marginal exceedance the 1/9 budget exists to
  absorb; split-half consistent (0.82) but the nested-window
  cross-check contradicts a real source: the 0.1 AU superset window
  (which contains every 2.5 R☉ cadence plus ten more days) has
  S = −19.4, and the inner 1.2 R☉ window does not exceed. A relay in
  the beam at b = 0.55 R☉ is in-beam for all three rungs
  simultaneously; the pattern is drift, and the margin is 2 % of the
  statistic's own scale.
- **teegarden 0.1 AU (margin +27.6): retained-ambiguous,
  non-promotable** — the frozen ladder does not cleanly veto it: the
  shape test is weakened by window truncation (only the first 7.5 d
  of the 11.9 d window are in-sector, so no egress exists to test),
  background correlation is moderate (−0.22), and the strict re-run
  persists. Annotations weigh heavily toward the same sector-end
  systematic that was *vetoed* in its sibling (which shares the same
  final-day cadences), but per the v2 discipline it is not dismissed
  by discretion. It cannot be promoted: the freeze requires
  independent recurrence, and teegarden's next covered window needs a
  future ecliptic sector — recorded as the designated follow-up,
  exactly as PS1's gj-1276 anomaly was handed to ZTF.

## Bottom line

0 candidates. The pulse-period cell (200 s → window length) closes
clean on both fully resolved grazing crossings — the first constraint
of its kind. The chord channel yields 1 vetoed exceedance, 1
budget-absorbed marginal, and 1 retained-ambiguous row awaiting a
future-sector recurrence test, atop the recorded lesson that the
undetrended chord statistic saturates its own error budget with
scattered-light drift. Remaining: injection completeness (both
temporal models) and `report/tess_crossings.md`.
