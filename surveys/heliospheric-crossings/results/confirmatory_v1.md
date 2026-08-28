# LASCO sunward-channels blind confirmatory run (2026-08-26)

**0 candidates.** The 14 blind units were fetched, measured, and
reduced ONCE under the frozen chain (threshold freeze v1.0 +
amendments v1.1/v1.2): **30 searched trials, 2 exceedances vs 3.3
expected control crossings — under budget** — with all recurrence
(S_stack) statistics null everywhere. Machine record
`confirmatory_search_v1.json`; measurements
`runs/heliospheric-crossings/dev/measurements_confirmatory_v1.jsonl`
(68,333 frames; fetch 68,330/68,522 ok, 92.7 GB new + 21,477 frames
shared from the dev cache; astrometric validity C2 ~80 % / C3
99.997 %).

## Unit table (S/T per statistic; exc = S > max(T, 0))

| unit | S_stack | S_event | S_pulse | disposition |
|---|---|---|---|---|
| van-maanen S1 2.5 R☉ (C2) | −1.84/0.49 | **4.51/1.19 exc** | 12.8/20.5 | searched; S_event adjudicated (below) |
| wolf-359 S1 2.5 R☉ (C2) | 0.08/0.31 | 0.11/1.20 | 6.6/14.2 | searched-null |
| teegarden S1 2.5 R☉ (C2) | 0.62/2.33 | 3.70/8.67 | 13.3/24.1 | searched-null |
| gj-1276 S1 2.5 R☉ (C2) | 0.25/0.60 | 0.76/1.70 | **18.7/17.4 exc** | searched; S_pulse vetoed (below) |
| van-maanen S1 0.1 AU (C3) | 18.3/35.7 | 13.1/47.3 | — | searched-null |
| wolf-359 S1 0.1 AU | 12.1/12.9 | 22.6/25.7 | — | searched-null |
| teegarden S1 0.1 AU | 2.0/21.0 | 5.9/54.2 | — | searched-null |
| gj-1276 S1 0.1 AU | 3.6/14.7 | 5.0/30.9 | — | searched-null |
| ross-128 S1 0.1 AU | — | — | — | **constraint-only: antipode permanently blended** — the fixed S1 sky position (356.94°, −0.79°) lies 121″ (2.2 px) from a VT 7.2 Tycho star, inside the frozen 3 px mask at every epoch of every window; C3's 56″ pixels cannot separate a ~100× brighter star. Ledger: nominal-covered / resolution-blended |
| van-maanen S2 0.1 AU | 8.1/22.1 | 10.8/44.7 | — | searched-null (stellar template applied) |
| wolf-359 S2 0.1 AU | 1.9/29.5 | 4.9/35.9 | — | searched-null |
| teegarden S2 0.1 AU | 4.5/19.3 | 5.5/24.1 | — | searched-null |
| gj-1276 S2 0.1 AU | 11.4/18.3 | 14.3/25.8 | — | searched-null |
| ross-128 S2 0.1 AU | −0.5/12.6 | 10.0/34.8 | — | searched-null (template applied) |

C3 wide-rung units show large |S| and larger T in both signs — the
coronal-background regime the ring-control construction absorbs by
design; the S > max(T, 0) rule holds everywhere.

## Exceedance adjudications (frozen veto ladder)

1. **gj-1276 S1 2.5 R☉ S_pulse = 18.7 vs T 17.4 — VETOED
   (single-frame / cosmic-ray class).** The pulse peak is one 12-min
   C2 frame (2014-09-04 07:12, differential +2078) with adjacent
   frames at −470/−264 (noise ~10²): the frozen persistence rule —
   single-frame excesses are never candidates — applies directly.
2. **van-maanen S1 2.5 R☉ S_event = 4.51 vs T 1.19 — adjudicated
   transient-corona (CME-period) systematic; retained in record,
   non-promotable.** The excess is one window (2020-10-05/06), all
   27 siblings ≤ 0.17 (S_stack −1.84: no recurrence — the frozen
   promotion rule fails regardless). Epoch structure: a ~3 h
   intermittent burst (21:36 → 00:24, spikes to +168 MAD) during
   which **both ±25° rings swing by hundreds in both signs** — a
   field-wide annulus disturbance, not a point source. CDAW lists a
   CME first seen in C2 at **20:48** (48 min before burst onset; CPA
   263°, with further CMEs 10/06 06:00 CPA 274° w83° and 15:12 CPA
   85° w48°). The catalogued CPAs differ from the source PA (122°) —
   the coupling is mechanical: a bright transient anywhere in the
   annulus shifts the azimuthal-median radial profile that detrends
   every position. Comet and planet hypotheses rejected (no
   catalogued planet within 25°; a comet is PA-localized, unlike the
   observed all-PA disturbance).

## Notes for the report

- Bright-planet in-FOV veto dropped 12,161 source epochs
  (window+baseline) across units; S2 stellar template applied to
  1,791 epochs (1,987 dropped for missing per-epoch ZP on
  template units).
- Construction blemish found at bookkeeping (no statistic re-run,
  blind preserved): baseline day-picks inside window days contribute
  ~5 % of baseline nights from in-window epochs — direction is
  conservative (dilutes a real signal slightly, inflates baseline
  scatter slightly); to be removed at the next freeze version if the
  survey is ever re-run.
- Remaining stages per the frozen order: completeness (stamp
  injections + Uranus positive control + level-1-era ZP validation)
  → `report/lasco_crossings.md`.
