# GALEX crossings coverage v1 (freeze v1.0 + amendment v1.1)

2026-08-26. Fresh snapshot-disciplined intersection of
`crossings/universal_v1` with GALEX gPhoton aspect coverage
(`scripts/coverage_intersect.py`; 47 query snapshots under
`runs/galex-crossings/`; per-event rows
`coverage_v1_events.ecsv`, machine summary
`coverage_v1_summary.json`). Gates exactly as frozen: flat-chord
windows, per-event per-z predicted positions (channel A: event star
position; channel B: apparent relay position on the 550–10,000 AU
5-point z-grid at the visit epoch), boresight ≤ 33′, aspect
flag % 2 == 0 (amendment v1.1), per-band live seconds. Coverage
counting only — no photon touched, no statistic formed.

## Population (era 2003-06-07 → 2013-05-01, side-filtered)

Channel A (`inbound`, axis side > 0): 70 events / 7 targets at
b < 0.1 AU. Channel B (`outbound`, side < 0): 69 events / 7 targets
(39 at b < 1.2 R☉, 49 at b < 2.5 R☉). The recon's scoping scan was
un-side-filtered (~140/channel); the side filter halves it without
changing any scoped unit.

## Covered units — the confirmatory five (freeze §9 confirmed)

| unit | t_ca | window | NUV s | FUV s | min boresight |
|---|---|---|---|---|---|
| gj-1276 B 0.1 AU | 2010-03-03 | 11.52 d | **1,637** | 0 | 23.7′ |
| gj-1276 B 0.1 AU | 2007-03-03 | 11.52 d | 109 | 109 | 32.4′ |
| gj-1276 A 0.1 AU | 2007-09-05 | 11.72 d | 92 | 92 | 24.0′ |
| wolf-359 A 0.1 AU | 2007-03-03 | 11.52 d | 97 | 97 | 22.8′ |
| ross-128 A 0.1 AU | 2007-03-17 | 11.53 d | 110 | 110 | 17.7′ |

Totals: 2,045 NUV s, 408 FUV s across 5 units / 5 events / 3
targets (gj-1276 B holds two events — the recurrence pair). Live
seconds are identical across the z-grid for both B units (the
~2–30″ z-position spread sits far inside the 33′ gate); z matters at
the search stage (aperture positions), not at coverage.

## Amendment v1.1 at this stage

The v1.0 aspect gate (`flag = 0`) zeroed gj-1276 B 2010: all 1,637
in-window seconds carry aspect flag 64, the pervasive state of that
field in 2009–2010. gPhoton's pipeline convention
(`PhotonPipe.py` L534–536: usable aspect ⇔ `flag % 2 == 0`; the
`flagDiv2` column) establishes bit 0 as the only bad-aspect bit;
the recon's 0.63″ astrometry check ran on flag-64 seconds. v1.1
(user-approved) adopts `flag % 2 == 0`; the v1.0-gate outputs are
preserved as `coverage_v1_flag0_superseded.*`. No other unit
changed; zero odd-flag seconds occur in any covered window. Dev
must still verify astrometry on flag-64 seconds formally
(off-window).

## Ledger entries (coverage-without-statistic)

- **B 1.2 R☉ (39 events) and B 2.5 R☉ (49 events): zero covered
  windows** — the ±0.3–0.4 d grazing windows never meet the visit
  cadence; incl. the van-maanen deep-graze family (6 in-era events
  per rung with visits at the field in 2004/2008, none in-window).
- **gj-908 A 0.1 AU 2007-09-21: rim-limited** — one in-window visit,
  110 s entirely at 36.3′ boresight distance (> 33′ usable gate);
  recorded, not searchable (freeze §9).
- **teegarden**: antipode has zero GALEX aspect coverage anywhere in
  the era; star has 2 visits, none in-window.
- **van-maanen**: both channels — visits at both fields (star 2006;
  antipode 2004 AIS + 2008 GII), none in-window at any rung.
- **ross-154, gj-908 B, wolf-359 B, ross-128 B**: no in-window
  seconds.

## Next

Threshold freeze: trials tally over the frozen statistic family
(S_rate / S_burst / S_period × band, minus the D1a `forced_dev`
S_rate lane on gj-1276 B 2010), control-ensemble sizes, and the
per-trial thresholds from dev pseudo-units (off-window visits) —
then the blind confirmatory run.
