# PTF crossings threshold freeze v1.0

Frozen 2026-08-26 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json` (`sha256:26e235ee…ef53ae`), bound
to content hashes of `hypotheses.md`, `results/coverage_v1_events.ecsv`
and `results/saturation_cut_v1.ecsv`. No pixel data touched. The
final ZTF/PS1 crossings constructions are adopted at freeze time
rather than rediscovered; the hypotheses-v1.0 substitutions (D3–D8)
are applied verbatim. Conventions inherited from the v2 engine:
8 designated controls, WEIGHT_CAP 20× effective-epoch floor,
fwer α = 0.05, KS α = 0.01.

## Substrate and calibration gate

Star-calibrated forced PSF photometry directly on level-1 science
cutouts (no public difference images). Per-frame zero point from
**PS1 DR2 mean stars** through the frozen band transform (g: PS1 g
direct; R: Jordi 2006 `R = r − 0.153(r−i) − 0.117`); gate: ≥ 5
unsaturated calibrators, robust ZP scatter ≤ 0.2 mag — a failing
frame is unusable, never fallback-calibrated (D3). Header MAGZPT is
a cross-check diagnostic only (`photcalflag ≈ 0` measured at
coverage). Exact dmask (fatal template 65533) defines usable pixels.
Per-epoch variances from the PSF fit + local background, rescaled by
k = median(r²/v)/0.4549 (floor 1) over clipped off-window epochs per
(unit, band).

## Statistics (freeze D4)

Unit = (target, radius, band); channel-B epochs per event are the max
over the nested z family (z multiplies nothing; controls take the
same max).

- **S_event** (primary): max over the unit's multi-epoch covered
  events of the weighted in-window mean excess vs the off-window
  baseline; requires ≥ 1 event with ≥ 2 in-window epochs.
- **S_stack** (secondary): recurrence stack over all covered events;
  applies only where ≥ 2 covered events.
- Single-epoch-only units: reportable single-epoch-exceedance class,
  zero trials; promotion requires recurrence.
- Channel A additionally gates on ≥ 8 usable off-window epochs and
  ≥ 8 valid pseudo-window offsets; the 1.0 AU rung is constraint-only
  (76 covered event rows reported, zero trials); saturation-excluded
  (target, band) units carry no statistic and no contrast claim.

## Controls and rule

- A: 8 temporal pseudo-windows at ±23/47/71/97 d (re-draws
  ±113/127 d), same-rung-only exclusion, mirror-gate validity,
  < 8 valid offsets → constraint-only (D7: assessed, not assumed).
- B: 8 spatial ring trajectories at 20/30/40″, the designated
  construction of ztf-v2/ps1-v2.
- Exceedance: S > max(T, 0), margin reported. **11 trials → expected
  control crossings 1.22** (dev 1, confirmatory 10 → 1.11).
- Veto ladder (calibrated): SkyBoT census; rate test (0.36–6.5″/day
  retrograde); **same-night repeat test** where pairs exist (11 of 14
  covered events; absence recorded per exceedance); flux-consistent
  catalogued-static test (`ptf_objects` snapshot, 2″); recurrence.
  Annotations only: cross-epoch chromatic consistency,
  campaign-cadence structure.
- Quality: no coverage-stage gate (D5); search on the exact dmask;
  strict re-run (seeing ≤ 2.5″, scatter ≤ 0.1) for any exceedance;
  the dev saturation gate verifies E_R 14.5 / E_g 15.0 against a
  bright-star frame (> 0.5 mag disagreement → amendment).

## Saturation cut (channel A, frozen rule §6)

86 era targets classified (`results/saturation_cut_v1.ecsv`): R —
80 excluded / 1 marginal / 5 ok; g — 73 / 4 / 9. Of the four covered
A 0.1 AU units, **only gj-1276 R survives** (R_est 15.47);
teegarden R (13.73), van-maanen g (12.40) and ross-128 R (9.86) are
excluded-class — the A channel thins to one dev unit, exactly the
outcome the frozen rule §6 anticipated. ross-128 A remains the
designated saturation-rule exercise (D8).

## Search units (from the frozen coverage × saturation intersection)

| unit | events (multi-epoch) | statistics | split |
|---|---|---|---|
| van-maanen B 2.5 R☉ g | 1 (1) — 2011 b = 0.28 R☉, same-night pair | S_event | confirmatory |
| van-maanen B 0.1 AU g | 2 (2) — 2009, 2011 | S_event + S_stack | confirmatory |
| van-maanen B 0.1 AU R | 3 (3) — 2010, 2011, 2013 | S_event + S_stack | confirmatory |
| ross-128 B 2.5 R☉ R | 2 (2) — 2010-09, 2012-09 | S_event + S_stack | confirmatory |
| ross-128 B 0.1 AU R | 2 (2) | S_event + S_stack | confirmatory |
| ross-128 B 0.1 AU g | 1 (1) — 45 epochs | S_event | confirmatory |
| gj-1276 A 0.1 AU R | 1 (1) — 3 epochs | S_event | dev |

Plus: wolf-359 B 0.1 AU R single-epoch class (dev, machinery);
3 saturation-excluded A units; B 1.2 R☉ and all remaining rows as
coverage-without-statistic ledger entries. Both headline B families
(the van-maanen April deep-graze recurrence and ross-128) are fully
confirmatory-blind; dev exposure is limited to the wolf-359 antipode
(1 epoch), the gj-1276 star and the ross-128 star.

## Completeness

Injections per searched unit (≥ 100 per event-window cell),
stamp-response (Moffat β = 3 at header SEEING) through the identical
star-calibrated chain (C1 rule), chord temporal profile per event,
line spectra 532 nm (g) / 658 nm (R, generic leakage), v2 recovery
window [2,1]. Positive control: one numbered main-belt asteroid
through the identical chain, recovered to ≤ 0.1 mag.
