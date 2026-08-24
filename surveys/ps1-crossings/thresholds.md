# PS1 crossings threshold freeze v1.0

Frozen 2026-08-24 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json`
(`sha256:b6f7c2d1…cf843`), bound to content hashes of `hypotheses.md`,
`results/coverage_v1_events.ecsv`, and
`results/saturation_cut_v1.ecsv`. No pixel data touched. The final ZTF
crossings constructions (freeze v1.0 + amendments v1.1/v1.2) are
adopted here at freeze time rather than rediscovered; PS1
substitutions are listed. Conventions inherited from the v2 engine:
8 designated controls, WEIGHT_CAP 20× effective-epoch floor, split
seed 20260824, dev fraction 0.30, fwer α = 0.05, KS α = 0.01.

## Substrate (both channels) — PS1 substitution

Star-calibrated matched-filter forced photometry directly on DR2 warp
cutouts (per-frame ZP from DR2 `mean` stars through the identical
filter, filter-median fallback where calibrator-sparse) — PS1
publishes no difference images. Exact warp mask (fatal template
16255) defines usable pixels. Per-epoch variances from the warp wt
map, rescaled by k = median(r²/v)/0.4549 over clipped off-window
epochs (floor 1) per (unit, band) — the ZTF v1.2 / v1
confusion-floor rule, frozen from the start for both channels.

## Statistics

**Channel A** (blended star, searchable bands only, 0.1 AU rung):
S_A = weighted stack over covered windows of the per-window mean
excess vs the same star's off-window baseline (robust mean, 3×3σ
clip); one-sided positive. Gates (any failure → constraint-only):
≥ 2 covered windows, ≥ 8 primary off-window epochs in band (resolved
at search time from the snapshotted listings), ≥ 8 valid pseudo-window
offsets. No parallax-factor systematics template: the star is measured
directly (no reference image), so the ZTF v1.1 PM-dipole construction
does not apply. **The 1.0 AU rung is constraint-only by frozen
declaration** (window ≈ observing season under the sidereal-year
phase-lock; ZTF dev lesson 1): coverage + per-window depth only, zero
trials (294 covered event rows reported). The ≤ 10″ station annulus
uses the channel-B construction minus the z-track.

**Channel B** (antipode track): S = max over the nested z family of
the weighted in-window stack along the (event, z) track; controls take
the same max — one trial per (event, radius, band). The static sky is
in the photometry: track nodes within 2″ of a DR2 catalogued static
source (≥ 3 detections) are annotated at search time; ring controls
cross static sources at the same areal rate, so the threshold absorbs
the confusion floor. Single-epoch units form the reportable
"single-epoch exceedance" class; promotion to candidate requires
recurrence at a second covered window.

## Controls and rule

- A: 8 temporal pseudo-window sets at frozen offsets ±23/47/71/97 d
  (re-draws ±113/127 d), same-rung-only overlap exclusion (ZTF v1.1);
  an offset is valid iff its pseudo-windows contain ≥ 2 windows with
  ≥ 1 primary epoch.
- B: 8 spatial ring trajectories, radii 20/30/40″, identical
  designated construction as ztf-v2/ps1-v2.
- Exceedance: S > max(T, 0), margin S − T reported (the R ratio was
  retired by ZTF v1.1). Expected control crossings across all 19
  units: **2.1**. Every exceedance individually adjudicated by the
  frozen veto ladder — no silent drops.
- Veto ladder (calibrated): MPC/known-object census; rate test
  (0.3–6.5″/day retrograde for the fitted z); **TTI-pair motion test**
  (same-night warp pairs ~15–40 min apart: ordinary movers displace
  arcsec-scale between pair members, the relay track < 0.2″ — new vs
  ZTF); flux-consistent catalogued-static test; recurrence.
  Annotations (never sole grounds): cross-epoch chromatic consistency
  (PS1 bands are not simultaneous), single-season/phase-lock
  structure.
- Quality: search on the exact-mask primary set; strict re-run
  (FWHM ≤ 2.5″) reported for any exceedance. **Saturation gate:**
  before any confirmatory channel-A search, the dev stage must verify
  the frozen per-band exclusion levels (E = g/r/i 14.0, z 13.0,
  y 12.0) against ≥ 1 bright-star warp (`CELL.SATURATION` +
  SAT/STARCORE bits); amendment required if off by > 0.5 mag.

## Search units (from the frozen coverage × saturation intersection)

| channel | units | constraint-only |
|---|---|---|
| A (target × band, 0.1 AU, ≥2 windows) | 2 (gj-1276 i ×3 windows, teegarden g ×3) | 2 single-window (gj-1276 g, r) + the whole 1.0 AU rung |
| B (event × radius × band, ≥2 epochs) | 17 | 1 single-epoch (gj-908 i) |

The 17 B units span all seven targets; the three van-maanen
`evt-4c4ea2c397` units (1.2 R☉ / 2.5 R☉ / 0.1 AU, i band, 2 epochs
with a TTI pair, b = 0.28 R☉) are the survey's distinctive grazing
constraint. gj-1276 i `evt-344c12d32a` is the deepest unit
(14 epochs, 13 pairs).

## Dev / confirmatory split (unit = target, seed 20260824)

- **A dev: empty** — every searchable A target is a singleton stratum
  (gj-1276 g+r+i+z+y, teegarden g+r, van-maanen y). A constructions
  stand as frozen (imported from the validated ZTF v1.0+1.1+1.2
  chain); the saturation gate still runs at dev on a bright-star warp.
- **B wide-rung dev (2): teegarden, wolf-359.** Only van-maanen has
  covered grazing events, so only it is dev-ineligible (its grazing
  windows are temporal subsets of its wide-rung windows); the grazing
  rung is fully confirmatory with procedures locked from B-wide dev.
- **B confirmatory:** gj-1276, gj-908, ross-128, ross-154, van-maanen
  (wide), plus the van-maanen grazing units.

## Completeness

Injections per search unit (≥ 100 per event-window cell),
stamp-response (Moffat β=3 at the warp header FWHM) on the warp
substrate through the identical star-calibrated filter, chord temporal
profile from the unit's (b, v⊥, radius), line spectra 532 nm (g) /
617/752/866/962 nm (r/i/z/y band centers, generic-leakage
interpretation), v2 recovery window [2,1].
