# LASCO sunward-channels threshold freeze v1.0

Frozen 2026-08-25, after the coverage intersect
(`results/coverage_v1.md`) and before any frame is fetched at a
survey position or any window-locked quantity is formed. Machine
form: `configs/threshold_freeze_v1.json` (sha256 recorded in
`results/threshold_freeze_v1_sha.txt`). Constructions inherit the
established chain (ZTF v1.0+1.1+1.2 rules, TESS shape tests, ATLAS
recurrence stack) with LASCO substitutions per `hypotheses.md`
(FROZEN v1.0, D1–D9).

## 1. Searchable-unit population (from coverage, frozen)

Unit = target × channel × rung × camera. 19 units: 5 × S1 2.5 R☉
(C2, full cadence — the wing regime, visible arc [2.2, 2.5] R☉),
7 × S1 0.1 AU + 7 × S2 0.1 AU (C3, 3-hourly subsample). Per-unit
covered+visible event counts in the config (29–31 everywhere).
S1 1.2 R☉ is `not_constrainable` (0 visible epochs; ledger rows
only). No z-grid exists in this survey: both channels' predicted
positions are z-independent (the axis, not the relay distance, sets
the geometry) — z enters only through the frozen beam-radius rungs.

## 2. Statistics (per unit; S > max(T, 0) everywhere)

All series are detrended (D3): per-frame azimuthal-median radial
profile subtraction (0.1 R☉ bins) → forced photometry at the
predicted position (empirical field-star PSF, per-frame star ZP) →
per-position 30 d running-median subtraction over off-window epochs
→ empirical variance rescale k = median(r²/v)/0.4549 (floor 1) from
off-window epochs.

1. **S_stack** (primary): inverse-variance-weighted mean of
   detrended in-window *visible* epochs, per event, standardized
   against the off-window baseline; stacked over the unit's events
   with per-event inverse-variance weights.
2. **S_event** (secondary): max over the unit's events of the
   per-event standardized excess.
3. **S_pulse** (tertiary, **full-cadence C2 units only** per D6):
   max single-frame standardized excess over visible in-window
   frames (the ≥ 1-frame / 12-min pulse cell).

Trials: 3 per C2 unit, 2 per C3 unit → **43 trials, expected
control crossings 43/9 ≈ 4.8** (dev 11 / 1.2, confirmatory 32 / 3.6).

## 3. Controls

- **Primary T: 8 same-radius PA-ring controls** — the full statistic
  chain on the source trajectory rotated about Sun center by
  ΔPA ∈ {±25°, ±50°, ±75°, ±100°}, same radii, same frames
  (radius-matched to the coronal background; D4). T = max over the 8.
- **Secondary (reported, never T): temporal pseudo-windows** at
  ±23/47/71/97 d on the 1/day off-window series, same-rung exclusion.

## 4. Unit and epoch gates

- Epoch usable: frame passes §7 masks of the hypotheses freeze
  (star-fit ≥ 5 stars, ZP scatter ≤ 0.2 mag, EXPTIME ≤ 2× synoptic
  norm) **and** the local-validity check: ≥ 90 % finite/non-fill
  pixels in the photometry aperture (catches pylon, edge, telemetry
  blocks without modeling them) **and** no major planet within 10 px
  of the aperture.
- Event included: ≥ 10 usable visible in-window epochs.
- S_stack searchable: ≥ 8 included events; else the unit is
  constraint-only (reported, not searched).
- Baseline: ≥ 200 usable off-window epochs (1 frame/day, ±110 d,
  all same-target-channel windows excluded).
- Exceedance adjudication: the frozen §6 veto ladder (known-object
  census incl. Sungrazer comets and planet-bleed annotation;
  multi-frame persistence — single-frame excesses are never
  candidates; TESS split-half ramp + background-anticorrelation
  shape tests; annulus double-passage discriminator; CDAW CME
  annotation; recurrence on the unit's other windows).

## 5. Split (D8; seed 20260825 for all randomized draws)

- **Dev** (5 units, machinery validation): ross-154 S1+S2 0.1 AU,
  gj-908 S1+S2 0.1 AU (wide-b, rung-insensitive), ross-128 S1
  2.5 R☉ (shallowest graze — validates the C2 wing machinery).
- **Confirmatory** (14 units, analysed once, after dev passes with
  any amendments frozen pre-pixel): van-maanen, wolf-359, teegarden,
  gj-1276 — S1 2.5 R☉ (C2) + S1/S2 0.1 AU (C3) each; plus ross-128
  S1/S2 0.1 AU (C3).

## 6. Completeness and positive control (per freeze D5/D6/§8)

Stamp-response injections into the real frames (empirical PSF,
C1-normalized at the calibration reference), spanning flux × window
phase × PA × duty cycle; the Uranus conjunction-passage positive
control (Neptune secondary) through the identical chain, ≤ 0.2 mag
recovery gate; the level-1 overlap era validates the level-0.5 ZP
chain (≤ 0.1 mag median) before any 2017-09+ depth is quoted.

## 7. Order of work

Dev search (5 units) → amendments if dev demands them (frozen before
any confirmatory pixel, versioned v1.x) → confirmatory (14 units,
once) → injections/completeness → adjudication → report
(`report/lasco_crossings.md`).

## Amendment v1.1 (2026-08-26, dev-driven; frozen before any confirmatory pixel)

Dev findings L1–L4 and the statistic-construction iterations
(`notes/dev_machinery_log_2026-08-25.md`) are consolidated here as
the amended frozen chain, machine form
`configs/threshold_freeze_v1_1.json`. The v1.0 construction taken
literally is numerically invalid (occulted-core fill-plateau epochs
carry zero formal error); v1.1 supersedes it as follows.

1. **Epoch validity** (was: frame ZP gate): astrometric solve
   required (roll resolved + translation; ≥ 2 matched stars on C2
   [L1], ≥ 5 on C3), aperture ≥ 90 % valid, finite non-degenerate
   error (err > 0 and ≥ 0.05 × series median — degenerate = occulted
   fill plateau), planet mask (frozen §4, now implemented), and the
   **Tycho-2 VT ≤ 11 star-proximity mask** (3 px; sha-pinned catalog)
   replacing nothing (new; Hipparcos V ≤ 8 was the interim).
2. **Statistic construction**: ring-differential — per epoch,
   source flux minus the *median* of the 8 ring-control fluxes on
   the same frame (controls: ring_i minus the median of the others);
   night-median binning; per-event z = (mean of window nights −
   median of baseline nights) / (SD of baseline nights / √n_win);
   S_stack = Σz/√n, S_event = max z, S_pulse = (max window epoch −
   baseline median)/SD_epoch (C2 full-cadence units only).
   Decision rule unchanged: S > max(T, 0), T = max over rings;
   S/ring-MAD reported alongside for interpretability (C2 z-units
   are under-dispersed; the ring comparison is unit-consistent).
3. **L2 — colour system** (`results/color_ensemble_v1.json`, 8,414
   calibrator records, seed 20260825): global colour coefficient
   c_C3 = 0.429 mag/(B−V) (c_C2 = 0.286, low-n), per-frame C3
   scatter 0.53 → 0.38 after correction. The frozen per-frame 0.2 ZP
   gate is **retired for search-epoch inclusion** (astrometry-only,
   item 1); photometric depths use per-window ensemble ZPs with the
   colour correction and a declared scale systematic (± ~0.1–0.2 mag,
   TESS-style), finalized at completeness.
4. **L3 — S2 stellar template**: for S2 units whose target has
   V ≤ 12.5 (catalog V/B−V table in the config), the predicted
   stellar flux 10^((zp_epoch − [V* + c(B−V* − 0.65)])/2.5)/EXPTIME
   is subtracted from source in-window epochs; epochs without a
   valid ZP are dropped for those units (counted). gj-908's
   window-locked stellar self-detection stands as the survey's
   in-situ positive control.
5. **L4 — low-galactic-latitude class**: fields at |b| ≲ 10° are
   star-transit dominated; ross-154's units carry a
   `high_background` flag — searchable under the deeper mask, but
   any exceedance there is adjudicated against the star-transit
   hypothesis first. All confirmatory targets are high-|b|.
6. **Baseline construction** (documented choice): corona-frame
   positions — the event's median visible radius at the current
   transverse azimuth, window differentials centred on the baseline
   differential median (cancels the stationary F-corona PA bias).

Dev trials under v1.1: 11 (ross-128 C2 ×3; gj-908/ross-154 × S1+S2
×2), expected control crossings 1.2.

## Amendment v1.2 (2026-08-26, dev adjudication round; pre-confirmatory)

Two rules found necessary while adjudicating the v1.1 dev exceedances
(machine form `configs/threshold_freeze_v1_2.json`):

7. **S2 source-star mask exemption**: the Tycho-2 proximity mask
   exempts the (S2, source, in-window) position — the source *is* the
   target star for bright targets (found on ross-154, whose S2 unit
   the mask deleted; gj-908 escaped only because Tycho-2 misses
   high-PM stars). The stellar template (item 4) models the star.
8. **Bright-planet in-FOV epoch veto**: epochs with Venus or Jupiter
   at elongation inside the camera FOV (< 8.5° C3 / < 1.7° C2) are
   dropped from window and baseline series — their stray-light halo
   and bleed columns are frame-wide, far beyond the 10 px proximity
   mask. Found via ross-154 S1: the four z = 20–32 events sit on the
   exact 8-year Venus synodic cycle (2000/2008/2016/2024 early-July
   windows, Venus elongation 5.9–8.0°, verified against the SOHO
   ephemeris).

**Dev adjudications recorded**: ross-154 S1 exceedance = Venus
stray light (vetoed, known-object class; rule 8 prevents it);
ross-128 S1 S_event = 1.45 marginal exceedance on the 2015-03-17
window — the documented St. Patrick's Day CME storm (2013-03-17, the
next-ranked event, is likewise a documented storm; a ring control
shows the same window at 1.23) → adjudicated control-crossing, CME
class, retained in record; CDAW verification tooling queued for the
confirmatory stage.
