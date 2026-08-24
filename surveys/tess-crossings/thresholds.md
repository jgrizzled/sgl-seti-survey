# TESS crossings threshold freeze v1.0

Frozen 2026-08-24 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json` (`sha256:af36042e…4d46a`), bound to
content hashes of `hypotheses.md`, the corrected coverage gate, the
refined per-cadence coverage, the TIC cut, and the TESScut recon. No
signal statistic has been formed. Chain to date: gate (corrected for
the HEASARC sector-date typo) → recon cube → hypotheses v1.0 →
coverage refinement (8 cubes, real cadences, spacecraft-UTC window
intersection via TIMECORR) → TIC cut (gj-908 excluded at T = 7.10;
gj-1276 / teegarden / van-maanen ok) → this freeze.

## Search units — channel B, 7 units × 2 statistics = 14 trials

| unit | b | sector (cadence) | in-window primary cadences | off-window |
|---|---|---|---|---|
| wolf-359 1.2 R☉ | 0.55 R☉ | s42 (600 s) | 81 | ≥ 1,956 |
| wolf-359 2.5 R☉ | 0.55 R☉ | s42 | 184 | — |
| wolf-359 0.1 AU | 0.55 R☉ | s42 | 1,116 | — |
| teegarden 1.2 R☉ | 1.04 R☉ | s91 (200 s) | 147 | ≥ 1,853 |
| teegarden 2.5 R☉ | 1.04 R☉ | s91 | 405 | — |
| teegarden 0.1 AU | 1.04 R☉ | s91 | 1,925 | — |
| ross-128 0.1 AU (partial window) | 2.07 R☉ | s42 | 48 | ≥ 921 |

Two frozen statistics per unit: **S_c**, the weighted least-squares
amplitude of the d = 1 chord top-hat (baseline and empirical variance
rescale from the ≥ 900-cadence off-window in-sector series;
WEIGHT_CAP 20×), and **S_p**, the per-cadence in-window S/N maximum
(pulses ≥ 1 cadence; grid-free — the frozen per-unit period grids
parameterize injections only). Controls: 8 spatial ring trajectories
at (±126″, ±189″, ±252″, 0/±189″) — 6/9/12 px — each carrying both
statistics; T per statistic = max over rings; exceedance
S > max(T, 0). **Expected control crossings: 1.6 over 14 trials.**

Unit gates (≥ 20 in-window, ≥ 200 off-window primary cadences): all
7 gate rows pass — the first crossings survey in the programme where
*no* covered narrow-rung unit is lost to cadence, masks, or control
geometry.

## Veto ladder

SkyBoT census mandatory at every pulse-statistic exceedance cadence
(21″ pixels in the near-ecliptic asteroid stream); chord-shape/timing
test against the predicted ingress/egress (flares and asteroid
transits both fail it); straylight adjacency (strict-mask re-run for
any exceedance: +0.05 d block margins, |POS_CORR| > 0.5 px); no
second covered window exists in-archive, so surviving exceedances are
retained-ambiguous with a light-curve morphology report — promotion
requires independent recurrence (a future ecliptic sector or another
archive), never in-window evidence alone.

## Split (seed 20260827, forced)

Dev = **ross-128** (the only non-grazing-target unit; partial window;
validates the full machinery). Confirmatory = **teegarden +
wolf-359** — every grazing unit, nested-window protected, the ZTF
precedent.

## Channel A — 4 reference rows, zero trials

gj-1276 (1,084 cadences, b = 0.97 R☉), van-maanen (1,475,
b = 0.52 R☉), teegarden s44 (460) + s71 (3,424): resolved on-star
crossing light curves and injection-calibrated *reference* depths
(threshold-free, labeled). gj-908 fell out twice over (TIC-excluded
and zero clean in-window cadences).

## Completeness

≥ 100 draws per unit per temporal model (d = 1 chord; d = 0.1 boxcar
on the frozen 6-period log grids, 2 × cadence → window length) via
the SPOC per-camera/CCD PRF through the identical kernel; magnitude
grid T 10–18; recovery against each statistic's own frozen threshold.
Flux scale: per-cutout TIC star calibration verified to ≤ 0.2 mag
before any depth is quoted.

Next: dev search on ross-128, then the fully-confirmatory grazing
units, completeness, `report/tess_crossings.md`.
