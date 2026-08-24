# TESS crossings hypothesis freeze v1.0 (Pipeline B)

Frozen 2026-08-24, after the corrected coverage gate
(`results/coverage_gate_v1.md`) and the TESScut recon
(`notes/tesscut_recon.md`), before any coverage refinement or signal
search. Input: `crossings/tess_v1` (`xng-a943f0f3dbe4`,
**TESS-spacecraft observer**, Horizons −95 SSB table sha-pinned) —
mandatory: the Earth-center list misplaces TESS-frame grazing events
by up to 0.33 R☉ and t_ca by ~3 h, comparable to the grazing windows
themselves. Era MJD 58303–61276 (2018-07-04 → 2026-08-24). Fourth
crossings survey; ZTF/PS1 constructions with TESS substitutions,
marked **[TESS]**.

## Scope (from the corrected gate — sector-table typo clamped)

Narrow-rung coverage exists only in the Ecliptic sectors (42–44, 71,
91; 600 s / 200 s FFI cadence):

| channel / rung | in-era events | gate-covered (full-window) |
|---|---|---|
| B 1.2 R☉ | 32 | 2 (2) — wolf-359 b = 0.55 R☉ s42; teegarden b = 1.04 R☉ s91 |
| B 2.5 R☉ | 40 | 2 (2) — same two events |
| B 0.1 AU | 56 | 3 (1) — + ross-128 partial |
| A 0.1 AU | 57 | 5 (2) — gj-1276 b = 0.97 R☉, van-maanen b = 0.52 R☉ full; teegarden ×2, gj-908 partial |
| A 1.0 AU | 818 | 292 events / 86 targets |

The era's two deepest grazes (van-maanen 0.078 R☉ and gj-1276
0.49 R☉, 2022 spring) fell between sectors and are declared **not
covered** — open cells, testable only by future ecliptic campaigns.

## Channels and the control-geometry declaration **[TESS]**

Channel definitions, beam-radius ladder, and window construction are
the ZTF freeze's, with b, v⊥, t_ca taken from the spacecraft-frame
events table. Two structural declarations, fixed at freeze:

1. **Channel A is constraint-only in v1.** TESS sectors last ~27 d
   and every A window is ≥ 9 d: at most one clean temporal
   pseudo-window fits beside the real one, so the 8-control standard
   is unreachable by geometry (the fourth appearance of this
   theorem). A-channel covered windows are still processed to
   light-curve products and injection-calibrated *reference* depths
   (threshold-free, labeled), because a fully resolved on-star
   crossing light curve at b ≈ 0.5–1 R☉ (gj-1276, van-maanen) exists
   in no other archive. The A 1.0 AU rung inherits its programme-wide
   constraint-only status.
2. **Channel B carries the discovery power.** Its controls are
   *spatial* — 8 ring trajectories at arcminute radii (2.1′ / 3.15′ /
   4.2′ = 6 / 9 / 12 px), valid for any window length — so the
   grazing and wide rungs are searchable with the full frozen
   discipline. The searchable-unit population is decided at the
   threshold freeze from the refined coverage; the gate bounds it at
   ~7 (event × rung) rows on 3 targets, expected control crossings
   ≈ 0.8.

## Band and wavelength **[TESS]**

One broad red band, ~600–1000 nm. Hypothesis flux is a monochromatic
line at the 786 nm effective wavelength converted via the band width
(generic in-band leakage). Notably the TESS band **contains 752 nm**
— the first crossings archive covering the PS1 gj-1276 anomaly's
line wavelength — but at T-band depths (~17–18 stacked) it cannot
test that anomaly's 22.8 AB amplitude; declared honestly. 1064 and
1550 nm remain outside and unconstrained.

## Duty cycle and temporal models **[TESS: the new axis]**

- **d = 1 chord — primary.** The frozen chord top-hat over
  [t_ca ± √(r²−b²)/v⊥] is, for the first time, *resolved*: 86–1,900
  cadences per window. The primary statistic is the weighted
  least-squares amplitude of the chord profile on the per-cadence
  forced-photometry series (baseline and empirical variance rescale
  k = median(r²/v)/0.4549, floor 1, from the off-window in-sector
  cadences; WEIGHT_CAP 20×).
- **Pulsed emission — secondary, first time constrainable.** A second
  frozen statistic per unit: the maximum per-cadence S/N inside the
  window (vs the same maximum over the ring controls) — sensitive to
  pulses ≥ 1 cadence. Injection completeness runs two temporal
  models: the d = 1 chord and a d = 0.1 boxcar-pulsed model on a
  frozen log grid of 6 periods from 2 × cadence to the window length.
  Pulse periods below one cadence integrate to ×d scaling as before.
- **Declared unconstrained:** transmitters scheduled to avoid
  Earth-crossing windows; sub-cadence pulse structure beyond the ×d
  scaling.

## Detection construction **[TESS: TESScut substrate]**

Per-cadence matched-filter forced photometry on TESScut FFI cutout
cubes (25 × 25 px per event; `PIXELS` table: FLUX e-/s with SPOC
FLUX_BKG subtracted per cadence, FLUX_ERR variances, TAN WCS from the
`APERTURE` extension): kernel = Gaussian at the effective PSF FWHM
(1.3–2 px, undersampled — sub-pixel phase carried by the
stamp-response machinery; SPOC per-camera/CCD PRF models are the
injection PSF, the SPHEREx lesson). Track positions per cadence from
the Earth-center apparent-relay computation (TESS-vs-geocenter
parallax at z ≥ 550 AU is ≤ 0.9″ = 0.04 px — budget term); windows
and b strictly from the spacecraft-frame list. Time: TESScut BTJD
(TDB, barycentre-corrected) vs the events table's UTC — conversion
exact; residual light-time terms ≤ 8 min ≪ the 0.3 d shortest window
(budget term).

**Veto ladder (B):** (1) SkyBoT known-object census — the 21″-pixel
near-ecliptic antipode apertures sit in the asteroid stream; an MBA
crosses the aperture in ~hours, so the days-long chord filter is
intrinsically robust, but every *pulse-statistic* exceedance is
census-checked at its cadence; (2) rate/shape test — chord symmetry
and the predicted ingress/egress times vs asteroid transit or flare
profiles (fast-rise-exponential-decay); (3) straylight/scattered-
light annotation — QUALITY flag blocks and background-correlation;
(4) recurrence where a second covered window exists. Channel-A
light-curve products additionally carry flare screening (chord
symmetry) as an annotation.

## Saturation rule (channel A) **[TESS]**

    excluded  T_est < 7.3   (single-cadence saturation ~6.8 + 0.5)
    marginal  7.3 ≤ T_est < 7.8
    ok        T_est ≥ 7.8

T_est from TIC `Tmag` at the cut stage. All seven narrow-rung
targets are T ≈ 11–14: the cut is expected to pass them all — frozen
for the record, verified against TIC before the threshold freeze.

## Quality masks **[TESS, from recon]**

Primary: `QUALITY == 0` (57.6 % of s42 cadences; dominant nonzero
flags are straylight bits in contiguous blocks). Strict (exceedance
re-runs): additionally drop cadences within 0.05 d of any
straylight-flagged block and cadences with |POS_CORR| outliers.
Flux scale: per-cutout star calibration against TIC magnitudes
through the identical kernel (the PS1 star-ZP rule); fallback
T = 20.44 − 2.5 log₁₀(e-/s) with the calibration verified on TIC
stars to ≤ 0.2 mag before any depth is quoted (the SPHEREx gate).

## Expected sensitivity, declared before search

Stacked-window depth T ≈ 17–18 → grazing-cone (2.5 R☉, z = 550)
transmitter-power limits of order kW — ~10× weaker in flux than the
optical surveys, but with the light curve resolved: the deliverable
is the first *temporal* constraint class (chord shape, pulse duty)
on grazing crossings, not a deeper flux limit. Thresholds must not
promise what the two-sector coverage cannot deliver: the searchable
population is ~7 unit-rows and the survey's statistical scale is set
accordingly.
