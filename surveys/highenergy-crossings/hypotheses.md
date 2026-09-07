# High-energy crossings hypothesis freeze v1.0 (Pipeline B; Fermi-LAT primary, Swift/BAT second family)

Drafted 2026-09-07 after the reachability recon
(`notes/highenergy_recon_2026-09-07.md`, results `results/recon_*_v0.json`),
before any further in-window data contact (pre-freeze contact declared
in §10 — one recon probe touched in-window photons; remedy D1).
**FROZEN 2026-09-07: the §11 decisions below are the
recommendations made concrete.** 
Eleventh crossings survey; the
construction transfers the GALEX freeze
(`surveys/galex-crossings/hypotheses.md` v1.0: photon-event substrate,
pulse cell) with the ATLAS freeze's recurrence-stack layer (per-window
events stacked over the era) and the TESS pulse-cell construction —
LAT-specific substitutions marked **[LAT]**.

Input: `crossings/universal_v1` (`xng-a09e2db7681d`, Earth-center,
1980→2028, both link directions, every b(t) minimum kept; target-state
uncertainty propagated 2026-09-06, ≤ 1e-4 R☉ on the grazing family).
Fermi and Swift are ~550 km LEO: Earth-center events valid under the
standing 0.010 R☉ / ~4 min budget. Era: **2008-08-04 15:43:36 →
2026-09-07 00:13:13 UTC** (LAT weekly files w009 → w953, the last
complete week at freeze).

## 1. Observable channels

Channel A (uplink interception, `inbound`, blended with the star,
star at solar opposition) and channel B (downlink pre-lens
interception, `outbound`, at the star antipode, z-parameterized
parallax-reflex track). Both sit at solar elongation 180° ± 6° inside
every window (the elongation gate that removes XMM-Newton and eROSITA
from this row; recon §"Geometry"). Sunward combinations out of scope
(§5.11/§5.16/§5.19/§5.23 own them). **[LAT]** The LAT PSF (68 %
containment ~5° at 100 MeV, ~0.8° at 1 GeV) is far larger than the
whole z-family locus (≤ 375″) and than the A/B channel offset from
the star (0 for A; the antipode is a fixed point for B): **one
aperture per channel position per event**, no z grid, no track
photometry — the SGL geometry enters only through the window times
and the two positions.

## 2. Beam-radius ladder

The standard ladder; all three rungs searched on both channels:

- **B 1.2 R☉ / A 1.2 R☉** (windows ±2.5–7.6 h; 4 targets: gj-1276,
  teegarden, van-maanen, wolf-359),
- **B 2.5 R☉ / A 2.5 R☉** (±11–16 h; + ross-128),
- **B 0.1 AU / A 0.1 AU** (±5.8 d; all 7 deep-family targets incl.
  gj-908, ross-154),
- **A 1.0 AU: out of scope** (programme-wide deferred status).

Flat-chord windows t_ca ± √(r² − b²)/v⊥ as in every crossings survey.
Rungs are nested hypotheses on the same photons (the standing
convention); each is its own unit.

## 3. Band **[LAT]**

Two lanes, frozen from the P8R3_SOURCE_V2 PSF (CALDB
`psf_P8R3_SOURCE_V2_FB.fits`, snapshotted under `runs/…/v1/caldb/`):

| lane | energy | aperture radius | role |
|---|---|---|---|
| L | 100 MeV – 1 GeV | 2.0° | the count-rich lane (λ ≈ 5–15 per grazing window) |
| H | 1 – 300 GeV | 0.7° | the near-zero-background lane (λ ≈ 0.5–2) |

A relay signal at any energy inside 0.1–300 GeV is the hypothesis; a
spectral model (power law Γ = 2) is used **only** for the exposure
weighting and the flux conversion of the injections, declared with
the ±0.3 in Γ sensitivity reported at completeness. Swift/BAT
(15–150 keV, §8) is the second family.

## 4. Duty cycle and temporal models

- **d = 1 persistent (window-integrated), stacked over the era —
  primary** (`S_stack`): total in-window counts over all searchable
  windows of the unit vs the total expected from the pseudo-window
  rate. The ATLAS recurrence cell: a relay that transmits at every
  crossing.
- **Single-window persistent** (`S_event`): the maximum over windows
  of the per-window excess — a relay active at one crossing.
- **Burst/pulse** (`S_burst`): the maximum over windows and over
  sliding 1 ks boxcars (100 s step) of the in-aperture count in the
  union band (0.1–300 GeV, r = 1.0°) — pulses shorter than a window.
- **Declared unconstrained:** structure below the 1 ks boxcar (photon
  numbers per window are ≲ 15 — no sub-ks test has power), coherent
  trains (no H-test: the LAT photon rate is too low for a period
  grid to be meaningful), and schedules that avoid Earth-crossing
  windows.

## 5. Detection constructions **[LAT: photon-event substrate]**

- **Photons.** Per (channel, target) one era-long LAT data-server
  query at the era-mean channel position (radius 3°, 100 MeV –
  300 GeV, "Photon" class), snapshotted; gates `EVENT_CLASS & 128`
  (P8R3 SOURCE), `EVENT_TYPE & 3` (front+back), zenith ≤ 100°, and
  the **good-time gate from the weekly spacecraft files**: a photon
  counts only if its 30-s interval has boresight angle θ ≤ 60° to the
  channel position, `DATA_QUAL > 0`, `LAT_CONFIG = 1`, and the
  interval is outside the Moon/Sun exclusion (below). The same
  intervals define the livetime and exposure of every window and
  pseudo-window — photons and exposure share one gate.
- **Exposure.** Per gated 30-s interval: livetime × A_eff(cosθ) with
  A_eff the Γ = 2 photon-weighted mean of the CALDB front+back
  effective area over the lane (`aeff_P8R3_SOURCE_V2_FB.fits`);
  φ-dependence and the efficiency (livetime-fraction) correction
  ignored — declared as a ≤ 10 % systematic on absolute fluxes,
  irrelevant to the relative statistics.
- **Moon / Sun exclusion.** The full Moon sits at the anti-Sun point
  — i.e. at the channel positions — once per synodic month, and the
  Moon is a ~10⁻⁶ ph cm⁻² s⁻¹ LAT source. Any 30-s interval with the
  geocentric Moon or Sun within **8°** of the channel position is
  excluded (from photons and exposure alike, real and pseudo windows
  alike; the ~1° LEO lunar parallax is inside the margin).
- **Null / rate model.** Per unit and window k, pseudo-windows of the
  same duration at offsets ±j·D from t_ca, D = window duration (min
  1 d), for all |offset| ≤ 60 d (grazing rungs) or ≤ 120 d (0.1 AU),
  excluding |offset| < 1.5 D; only pseudo-windows with exposure
  ≥ 50 % of the real window's are kept (else the exposure match is
  poor). Expected count λ_k = (Σ_pw n / Σ_pw exposure) × exposure_k,
  per lane. `S_stack` = −log10 P(N ≥ Σ n_k | Σ λ_k) (Poisson, upper
  tail), `S_event` = max_k −log10 P(N ≥ n_k | λ_k), `S_burst` =
  max over windows and boxcars of −log10[1 − (1 − p_box)^{N_box}]
  with p_box the Poisson upper tail of the boxcar count at λ_box =
  rate × boxcar exposure and N_box = the number of independent
  boxcars (window exposure / boxcar exposure), i.e. a per-window
  look-elsewhere correction inside the statistic.
- **Coverage gate.** A window is searchable when its lane-L in-FoV
  livetime is ≥ 1 ks (post-2018 windows with the boresight held
  away from the anti-Sun point are coverage-without-statistic; the
  recon table gives every value). A unit is searchable when ≥ 1
  window is.
- **Veto ladder (adjudication of any exceedance, frozen order):**
  (1) **4FGL-DR4 association** — the excess photons' centroid within
  the lane's 68 % PSF radius of a 4FGL-DR4 source, or a 4FGL source
  within the aperture whose 4FGL variability index > 18.48 (the
  catalogue's 99 % threshold) → `known_source`; (2) **Moon/Sun/Earth
  limb** — re-run with the exclusion radius doubled and zenith ≤ 90°;
  (3) **solar-system census** — planets (Jupiter is not a LAT source;
  no planet veto beyond the Moon) — recorded as not applicable;
  (4) **detector/attitude** — the interval-level `DATA_QUAL`,
  `LAT_MODE` and rocking-angle transitions within the window;
  (5) **recurrence** — the same unit's other windows (`S_stack`
  already tests it) and, for channel A, the star's own X-ray/UV
  flare record (M dwarfs do not flare at > 100 MeV: an A-channel LAT
  excess has no stellar-flare interpretation, unlike GALEX/XRT).
- **Dispersion gate.** Per trial, the pseudo-window p-values must be
  consistent with uniform: if > 10 % of the pseudo-windows have
  p < 0.05 (over-dispersion from a variable 4FGL source in the
  aperture or exposure mis-modelling), the trial is `constraint_only`
  — its threshold is not calibrated.

## 6. Saturation / pile-up rule (channel A)

Not applicable **[LAT]** — no target star is a γ-ray source; the
A-channel aperture is sky like the B-channel aperture. Recorded so
the standing rule has an entry.

## 7. Quality gates

- Photon: SOURCE class, front+back, zenith ≤ 100°, inside a gated
  interval (§5).
- Interval: θ ≤ 60°, `DATA_QUAL > 0`, `LAT_CONFIG = 1`, Moon/Sun
  > 8°.
- Strict (exceedance re-runs): θ ≤ 50°, zenith ≤ 90°, Moon/Sun
  > 16°.
- Window: lane-L livetime ≥ 1 ks (§5).
- Depths from the injection chain only (C1). Scoping expectation,
  declared before search: λ_L ≈ 8 per 1.2 R☉ window, ≈ 20 per 2.5 R☉,
  ≈ 200 per 0.1 AU; a 5σ persistent excess ≈ 3× those; single-window
  lane-L flux floor ~10⁻⁷ ph cm⁻² s⁻¹ (> 100 MeV). Expectations, not
  constraints.

## 8. Swift/BAT second family **[BAT]**

The recon found ≥ 1 ks of BAT survey exposure within 20° of the
boresight in 7–11 of ~22 grazing windows per target-channel and in
~90 % of the 0.1 AU windows. The BAT arm is **constraint-only at
freeze** (its 15–150 keV floor is mCrab-class, its arbitrary-position
light curves need the HEASoft `batsurvey` chain, and its systematics
are not characterised on this machine): it runs after the LAT chain,
on the LAT-searchable windows with ≥ 1 ks inside 20°, as
`batsurvey` per-pointing rates at the 14 positions vs the same
position's out-of-window survey points; no confirmatory trial, no
threshold, a per-window rate limit only. If HEASoft cannot be run
(no native install; a container is the route), the arm is recorded
`blocked` with the reason and the recon coverage table stands as its
ledger.

## 9. The pointed unit and the corridor screen

- **Swift XRT unit** (wolf-359 A 0.1 AU, obsid 00010119025,
  2018-03-06, PC 4.2 ks, +2.9 d from t_ca): a pre-registered
  **look**, not a trial — the on-star 0.3–10 keV light curve of the
  visit in 100-s bins against the star's 2017 Swift campaign (113
  visits, off-window) as the null; any flare is reported with the
  §5.12 flare-morphology description; no threshold, no constraint
  claimed. The UVOT UV-grism spectrum of the same visit is handed to
  §5.17 (recorded, not analysed here).
- **Corridor screen** (Pipeline A step 1, separate hypothesis file
  `corridor_screen.md`): eRODat upper limits at the 88 anti-star
  positions (DR1 eRASS1 + DR2 eRASS:3), catalogue cones (5XMM-DR15
  per-detection, CSC 2.1, 2SXPS/LSXPS, eRASS1/eRASS:3, BAT-157m,
  4FGL-DR4) and an **SGL-track test** on every hit within 7′: the
  relay grid 550–10,000 AU has an annual parallax of 20–375″, so any
  source whose positions at epochs ≥ 3 months apart agree to < 10″
  (or that is a compact detection in a multi-scan stack) is
  `fixed_sky`, i.e. excluded as a relay at every z on the grid;
  single-epoch hits are `single_epoch_unresolved`. Ledger rows
  `catalogue_screen_only`; no flux search.

## 10. Pre-freeze data-contact declaration

1. **Recon photon pull at a confirmatory locus**: the LAT data-server
   test query (van-maanen B, 5°, 2020-04-01 → 04, ≥ 100 MeV; 150
   photons) examined in-window counts within 0.5/1/2/5° and > 1 GeV
   — rate-class contact on the van-maanen B 2020-04-02 event at all
   three rungs (the 3-day span contains the 0.1 AU window's centre
   and the whole grazing windows). No background, exposure or
   pseudo-window comparison was formed. **Remedy D1.**
2. **Era-long photon pull at the van-maanen B position** (the
   query-timing test): only the per-file photon totals and time
   spans were read (six ~3-year chunks); no window association.
   Declared, no remedy needed beyond D1.
3. **Spacecraft-file livetime** at all 14 positions for every
   grazing window and every 0.1 AU window (the recon tables) — the
   class every coverage stage touches pre-threshold-freeze.
4. **Master-table and catalogue metadata** (HEASARC, XSA, CSC, eRODat
   cones and upper limits) at the 14 channel positions and 88
   corridors — annotation-class.
5. The eRASS1 tile 194096 event list was opened to read scan times
   near the van-maanen antipode (100 d off-window) — off-window.

## 11. Freeze decisions (adopted 2026-09-07)

- **D1 — pre-freeze contact remedy**: the van-maanen B **2020-04-02
  event is `forced_dev`** — excluded from the van-maanen B
  confirmatory statistics at all three rungs (it becomes a dev
  window for machinery checks); the other 17 van-maanen B windows
  stay blind. All recon photon files stay under `runs/…/recon/`,
  the survey re-pulls fresh.
- **D2 — era and input**: 2008-08-04 → 2026-09-07 (w009–w953);
  `crossings/universal_v1` Earth-center, standing LEO budget; no LAT
  observer list.
- **D3 — substrate**: §5 (data-server era-long photon pulls; SOURCE
  class; shared interval gate for photons and exposure; CALDB
  effective area; Moon/Sun 8° exclusion).
- **D4 — units, statistics, trials**: unit = (target, channel,
  rung): 32 units (8 at 1.2 R☉, 10 at 2.5 R☉, 14 at 0.1 AU). Five
  frozen statistics per unit — `S_stack` and `S_event` in lanes L
  and H, `S_burst` in the union band. Dev units (D8) excluded:
  **30 confirmatory units × 5 = 150 trials**, FWER α = 0.05 (Šidák
  per-trial α = 3.42 × 10⁻⁴, T = 3.466 in −log10 p), as in every
  survey.
- **D5 — gates**: §7.
- **D6 — thresholds**: analytic Poisson thresholds at T (D4),
  validated per trial on the pseudo-window ensemble by the §5
  dispersion gate (> 10 % of members with p < 0.05 →
  `constraint_only`); the ensemble is also the empirical check that
  the exceedance rate over all pseudo-window members is ≤ 2× the
  nominal at T. Seed 20260907 for every random draw (injections).
- **D7 — controls**: pseudo-windows per §5 as the null family;
  **positive controls** — (i) GRB 130427A (RA 173.136, Dec +27.699,
  T0 2013-04-27 07:47:06 UTC) run through the identical chain as a
  fake unit with a 1.2 R☉-length window centred on T0: `S_event` and
  `S_burst` must exceed T; (ii) the 3C 454.3 November-2010 flare
  (RA 343.491, Dec +16.148; window 2010-11-17 → 2010-11-22) as the
  persistent-excess control: `S_event` lane L must exceed T. Both
  recorded whatever the outcome.
- **D8 — dev/confirmatory split**: **dev = gj-908 A and B (0.1 AU
  only — its b_min is 12 R☉, no grazing rungs) plus the D1 window
  plus the positive controls**; **confirmatory = the other 30 units,
  blind**. The dev units also carry the machinery checks: gate
  census, exposure vs data-server GTIs, pseudo-window uniformity,
  Moon-exclusion effect.
- **D9 — completeness**: photon-level injections into the real
  gated photon lists: persistent in-window excess at flux F (Γ = 2,
  PSF-scattered positions from the CALDB PSF, Poisson counts from the
  window exposure) for `S_stack`/`S_event`, and 1-ks pulses of
  fluence Φ for `S_burst`; 90 % recovery at the frozen thresholds per
  unit and lane; flux → power through the rung cone with the §5.17
  convention (flat-top footprint of radius b_rung, P = F_E · π b²,
  F_E the energy flux over the lane) and reported with the Γ ± 0.3
  and ±10 % exposure systematics.
- **D10 — BAT**: §8, constraint-only, after the LAT chain.
- **D11 — the XRT look and the corridor screen**: §9.

## Amendment v1.1 (2026-09-07, at the threshold freeze, before any confirmatory statistic was formed except as declared here)

1. **Jeffreys floor on the pseudo-window rate.** The first machinery
   test (A gj-1276 1.2 R☉) showed lane-H pseudo-window pools with a
   single photon: the leave-one-out rate is then exactly zero and any
   photon gives p = 0. Rule: every rate estimate — real-window λ_k and
   the leave-one-out validation λ — is (Σn + ½)/Σexposure (the
   Jeffreys prior), applied identically to units, pseudo-windows,
   controls and injections. Thresholds untouched.
2. **A gj-1276 1.2 R☉ demoted to `forced_dev`.** Its five statistics
   were printed during that machinery test (all below T: S_stack_L
   0.04, S_event_L 0.80, S_stack_H 0.69, S_event_H 1.86, S_burst 2.10),
   before the threshold freeze. The unit leaves the confirmatory
   family (constraint-only, labelled) and joins the dev set; the
   family is now **29 units × 5 = 145 trials**, Šidák per-trial
   α = 3.54 × 10⁻⁴, T = 3.451. No other real-window quantity was seen.

## Amendment v1.2 (2026-09-07, at the threshold freeze, off-window data only)

**Per-trial ensemble exceedance gate.** D6 stated the exceedance check
over all pseudo-window members globally. The ensembles showed why it
must be per trial: the B van-maanen lane-L apertures contain the
blazar 3C 279 (1.8° from the antipode), whose flares put 17 / 21 / 33
pseudo-windows above T at the three rungs (0.5 expected each) and
pseudo-stack values up to 76, while the 10 % dispersion gate caught
only the 2.5 R☉ and 0.1 AU rungs. Rule: a trial is `constraint_only`
when its pseudo-window exceedances at T exceed max(3, 3 × expected).
Applied before any confirmatory statistic was formed; T and the trial
count (145) unchanged. Result: 9 constraint-only trials, all on B
van-maanen (S_stack_L, S_event_L, S_burst at the three rungs); the
remaining 52,686 pseudo-windows show 15 exceedances vs 18.6 expected.
