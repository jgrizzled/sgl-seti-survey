# S2 1 AU outer-skirt blended search — hypothesis freeze v1.0 (Pipeline B, plan §5.25)

Drafted 2026-09-07 after the recon
(`notes/skirt_recon_2026-09-07.md`), before any archive data is
touched at a survey position (two non-target probe stars only, §10).
**FROZEN 2026-09-07**: the user approved every §11 decision as recommended (D1–D8) before any survey-position data was pulled. Constructions inherited from
the ATLAS crossings freeze (`surveys/atlas-asassn-crossings/
hypotheses.md`) are marked **[ATLAS]**; the new cell's constructions
are marked **[SKIRT]**.

Input: `crossings/universal_v1` (Earth-center, 1980→2028). Earth-
center is valid for ground observers (site baseline ≤ R⊕, §1 of the
ATLAS freeze) and the geometry of this cell is a function of solar
elongation alone (§3), so no derivative observer list is needed.

## 1. The cell and why it is searchable at all **[SKIRT]**

The 2026-08-25 sunward geometry study left one sunward combination
outside the coronagraph substrate: **S2** (uplink `inbound`, Earth on
the anti-target side), apparent source = the target star at solar
elongation ε near conjunction, Earth's distance from the Sun–star axis
b_e = r_E sin ε (exact; r_E the Earth–Sun distance). The frozen 1 AU
uplink rung ("~1 m-class transmitter, effectively any minimum") spans
ε up to ~89°, and its outer skirt at ε ≳ 30–40° is ordinary
evening/morning sky — night-sky visible, hence not exclusive to the
heliospheric imagers.

Two structural facts (recon §2) shape the search:

1. **The 1 AU rung itself has no temporal signature.** b_e ≤ 1 AU for
   every ε, so Earth is inside a 1 AU-radius beam all year; a
   persistent transmitter is a constant excess blended with the star,
   indistinguishable from the star. The 1 AU rung is therefore a
   **coverage ledger** (epochs per S2 window at ε ≥ the ATLAS floor),
   not a statistic — the same status the ATLAS survey gave the
   channel-A 1.0 AU rung (window ≈ season).
2. **Sub-1 AU rungs have an elongation-locked edge that IS a
   signature.** A top-hat beam of radius r_b < 1 AU is "on" only while
   r_E sin ε < r_b, i.e. ε < ε_b = arcsin(r_b/r_E) on the S2 side and,
   by the same inequality, ε > 180° − ε_b on the channel-A side — the
   beam passes the Sun at b ≥ 0.5 AU unlensed and unocculted
   (deflection 0.016″), so it is the same beam on both sides. A
   sub-1 AU uplink therefore predicts an annually recurring,
   **elongation-symmetric step**: the star is brighter by a constant
   F_tx inside ε < ε_b (conjunction skirt) and ε > 180° − ε_b
   (opposition) than in the quadrature band between them. Rungs whose
   edge ε_b lies between the ATLAS solar-elongation floor and ~75°
   have both in-beam and out-of-beam epochs in every year: r_b =
   0.85 / 0.90 / 0.95 AU ↔ ε_b = 58° / 64° / 72° (r_E = 1).

The plan's §5.25 cell — the S2 skirt — is the sunward half of that
step. The A-side half is the same beam seen from the other side of
the Sun; including it is a geometric identity, not a design choice,
and it closes the ATLAS report's deferred "A 1.0 AU ledger addendum"
with the same pull. Per-side statistics are reported separately so
the skirt's own contribution is visible (§6).

## 2. Archive and roles **[ATLAS]**

- **ATLAS forced photometry, reduced mode** (`use_reduced`, tphot on
  the target images): the substrate for every star — targets and the
  null-ensemble controls alike — full history MJD 57227 → coverage
  start, o primary / c annotation, all-sky. Reduced mode because the
  ATLAS survey established that difference-mode channel A on a
  high-PM star is the PM dipole (k 10²–10⁷, no constraint); the
  learnings' remedy is *reduced-mode photometry at a per-epoch
  PM-propagated position*.
- **Proper motion is applied server-side**: the queue accepts
  `propermotion_ra`, `propermotion_dec` (mas/yr) and
  `radec_epoch_year` (OpenAPI snapshot; recon §4). One full-history
  task per star; the output RA/Dec columns (the forced position per
  epoch) verify the application against the Gaia-propagated track
  (dev gate G1, §7). Parallax (≤ 0.3″) is below the response
  threshold (R ≥ 0.98) and is a declared budget term.
- ASAS-SN Sky Patrol v2: **not used** (different instrument; controls
  must share the target's exposures). Sky Patrol v1: manual-only, as
  before, for enumerated follow-ups.

## 3. Eras, geometry, in-era scope **[SKIRT]**

| item | value |
|---|---|
| ATLAS era | MJD 57227 → coverage-start date (re-measured then; the probes end MJD 61252); δ < −50 targets 2022+ only |
| Geometry | per target, daily ε(t) and b_e(t) = r_E sin ε over the era (`scripts/skirt_geometry.py`, `results/skirt_geometry_v1.json`); yearly conjunction minima reproduce the frozen list's S2 `b_min` to a median 0.0003 AU (max 0.025 AU, daily-grid quantisation for β ≈ 0) |
| S2 events in era | 942 / 85 targets (one conjunction per target per year); A events the mirror set |
| ATLAS solar-elongation floor | **50° declared** (recon §2: two full-history probes at ecliptic field stars — min 40.6°, 1st percentile 52–54°; night efficiency 5–8 % at 45–55°, 18–21 % at 60–70°, 33–45 % beyond 100°; southern units 55–58°) |
| Rung ladder (D2) | r_b = 0.90, 0.95 AU searched (ε_b 64°, 72°); 0.85 AU (58°) not proposed (2–5 skirt nights/yr at the measured efficiency); 1.0 AU ledger-only |

**Expected in-beam nights per year** (geometry × measured
efficiency, recon §3): skirt side 4–8 (r_b 0.90) / 7–13 (0.95) for
|β| < 45°, 18 at 0.95 for |β| 64–66°; opposition side 45–58;
quadrature (out-of-beam) band 21–41. High-|β| targets never reach
b_e < r_b (|β| ≥ ε_b) or ride the edge (|β| within 2° of ε_b, where
membership flips with r_E): **ledger-only** — wolf-1069 and gj-293 at
both rungs, gj-13157 and gj-3112 at 0.90.

## 4. Target population **[SKIRT]**

Photometric scope (`results/target_scope_v1.json`): the blended cell
needs the star unsaturated on the ATLAS PSF. Screening estimate
o ≈ G − 0.22 + 0.13 (BP−RP) (four substrate-measured unsaturated
anchors, rms 0.13 mag); the **D2 rule of the ATLAS freeze decides on
the measured reduced-mode value** (excluded < 12.5, marginal
12.5–13.0, ok ≥ 13.0), exclusions reported per target. Estimated
scope: 59 targets excluded (o_est < 12.5), 7 bright non-Gaia stars
excluded (Sirius, α Cen, Procyon, Fomalhaut), WISE 0855 no optical
counterpart, van-maanen excluded on its measured o = 12.30,
**luhman16-a/b excluded** (no verified Gaia counterpart; two
comparable background stars within 25″ swept by the pair's 2.8″/yr
motion — recon §4); **18 targets** — marginal gj-1111, wolf-1069,
teegarden (measured 12.59), gj-915; ok gj-3512, gj-3306, gj-3112,
gj-293, gj-1087, gj-9193, gj-2012, gj-13157 (G 17.2 neighbour at
4.7″, 7 % constant dilution, flagged), gj-518, gj-12724, gj-11547,
gj-1276 (measured 15.19), gj-11068, eps-ind-b (o_est 18.5, near the
per-exposure limit — kept: absolute contrast sensitivity improves
for fainter stars). Searchable at r_b 0.90: 14; at 0.95: 16; the
rest ledger-only (§3).

**Sensitivity scaling (declared before search).** A blended step is
contrast-limited: the detectable fractional excess is the
systematics floor of monthly-mean reduced photometry, ~1 % for
o 13–15 (recon §3 precision), so the flux limit is ~5 mag below the
star (o 13.5 → ≈ 18.5 AB; gj-1276 15.2 → ≈ 20.2) — shallower than the
antipode channels' 19.9–20.5 but through a beam 81–90× the 0.1 AU
rung's area. Order of magnitude: P = F_line π r_b² ≈ 3–15 MW for a
1 % excess on an o 13.5–15 star (o-band line, 690 nm, 260 nm width),
≈ 0.1–0.5 MW on eps-ind-b — a sub-metre-class-aperture uplink cell
at MW power, the same class the heliospheric surveys deliver for the
sunward downlink. Depths are not promised past what the measured
coverage and control ensembles deliver.

## 5. Duty cycle and temporal models

- **d = 1 while geometry holds — primary**: the in-beam step recurs
  every year on both sides; the stack over years is the primary
  statistic.
- **Single-year step — secondary** (S_event analogue): max over years
  of the per-year step, for non-recurring transmissions.
- **Declared unconstrained**: the 1 AU rung's constant excess; sub-
  day pulses (the pulse cell belongs to the antipode surveys —
  on-star per-exposure outliers are stellar flares on M dwarfs);
  schedules avoiding Earth's in-beam months; Gaussian-profile beams
  (top-hat per programme; a Gaussian edge softens the step and lowers
  recovery, quantified only if a v2 adopts it).

## 6. Detection construction **[SKIRT]**

**Series.** Per star (target or control) and band: reduced-mode
full-history task with server-side PM → FAQ primary mask
(`atlas_api.faq_quality_mask`) → per-epoch response factor
R_i = exp(−d_i²/2s_i²) from the output forced position vs the
Gaia-propagated position (d_i ≈ 0 after PM; epochs with R_i < 0.9
dropped — a PM-application failure shows here, gate G1) →
**nightly-unit means** (one epoch per unit-night: median of the quad,
≥ 2 exposures) → fractional flux f = F/median(F).

**Systematics layer — airmass regression (D5).** Per star and band,
f = a + b (X − 1) fitted over all nightly epochs (X from the unit's
site and the epoch's MJD), residual r = f − fit. Rationale: the S2
skirt epochs sit at X 1.5–3 in twilight while quadrature/opposition
epochs sit near X 1–1.5; the exposure zero-point is set by field
stars of mixed colour, so a red target carries a colour-dependent
extinction residual of order 0.02–0.05 mag per airmass — larger than
the ~1 % signal — that colour-matched controls only partly share
(the reddest targets have no colour-matched field stars, recon §5).
The same layer runs on every control and inside every injection, so
whatever signal it absorbs is charged to completeness, not to the
false-alarm rate.

**Per-rung sets** (r_b): IN = nightly epochs with r_E sin ε < r_b
(skirt side ε < 90° **and** opposition side ε > 90°), OUT = the rest
(quadrature band). Per year y (the S2 event's conjunction year, from
the mid-year of the quadrature band to the next):
Δ_y = mean(r ∈ IN_y) − mean(r ∈ OUT_y),
v_y = s²(IN_y)/n_IN + s²(OUT_y)/n_OUT (in-year sample variances).

**Statistics per unit** (unit = target × rung × band; **D4**):

1. **S_sym** = Σ_y Δ_y/v_y / √(Σ_y 1/v_y) — recurrence-stacked
   symmetric step (primary).
2. **S_skirt** — the same construction with IN_y restricted to the
   S2 side (ε < 90°; the plan's cell on its own) and OUT_y to the
   S2-side quadrature band.
3. **S_year** = max_y Δ_y/√v_y (secondary; single-year step).
   S_opp (opposition side only) is computed and reported as an
   **annotation**, not a trial — the symmetry check for any
   exceedance (a beam gives equal Δ on both sides; airmass does not).

**Controls — the null ensemble (D6) [SKIRT].** Temporal pseudo-
windows are unavailable (the in-beam bands are months long against
an annual recurrence — the standing window ≈ season theorem, four
confirmations). Per target: **8 field stars** drawn (declared seed)
from the recon pool of 16 within 1.0° of the target (same ATLAS
exposures ⇒ identical ε/airmass/twilight history), matched in G
(±0.3, widened to ≤ ±1.5 where the pool is thin) and BP−RP (±0.4,
widened to ≤ ±4.5 for the reddest targets), Gaia DR3 non-variable,
RUWE < 1.4, PM < 0.1″/yr, isolated (no neighbour within 10″ brighter
than G + 3). Each control runs the identical chain and yields the
identical statistics; **T = max over valid controls; exceedance
S > max(T, 0)**; per-trial crossing probability 1/(n_valid + 1)
(mirror-gated validity: a control is valid for a statistic when its
series passes the same coverage gate the target's does, §7).

**Veto / adjudication ladder** for an exceedance: (1) S_opp vs
S_skirt symmetry; (2) the c-band series (a line in o is absent or
different in c; a colour-dependent systematic tracks BP−RP); (3)
per-year consistency (a recurring step vs one anomalous year — flare
seasons on active M dwarfs); (4) the epoch-level residuals against
the control ensemble's mean residual curve (an exposure-set
systematic shows in all nine stars); (5) SkyBoT / Gaia neighbour
census for a transient blend; (6) an enumerated Sky Patrol v1 manual
follow-up if (1)–(5) leave it retained-ambiguous.

## 7. Coverage gate, trials, budget **[SKIRT]**

A statistic of a unit is a **searched trial** when the target's
series has ≥ 3 years each with ≥ 5 IN nights and ≥ 5 OUT nights
(S_skirt: ≥ 3 skirt nights on the S2 side count as IN), **and** ≥ 4
controls pass the same gate; S_year needs ≥ 1 such year. Otherwise
constraint-only. Expected control crossings = Σ_trials 1/(n_valid+1),
tallied at the threshold freeze from the measured series (nothing
in-window is formed for that tally: the gate reads epoch counts
only). Trials per fully-covered target: 2 rungs × 3 statistics = 6
(o band); c annotation only.

**Dev gates** (run on the dev split before the confirmatory run):
G1 PM application — forced RA/Dec at the first and last epochs within
0.3″ of the Gaia-propagated position for every dev star, else the
fallback is per-interval tasks at PM-propagated positions
(≤ 0.5″ drift per interval); G2 airmass layer — the fitted slope b of
the controls is consistent with zero-mean scatter across the ensemble
(|median b| < 2 MAD) so the layer removes rather than injects
structure; G3 scale — inherited from the ATLAS D6 gate ((15000) CCD,
±0.1 mag) with the mode caveat recorded (reduced vs difference
photometry share tphot's calibration).

**Task budget.** One full-history reduced-mode task per star:
18 targets + 8 × 18 controls = **162 tasks**; measured full-history
reduced-mode execution 12–15 min → **~32–40 h of serial queue
drain**, resumable (`atlas_drain.py` sidecars). Ledger-only targets
(wolf-1069, gj-293) pulled without controls: 146 tasks — decision D7.

## 8. Flux scale, injections, completeness

Response-model injections into the real target series: a top-hat
step of flux f (o-band AB grid) added to every IN epoch of the rung
(both sides for S_sym/S_year, the skirt side for S_skirt) before the
airmass layer, then the frozen statistic against the unit's frozen
T; 100 draws per magnitude (the draw randomises nothing but the
per-epoch photon noise realisation — the step is deterministic — so
m90 is the faintest grid magnitude with S > max(T, 0) in ≥ 90 % of
noise draws), seed declared at threshold freeze. Power
P = F_line π r_b² (line through o, 690 nm, 260 nm effective width),
per rung. ±0.1 mag scale systematic inherited (G3).

## 9. Expected sensitivity, declared before search

Recon §2: reduced-mode per-exposure precision for a G 14 field star
1.6 % (fractional MAD, o; photon term 0.4 %), 1.8 % at ε < 70°. A
skirt band of 4–13 nights per year (r_b 0.90/0.95, floor 50°) over
≤ 11 years, plus the opposition band's ~50 nights/yr, gives
statistical precision ≪ 1 % on the stacked Δ; the systematics floor
of the control ensemble sets T. Declared expectation: **m90 ≈ o_star + 4 to
+5.5** per unit (1–0.6 % excess), i.e. o ≈ 17.5–19 for the
o 13–14 targets and ≈ 19–20.5 for gj-1276 / gj-11068 / luhman16 —
3–15 MW through the 0.90–0.95 AU beams. Nothing deeper is promised.

## 10. Pre-freeze data-contact declaration

Two full-history reduced-mode tasks were pulled on 2026-09-07 at
**non-target Gaia DR3 field stars** on the ecliptic (3698956337898276352
at RA 180.05, Dec +0.00 — 3.2° from ross-128, a saturated target
outside this survey's population; 2546036791796380672 at RA 0.23,
Dec +0.10), to measure the solar-elongation floor, cadence and
precision. No survey position was touched; no target-side or
elongation-locked quantity of any target was formed. No remedy is
needed. The control pools were selected from the Gaia catalogue only.

## 11. Freeze decisions — all adopted as recommended (user, 2026-09-07)

- **D1 — cell definition.** Search the sub-1 AU elongation-locked
  step on both sides of the Sun (S2 skirt + opposition mirror, the
  same beam by geometric identity), with per-side statistics; the
  1 AU rung is a coverage ledger (no temporal signature). Alternative:
  S2-side-only statistics (S_skirt alone) — thinner (4–13 nights/yr)
  and blind to the symmetry veto.
- **D2 — rung ladder.** r_b = 0.90 and 0.95 AU (edges 64°, 72°);
  0.85 AU not searched; 1.0 AU ledger. Alternative: 0.95 only (halves
  the trials).
- **D3 — population.** The 18 targets of §4; the ATLAS D2 saturation
  rule on the measured reduced-mode value at the cut stage;
  luhman16 excluded; gj-13157 kept with its dilution flag; ledger-only
  where |β| ≥ ε_b − 2°.
- **D4 — statistics and trials.** S_sym (primary), S_skirt, S_year;
  S_opp annotation; o searched, c annotation; trials = 3 per unit,
  2 rungs → ≤ 6 per target, expected crossings Σ 1/(n_valid + 1).
- **D5 — airmass regression layer** on every star (fit over all
  nightly epochs), inside injections too; nightly-unit means as
  epochs; FAQ mask primary; R ≥ 0.9 gate on the forced position.
- **D6 — controls.** Null ensemble of 8 field stars per target drawn
  by seed from the recon pool of 16; T = max valid control; validity
  = the target's coverage gate (≥ 3 years with ≥ 5 IN and ≥ 5 OUT
  nights; S_skirt ≥ 3 skirt nights); ≥ 4 valid controls for a
  searched trial.
- **D7 — task scope.** All 18 targets with controls (162 tasks,
  ~32–40 h drain) — the ledger-only pair included so the 1 AU ledger
  is measured, not modelled. Alternative: skip controls for the
  ledger-only pair (146 tasks).
- **D8 — split and seed.** Dev = teegarden (highest PM, marginal
  saturation — stresses G1 and D2) + gj-2012 (ordinary, colour-
  matched controls); all other targets confirmatory; seed 20260907
  (control draw and injection noise). Dev gates G1–G3 of §7 must
  pass before the confirmatory run; a G1 failure triggers the
  declared per-interval fallback, not a redesign.

## 12. Pre-confirmatory amendments (recorded before any confirmatory-target statistic was formed)

- **A1 — G2 wording (2026-09-07, at the dev stage).** §7 tested the
  airmass layer by "|median b| < 2 MAD across the ensemble", a
  zero-centred test that is wrong in expectation: the fitted slope
  carries the colour-dependent extinction term the layer exists to
  remove, so an ensemble of red controls legitimately has a non-zero
  median slope. G2 is replaced by the layer's actual purpose — it
  must not inject structure: MAD(r) ≤ MAD(f − median) for ≥ 75 % of
  the dev stars, and each dev target's slope within its control
  ensemble's range extended by 3 MAD (`scripts/dev_gates.py`).
- **G1 result (dev, 2026-09-07).** The queue applies
  `propermotion_ra` as μα cos δ: forced positions track the Gaia
  propagation to 0.014″ (gj-2012, 0.6″/yr) and 0.03″ (teegarden,
  5.1″/yr; span 55.6″ over the era vs 56.7″ expected); the raw-μα
  reading drifts to 0.4″ / 1.6″. No per-interval fallback needed.
- **A2 — G2 outcome (dev, 2026-09-07; `results/dev_gates_v1.json`,
  `results/dev_layer_test_v1.json`).** The A1 scatter-reduction form
  is insensitive (nightly MAD changes by < 1 % either way — the
  airmass span of most nights is too small for a linear term to move
  a robust scatter). The layer was therefore tested on its purpose:
  the null-ensemble floor with vs without it (RMS of the 8 controls'
  statistics, both dev fields, both rungs). Result: neutral —
  teegarden field 0.99 vs 1.06 (S_sym, 0.90 AU) and 1.57 vs 1.59
  (0.95); gj-2012 field S_skirt 0.64 vs 1.71 (better) against S_sym
  1.89 vs 0.72 (worse; the field's common −2.6 %/airmass slope,
  shared by target and controls, is re-expressed by the fit). The
  ensembles show coherent field-wide offsets (all 8 teegarden
  controls at S_skirt ≈ −1.9) that the max-rule threshold absorbs.
  What the layer uniquely removes is the target-specific colour term:
  teegarden's slope is +3.3 %/airmass against ≤ 0.7 % for its
  (bluer) controls, a +1.5 %-of-flux bias on the high-airmass skirt
  nights that no control shares — the false-positive mechanism D5 was
  frozen against. **Disposition: the frozen D5 layer stands
  unchanged; G2 is recorded as neutral-on-floor / pass-on-purpose.**
  No construction was altered on dev data.
- **A3 — parallax response correction (dev, 2026-09-07).** The dev
  search (12 trials, 8 exceedances vs 1.3 expected; all six teegarden
  trials, S_sym 3.6 against controls all ≤ −0.3) exposed a mechanism
  §2 had waved through as a "declared budget term": both dev targets
  sit above every one of their controls, and the one property every
  target has and no control shares is a parallax. The parallactic
  displacement peaks at quadrature — the out-of-beam band — so a
  fixed-position PSF fit loses flux there and nowhere else: an
  elongation-locked step of the signal's sign and size (teegarden
  π = 0.26″, forced-vs-apparent offset 0.26″ at ε 80–100° vs 0.07″
  at opposition). Amendment: the response factor now uses the
  apparent position (PM + parallax from the registry, Earth's
  barycentric position, no aberration — the WCS is ICRS) and the
  measured flux is corrected by 1/R_i; injections add f to the
  corrected flux (a transmitter at the star shares the star's
  offset). The response form for a fixed-position fit is not the
  peak-value factor exp(−d²/2σ²) of the ATLAS injection model but
  the cross-correlation exp(−d²/4σ²) for Gaussian PSFs; because the
  two differ by the full size of the effect (0.5 % vs 1.0 % at
  teegarden's quadrature), the exponent is **measured**, not assumed:
  six forced-photometry tasks on two quiet dev controls at Dec
  offsets 0.5″, 1.0″, 2.0″ (`scripts/response_calib.py`,
  `results/response_calib_v1.json`) fit R = exp(−d²/(k σ²)); the
  fitted k enters the correction. The task list gains
  `parallax_mas` (registry for targets; 0 for the field-star
  controls — G 12–15 field stars have π < 3 mas, R > 0.99999).
  With the analytic k = 4 the dev search gives 7 exceedances / 12
  (teegarden S_sym 2.6 / 2.5; gj-2012 S_sym 0.16 vs T −1.02).
- **A4 — colour-unmatched ensembles (pre-confirmatory rule, 2026-09-07).**
  Five targets (teegarden, gj-1111, gj-3512, gj-12724, gj-11547;
  BP−RP 3.7–5.2) have no field star within 1.0 mag of their colour in
  their ensemble (recon §4). A colour-dependent systematic (telluric
  water bands at 720/820 nm inside o, differential chromatic
  refraction, the calibrators' colour term) is elongation-locked
  through season and airmass and is not shared by bluer controls;
  the ladder's rule (2) colour-trend test cannot extrapolate over
  3 mag. Rule, fixed before the blind run: an exceedance on such a
  target is classed **`retained_ambiguous_colour_unmatched`** — not
  promotable, reported with the control colour-trend extrapolation
  as annotation — and the unit's depth is reported as
  systematics-limited. Targets with a control within 1.0 mag of
  their colour are adjudicated by the frozen ladder unchanged.
- **A3 result (dev, 2026-09-07 late).** Measured response
  (`results/response_calib_v1.json`, 6 offset tasks, 9,219 matched
  exposures): R = exp(−d²/(k σ²)) with **k = 2.42** in the small-offset
  regime (0.5″/1.0″; 2.47 over all offsets; per-star medians 2.1–3.6
  at 0.5″), between the two analytic forms; the median PSF σ on these
  fields is 1.80″, so teegarden's quadrature loss is 0.6–0.9 % (k
  3.6–2.1). The frozen correction uses k = 2.42; the star-to-star
  spread is a ±0.2 % systematic on Δ for π ≈ 0.26″ targets
  (teegarden, eps-ind-b, gj-1111), ≤ 0.05 % for π ≤ 0.15″. Dev search
  after A3: 7 exceedances / 12 trials — teegarden S_sym 2.0 / 2.0
  (Δ̄ +0.52 / +0.58 %, both sides positive: skirt +0.70 / +0.29 %,
  opposition +0.43 / +0.47 %; controls −0.05 to −0.45 %), S_year 4.0 /
  4.5 (the 2020-21 cycle, +4 %); gj-2012 S_sym 0.001 vs T −1.02
  (vetoed asymmetric). Adjudication: teegarden
  `retained_ambiguous_colour_unmatched` (A4; nearest control colour
  2.7 mag bluer; a +0.5 % symmetric annual excess on an M7 flare star
  with a +4 % season is stellar activity or a colour systematic before
  it is a transmitter), one S_skirt trial `non_promotable_single_cycle`.
  The construction is not changed further; the confirmatory run uses
  A3 (k = 2.42) and A4 as frozen here.
- **A5 — depth reporting for exceedance trials (recorded at the
  completeness stage, 2026-09-10, after the blind run; a reporting
  rule, not a search change).** The §8 injection adds the step to the
  target's real in-beam set and tests S > max(T, 0); for a trial
  whose real statistic already exceeds the threshold every injection
  "recovers" and m90 lands on the grid ceiling. Such trials (24 of
  81) carry no depth in the report; the 57 non-exceedance trials do.
