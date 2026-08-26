# GALEX crossings hypothesis freeze v1.0 (Pipeline B)

Drafted 2026-08-26 after the reachability recon
(`notes/galex_recon_2026-08-26.md`) and the era scoping
(`results/era_scope_v0.json`), before any further in-window data
contact (pre-freeze contact declared in §10 — one recon probe touched
in-window photons; remedy D1). **FROZEN 2026-08-26: the user approved
every §11 decision as recommended (D1–D8).** Eighth crossings survey; the
construction transfers the PTF freeze
(`surveys/ptf-crossings/hypotheses.md` v1.0) with the TESS freeze's
temporal-statistic layer, adapted to a photon-event substrate —
GALEX-specific substitutions marked **[GALEX]**.

Input: `crossings/universal_v1` (`xng-a09e2db7681d`, Earth-center,
1980→2028, both link directions, every b(t) minimum kept). GALEX is
~700 km LEO: Earth-center events valid under the standing
0.010 R☉ / ~4 min budget (the WISE/SPHEREx convention). Era:
**2003-06-07 → 2013-05-01 UTC** (the measured global span of the
photon DB `aspect` table). Times below are "GALEX time" =
UNIX − 315,964,800 s; DB columns are ms; photon stamps are 5 ms ticks
(verified).

## Amendment v1.1 (2026-08-26, at coverage; user approved)

**Aspect-flag gate: `flag = 0` → `flag % 2 == 0`.** The v1.0 gate
was written on the recon's single aspect-row observation (flag 0, a
2007 GI visit). The coverage run found the gj-1276 B 2010 visit —
and the 2009–2010 era at that field generally — carries aspect
flag 64 on every second; gPhoton's own pipeline
(`PhotonPipe.py` lines 534–536: photons get bad-aspect flag 12
unless the neighbouring aspect seconds satisfy `aspflags % 2 == 0`)
and the aspect table's `flagDiv2` column both establish the
mission/gPhoton convention: **only bit 0 marks a failed aspect
solution; even flag values are usable.** The recon's astrometric
check (0.63″ centroid on an MCAT source) ran on flag-64 seconds.
Statistics, thresholds, and every other gate untouched; the odd-flag
cut stays. Consequence: gj-1276 B 2010 is restored (1,637 usable
NUV s; its 3 odd-flag-81 seconds drop), no other unit changes. The
v1.0-gate coverage outputs are preserved as
`results/coverage_v1_flag0_superseded.*`. Dev-stage item: the formal
astrometric verification on flag-64 seconds (off-window).

## Amendment v1.2 (2026-08-26, at dev; user approved)

**S_period de-quantization jitter.** Dev step v found the H-test
grid corrupted by the 5 ms photon-tick quantization: every grid
frequency with a harmonic commensurate with the 200 Hz tick rate
accumulates coherent quantization power (null H ≈ 2,700–3,800 vs
the ~32 χ² expectation; P = 10 s injections unrecoverable). v1.2:
**S_period is computed on tick-jittered times** — U(0, 5 ms) added
to each photon time before phase folding, RNG seed 20260826 + a
per-series offset, applied identically to units, controls, and
injections. Validated on the off-window B pseudo-unit: null H
29–38, unit clean (31.2 vs T 37.8), P = 10 s recovery restored
(`results/jitter_validation_v0.json`). Grid, harmonic count,
photon gate, and the other two statistics untouched.

## 1. Observable channels

Channel A (uplink interception, `inbound`, blended with the star,
near solar opposition) and channel B (downlink pre-lens interception,
`outbound`, at the star antipode, z-parameterized parallax-reflex
track 0.36–6.5″/day). Sunward combinations out of scope (§5.11 owns
them). **[GALEX]** No Pipeline-A corridor infrastructure exists
(GALEX was never a Pipeline-A archive); positions come fresh from the
events table per unit.

## 2. Beam-radius ladder

The standard ladder with the era-scoped structural outcomes declared
at freeze:

- **B 1.2 R☉, B 2.5 R☉ (and A at the grazing radii):
  coverage-without-statistic.** The era scoping found **zero
  in-window seconds** on the ±0.3–0.4 d grazing windows at all 10
  grazing-family positions (the van-maanen deep-graze family has
  visits in 2004/2008, none in-window). These rungs enter the
  covered-window ledger as structurally uncovered.
- **B 0.1 AU and A 0.1 AU: the searched rungs.** Flat-chord windows
  t_ca ± √(r²−b²)/v⊥ ≈ ±5.8 d.
- **A 1.0 AU: out of scope** — the d = 1 wide-beam rung keeps its
  programme-wide deferred status (§5.7: v2-style null-ensemble
  redesign required before it is ever searched).

## 3. Wavelength **[GALEX]**

NUV (177–283 nm, λ_eff 231 nm) and FUV (134–179 nm, λ_eff 154 nm).
Hypothesis flux = monochromatic line converted via band effective
width. **The NUV band contains 266 nm — the frequency-quadrupled
Nd:YAG line — the programme's first band containing a harmonic of
the 1064 nm hypothesis family**; FUV contains no declared line and
carries the generic-leakage interpretation. 532/1064/1550 nm remain
outside and unconstrained. A monochromatic line lands in exactly one
band: single-band excess with a flat simultaneous other band is the
laser discriminator (the ZTF chromatic-anomaly rule), directly
measurable on the 4 units with both detectors live.

## 4. Duty cycle and temporal models **[GALEX: the 5 ms axis]**

The survey's cell is sub-exposure time structure inaccessible to
every prior archive (TESS floor: 1 FFI cadence = 200/600 s):

- **d = 1 persistent — primary.** Window-locked in-aperture rate
  excess over the visit (the chord is unresolved: visit ≪ window;
  this is the standard d = 1 statistic sampled at one point of the
  chord).
- **Burst/pulse — the new cell.** Sliding-boxcar count excess at
  frozen timescales **{0.05, 0.5, 5, 50 s}** — sensitive to single
  pulses and to any train's brightest pulse.
- **Coherent trains.** H-test (de Jager) on photon arrival phases
  over a frozen log period grid **P ∈ [20 ms, T_visit/3]**,
  oversampling 5× the independent-period count. Spacecraft orbital
  motion smears topocentric phase by up to ~40 ms across the 1.6 ks
  visit (~2.5 ms across 100 s visits): **no orbit correction is
  applied; the smear is carried by the injections** (worst-case
  linear drift model), so completeness curves state the short-period
  sensitivity loss honestly (C1 rule) rather than a correction claim.
- **Declared unconstrained:** pulse periods > T_visit/3 (window ≫
  visit), schedules avoiding Earth-crossing windows, and sub-5 ms
  structure (tick floor).

## 5. Detection constructions **[GALEX: photon-event substrate]**

No images. Both channels use **direct photon-event retrieval** from
the MAST gPhoton DB (`NUVPhotonsV`/`FUVPhotonsV`, flag = 0 photons;
usable-detector gate boresight ≤ 33′ from `fGetNearbyAspectEq`;
live time = distinct aspect seconds passing the gate).

- **Apertures.** r_ap = 8″ (NUV PSF FWHM ≈ 5.3″; clustering probe
  RMS 4.2″). **B (track):** one aperture per z-grid point
  (550/1000/2500/5500/10000 AU — the frozen programme grid) at the
  apparent-relay position computed at the visit epoch (within-visit
  track motion ≤ 0.12″, static per visit; the z family spans a
  ~2–30″ locus segment at the scoped 3.9–5.3 d offsets from t_ca);
  overlapping apertures deduplicated; statistic = max over the
  deduplicated set, mirrored exactly in the controls. **A (blended):**
  one aperture at the propagated star position; the star's own
  photons are in-aperture by construction — the statistic is excess
  over the star+background rate, not over sky.
- **Background/null.** Per-unit empirical null from **pseudo-position
  controls**: apertures (full z-family construction for B) at ≥ 8
  positions on the same visit, same boresight-distance annulus,
  ≥ 30″ from any MCAT source above the confusion gate and from the
  real locus. Window ≫ visit makes temporal pseudo-windows
  structurally unavailable **within** a visit; **off-window visits at
  the same position** are the temporal control family (rich: e.g. 115
  visits at the gj-1276 antipode) — both families frozen here,
  thresholds set at the threshold freeze.
- **Calibration.** Rate → flux via in-visit MCAT-star regression
  through the identical aperture/live-time chain (the programme's
  star-ZP rule), verified to ≤ 0.2 mag before any depth is quoted
  (the SPHEREx gate). The recon's naive ZP conversion disagreed with
  a raw rate by ×2.6 — live-time/aperture/dead-time corrections are
  a dev-stage gate, not assumed. Global dead time and local
  nonlinearity carried as declared budgets; relative statistics
  (locus vs same-visit controls) cancel them to first order.
- **Veto ladder:** (1) stellar-flare veto (A units — all three
  scoped A stars are M dwarfs): flare morphology test
  (fast-rise-exponential-decay template vs window-persistent or
  strictly periodic structure) **and** the two-band discriminator —
  flares brighten FUV+NUV together (FUV-loud), a line is
  single-band; FUV is live on all three 2007 A units. (2) SkyBoT
  known-object census per exceedance (B antipodes sit in the
  opposition asteroid stream; a crossing MBA transits the 8″
  aperture in ~minutes — timescale-distinct from the persistent
  statistic but exactly the burst statistic's band: every burst
  exceedance is census-checked). (3) hotspot/artifact check —
  detector-fixed (xi/eta) clustering vs sky-fixed; known NUV hotspot
  list. (4) recurrence at any second covered window (gj-1276 B has
  two: 2007 + 2010).

## 6. Saturation / nonlinearity rule (channel A) **[GALEX]**

NUV local nonlinearity reaches ~10 % near 100–200 ct/s
(m_NUV ≈ 14.5–15). Rule on the MCAT NUV mag at the cut stage:

    excluded  m_NUV < 14.5
    marginal  14.5 ≤ m_NUV < 15.5
    ok        m_NUV ≥ 15.5

Expected consequence, stated for the record: all three scoped A
stars (wolf-359, gj-1276, ross-128) are UV-faint M dwarfs
(NUV ≳ 19) — expected ok-class; the rule is frozen so the cut stage
decides, not the expectation.

## 7. Quality gates **[GALEX, from recon]**

- **Photon:** flag = 0 (census: 100 % of 8,621 off-window photons at
  the scoped field; the rule also guards any flagged class the DB
  serves elsewhere).
- **Aspect/live-time:** aspect seconds with boresight distance ≤ 33′
  (the usable-detector cut; beyond sits the rim — the gj-908 A case,
  §9) and aspect `flag` = 0 (schema pinned at recon; distribution
  checked at dev).
- **Strict (exceedance re-runs):** boresight ≤ 25′, and the
  detector-fixed clustering check from the veto ladder.
- **Depth accounting:** per-unit sensitivity comes from the
  injection chain only (C1). Scoping expectation, declared before
  search: background ≈ 0.2 ct/s in the 8″ aperture → persistent-rate
  5σ ≈ NUV 23 (1.6 ks unit) / ≈ 21.6 (100 s units); single-pulse
  fluence ~4–5 photons at the 1 s timescale. These are expectations,
  not constraints.

## 8. Injections, completeness, positive control

Photon-level injections: synthetic arrival-time sets (persistent
rate, boxcar pulses on the frozen timescales, periodic trains on the
frozen P grid **including the orbital phase-smear model**) with
PSF-scattered positions, superposed on the real photon lists,
through the identical retrieval→gate→statistic chain; calibrated per
unit (C1). Positive controls: (i) MCAT-star rate recovery ≤ 0.2 mag
through the full chain (per unit); (ii) if the dev-stage off-window
flare census at the A-unit stars yields a flare, the burst statistic
must recover it and the flare veto must classify it — recorded
either way.

## 9. Era scope (from the frozen universal list + recon aspect scan)

`results/era_scope_v0.json`, strict ≤ 33′ cut. The searched-unit
candidates (all 0.1 AU rung):

| unit | visit (UTC) | offset from t_ca | NUV s | FUV s |
|---|---|---|---|---|
| gj-1276 B 2010-03-03 | 2010-03-06 21:57 | +3.9 d | **1,637** | 0 |
| gj-1276 B 2007-03-03 | 2007-02-26 11:19 | −4.8 d | 109 | 109 |
| gj-1276 A 2007-09-05 | 2007-09-01 01:52 | −4.1 d | 92 | 92 |
| wolf-359 A 2007-03-03 | 2007-02-26 11:21 | −5.3 d | 97 | 97 |
| ross-128 A 2007-03-17 | 2007-03-22 15:40 | +4.7 d | 110 | 110 |

**gj-908 A 2007-09-21 is excluded**: its only in-window visit sits
at 36.3′ minimum boresight distance — outside the 33′ usable cut
(detector rim; no calibrated rim response exists to admit it) —
recorded as coverage-without-statistic, rim-limited. Zero-coverage
ledger entries: both grazing rungs everywhere (§2); the teegarden
antipode (no GALEX coverage at all); van-maanen both channels
(visits, none in-window); ross-154 both channels; all other
unit-events without in-window seconds. The coverage stage re-derives
this table with fresh snapshot-disciplined queries and exact
per-event, per-epoch positions; §9 is scoping, not the frozen
coverage record.

## 10. Pre-freeze data-contact declaration

Probe files live only in the session scratchpad and are re-pulled
fresh at the coverage stage. Contact at survey positions before this
freeze, in full:

1. **Aspect metadata (timestamps, band, boresight distance)** at all
   14 positions, including in-window second counts (the §9 scan) —
   the class every coverage stage touches pre-threshold-freeze.
2. **In-window photon contact at the gj-1276 B 2010 unit** (the
   reachability verification): one aperture-scale box **total count**
   (1,695 photons / 30″×30″ / full visit) and **10 photon rows**
   (times, positions, flags; 6 examined) — retrieved before the
   off-window discipline was articulated for this archive. No
   background comparison, live-time normalization, control, or
   per-time-bin structure was formed; the count is consistent with
   a priori sky background. This is nonetheless rate-class contact
   at a confirmatory locus → **remedy D1**.
3. **Off-window photon contact**: box counts and row samples at the
   van-maanen antipode (2004-04-08 visit — 3 d outside the widest
   window; 2010 era-probe pulls), and the §7 flag census + §5
   clustering/rate check at the gj-1276 antipode and an MCAT field
   star (2009-02-16, 15 d off-window).
4. **Catalog rows**: 6 `visitphotoobjall` rows in a 0.01° box at the
   ross-128 star (unattributed visits; one row lies 6″ from the star
   — possibly the star itself, possibly from the in-window visit;
   magnitudes were read, no window association was made) and
   `photoobjall` cones at the gj-1276 antipode (off-window
   calibration probes). Era-integrated `photoobjall` means are
   annotation-class; the `visitphotoobjall` pull is declared as
   potentially window-adjacent (remedy option in D1).

## 11. Freeze decisions — all adopted as recommended (user, 2026-08-26)

- **D1 — pre-freeze contact remedies.** (a) gj-1276 B 2010: the
  contacted quantity class (aperture rate) is demoted — **S_rate on
  this unit is `forced_dev`** (constraint-only, labeled, excluded
  from the confirmatory family); the burst and period statistics,
  whose structure was never examined (10 of 1,695 timestamps seen,
  for tick verification), **remain blind confirmatory**. (b)
  ross-128 A: retain fully blind — the catalog pull made no window
  association and the A statistic is rate-vs-own-star, which no
  static magnitude prefigures; declared here so the adjudication of
  any ross-128 A exceedance must cite it. (c) All scratchpad probe
  files stay out of the repo; coverage re-pulls fresh.
- **D2 — era and input**: era 2003-06-07 → 2013-05-01 (measured
  aspect span); `crossings/universal_v1` Earth-center with the
  standing LEO budget (0.010 R☉, ~4 min); no GALEX observer list.
- **D3 — substrate**: direct photon-event retrieval per §5 (flag 0,
  boresight ≤ 33′, live time from gated aspect seconds); MCAT-star
  in-visit rate calibration primary, ≤ 0.2 mag gate; no image
  products.
- **D4 — units and statistics**: unit = (target, channel-rung,
  band). NUV primary on all 5 units; FUV searched with the identical
  constructions on the 4 FUV-live units (a line lands in one band —
  FUV is a search band, not only a veto). Three frozen statistics
  per unit — S_rate (d = 1 persistent), S_burst (4 timescales),
  S_period (H-test, frozen grid) — trials tallied at the threshold
  freeze (≤ 27 = 15 NUV + 12 FUV, minus the D1a demotion), FWER
  α = 0.05 across the family as in every survey.
- **D5 — gates**: §7 (photon flag 0; ≤ 33′; strict ≤ 25′ +
  detector-fixed clustering check; depths from injections only).
- **D6 — channel-A nonlinearity rule**: §6 (E = 14.5/15.5 marginal
  band; MCAT-based, cut-stage verified).
- **D7 — controls**: ≥ 8 pseudo-position controls per unit on the
  same visit (full z-family construction mirrored for B) primary;
  off-window same-position visits as the temporal control family;
  both ensembles sized at the threshold freeze.
- **D8 — dev/confirmatory split**: **dev = pseudo-units only** —
  off-window visits at the scoped positions (incl. the gj-1276
  antipode long visits, the wolf-359/ross-128 star off-window
  visits for the flare census, and the D1a forced_dev rate lane);
  **confirmatory = all 5 §9 units, blind** (with S_rate on gj-1276
  B 2010 carried as forced_dev constraint-only). No searched unit
  is spent on machinery development — the off-window archive is
  rich enough to exercise every construction including the A-channel
  flare chain.
