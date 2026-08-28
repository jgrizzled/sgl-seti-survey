# DASCH crossings hypothesis freeze v1.0 (Pipeline B)

Drafted 2026-08-26 after the reachability recon
(`notes/dasch_recon_2026-08-26.md`), the crossings backward
extension (`crossings/universal_1885_v1`), and the flat-chord era
scoping (`results/flat_chord_scan_v0.json`). **FROZEN 2026-08-26:
the user approved every §11 decision as recommended (D1–D8).**
Eighth crossings survey; the construction adapts the PS1/PTF
freeze line to the programme's first **catalogue-level-first
substrate** (per-exposure calibrated archive photometry, not
survey-side pixel photometry), each departure marked **[DASCH]**.

Input: `crossings/universal_1885_v1` (`xng-e2d1063af9d0`,
Earth-center, 1885→1993, both link directions, every b(t) minimum
kept; shared-era rows verified against `universal_v1` to 23 s /
3.4×10⁻⁷ R☉). DASCH era: **~1885 → ~1990** (measured at the scoped
positions: first exposures 1885, last calibrated exposures 1989;
per-unit era ends measured at the coverage stage). Harvard plus the
Arequipa/Bloemfontein southern stations: the topocentric offset is
bounded by the geocentric radius = 0.0092 R☉ — inside the standing
0.010 R☉ observer budget used programme-wide. The pre-1941 events
carry sglseti's `long_propagation_span` annotation; the measured
old-epoch budget (recon: ephemeris ≤ 1×10⁻⁵ R☉ vs DE440S,
astrometric sensitivity ≤ 2×10⁻⁴ R☉, |Δt_ca| ≲ 5 s) is carried as
a declared term and the annotation is not a validity gate.

## 0. Scope **[DASCH]**

The five grazing-family targets only — van-maanen, wolf-359,
teegarden, gj-1276, ross-128 (the complete set with any era event
at b ≤ 2.5 R☉; census in the recon addendum). The other 83
registry endpoints' 0.1 AU rungs (~10⁴ windows) are declared out of
scope for this survey — a future queue item, not a
coverage-without-statistic ledger entry.

## 1. Observable channels

Channel A (uplink interception, `inbound`, at the star, near solar
opposition) and channel B (downlink pre-lens interception,
`outbound`, at the star antipode). Sunward combinations stay with
§5.11. **[DASCH]** Channel B searches the archive's own calibrated
detections in the ~13′ platephot subregion around the per-event
relay locus (`relay_icrs_*` from the events table); channel A
searches the target star's DASCH lightcurve (calibrated mags +
per-epoch `limiting_mag_local`). No Pipeline-A corridor
infrastructure exists for DASCH; the adapter is REST-level
(queryexps / querycat / lightcurve / platephot / cutout) with
snapshot discipline.

## 2. Beam-radius ladder

The ZTF/PS1/PTF ladder verbatim: B 1.2 R☉ / B 2.5 R☉ / B 0.1 AU /
A 0.1 AU / A 1.0 AU, flat-chord window t_ca ± √(r²−b²)/v⊥ per
event.

- **A 1.0 AU: constraint-only** (standing window≈season theorem).
- **[DASCH]** All grazing B rungs are *searchable* here — the
  century cadence covers even the 4–7 h chords (§9): B 1.2 R☉
  holds 12/7/4 covered windows at van-maanen/gj-1276/wolf-359.
  teegarden B 1.2 R☉ (38 era events, 0 covered) enters the ledger
  as coverage-without-statistic; ross-128 has no era events
  ≤ 1.2 R☉ (rung not defined).

## 3. Wavelength **[DASCH]**

Photographic blue calibrated to **APASS DR8 B** (the DR7 science
refcat; ATLAS-refcat2 g has documented false long-term trends and
is annotation/astrometry-only — recon). Line hypotheses: 532 nm
(doubled Nd:YAG) is outside B (~390–490 nm) — **declared
unconstrained**, unlike every CCD-era optical survey; the
**355 nm tripled-Nd:YAG line is conditionally in-band**
(blue-sensitive emulsions typically respond ~340–500 nm; the UV
edge is emulsion/series-dependent) — constraints at 355 nm carry a
declared per-series caveat, resolved via colorterm/emulsion
metadata only if an exceedance demands it; 1064/1550 nm out of
band. Hypothesis flux = monochromatic line converted via band
effective width; broadband-leakage interpretation primary.

## 4. Duty cycle

d = 1 while geometry holds (primary; channel-B geometry enforces
window-locked transience); d = 0.1 as a flat ×d flux scaling
(secondary, no re-search). **[DASCH]** No pulse statistic and no
sub-exposure cell: the median exposure is 60 min and integrates
the chord — pulse periods below the exposure time are declared
unconstrained. The survey's new cell is the **century recurrence
stack** (S_stack, §5): a schedule transmitting during
Earth-crossing windows across ~10–60 covered annual recurrences
per unit.

## 5. Detection constructions **[DASCH: catalogue-level substrate]**

DASCH photometry is plate-level SExtractor + multi-stage
calibration served per exposure; there is no survey-side pixel
photometry in the primary chain (cutouts are resampled 1.44″/px
products — vetting substrate only).

- **A (on-star):** the star's lightcurve (querycat → lightcurve,
  APASS refcat). Per covered window: window-locked brightening =
  calibrated in-window mag(s) vs the star's off-window baseline
  (era-local, per-series where sample permits). Where the star is
  fainter than the plate limit (teegarden, B ≈ 18: limits-only
  regime), the statistic is a window-locked *detection* at the
  star position above `limiting_mag_local` — the archive's
  non-detection rows carry the per-epoch limit. Source-splitting
  (known issue) handled by merging catalog entries within a frozen
  radius at the coverage stage.
- **B (locus track):** per covered window, platephot on every
  in-window exposure at the event's relay locus; candidate = a
  detection within the frozen match radius of the locus
  (astrometric tolerance + locus motion over the window; radius
  set at the threshold freeze), **uncatalogued** (no refcat match)
  or refcat-matched-but-flux-anomalous, passing the fatal-bit
  screen (§7).
- **Statistics (per unit = target × channel × rung):**
  **S_event** = the most significant single-window excess /
  locus-consistent detection; **S_stack** = the recurrence
  statistic across all covered windows of the unit (stacked excess
  for A; window-locked detection count vs control expectation for
  B). Both primary; trials tallied at the threshold freeze.
- **Null and thresholds:** 8 spatial ring-control loci per unit
  (offset positions at matched plate-center distance inside the
  same platephot subregions — same plates, same windows, free with
  the primary pull) give the primary error rate; off-window
  temporal draws at the true locus (depth-matched per decade) are
  computed as annotation (the standing theorem: temporal controls
  are structurally weaker for annual phase-locked windows).
- **Veto ladder (B exceedances):** (1) fatal AFLAGS screen (§7);
  (2) morphology vs plate PSF (fwhm/ellipticity from the platephot
  row — a 60-min exposure trails an ordinary mover ≳ 15″/hr into
  an elongated track while the relay track moves < 1″); (3) known
  Solar-System object census at the exceedance epoch (Horizons/
  SkyBoT cover the plate era); (4) same-plate multiple-exposure
  ghost check (expnum family + PlateQualityFlags); (5) cutout
  visual/defect inspection (documented DASCH defect classes);
  (6) recurrence on the recomputed locus at the unit's other
  covered windows — with ~10–60 windows this is the decisive
  discriminator. Only vetoes with measured selection functions
  reject (v2 rule); the rest annotate.

## 6. Brightness rule (channel A) **[DASCH]**

Plate saturation thresholds are series-dependent and bright
(B ≲ 9–10 on the common series); every scoped target is fainter
(van Maanen measured B ≈ 12.6; the M dwarfs B ≈ 13–18). Rule: a
unit's in-window epochs with `TOO_BRIGHT` or `SATURATED`
(AFLAGS/BFLAGS) are unusable rows; no per-target exclusion class
is expected — stated for the record and verified at dev from the
stars' own lightcurves, ± amendment if a series proves saturated
at B ~ 13.

## 7. Quality gates and timing **[DASCH]**

- **Listing-level (coverage):** an exposure is *covered* iff it has
  an APASS photometric calibration (`limMagApass` present in
  queryexps ⇒ astrometric + photometric solution) and
  `wcssource = imwcs`/`catalog` (never `logbook`). No series
  exclusions at coverage; the per-series census is recorded.
- **Fatal AFLAGS template (search-stage, candidate detections):**
  `SUSPECTED_DEFECT | PICKERING_WEDGE | MULT_EXP_UNMATCHED |
  MULT_EXP_BLEND | REJECTED_BLEND | UNCERTAIN_DATE` (bit source:
  daschlab `photometry.py`, snapshot
  `snapshots/daschlab_photometry_2026-08-26.py`). Blend-class bits
  (`CASE_B/C/BC_BLEND`, `SXT_BLEND`) and `BAD_PLATE_QUALITY` are
  **strict-template** additions — they can lose in-scope sources,
  so their cost is measured (§8) and charged, never silent.
  `UNCERTAIN_CATALOG_MAG` is never fatal (it flags the refcat
  entry, not the plate datum — set on 100 % of van Maanen rows).
- **Timing gate:** an exposure enters a rung's searched set only if
  its timing uncertainty ≤ 0.5 × the window half-width — from
  `time_accuracy_days` where a lightcurve/platephot row provides
  it, else the date-only heuristic (timestamp exactly 00:00:00 ⇒
  1.0 d). Date-only exposures are therefore searchable at the
  0.1 AU rungs (±5.8 d) but excluded from grazing chords
  (coverage-without-statistic, listed per window). Window overlap
  is **interval-based**: [t_start, t_start + exptime] (exptime in
  minutes — recon) vs the flat-chord window; in-window weight =
  overlap fraction.
- **Depth accounting:** per-exposure `limMagApass` (coverage) and
  per-row `limiting_mag_local` (search) are the depth ladder,
  validated by the §8 efficiency curves; no header-derived depths
  exist to mistrust.

## 8. Completeness and positive control **[DASCH]**

No pixel injections (the substrate is the archive's own catalogue
chain). Completeness is **measured, per covered window**, from
refcat field-star recovery: platephot detections vs the querycat
truth set in the same ~13′ subregion give a per-plate detection-
efficiency-vs-magnitude curve (Wilson intervals); the fatal- and
strict-template costs are measured on the same population; the
per-window m90 comes from that curve, cross-checked against
`limiting_mag_local`. The C1-rule analogue: **no depth statement
from `limMag*` alone** — only the measured efficiency curve
qualifies a constraint. Positive control (dev stage, criteria
frozen now): one documented high-amplitude variable (eruptive or
large-amplitude periodic, B range crossing plate depths, within or
near a scoped subregion) recovered end-to-end through the
identical querycat → lightcurve → statistic chain, with its
in-archive behavior reproduced (mags to ≤ 0.3 mag; the DASCH RY
Cnc tutorial object is the fallback if no in-field candidate
exists).

## 9. Era scope (flat-chord scan, calibrated exposures, recon epoch lists)

`results/flat_chord_scan_v0.json` (definitive flat-chord windows,
interval overlap, per-event loci at the probed positions). Covered
windows / in-window exposures per searched-unit candidate:

| Unit | B 1.2 R☉ | B 2.5 R☉ | B 0.1 AU | A 0.1 AU |
|---|---|---|---|---|
| van-maanen | 12 / 17 | 21 / 36 | 62 / 238 | 63 / 219 |
| gj-1276 | 7 / 12 | 16 / 31 | 60 / 215 | 58 / 176 |
| wolf-359 | 4 / 6 | 15 / 23 | 59 / 177 | 60 / 217 |
| teegarden | 0 (ledger) | 9 / 15 | 56 / 194 | 61 / 236 |
| ross-128 | — (no era events) | 10 / 12 | 57 / 208 | 56 / 197 |

**18 searched-unit candidates** (3 × B 1.2 R☉, 5 × B 2.5 R☉,
5 × B 0.1 AU, 5 × A 0.1 AU); median chord half-widths 4–7 h
(1.2 R☉), 11–16 h (2.5 R☉), ±5.8 d (0.1 AU). The coverage stage
re-derives this table with fresh snapshot-disciplined queries, the
exact per-event loci, the timing gate, and the §7 listing gate;
§9 numbers are scoping, not the frozen coverage record.

## 10. Pre-freeze data-contact declaration

Probe files live only in the session scratchpad (summaries under
`results/` are derived counts, not photometry); everything is
re-pulled fresh at the coverage stage. Contact at survey positions
before this freeze, in full:

1. **Metadata** (queryexps exposure lists: timestamps, exptime,
   limiting mags, series, WCS source) at all 10 scoped
   star/antipode positions — the standard pre-freeze class;
   in-window timestamps were counted (§9) but no flux at a locus
   was formed into a window-locked statistic.
2. **Photometry, van-maanen A:** the star's **full APASS
   lightcurve** (3,544 rows, 1890–1989, necessarily including
   in-window epochs) was pulled and examined at recon (global mag/
   flag/timing distributions only; no window-locked quantity
   formed). This is substantive on-star contact for a
   would-be-confirmatory unit → **remedy: `forced_dev` for
   van-maanen A** (D1, the ATLAS wolf-359 precedent).
3. **Pixel/subregion, off-window only:** cutouts and platephot at
   the van-maanen *star* position on plates mf10951
   (1926-10-13, 9.6 d from the nearest t_ca — outside every rung
   window incl. ±5.8 d) and a26884 (1949-07-06, 92 d off);
   platephot repeat at a 0.06° offset centre (subregion-size
   probe). No antipode photometry anywhere.
4. **querycat** at the van-maanen star (2 rows, catalog metadata)
   — static-catalog class, non-photometric.

## 11. Freeze decisions — recommended, awaiting approval

- **D1 — data contact & remedy:** accept §10; **van-maanen
  A 0.1 AU is `forced_dev`** (its only searchable A rung); all
  other units blind-eligible. Scratchpad probes never enter the
  repo; coverage re-pulls under snapshot discipline.
- **D2 — input & era:** `crossings/universal_1885_v1`
  (Earth-center; 0.0092 R☉ topocentric bound inside the standing
  0.010 R☉ budget; `long_propagation_span` carried as annotation
  with the recon-measured budget, not a gate). Scope = the
  5-target grazing family (§0).
- **D3 — substrate:** catalogue-level primary (lightcurve for A,
  platephot for B, APASS refcat); resampled cutouts are vetting
  substrate only; no survey-side pixel photometry in the primary
  chain; ATLAS-refcat annotation-only.
- **D4 — units & statistics:** unit = (target, channel-rung),
  band B only; S_event + S_stack both primary per searched unit
  (≤ 36 trials over 18 candidate units, tallied exactly at the
  threshold freeze); d = 0.1 as ×d scaling; no pulse statistic
  (sub-exposure cell declared unconstrained).
- **D5 — quality & timing:** §7 verbatim — coverage gate
  (`limMagApass` present + non-logbook WCS), fatal template
  (defect/wedge/mult-exp/rejected-blend/uncertain-date), strict
  template (+ blend-class + plate-quality bits, cost measured),
  timing gate at 0.5 × half-width with the 00:00:00 date-only
  heuristic, interval-based window overlap.
- **D6 — brightness rule:** §6 — row-level `TOO_BRIGHT`/
  `SATURATED` exclusion only; no target-level exclusion class
  expected; dev verification with amendment trigger.
- **D7 — controls:** 8 spatial ring-control loci per unit from the
  same platephot subregions (primary error rate); off-window
  temporal draws as annotation; per-cell normalisation and
  heavy-tail exclusion decided at the threshold freeze on the
  ring ensemble (v2 rule).
- **D8 — dev/confirmatory split:** **dev** = van-maanen A 0.1 AU
  (forced, D1), teegarden A 0.1 AU (limits-only-regime
  machinery), wolf-359 B 2.5 R☉ (B-chain machinery on a
  mid-sized unit); **confirmatory (15 units)** = all three
  B 1.2 R☉ units (van-maanen, gj-1276, wolf-359), the remaining
  B 2.5 R☉ units (van-maanen, gj-1276, teegarden, ross-128), all
  five B 0.1 AU units, and A 0.1 AU for gj-1276, wolf-359,
  ross-128. The headline century families (van-maanen B,
  gj-1276 B) stay blind end-to-end. Zero-coverage rungs enter the
  ledger without a split.

## 12. Amendment v1.1 (2026-08-26, dev stage; before any confirmatory contact)

**Trigger (the forced_dev unit doing its job).** The frozen §5
channel-A construction assumed the DR7 APASS lightcurve tracks the
star. Dev found the APASS refcat carries **no proper motion for two
of the five A stars** (van-maanen, wolf-359: dummy entries, pm = 0,
stdmag 99.9), so DR7's matching cannot follow a 1–5″/yr star across
the century — the van-maanen APASS "lightcurve" holds 9 spurious
detections at B ≈ 10.6 vs the true 1,552-detection ATLAS lightcurve
at B ≈ 12.6. Probes (querycat metadata, all five stars, both
refcats) established: APASS PM valid for gj-1276 and ross-128;
ATLAS PM valid for van-maanen (matches Gaia to 0.1%); ATLAS entry
for wolf-359 is source-split with wrong PM.

**Amended rule (applies identically to controls and injections):**

1. **Lightcurve routing.** A detected-regime A unit uses the refcat
   entry whose catalogued PM matches the registry µ within
   max(10 %, 50 mas/yr) (vector) and has a valid stdmag; preference
   APASS, else ATLAS with the documented false-long-term-trend
   caveat annotated. If neither refcat tracks the star, the unit
   falls to the **limits-only regime** (the hit statistic needs no
   refcat entry). Routing outcome: van-maanen → ATLAS
   (`ref 9717808673`); ross-128 → APASS (`411474441004816`);
   gj-1276 → limits-only (APASS PM valid but stdmag 18.0 is below
   every plate limit); wolf-359, teegarden → limits-only.
2. **Era-local baseline.** Detected-regime baselines are ±5 yr
   around each window (padded windows excluded), gate ≥ 20 rows —
   this defuses ATLAS's slow false trends and century emulsion
   drift alike. **S_stack** becomes the Stouffer combination
   Σz_w/√N_w over windows with an in-window row (the global-pooled
   variant is dropped before ever being used on data).
3. **Entry matching** uses minimum separation to the Gaia-propagated
   track over 1885–2020 (the frozen 10″ radius unchanged) — the
   refcat `pos_epoch` column is a placeholder (2000.0 everywhere)
   and APASS positions sit at the ~2012 observation epoch.

Statistics untouched for channel B and the limits-only regime; the
searched-unit list, trials tally (36), and dev/confirmatory split
are unchanged.

## 13. Amendment v1.2 (2026-08-26, dev stage; before any confirmatory contact)

**Trigger (B-chain positive control).** (7) Iris (V 9.9, ~7″ trail)
on plate me00943 (1911-03-01, 30 min, lim 13.2) — a real moving
point-source transient run through the frozen hit chain — was found
on the plate as an uncatalogued mag-10.3 detection but **rejected by
the frozen fatal template**: DASCH's classifier stamps it
`SUSPECTED_DEFECT` (+ `BAD_PLATE_QUALITY`). A single-plate
uncatalogued transient is precisely what the defect classifier is
trained to eat — and precisely the survey's signal class. Under the
v2 calibrated-veto rule the flag cannot remain a silent fatal bit:
its selection function against in-scope sources is now measured
(1/1 real transients killed).

**Amended rule (identical for locus and ring controls):**

1. **Fatal template** for candidate hits becomes `PICKERING_WEDGE |
   MULT_EXP_UNMATCHED | MULT_EXP_BLEND | REJECTED_BLEND |
   UNCERTAIN_DATE` (mask 152045824). `SUSPECTED_DEFECT` is demoted
   to an **annotation**: a hit carrying it enters the statistic and,
   if it contributes to an exceedance, the §5 veto ladder's cutout
   inspection (rung 5) decides defect-vs-source with the plate
   image — the discrimination moves from a silent bit to a
   documented adjudication step.
2. The Iris control's 22.6″ offset from the Horizons prediction is
   dominated by plate-timing error against the asteroid's 23″/hr
   motion (offset anti-parallel to the motion vector) — a term that
   scales with apparent rate and is ≤ 0.3″ for relay loci
   (≤ 6.5″/day): **r_match = 10″ stands** for the survey; positive
   controls with fast movers are matched against the trail segment
   extended by ±1 h × rate (the measured plate-timing class).

Dev hit units are re-scored from the cached snapshots under the
amended template before the confirmatory freeze-chain runs.

## 14. Amendment v1.3 (2026-08-26, dev stage; before any confirmatory contact)

**Trigger (teegarden dev exceedance adjudication).** Under v1.2 the
teegarden A unit produced the dev programme's one exceedance
(S_stack 2 vs T_stack 1; 1 vs 0.67 expected). Adjudication: the
1945 hit (mag 16.83, 2.2″, plate limit 16.9) is **the quiescent
star itself** marginally extracted at the plate limit — the
off-window same-locus measurement (12 deep off-window plates)
found the identical signature at the same flux (mag 16.89, 1.4″,
1/12 plates); the 1941 hit (mag 14.64 at limit 14.7) is
**defect-class by cutout inspection** (a ≤ 3.4σ diffuse background
bump at the position, no point source; the 21σ field star on the
same cutout shows what mag 14.6 looks like). Structural finding:
**for a limits-only unit the star supplies its own background of
marginal at-limit extractions, invisible to spatial ring controls.**

**Amended rule (limits-only regime only; before the confirmatory
units wolf-359 A and gj-1276 A are touched):**

1. **Quiescent-level measurement (off-window, pre-search).** Per
   limits-only unit, the star's DASCH quiescent level and marginal-
   detection rate are measured from uncatalogued detections within
   r_match of the propagated star position on **off-window** usable
   plates (outside every padded window), stratified by plate limit.
   This is declared off-window data contact of the PTF-§10 class.
2. **Quiescent-star gate (calibrated veto).** A locus hit that is
   not ≥ 1.0 mag brighter than the measured quiescent level is
   classified `quiescent_star` — annotated, excluded from S. The
   veto can only lose a relay emitting at the star's own quiescent
   flux; its selection function is the measured off-window rate,
   and its completeness cost is charged to the constraint. Ring
   hits are never gated (no star there) — the asymmetry is
   conservative (raises T only).
3. Where no off-window uncatalogued detection exists at the star
   across the deep off-window plates, the measured level is the
   stratified plate-limit distribution and the gate reduces to
   requiring the hit ≥ 1.0 mag brighter than the plate limit of
   its own exposure.

Dev is re-scored under the full v1.0+v1.1+v1.2+v1.3 chain; the
adjudication record above stands regardless of the re-score.
