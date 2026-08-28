# DASCH beam-crossings survey (Pipeline B) — report

Eighth Pipeline-B (§3.5) survey and the programme's first reach into
the photographic century: a search of the DASCH DR7 scanned Harvard
plates (~1885–1990) for optical emission during Earth's crossings of
hypothesized Sun–star relay beam axes — the first coverage of any
kind before the 1980 universal-list start, and the first survey on a
**catalogue-level substrate** (the archive's own calibrated
per-exposure photometry; no survey-side pixel photometry in the
primary chain). Executed end-to-end 2026-08-26. **No candidates: 30
blind confirmatory trials, 2 exceedances vs 3.33 expected control
crossings, both adjudicated defect-class in the plate pixels.**

## 1. Construction and provenance chain

1. **Crossings backward extension** `crossings/universal_1885_v1`
   (`xng-e2d1063af9d0`): every b(t) minimum for 88 endpoints × 2
   link directions, Earth-center, 1885→1993 — built as the §5.8
   item-5 prerequisite. Regression against `universal_v1` over the
   shared 1980–92 era: all 4,104 events match 1:1, max |Δt_ca| 23 s,
   max |Δb| 3.4×10⁻⁷ R☉. The §7 old-epoch model budget was measured
   at recon: astropy-builtin vs DE440S heliocentric Earth ≤ 6 km
   (1×10⁻⁵ R☉) at 1885; grazing-family astrometric sensitivity (5
   km/s RV + 0.1 mas/yr µ perturbations through rigorous space
   motion) ≤ 2×10⁻⁴ R☉ and ≲ 5 s in t_ca. The pre-1941
   `long_propagation_span` annotations are conservative bookkeeping
   against the declared 75-yr linear bound, not physical limits.
   **Census finding: van-maanen, gj-1276 and wolf-359 are
   photosphere-grazing (b ≲ 1.2 R☉) on essentially every annual
   crossing of the DASCH century** (van-maanen b 0.59–0.69 R☉
   throughout — the modern deep-graze family is the tail of a
   secular deepening); teegarden B holds 38 grazing events;
   ≤ 2.5 R☉ adds ross-128.
2. **Recon** (`surveys/dasch-crossings/notes/dasch_recon_2026-08-26.md`):
   DR7 rides the anonymous Starglass REST API (queryexps /
   querycat / lightcurve / platephot / cutout / mosaic_package);
   ~9,000–15,000 exposures overlap each grazing-family position,
   ~2,300–3,100 with APASS photometric calibration, 150–340 at
   B ≥ 15; Menzel gap 1953–68 at every position; APASS DR8 B is the
   science refcat (ATLAS-refcat2 g has documented false long-term
   trends).
3. **Hypothesis freeze v1.0** (`hypotheses.md`, D1–D8 approved
   2026-08-26): scope = the 5-target grazing family; ladder
   B 1.2 R☉ / B 2.5 R☉ / B 0.1 AU / A 0.1 AU (+ A 1.0 AU
   constraint-only); catalogue-level substrate (lightcurve for
   detected-regime A, platephot for B and limits-only A; resampled
   1.44″/px cutouts vetting-only); 532 nm outside the blue plates'
   band (declared unconstrained — unlike every CCD-era survey),
   355 nm conditionally in-band, broadband leakage primary; no
   pulse statistic (60-min median exposures integrate the chord);
   the century **recurrence stack** is the survey's new cell;
   van-maanen A `forced_dev` after recon lightcurve contact (D1).
4. **Coverage** (`results/coverage_v1_windows.ecsv`, fresh
   snapshot-disciplined queryexps pulls; listing gate = APASS
   calibration + non-logbook WCS; interval-based window overlap —
   `exptime` is in minutes; timing gate σ_t ≤ 0.5 × half-width with
   the 00:00:00 date-only heuristic): **18 searched units**.
   Headline rows: van-maanen B 1.2 R☉ 12 covered windows / 17
   exposures (median chord half-width 7.1 h), gj-1276 B 1.2 R☉
   7/12, wolf-359 B 1.2 R☉ 4/6; the 0.1 AU rungs hold 56–63
   covered windows per unit (~177–240 exposures each) — an
   order-of-magnitude larger grazing recurrence sample than all
   prior surveys combined. Ledger (coverage-without-statistic):
   teegarden B 1.2 R☉ (38 events, 0 covered), ross-128 B 1.2 R☉
   (no era events).
5. **Threshold freeze v1.0** (`thresholds.md`,
   `configs/threshold_freeze_v1.json` sha256:83dae35e…, bound to
   hypotheses + coverage hashes): B hit statistic (uncatalogued
   platephot detections within r_match = 10″ of the z-family track,
   z ∈ {550, 1000, 2500, 5500, 10000} AU) with 8 same-pull ring
   controls (±60/90/120″); A detected-regime robust-z lightcurve
   statistic with 8 temporal pseudo-windows (±23–97 d); S_event +
   S_stack per unit (36 trials; dev 6 / confirmatory 30); exceedance
   = S > max(T, 0), 1/9 budget.
6. **Dev stage** (`results/dev_v1.md`): 0 exceedances / 6 trials
   under the final chain; both positive controls passed — RY Cnc
   (A-chain: 1,008 usable rows, 85 faint excursions 75 %
   phase-locked at the eclipse phase vs 30 % uniform null, depth
   0.79 mag, scatter 0.18 mag) and **(7) Iris 1911** (B-chain:
   recovered as an uncatalogued mag-10.3 detection 0.2″ from the
   timing-extended Horizons trail, correctly trail-elongated — the
   morphology veto's positive exemplar). Three dev findings became
   amendments, all before any confirmatory contact:
   - **v1.1** — the DR7 APASS refcat has **no proper motion for
     van-maanen or wolf-359** (dummy entries): DASCH matching cannot
     track a 1–5″/yr star, so the APASS "lightcurve" of van Maanen
     holds 9 spurious detections vs the true 1,552-detection ATLAS
     lightcurve. PM-match routing rule (APASS-preferred, ATLAS with
     the false-trend caveat, else limits-only), era-local ±5 yr
     baselines, Stouffer S_stack.
   - **v1.2** — `SUSPECTED_DEFECT` demoted from fatal to
     annotation: DASCH's defect classifier stamps real single-plate
     moving transients (the Iris control was killed by it, 1/1) —
     exactly the survey's signal class. Defect discrimination moves
     to the cutout-inspection rung.
   - **v1.3** — limits-only quiescent-star gate: a sub-limit star
     supplies its own background of marginal at-limit extractions,
     invisible to spatial ring controls (teegarden dev exceedance:
     the 1945 "hit" was the quiescent star at its measured
     off-window rate — 2/12 deep plates, median B 17.11 vs expected
     ≈ 17.3; the 1941 hit was defect-class by cutout). Gate: hit
     ≥ 1.0 mag brighter than the measured off-window quiescent
     level; ring hits never gated (conservative).

## 2. Results

**Blind confirmatory: 0 candidates.** 15 units / 30 trials under
the frozen v1.0+v1.1+v1.2+v1.3 chain: **2 exceedances vs 3.33
expected control crossings**, both in wolf-359 A 0.1 AU, both
adjudicated (`results/adjudication_wolf359A_v1.json`):

- The unit's quiescent measurement worked exactly as designed:
  12/12 deep off-window plates detect the PM-less-refcat star as an
  uncatalogued source, measured quiescent B = 15.66 (15.45–16.07) —
  the gate passed only two of the unit's four in-window locus hits.
- **1943 / ac37753** (hit B 14.58, 0.35 mag above the plate limit,
  7.4″ from the star, fwhm 10 px, ellip 0.53): the plate pixels at
  the position show a ≤ 3.8σ diffuse bump; field stars on the same
  cutout reach 19.6σ. Vetoed defect-class.
- **1982 / dny00425** (hit B 13.02, 0.5 mag above the limit, 7.1″,
  ellip 0.40): pixels show ≤ 4.2σ diffuse structure vs 16.3σ field
  stars; the star itself (B 15.66) is correctly invisible on this
  lim-13.5 plate. Vetoed defect-class.
- Supporting: genuine star extractions in the unit sit at
  1.4–2.4″ (the hits at ~7″); no recurrence across the 56 covered
  windows; the cutout-inspection selection function is anchored by
  the dev exemplars (real transient → compact 15–20σ; defect →
  ≤ 4σ diffuse). The d = 1 persistent hypothesis is independently
  refuted by the unit's own 56-window coverage regardless of the
  hits' nature.

Annotations (below threshold, recorded): van-maanen B 0.1 AU holds
one locus hit (1905, B 14.05, 5.7″, S_stack 1 = T_stack 1 — not an
exceedance); teegarden and wolf-359 B rings caught single stray
detections (T = 1), locus clean; ross-128 A (the one detected-regime
confirmatory lightcurve unit, APASS-routed, era-local baselines)
S_event 4.06 vs T 7.09 — the photographic A-chain nulls are
heavy-tailed (defect-class rows inside control windows; dev saw
z = 9.65 in a van-maanen control) and the max-of-8-controls
threshold absorbs them at depth cost, the standing v2 heavy-tail
behaviour.

Dev results (dev units, not blind): wolf-359 B 2.5 R☉ clean null
(S = T = 0); teegarden A 0.1 AU S = T = 1 (no exceedance;
adjudication above); van-maanen A 0.1 AU (forced_dev) S_event 1.93
vs T 9.65, S_stack −1.09 vs T 1.20.

## 3. Constraints

Completeness is measured per covered window from refcat field-star
recovery in the same platephot subregion (Wilson intervals; the C1
analogue — `limMag*` columns alone never qualify a constraint):
`results/completeness_v1.json`. Median / best per-window m90 per
unit and the derived power statements — isotropic-equivalent relay
power through the rung cone, P = F(m90) · π r² (distance-independent
for the relay at z ≫ 1 AU; Johnson-B Δν = 1.42×10¹⁴ Hz, AB−Vega
−0.09):

| Unit | Covered win / exp | m90 (lim<13 / 13–15 / ≥15) | Power at deepest resolved stratum |
|---|---|---|---|
| gj-1276/A_0.1AU | 58 / 176 | — / 12.5 / 13.5 | ≳ 15.7 MW |
| gj-1276/B_0.1AU | 60 / 215 | 10.5 / 12.5 / 13.5 | ≳ 15.7 MW |
| gj-1276/B_1.2Rs | 7 / 12 | 10.5* / 12.5* / 13.5 | ≳ 49 kW |
| gj-1276/B_2.5Rs | 16 / 31 | 10.5 / 12.5* / 13.5 | ≳ 213 kW |
| ross-128/B_0.1AU | 57 / 208 | 8.5 / 8.5 / 8.5 | ≳ 1.6 GW† |
| ross-128/B_2.5Rs | 10 / 12 | 8.5* / 8.5* / 8.5* | ≳ 21 MW† |
| teegarden/A_0.1AU | 61 / 240 | — / 12.5 / 13.5 | ≳ 15.7 MW |
| teegarden/B_0.1AU | 56 / 194 | 10.5 / 12.5 / 13.5 | ≳ 15.7 MW |
| teegarden/B_2.5Rs | 9 / 15 | 10.5* / 12.5 / 13.5* | ≳ 213 kW |
| van-maanen/B_0.1AU | 62 / 238 | 10.5 / 12.5 / 14.5 | ≳ 6.3 MW |
| van-maanen/B_1.2Rs | 12 / 17 | 10.5 / 12.5* / 14.5* | **≳ 20 kW** |
| van-maanen/B_2.5Rs | 21 / 36 | 10.5 / 12.5 / 14.5 | ≳ 85 kW |
| wolf-359/A_0.1AU | 60 / 217 | — / 12.5 / 13.5 | ≳ 15.7 MW |
| wolf-359/B_0.1AU | 59 / 177 | — / 12.5 / 14.5 | ≳ 6.3 MW |
| wolf-359/B_1.2Rs | 4 / 6 | — / 12.5 / 14.5* | **≳ 20 kW** |
| wolf-359/B_2.5Rs | 15 / 23 | — / 12.5 / 14.5* | ≳ 85 kW |

Notes. `*` = field-curve fallback: the stratum was too thin in the
unit's own windows (< 10 pooled truth stars per bin), so the same
target-field's 0.1 AU-unit pooled curve — the identical plate
population and subregion — supplies it. `†` = the ross-128 antipode
field's measured recovery never reaches 90 % beyond B 8.5 in any
stratum (an anomalously poor-recovery field; the raw efficiencies
sit at 60–85 % over B 9–14, so the field is far from blind — the
90 %-completeness convention is simply not met; flagged for a
follow-up astrometric/source-splitting diagnosis). The deepest-
stratum power statement applies to the subset of covered windows
holding lim ≥ 15 exposures; the shallow-stratum bulk of each unit
carries the corresponding brighter limit (e.g. van-maanen B 1.2 R☉:
m90 10.5 → ≳ 780 kW over the lim < 13 majority of its windows).
**Measured, not assumed:** 90 % recovery sits **1.5–2.5 mag above
the archive's `limMag` columns** — the C1 lesson at catalogue level;
every number above comes from the pooled field-star recovery
curves (`completeness_pooled_v1.json`), never from `limMag`.
Template costs on the same population
(`template_costs_v1.json`, 8,714 matched rows): fatal 1.2 %,
strict 22 %, `SUSPECTED_DEFECT` prevalence on genuine stars 13.9 %.
Lightcurve units (van-maanen A forced-dev, ross-128 A): the
detected-regime statistic bounds window-locked brightenings of the
star itself — excursions ≳ 0.2–0.3 mag (≈ 25 % of the stellar B
flux) sustained over a window are excluded at the unit thresholds;
no platephot depth applies.

These are the programme's **first power constraints of any kind
before 1980**: tens-of-kW-class through the photosphere-grazing
cones on the deep-plate windows (van-maanen and wolf-359 B 1.2 R☉
≳ 20 kW), sub-MW across the shallow-plate century bulk reaching
back to the 1890s, and MW-class through the 0.1 AU cones — across
~10–60 independent annual recurrences per unit; the century
recurrence-stack cell (S_stack) closes clean on every searched
unit. Channel-A flux limits carry the
remote-uplink interpretation (10-m aperture at the star's distance,
η = 0.5): at B = 14, P_tx ≳ 0.4 MW (wolf-359), ~1 MW (teegarden,
van-maanen, ross-128), ~5 MW (gj-1276). The 532/1064/1550 nm line
hypotheses are declared unconstrained (out of band); 355 nm
statements carry the per-series emulsion-response caveat
(hypotheses §3).

Structural cells left open: teegarden B 1.2 R☉ and ross-128
B 1.2 R☉ (ledger rows); pulse periods below the exposure time;
duty-cycle schedules avoiding Earth-crossing windows; d = 0.1
statements are the frozen ×d flux scaling of the same limits.

## 4. Lessons (fed to `notes/learnings.md`)

1. **A catalogue-level substrate inverts the veto problem.** The
   archive's own quality flags are tuned to *suppress* single-plate
   transients (`SUSPECTED_DEFECT` killed the real-asteroid control
   1/1): on such a substrate, archive quality bits may gate the
   statistic only after their selection function against in-scope
   sources is measured — otherwise they are annotations feeding
   pixel-level adjudication.
2. **High-PM stars break refcat-keyed archives silently.** DASCH's
   APASS refcat carries pm = 0 dummy entries for two of five
   grazing-family stars; the resulting "lightcurves" are spurious
   without any error signal. The PM-match routing test (catalogued
   PM vs registry µ) is cheap and should be a standing check for any
   catalogue-level archive.
3. **A sub-limit star is its own background.** Marginal at-limit
   extractions of the quiescent star appear in-window and
   off-window alike and are invisible to spatial ring controls —
   the off-window same-locus rate measurement (v1.3) is the
   calibrated veto, and it doubles as a free DASCH-photometry
   validation (teegarden measured B 17.11 vs expected ≈ 17.3;
   wolf-359 15.66).
4. **At-limit + non-PSF + ~7″ offsets is the plate-defect
   signature.** All four adjudicated hits across dev + confirmatory
   sat ≤ 0.5 mag above their plate limits with fwhm 2–4× the field
   PSF; cutout pixels settled every one in minutes. Cheap, decisive,
   and anchored by the Iris positive exemplar on the same footing.
5. **The 1/9 budget behaved again**: 2 exceedances vs 3.33 expected
   (confirmatory), 1 vs 0.67 (dev, closed by amendment) — fifth
   consecutive survey on budget.

## 5. Ledger and follow-ups

- Covered-window ledger rows (for the §5.7 refresh): 18 searched
  units with their per-window states in `coverage_v1_windows.ecsv`;
  teegarden B 1.2 R☉ and ross-128 B 1.2 R☉ as
  coverage-without-statistic; date-only exposures excluded at
  grazing rungs are enumerated in the coverage table.
- The pre-1941 windows carry the `long_propagation_span` annotation
  with the measured budget (recon §"model-accuracy") — any future
  per-event follow-up can invoke `sglseti.crossing_uncertainty`.
- No retained-ambiguous rows: both confirmatory exceedances closed
  defect-class with pixel evidence.
- Follow-up: diagnose the ross-128 antipode field's anomalously poor
  platephot recovery (60–85 % over B 9–14 in every stratum —
  astrometric solution or source-splitting suspected); its power
  rows are honest but far shallower than the field's plates justify.
- Standing maintenance: fold the DASCH rows into the unified ledger
  at the next crossing-list refresh; the 83 out-of-scope targets'
  0.1 AU rungs (~10⁴ windows) remain a declared future queue item.
