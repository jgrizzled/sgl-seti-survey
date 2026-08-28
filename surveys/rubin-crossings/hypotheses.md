# Rubin DP2 crossings hypothesis freeze v1.0 (Pipeline B)

Drafted and **FROZEN 2026-08-26**, after the Rubin reachability recon
(`surveys/rubin/notes/rubin_recon_2026-08-26.md`) and its geometry
intersect (`surveys/rubin/results/recon_intersect_v0.json`), before
any further data contact at survey positions (pre-freeze contact
declared in §10). Adopted under the session directive to run the
crossings survey; decisions D1–D8 recorded in §11 — any user
amendment follows the standing amendment process (new version,
declared before the affected stage). Ninth crossings survey; the
construction transfers the GALEX freeze
(`surveys/galex-crossings/hypotheses.md` v1.0–v1.2) adapted to a
**catalog-level substrate** — the programme's first survey whose
detection layer is another pipeline's difference-image catalog rather
than pixels we process ourselves. Rubin-specific substitutions marked
**[RUBIN]**.

Input: `crossings/universal_v1` (`xng-a09e2db7681d`, Earth-center,
1980→2028, both link directions, every b(t) minimum kept). Rubin is
ground-based (Cerro Pachón): the geocenter→site offset at z ≥ 550 AU
is ≤ ~16 mas in angle and ≤ R⊕ = 0.009 R☉ in impact parameter — inside
the standing 0.010 R☉ Earth-center budget (the DECam observer
validation measured the site effect at 7–16 mas). Era: **MJD
60790.117 → 61047.354** (2025-03-26 → 2025-12-08 UTC, the measured
`dp2.Visit` span; the release notes say "to January 2026" — the
coverage stage re-measures the span and the larger value governs).
Processing pinned: DP2 (Early), Rubin Science Pipelines **v30**,
processing run **DM-55060**, TAP schema `dp2` at `data.lsst.cloud`.
If the archive republishes DP2 under a different processing run, that
is a new survey version, not a silent re-query.

## Amendment v1.1 (2026-08-26, at dev)

**Control locus-avoidance: ≥ 10″ → ≥ 2.5″ (2.5 × r_assoc) per
position.** The v1.0 D7 rule required every control position ≥ 10″
from the real locus segment; the dev stage found this geometrically
impossible for the inner-z pattern positions (the z = 10,000 AU
offset from the axis point is ~2″ at the unit epoch — no rotation
about the axis point can move it 10″ from itself), leaving several
pseudo-units — and, by the same geometry, the confirmatory unit —
with no valid control patterns. The avoidance rule's purpose is
non-overlap of the 1″ association apertures; **2.5″ guarantees
non-overlap with margin**. The +5° resolution rule, the 8-pattern
count, the statistic, and every other gate are untouched. Patterns
still unresolvable after a full turn are recorded invalid (census).
The v1.0-rule dev outputs are preserved as
`results/dev_v1_rule10as_superseded.json`.

## 1. Observable channels

Channel A (uplink interception, `inbound`, blended with the star,
near solar opposition) and channel B (downlink pre-lens interception,
`outbound`, at the star antipode, z-parameterized parallax-reflex
track). Sunward combinations out of scope (§5.11 owns them).
**[RUBIN]** The Pipeline A corridor machinery is deliberately not
used (that survey is tabled, plan §6); positions come fresh from the
events table per unit.

## 2. Beam-radius ladder

The standard ladder with the recon-scoped structural outcomes
declared at freeze (the coverage stage re-derives all of them with
the exact per-event windows, side-of-axis and validity cuts, and the
detector-footprint gate — §9 scoping is not the coverage record):

- **B 1.2 R☉, B 2.5 R☉ (and A at the grazing radii):
  coverage-without-statistic expected.** The recon found zero
  in-window boresight-cone visits on the ±0.35 d grazing windows
  (5 + 6 in-era events). These rungs enter the covered-window ledger.
- **B 0.1 AU and A 0.1 AU: the searched rungs.** Flat-chord windows
  t_ca ± √(r²−b²)/v⊥. Scoped in-window content: ross-128 B
  (b = 1.86 R☉ — a grazing-family event searched through the 0.1 AU
  rung window, sampling b(t_visit) = √(b_min² + (v⊥Δt)²)) and
  ross-154 A (expected lost to the §6 saturation rule).
- **A 1.0 AU: out of scope** — the d = 1 wide-beam rung keeps its
  programme-wide deferred status (§5.7: v2-style null-ensemble
  redesign required before it is ever searched). The recon's 36
  in-window targets on this rung are recorded in the ledger as
  geometry-only.

## 3. Wavelength **[RUBIN]**

LSST ugrizy (SDSS-like; r ≈ 552–691 nm, i ≈ 691–818 nm). The scoped
searched unit is r-band. No declared line of the 1064 nm hypothesis
family falls in r or i (532 nm sits in g — no in-window g visit was
scoped); detections carry the broadband reflected-light /
broadband-beacon interpretation, and any constraint is quoted for the
covered band only. Colour discrimination is available only where a
DiaObject anchor with multi-band ForcedSourceOnDiaObject epochs
exists — an annotation, not a frozen statistic.

## 4. Duty cycle and temporal models **[RUBIN]**

One 30 s visit per searched unit (scoped): the accessible cell is
**d = 1 persistent-during-window** (and any duty cycle covering the
visit — a 30 s exposure inside a ±5.8 d window samples duty cycles
≥ ~6 × 10⁻⁵ with probability ≈ duty cycle; stated in the report, not
searched as separate trials). Declared unconstrained: pulse structure
shorter than the exposure (TESS/GALEX own those cells), schedules
avoiding the covered 30 s, and periods > the visit span. No temporal
statistic is frozen — the substrate delivers one epoch per unit.

## 5. Detection constructions **[RUBIN: catalog-level substrate]**

No images. The detection layer is the DP2 difference-image pipeline
itself; our statistic is built on **DiaSource association at
predicted positions**, with all queries snapshot-disciplined (SQL +
raw response under `runs/rubin-crossings/`).

- **Predicted positions.** **B (track):** one position per z-grid
  point (550/1000/2500/5500/10000 AU — the frozen programme grid) at
  the apparent-relay position computed at the visit epoch (the GALEX
  `relay_apparent` construction: Sun − z·û_star viewed from Earth,
  astropy barycentric ephemerides), deduplicated at the association
  scale; statistic = max over the deduplicated set, mirrored exactly
  in the controls. **A (blended):** the propagated star position at
  the visit epoch (annotation only if the §6 rule excludes the unit).
- **Association.** r_assoc = **1.0″** (≈ 2× the astrometric scatter
  expected at the faint end; LSST single-epoch astrometry ≲ 50 mas at
  SNR ≥ 10, seeing ~1.2″ — the generous radius costs background
  association rate, which the controls carry). A DiaSource associates
  if within r_assoc of any deduplicated predicted position, in the
  in-window visit.
- **Statistic (single, frozen).** **S_det** = max over associated
  gate-passing DiaSources of SNR = `psfFlux`/`psfFluxErr`; S_det = 0
  if no association survives. Threshold = max over the 8
  pseudo-position controls (§7 of the threshold freeze); exceedance
  S > max(T, 0) — the standing crossings convention.
- **Real/bogus handling.** `reliability` is **not** given an absolute
  cut (its distribution at the unit cone was aggregate-touched
  pre-freeze, §10): it enters only through the veto ladder
  (annotation at adjudication) and through gates frozen on dev
  off-window fields *away from the unit cone* if dev shows the
  flag-gates alone leave a pathological association-rate null.
  Any such gate becomes a numbered amendment before confirmatory.
- **Static-sky exclusion.** A predicted position within **2″** of a
  `dp2.Object` (deep-coadd static source) with any-band
  `psfMag < magLim(visit band) + 0.5` is per-position excluded
  (confusion/subtraction-residual guard — the flux-consistent-static
  veto class); if every z-position of a unit-visit is excluded the
  unit is `not_constrainable` for that visit. Mirrored in controls.
- **Veto ladder (frozen order):** (1) known-object census — SsSource/
  SSObject association in DP2 plus a SkyBoT cone check at the visit
  epoch (the antipode sits in the opposition asteroid stream);
  (2) static/residual check — deep-coadd Object and
  ForcedSourceOnDiaObject history at the association position (a
  source persistent off-window is not a crossing candidate);
  (3) artifact flags — the §5-gate pixelFlags plus `reliability` and
  trail fits as annotation; (4) recurrence — any second covered
  window (none scoped in-era for ross-128 B; recurrence then falls to
  future data releases, recorded retained-ambiguous per the standing
  rule: promotion never on in-window evidence alone).

## 6. Saturation / bright-star rule (channel A) **[RUBIN]**

LSSTCam saturates near **r ≈ 16 (± ~0.5, seeing-dependent) in 30 s**.
Rule, applied at the cut stage from catalog annotation (Gaia DR3 G as
proxy, snapshot-disciplined):

    excluded  G < 16.5
    marginal  16.5 ≤ G < 17.5   (admit only if the dev-stage
                                 bright-star DiaSource census shows a
                                 usable association null)
    ok        G ≥ 17.5

Expected consequence, stated for the record: ross-154 (V ≈ 10.4,
G ≈ 9.6) → **excluded**; the A 0.1 AU rung then closes as
coverage-without-statistic (saturation-limited). The rule is frozen
so the cut stage decides, not the expectation.

## 7. Quality gates **[RUBIN, from recon]**

- **Visit/detector:** the predicted position must fall on an active
  detector of the visit — point-in-quadrilateral against the
  `dp2.VisitDetector` corner coordinates (`llcra…urcdec`), with that
  detector's `magLim` non-null. Boresight cones are discovery only.
- **DiaSource:** exclude detections with any of
  `pixelFlags_bad`, `pixelFlags_saturatedCenter`, `pixelFlags_crCenter`,
  `pixelFlags_edge`, `pixelFlags_nodataCenter`,
  `pixelFlags_interpolatedCenter`, `pixelFlags_suspectCenter`,
  `pixelFlags_streakCenter`, `psfFlux_flag`, `centroid_flag`; and
  **exclude archive injections**: `pixelFlags_injected` or
  `pixelFlags_injectedCenter` or `pixelFlags_injected_template` or
  `pixelFlags_injected_templateCenter`.
- **Depth accounting:** per-unit sensitivity statements come from the
  measured archive-injection detection efficiency where recoverable
  (§8), else from `magLim` **explicitly labeled not
  injection-calibrated** (C1 rule: no exclusion claim on an
  uncalibrated threshold — the WISE W3/W4 precedent).

## 8. Injections, completeness, positive control **[RUBIN]**

We cannot inject into pixels we do not have. Frozen completeness
route, in order: (i) **archive-side injections** — DP2 flags injected
sources (`pixelFlags_injected*`); if the injected population near the
unit fields (same band, comparable magLim) is recoverable with
sufficient statistics, the detection efficiency vs magnitude measured
from it calibrates the unit constraint (measured, not assumed —
whether a truth table or only the flagged detections are available
decides how far this goes; dev-stage gate); (ii) otherwise the
constraint is a `magLim`-referenced threshold statement labeled not
injection-calibrated. **Positive control:** known solar-system
objects — SsSource-linked DiaSources crossing the control fields (and
the unit detector where available) must associate and pass the gates
when their predicted MPC positions are run through the identical
chain (the DECam/PS1 asteroid-control pattern at catalog level);
control fails → confirmatory does not run.

## 9. Era scope (recon, boresight-cone only — superseded by coverage)

From `recon_intersect_v0.json` (1.5° boresight cones, no footprint or
validity/side cuts):

| candidate unit | t_ca (MJD) | b_min | in-window visits |
|---|---|---|---|
| ross-128 B 0.1 AU | 60937.94 | 1.86 R☉ (0.0086 AU) | 1 (r) |
| ross-154 A 0.1 AU | 60859.41 | 3.41 R☉ (0.0158 AU) | 9 (i r y z) — §6-excluded expected |

Grazing rungs: 0 in-window visits (5 + 6 in-era events). A 1.0 AU:
36 in-window targets, out of scope (§2). The ross-128 antipode field
has ~45 boresight-cone visits total → ~44 off-window: the dev-stage
pseudo-unit family.

## 10. Pre-freeze data-contact declaration

All recon probes are recorded in `surveys/rubin/results/recon_*.json`
(queries included). Contact at survey positions before this freeze,
in full — all at the **ross-128 antipode** (the B unit locus):

1. **Aggregate DiaSource counts** in a 6′-radius cone (era-wide,
   includes the in-window visit in aggregate): per-band counts, mean
   `reliability`, min/max MJD. No per-row position, flux, time, or
   score was read. 6′ = 360× the association radius; the searched
   statistic (per-visit SNR at five ~arcsec-scale positions) is not
   prefigured by band-level cone counts, but the **reliability means
   were seen** → the §5 rule that `reliability` gets no absolute
   frozen cut, and remedy D1.
2. **Aggregate counts**: DiaObjects (436) and Objects (6,786) in the
   same cone; 5 Object rows (id + coordinates only) in a 36″ cone;
   one sample Object's ForcedSource **per-band epoch counts** (no
   fluxes, no epochs read).
3. **Deep-coadd cutouts** (2′, r-band incl. mask/variance/PSF) at the
   antipode — era-integrated static sky, annotation-class (the GALEX
   `photoobjall` precedent); no per-visit quantity.
4. **Global aggregates** everywhere (visit table, per-band magLim
   averages) — coverage-class.

No per-visit, per-position, or per-row quantity at any predicted
locus position was formed. ross-154 and every other event position:
visit metadata only.

## 11. Freeze decisions (D1–D8) — adopted 2026-08-26

- **D1 — pre-freeze contact remedies.** (a) ross-128 B **remains
  blind confirmatory**: the contacted classes (band-level cone
  aggregates, static coadd) sit far above the statistic's resolution;
  this declaration must be cited at the adjudication of any ross-128
  exceedance. (b) Because reliability means at the unit cone were
  seen, `reliability` is structurally barred from absolute frozen
  cuts (§5); any dev-derived gate is an amendment frozen off-unit
  before confirmatory. (c) Recon result files stay; the coverage
  stage re-pulls everything fresh with snapshots.
- **D2 — era and input**: era MJD 60790.117 → 61047.354 (re-measured
  at coverage; larger span governs); `crossings/universal_v1`
  Earth-center under the standing 0.010 R☉ ground-observer budget;
  processing pinned v30/DM-55060.
- **D3 — substrate**: catalog-level per §5 — DiaSource association
  primary; VisitDetector footprint gate; deep-coadd Object for the
  static exclusion; ForcedSourceOnDiaObject as veto/annotation only;
  no image-level statistic; every query snapshotted.
- **D4 — units and statistics**: unit = (target, channel-rung, band);
  one frozen statistic **S_det** per unit (§5); expected family from
  scoping: **1 unit × 1 statistic = 1 trial** (ross-128 B 0.1 AU r),
  final tally at the threshold freeze; FWER α = 0.05.
- **D5 — gates**: §7 exactly (footprint gate; pixelFlags exclusion
  set incl. archive injections; magLim non-null).
- **D6 — channel-A saturation rule**: §6 (G 16.5/17.5 bands,
  Gaia DR3 annotation, cut-stage verified).
- **D7 — controls**: 8 pseudo-position controls per unit-visit on the
  same visit and detector: the full deduplicated z-family pattern
  translated to positions at the locus offset rotated k·40°
  (k = 1…8) about the field position, each ≥ 10″ from the real locus
  segment; a control position failing the §5 static-sky exclusion or
  the §7 footprint gate rotates +5° until valid (frozen resolution
  rule). Off-window visits at the unit position are the dev-stage
  pseudo-unit family (D8), not confirmatory controls.
- **D8 — dev/confirmatory split**: **dev = pseudo-units only** — the
  ~44 off-window ross-128-antipode visits (association-rate null,
  control-rule census, positive-control chain, archive-injection
  completeness probe, §6 marginal-band census if needed).
  **Confirmatory = every §9-scoped searched unit that survives the
  coverage gates, blind.** No searched unit is spent on machinery.
