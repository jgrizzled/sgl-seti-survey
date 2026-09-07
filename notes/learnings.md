---
title: "sgl-seti-survey — consolidated learnings"
date: 2026-08-24
status: "living reference — add to it as surveys close; keyed by theme, not chronology"
---

# Consolidated learnings

Durable design lessons from the v1 surveys, the WISE scientific review,
the v2 programme, and the crossings programme — the material future
adapters and freezes must not rediscover. Sources: the v1 reports and
`surveys/wise/{scientific_review,v2_plan}.md` were deleted at v1
retirement / reorganization and live in git history (recover with
`git show 8fb226f:surveys/wise/v2_plan.md` and
`…:surveys/wise/scientific_review.md`); the frozen v2 decision rule
itself is live in `surveys/wise/hypotheses.md` v2.1 (and each survey's
hypotheses doc); crossings lessons are detailed in the live
`report/*_crossings.md`.

## 1. What a v1-style search can and cannot claim (WISE scientific review, 2026-08-21)

The review's verdict, which applied to every v1 survey (they copied the
same constructions): the infrastructure was sound and transparent, but
the statistical experiment supported only *"no compelling candidate
remained after the survey's heuristic review rules"* — not a
90%-complete exclusion, not a controlled global false-alarm rate, not
the W3/W4 warm-structure limits, not any population/Picky-Network
inference. The nine findings, condensed:

1. **"T = max of N offset controls" is a rank statement, not a false-alarm
   rate.** It gives a 1/(N+1) per-search crossing probability (WISE: ≈78
   expected crossings across 704 searches vs 70 seen — exactly the null);
   no family-wise error rate existed across endpoints × roles × bands ×
   grid. Vetoes tuned after looking at candidates cannot be a confirmatory
   test.
2. **Completeness that stops at threshold crossing is not survey
   completeness.** v1 injections never passed through the vetoes; rules
   like the W3/W4 "single-visit cryo cell" veto auto-rejected 326/352
   cells (the *expected* cryo cadence), silently gutting the claimed
   warm-structure exclusion.
3. **Declared geometry must be computed geometry.** The "99% confidence
   locus" was stored as metadata while the search evaluated the nominal
   locus; covariance was never propagated.
4. **Tensor-level analytic injections measure the model, not the
   instrument.** Gaussian-on-grid injections cannot see PRF throughput,
   masking, background, or calibration failures; 32 reps/point gives no
   CI on a 90th percentile.
5. **Implementation must match the frozen hypothesis exactly**: prior vs
   grid distinguished; L∞ vs L2 motion bound declared; constraint
   intervals must tile the distance domain (v1 left 9.2% of prior mass
   in gaps); one duty-cycle number hides distinct temporal models.
6. **Uncalibrated vetoes are annotations.** Phase balance, proximity,
   significance ratios, "star-like" flags can all reject an in-scope
   source; without a measured selection function they may only flag,
   never reject.
7. **Quality masking needs a pre-specified primary mask plus a
   strict/loose sensitivity analysis** reported as a systematic.
8. **A frozen targeted portfolio is not a population sample** — no
   occurrence-rate inference without a generative selection model.
9. **Provenance claims must match the implementation**: hash contents,
   not filenames; stale-product detection mandatory; report tables built
   from the ledger with the build failing on mismatch.

## 2. The v2 statistical design — and what the data changed

Principle: **one frozen decision rule, applied blind to one end-to-end
injection set, one exchangeable null ensemble, and the real data — in
that order.** Development/confirmatory endpoint split stratified by
confusion class, seed in the freeze; any rule change after touching the
confirmatory set voids the claim.

Amendments forced by evidence on development sets (the retrospectives):

- **Ring-only null.** Of the three planned null constructions, only the
  48-offset spatial ring is exchangeable. Per-epoch time scrambling
  destroys the static-sky coherence of the annual parallax return
  (scramble maxima far below spatial controls); trajectory
  randomisation samples other parts of the corridor (KS vs ring fails
  in ~55–64% of cells). Both are kept as per-cell annotations
  (p_phase, p_trajectory), never pooled into the error rate.
- **Per-cell normalisation with heavy-tail exclusion.** Family statistic
  max R̃ = R/q95(ring); cells void if the ring max exceeds 2.5×q95 or
  inner/outer ring members differ (KS α=0.01). The family-wise α=0.05
  price scales with the heaviness of the local nulls: R̃_FWER 1.47
  (PS1, clean) → 1.66 (ZTF) → 1.78 (joint) → 2.35 (SPHEREx dev) → 2.44
  (WISE); depth cost vs the uncontrolled v1 numbers 0.3–0.6 mag
  (optical) to 1–1.5 mag (WISE).
- **Calibrated veto vs annotation.** A rule may reject a candidate only
  if it (i) rests on independent evidence of a non-SGL origin, (ii) can
  only *lose* a genuinely swamped in-scope source (a completeness cost,
  charged to the final-candidate curve), and (iii) has its selection
  function measured by injections. Everything else is an annotation on
  the Candidate record; survivors are `retained-ambiguous`, never
  discretionarily dismissed. Only the flux-consistent catalogued-static
  / halo test (and, for WISE, the W3/W4 cross-band confirmation
  procedure) qualified.
- **The held-out-epoch prediction test failed as a veto** (passes
  59–63% of null trajectories, rejects 74–76% of in-scope long-block
  sources) → annotation. Epoch hold-outs: pre-register a true hold-out
  where epochs keep arriving (ZTF, SPHEREx QR3); where the mission is
  over (WISE, PS1) the post-hoc test is an annotation with measured
  pass rates.
- **Static-source double counting.** Wherever the search image already
  removes the static sky, a catalogue flux-consistency veto is wrong:
  in ZTF's difference regime scale the predicted static flux by
  (1 − dfrac) per epoch (without it the veto "explained" injections and
  m90 collapsed 3 mag); with SPHEREx's static template, no catalogue
  veto at all. In a joint stack, sum the per-archive catalogue
  predictions against the joint peak — OR-ing per-archive flags
  over-rejects.
- **The single-epoch clip is a bright completeness limit.** The layered
  search clips any epoch with a single-frame ~5σ detection, so stack
  recovery is non-monotonic near that limit; fit completeness only
  fainter than each cell's single-epoch limit and carry the bright
  bound on the constraint (the catalogue layer owns brighter sources).
  Corollary: a stack-regime positive control must be *fainter* than the
  single-frame limit (the ZTF/DECam asteroid controls are
  catalogue-layer objects under the frozen clip).

## 3. Injection and completeness

- Inject into the calibrated image (DN, before background/masking/
  matched filter), masks untouched; the stamp-response linearity trick
  (cache the unmodified flux map, re-convolve only a stamp) makes
  ~10⁵ injections affordable.
- Empirical PRF per archive (WISE PRF grid, ZTF/PS1 per-quadrant/
  per-skycell seeing Moffat, SPHEREx PSF cube), subpixel phase drawn
  uniformly; the matched filter may stay Gaussian — the
  Gaussian-on-PRF throughput is then *measured* (WISE: W1 0.80 / W2
  0.75 / W3 0.46 / W4 0.58), not assumed.
- Inject a physical spectrum with colour corrections, not a magnitude —
  and for a per-pixel-wavelength instrument (SPHEREx), an SED.
- Continuous placement: z log-uniform from the prior, µ uniform on the
  declared L∞ box, cross-track from the propagated envelope, four
  temporal models (persistent / exposure flicker / visit-scale /
  long-block); coverage claims use the worst of the four. No grid
  nodes, no analytic mismatch factor.
- Report threshold and final-candidate completeness separately, with
  Wilson intervals; bin edges at reciprocal-distance midpoints so
  intervals tile the domain gap-free.
- Measure end-to-end throughput against the archive pipeline's own
  photometry of a known mover (the WISE asteroid control: −0.42 mag vs
  NEOWISE photometry, consistent with PRF + sampling), not against
  H/G predictions.
- Run at least one positive control blind through the full chain before
  any freeze (asteroid via SSOIS/Horizons; every archive got one).

## 4. Geometry and astrometry

- Propagate the endpoint covariance (seeded MC, ~2,000 draws); where
  the 99% cross-track envelope ≤ 0.5×PSF FWHM the nominal locus
  suffices; otherwise add a cross-track tensor dimension counted in the
  null. Optical PSFs make the test ~5× stricter than WISE: 16/138 (ZTF)
  and 18/138 (PS1) cells needed the dimension vs 3/176 (WISE).
- The envelope is observer-independent to 0.02″ — run the covariance MC
  with the Earth-centre observer (the terrestrial-site observer is 10×
  slower per locus evaluation).
- sglseti's Rx locus differs from a naive anti-star construction by
  2µz/c (0.18″ at 550 AU → 3.3″ at 10,000 AU for Barnard's Star); Tx
  agrees ≲0.05″ once the 2d/c light time is included. Spacecraft vs
  geocentre moves the WISE locus ≤ 17 mas; for TESS's HEO the
  correction is real (|Δb| up to 0.33 R☉, |Δt_ca| up to 3.2 h).
- Distance grids uniform in 1/z, never log z (a log grid left 17″ gaps
  at small z).
- Joint/multi-archive µ grids need a common µ reference epoch chosen
  mid-baseline: at T0 far from an archive's epochs a 0.5″/yr µ step
  under-samples the family (PS1 at T0=59800: off-grid i depth
  collapsed; g/r lost 0.45 mag) — use ~0.1″/yr sampling or
  interpolation.
- Endpoint curation dominates effort: Hipparcos/hip2 solutions for
  Sirius/Procyon are already barycentric; Gaia component solutions for
  tight binaries (GJ 65, RUWE ≈ 11) are orbit-corrupted; propagate CNS5
  epochs before cross-matching; always Kepler-check extracted orbits.
- Use the full (SIP) WCS when a cutout spans arcminutes of distortion
  (SPHEREx: local-linear sampling defeated the template, 30×
  contamination); read the pixel scale from the header everywhere (v1
  hardcoded 2.75″ for W4's 5.52″ pixels — kernel 2× too wide).

## 5. Photometry and flux calibration

- The shared matched filter must return total flux, not PSF peak
  amplitude (v1 depth labels were optimistic by 1.9–3.3 mag); and the
  AnalysisRun config hash must include the photometry convention — the
  first recalibration silently reproduced the old run id.
- Per-frame star-calibrated zero points through the identical matched
  filter are the default wherever a per-epoch star catalogue exists
  (PS1 lesson; confirmed at DECam, where header MAGZERO is a counts
  convention with +3.3 mag outliers).
- Confusion, not pixel noise, sets depth at ≥6″ resolution; thresholds
  come from local controls, per-pixel uncertainties understate the
  floor.
- Cap per-epoch (v2: per-frame, 20× band median) weight so one frame
  cannot dominate a stack (v1 saw S=398/399 from one frame); report
  N_eff.
- Filter sentinel magnitudes (PS1 −999) before any 10^(−0.4m).
- Verify the flux scale against catalogued stars (≤0.2 mag) before
  quoting any depth.

## 6. Screening and vetting

- The parallax-phase test is the decisive cheap veto (a static source
  recurs at one day-of-year window; a relay must appear at both
  phases), but it loses power in dense corridors and cannot run at all
  on single-phase cadence (PS1 3π revisits at 97:3) — cross-archive
  stacks supply the missing phase.
- The catalogued-static test is a calibrated veto only with the
  flux-consistency check (catalogue mag → matched-filter S within ~2×)
  and the phase signature; bare proximity is an annotation. Require ≥3
  catalogue detections (an SGL source at 20″–375″/yr parallax cannot
  form a compact multi-epoch catalogue entry).
- ZTF i-band is too sparse/clustered to carry the phase test; carry it
  in g/r.
- Same-night TTI pairs are a free ordinary-mover veto (PS1).
- A real source may be absent from static catalogues, split across
  entries, or live only in reject tables — catalogue absence is neither
  a detection nor a null.

## 7. Per-archive facts worth keeping (v1 reports deleted)

- **WISE:** the ~130 KB `-msk` products carry exact WCS + usable pixels
  without the image — precise passes are nearly free. W3/W4 exist only
  in the cryo mission (~12 exposures, one visit): no long-block
  constraint is possible there, and any W3/W4 exceedance needs the
  cross-band (W1/W2 flux-ratio + post-cryo persistence + catalogue)
  confirmation procedure, never a cadence veto.
- **ZTF:** fixed-grid CCD gaps hit ~13% of sky (all three pilot
  antipodes); reference-image edge strips need a hybrid search image /
  per-corridor references; difference-image injections must go into
  the *science* image before differencing.
- **PS1:** full skycell masks (3.3 MB fpack) give exact WCS +
  usability; `CONV.BAD` marks resampled OTA-gap bands that hold finite
  image values — masks, not NaNs, define usability.
  `detection.obsTime` trails warp `MJD-OBS` by ~50–60 s and ~15% of
  warp epochs have no catalogued detections. Cutouts are purged after
  tensoring — image-level work re-fetches per batch.
- **SPHEREx:** sub-pixel-phase matched filter for the undersampled PSF
  (+0.4 mag); a static-sky template substitutes for difference images
  (~3 mag) and is the static treatment (no catalogue veto on top);
  deep-field observing seasons are detector-dependent; an Rx/Tx
  role-coincidence flag catches template artefacts; template
  absorption of slow real sources is not yet injection-modelled.
- **DECam:** per-exposure EXTNAME→HDU maps with fetch-time asserts;
  NSC DR2 catalogue epochs end 2017–2019 (time-partial screening);
  ~24% of arc centres fall in chip gaps and some exposures span 2–3
  CCDs — select CCDs by the covered locus, not the arc centre. Access
  details: `surveys/decam/notes/decam_recon_2026-08-24.md`.

## 8. Crossings (Pipeline B)

- The crossing geometry phase-locks to the sidereal year: wide-beam
  rungs have window ≈ observing season and are constraint-only by
  construction — design narrow-rung-first for optical archives.
- Two structural theorems, confirmed independently four times: an
  elongation-90° surveyor cannot see channels B / A-0.1 at all
  (elongation gate — Kepler/K2 2026-09-04 is the spacecraft case:
  a solar-array pointing constraint that keeps the boresight 60°–145°
  from the Sun makes a whole mission blind to every rung ≤ 0.1 AU,
  and a footprint intersect can be closed by geometry before any
  pixel is touched); and ~annual multi-month windows defeat the
  8-offset pseudo-window temporal control family (ZTF, PS1, WISE,
  TESS channel A alike) — a v2-style null ensemble is the designated
  route if the wide rung is ever searched with calibrated error rates.
- Blended on-star channels are systematics-dominated (PM dipole against
  multi-year references, k factors 10²–10⁴); the parallax-factor
  systematics template + empirical variance rescale recipe is reusable.
  The antipode (downlink pre-lens) channel is the workhorse.
- "Archivally virgin" is rung- and mode-specific (radio, 2026-09-04):
  the antipode channel has zero *targeted* radio pointings, but
  wide-field ASKAP continuum epochs (VAST/RACS/FLASH) do fall inside
  the ±6 d 0.1 AU windows — 6 validated epochs on 5 targets, none on
  the ≤ 1 d grazing windows. State which rung and which archive class a
  "nobody has looked" claim covers.
- A spacecraft substrate needs its own crossing list before any
  intersect, and the cost is now ~15 min (`--fetch-observer` +
  `--observer` in `sglsurvey/crossings.py`; TESS, SOHO, Kepler). The
  Earth-center list is wrong for a heliocentric-orbit observer at
  *every* rung (Kepler: |Δt_ca| median 33 d), not just at grazing b as
  for HEO/L1 spacecraft.
- MAST CAOM TAP: the sync endpoint 504s at 60 s, the async UWS
  endpoint completes the same query; `GROUP BY sequence_number`
  campaign ranges must be clipped to planned dates (a K2 product is
  mis-tagged C12 → 14); EPIC/2MASS-epoch catalog positions put a
  4.7″/yr star 80″ from its own K2 target row by 2017 — search by
  epoch-propagated position with a ≥ 2′ box.
- CASDA ObsCore is the cheapest revisiting-radio substrate: `cube` +
  `cont.restored.t0` rows carry `s_region` and per-SBID intervals for
  every ASKAP project at once (one TAP cone per position, ~7 s); the
  VAST-pilot cubes lack `t_min` (join `casda.observation` by SBID),
  several cubes per SBID need a quality-ranked dedupe, and `s_fov` is
  the image bound, not the footprint — declare your own. LoTSS DR3
  `pointings.dateallobs` gives every 8-h run's mid-MJD; DR2 has only
  `dateobs`.
- Radio quick-look access map (2026-09-04): CASDA *catalogues* are
  anonymous TAP — VAST full-survey SBIDs each have an epoch-resolved
  component table (`AS207.vast_extragal_dr1_<field>_sb<sbid>_components_v01`),
  RACS has release tables (only `racs_mid_*` and `racs_low2_v_*` carry
  `sbid`), and the level-5 Selavy files of pilot/guest/FLASH SBIDs are
  not exposed to TAP but downloadable; CASDA *images* and level-5
  files are HTTP 401 without an OPAL login. With one (`.env`
  `OPAL_USERNAME`/`OPAL_PASSWORD`): basic auth on the DataLink VOTable
  → per-product `cutout_service` token → POST `ID=<token>` to
  `casda_data_access/data/async`, then `/parameters` (CIRCLE) and
  `/phase` RUN, no auth; ~25 s per 3′ cutout. Permissions are per
  product (pilot v1 and REJECTED images closed). VLASS is the opposite: CADC CAOM2 has exact quick-look
  tile time bounds (VLASS4.1 planes lack them — read the tile header)
  and SODA cutouts (`minoc/files/<uri>?CIRCLE=ra+dec+r`) are anonymous.
  ASTRON `lotss_dr3.main_sources` is a mosaic of all runs, never
  epoch-resolved.
- The 1/9 exceedance budget behaved exactly as designed in both
  searched surveys (ZTF 3 vs 4.0 expected; PS1 2 vs 2.1).
- **Correlated exact-mask attrition** is the single-phase-cadence
  failure mode: a window's epochs are 1–2 same-pointing nights, so
  ~25% per-epoch mask attrition becomes ~40% whole-unit loss (it took
  the era's best grazing event, van-maanen b = 0.28 R☉). Exact-mask-test
  the narrow-rung nodes at the *coverage* stage, before the freeze
  counts a window as covered.
- The warp-direct (no-differencing) substrate works: the confusion
  floor lands in the ring-control threshold, and a deep static *stack*
  catalogue replaces difference imaging as the static-vs-transient
  discriminator.
- TESScut serves collateral (off-science) CCD pixels — a science-array
  usability screen belongs in the coverage stage (it removed 3/8
  cutouts, including a dev unit).
- The pre-registration ordering discipline works: two WISE-crossings
  freeze defects (nonexistent qa_status value; mis-scaled control pool)
  were caught and amended before any pixel was touched.
- **An undetrended chord/stack statistic saturates its error budget
  with low-frequency drift** (TESS sector-scale scattered light: S and
  T in the tens, 3 exceedances vs 1.3 expected, all
  systematics-adjudicated): the empirical variance rescale k fixes the
  *variance*, not low-frequency structure, so the 1/9 control-crossing
  budget is only approximate for such statistics. Per-cadence /
  differential statistics (the pulse max) behave exactly to budget —
  any high-cadence crossings v2 needs a detrending layer under the
  chord filter.
- **Injection flux calibration through a star-measured ZP must
  normalize the stamp response at the calibration reference** (TESS
  finding C1): the chain-measured ZP already absorbs the kernel↔PRF
  throughput, so applying the raw stamp response to a ZP-converted
  flux double-counts it (×R error — 1.3× optimistic on one cube, 2×
  conservative on another). Also measure the ZP through the *unit's*
  fitted kernel, not a default one.
- Sector/tile boundary truncation is a structural weakness for
  short-window archives: both TESS chord exceedances needing
  adjudication and the retained-ambiguous row sat in a sector's final
  day (no egress to shape-test), and the era's two deepest grazes fell
  between sectors entirely.

- **Windows shorter than the revisit interval are covered by luck**
  (ATLAS crossings, 2026-09-04): nightly cadence delivered 0.1 AU
  windows (10–12 d) at 95 % but 0.3–1.3 d grazing windows at 18–45 %,
  and the pre-declared recurrence stack (100–180 exposures) was 14–30
  in reality with ≤ 3 valid pseudo-stacks — the recurrence cell needs
  a substrate whose cadence is well inside the window. Also: count
  crossing windows with the axis-side filter — link direction alone
  double-counts the sunward crossing of each year (caught pre-data by
  an elongation check).
- **Mirror-gated control validity with 1/(n+1) per-trial accounting**
  keeps the exceedance budget honest when the nominal 8 controls do
  not exist; the frozen max rule then makes a trial *artefact-set*
  (insensitive, still valid) whenever one pseudo-window holds a
  single-frame outlier — a chi/N mask term and a ≥ 2-epoch S_event
  gate are the cheap fixes for a v2.
- **Channel A on a high-PM star in difference imaging is the PM
  dipole**: at 5″/yr against a multi-year template the "difference
  flux" is the whole star (k 10²–10⁷). Either carry the parallax/PM
  systematics template or use reduced-mode photometry at a per-epoch
  PM-propagated position.

## 9. Reproducibility and engineering

- Content-hash (bytes + schema version) every input and derived
  product; `--only-missing` must refuse products whose recorded input
  hash changed; report tables are generated from the ledger and the
  build fails on count mismatch.
- Freeze files stay byte-identical to the pre-registered versions; a
  freeze's `hypotheses_hash` refers to the document *as of freeze time*
  (git history) — future freezes hash the complete, self-contained
  active hypothesis document so this ambiguity cannot recur.
- Set OMP/OPENBLAS threads to 1 under multiprocessing; `set -o
  pipefail` before any stage whose failure would let a purge step
  delete its inputs (a hidden exit status behind `| grep` cost one
  rebuild, a locus-cache eviction bug another).
- v2-scale costs for planning: WISE 176 pairs ≈ 1.3 h tensor build /
  1.5 h injections per hold-out set on 7 workers; ZTF ≈ 5 h + 6 h;
  PS1 16 fetch/build/purge batches ≈ 7.5 h; plan-stage compute
  estimates ran ~3× pessimistic, analysis-time estimates optimistic.

## 10. Open items carried forward (from the v2 retrospectives)

- ~~A stack-regime (fainter-than-single-frame) asteroid control for ZTF
  and PS1.~~ Done 2026-08-25: asteroid (220000) at V 21.9–23.05 (every
  ZTF frame below the clip; PS1's one >5σ frame clipped by the frozen
  rule) recovered by the v2 rule — ZTF zr R̃ 2.90 / zg 1.95, PS1 i
  1.77 (z an honest non-detection at its depth), and all three joint
  g/r/i families above the family threshold (R̃ 2.21/3.04/1.77 vs
  1.54); throughput −0.14…+0.29 mag vs Horizons+solar colours.
  Fetch: `surveys/{ztf,panstarrs}/scripts/asteroid_stack_fetch.py`;
  scoring `surveys/joint/scripts/asteroid_control_stack.py`; products
  `runs/{ztf/v2,panstarrs/v2,joint/v3}/control/220000_stack/`.
- ~~SPHEREx template refit with injections (slow-source absorption
  unmodelled); a six-detector SPHEREx joint cell from the stored
  accumulators.~~ Done 2026-09-04 (`report/spherex_joint6.md`; §13
  below): joint cell 0 candidates (122 confirmatory cells, persistent
  m90 20.93); the template absorbs a median 49 % of a slow source's
  flux — the v2 depths are overstated by ~0.7 mag for slow sources.
- ~~The WISE-joint (three-archive) stage on the common µ reference
  epoch, with Vega→AB and surface-brightness conventions reconciled —
  the 0.5–22 µm colour-consistency test then becomes a calibrated
  veto.~~ Done 2026-08-24 → 25 as joint v3.0 (`report/joint_ps1_ztf.md`;
  history §8): 0 candidates blind, optical bit-identical to v2, colour
  veto's measured false-veto rate 0 in the fitted population. Lessons:
  a µ-reference-epoch change is a rebuild, never a relabel; report the
  frame zero point in the target system (Vega MAGZP + AB offset)
  and the engine's zp_ref scaling does the rest; floor a
  confusion-limited archive's per-node σ with the empirical ring
  scatter or the veto over-fires; the colour veto's only firings live
  above the optical bright limit (PRF throughput × temporal-window
  overlap), so charge it to the final-candidate curve and it costs the
  fitted constraints nothing; 0.5–22 µm remains 0.5–4.6 µm in practice
  (W3/W4 stay excluded by the review rule).

## 11. Catalogue-level substrates (DASCH; cf. Rubin DP2)

The 2026-08-26 surveys added a substrate class where the primary
chain never touches pixels — the archive's own calibrated detection
catalogs. Lessons from DASCH (`report/dasch_crossings.md` §4):

- **Archive quality flags are tuned against the signal.** A
  catalogue pipeline's defect classifier is trained to suppress
  single-plate/single-epoch transients — DASCH's `SUSPECTED_DEFECT`
  killed the real-asteroid positive control 1/1. On this substrate
  an archive quality bit may gate the statistic only after its
  selection function against in-scope sources is measured;
  otherwise it is an annotation feeding pixel-level adjudication.
  (The v2 calibrated-veto rule, rediscovered at the flag level.)
- **High-PM stars break refcat-keyed archives silently**: DASCH's
  APASS refcat carries pm = 0 dummy entries for 2 of 5 grazing
  stars, yielding confidently wrong lightcurves with no error
  signal. Standing check: catalogued PM vs registry µ
  (vector match within max(10 %, 50 mas/yr)) before trusting any
  archive lightcurve of a nearby star; route to an alternate refcat
  or a positional (limits-only) search on failure.
- **A sub-limit star is its own hit background**: marginal at-limit
  extractions of the quiescent star appear window-independently and
  are invisible to spatial ring controls — measure the off-window
  same-locus rate and gate hits to ≥ 1 mag brighter than measured
  quiescence (cost charged). The measurement doubles as a free
  photometry validation (teegarden 17.11 vs expected ≈ 17.3 B).
- **At-limit + non-PSF morphology + ~3× the genuine astrometric
  scatter is the plate-defect signature**; cutout pixels adjudicate
  it in minutes, with the positive-control transient (compact,
  15–20σ) anchoring the discrimination. All 4 adjudicated DASCH
  hits fit it.
- Old-epoch geometry is not the hard part: the 1885 model budget
  measured ≤ 2×10⁻⁴ R☉ (astrometry) + 1×10⁻⁵ R☉ (ephemeris) —
  plate *timing* (logbook dates, timezone-era conventions; the
  Iris control implies ~1 h errors) dominates and belongs in
  coverage gates and control-match windows, not geometry budgets.
- JPL `sb_ident` 500s on pre-1950 epochs; the workaround is direct
  Horizons per-asteroid ephemerides × plate-date coincidence
  scanning (asteroids 1–100, 1-day steps, minutes of wall time).

- **Catalog-level light-curve releases can be un-searchable for
  exactly our targets** (PGIR DR1, 2026-09-04): forced photometry at
  the reference-catalog epoch (2MASS, 1998–2001) puts every nearby
  high-proper-motion star 15–120″ from its own entry two decades
  later — the entry reads sky, the flux lands in drifting neighbour
  entries — and empty-sky antipodes have no entry at all. Check
  three things before trusting any such catalog: is the position
  PM-propagated, are non-detections kept (they were; the visit list
  survives), and what is the time column's storage type (float32 JD
  = 0.25 d bins at JD 2.46×10⁶ — invisible in a schema description
  that says "Julian date"). A visit list is still worth recording:
  it turns a future data ask into a 56-window request.

## 12. Heliospheric imagers with a drifting observer (STEREO-A HI-1; cf. LASCO §8)

- **Read the per-day pointing history before the census.** STEREO-A
  flew rolled 180° for eight years (2015-11 → 2023-08): HI-1 looked
  west instead of east, and every sunward arc moved from before t_ca
  to after it. A recon frame from one era is not the mission; one
  8-KB Range request per day (6,675 headers in 10 min) is.
- **Level-1 heliospheric frames need the instrument team's
  background removal, not a local one.** A 31-px median high-pass
  fails within its half-width of the CCD edge on a 2 DN/s/px gradient
  (+16 DN/s ramps) and a 5–9 px annulus is curvature-biased on the
  F-corona ridge — and the ridge is exactly where every ecliptic
  source sits while HPLT-offset controls sit on the flanks, so no
  same-frame differential cancels it. The level-2 per-pixel running
  lowest-quartile background removes the static F-corona exactly
  because the sources drift 54 px/day. Random-position bias tests
  cannot diagnose this (avoiding detected peaks selects troughs);
  pixel stamps at the actual patch can.
- **Frame-wide bright-body vetoes do not scale to 20° fields:** Venus,
  Jupiter and Earth are inside HI-1 for weeks to years; the measured
  scatter is unaffected beyond 2° of the patch. Veto by proximity,
  and measure the radius before freezing it.
- **ICRS-fixed sky patches + the patch's own transit baseline** are
  the right construction for a moving-sky imager: no star transits
  through the aperture during a window, and the static field
  subtracts itself. The residual is the regional level-2 background
  offset between the baseline and arc elongations (±0.3 units) — the
  plateau class of both adjudicated S_event exceedances — a v2 should
  take the baseline from the same elongation band in other years.
- **Extended fronts masquerade as pulses.** Both pulse exceedances
  were CME/streamer fronts sweeping over the patch (a 10-px band at
  2–3 px/frame; a whole-stamp lift for two frames); the persistence
  rule alone passes them — the point-source morphology test in the
  pixels is what vetoes them. Pulse thresholds from control patches
  that contain bright static stars are 5–100 (their own excursions);
  match control patches to the source's static content in a v2.
- Horizons accepts `CENTER='@-234'` for any spacecraft-observer
  ephemeris; SkyBoT does not know MPC code C49 — keep a bright-
  asteroid list for known-object censuses from spacecraft. A bright
  asteroid arc passage is a clean positive control for a moving
  observer (Pallas: 0.01 mag).
- A 0.60 mag/(B−V) colour term turns HI-1's 0.27-mag star scatter into
  0.08 (0.047 with a 2-D ZP); forgetting the sign doubles it — the
  instrumental magnitude is V − c(B−V − 0.65), red stars brighter.
- Crowded low-latitude fields (ross-154, b = −13°) blow a per-star ZP
  MAD gate set on a clean field (0.05 → 0.12–0.17 from blending) even
  though the ensemble ZP stays good to 0.006 mag; gate on the ensemble
  precision, or class the target, before the freeze.
- Fork-safety: `np.load` of a compressed npz is lazy; forked workers
  sharing the handle corrupt each other's reads ("Error -3 while
  decompressing"). Copy arrays eagerly at import.

## 13. Static-sky templates with few visits (SPHEREx deferred controls, 2026-09-04)

- **A static template fitted to data that contain the source absorbs it
  in proportion to the source's share of the node's weight — and with
  2–3 visits per node that share is ~half.** SPHEREx QR2: median 49 %
  of a persistent source's stacked flux is removed by the (a, b·λ)
  fit, 80–90 % in single-visit corridors (a source present in the
  node's only visit *is* a star to the fit), 5–10 % in the deep field.
  The effect is set by the visit count, not by z (a 550 AU source
  moves 375″/yr but is static within a visit), not by residual motion,
  and it is linear in magnitude below the clip. Injections added after
  the template subtraction cannot see it: measure it by refitting the
  template with the injected per-epoch flux added, and anchor the fast
  refit against a full image-level injection (agreement < 0.01 here).
  Visit- and block-scale sources are absorbed like persistent ones;
  only exposure-flicker is absorbed proportionally less. Forward fix: a
  source-excluded template (per node, drop the epochs during which the
  hypothesised track is within ~2 FWHM) — exact by construction.
- **A ring normaliser must be a positive scale.** In an over-subtracted
  field every trajectory's S_max can be negative; T becomes a tiny
  positive number and q95 negative, and R̃ = R/q95 turns the most
  negative ring values into the family's largest statistics (joint dev
  R̃_FWER 5.26 → 1.59 after voiding q95 ≤ 0). The per-detector engine
  lacks the guard and was saved by the heavy-tail rule's sign; add it
  at the QR3 re-run.
- **A joint cell from stored accumulators is cheap and worth ~√N.**
  Σ_b A_b / √Σ_b B_b over six detectors with the same grid needs no
  image, only a re-run of the injection chain with a common magnitude
  window so the j-th injection is one source everywhere (the profile's
  `inj_subdir` is now generic). Gain 0.6–0.9 mag over the best
  detector; the joint injections still draw visit/block on-patterns
  per detector although dichroic pairs share exposure times — fix in
  the chain, not the analysis.
- **Kept products must be real files, not links into directories
  scheduled for deletion.** The v3 templates survived the v1 retirement
  only as a dangling symlink; the fix was regeneration from the
  retained cutouts plus a rebuilt-tensor comparison. When comparing
  rebuilt tensors, match epochs by observation id — exposures of a
  dichroic pair share an MJD and `argsort` orders them arbitrarily.
- **Two memory-hungry pools on one 60 GB box will OOM each other, and a
  `multiprocessing.Pool` with a killed worker hangs rather than fails.**
  Size the per-worker peak first (template epoch cache = nodes × epochs
  × 8 B; injection maps ≈ 5 MB × exposures per band), sample only the
  nodes you need, and watch memory in the monitor loop.

## 14. Archived spectra as a crossings substrate (spectral-archive family, 2026-09-05)

The first non-imaging survey. Reduced 1D spectra with sub-minute time
stamps are the cleanest substrate the programme has met (no WCS, PSF,
mask or ZP machinery), and the same star's out-of-window spectra are a
free null ensemble. The lessons are about what a same-star ensemble on a
barycentric grid does *not* see.

- **Observer-frame emission smeared across a barycentric ensemble is
  invisible to the per-pixel σ and lands on the threshold instead.**
  Airglow (OH Meinel bands, [O I] 5577/6300), lamp lines and
  telluric-correction residuals are fixed in the observer frame; 60
  nulls spanning the year put each such line at 60 different
  barycentric pixels (±30 km/s), so no pixel's MAD is inflated, while
  the in-window spectra of one week put it at one pixel — 12 of the 18
  confirmatory exceedances were exactly this. v2 must mask sky
  emission in the observer frame *before* resampling, and treat
  telluric transmission < 0.9 as a mask for the NIR line cells (the
  wolf-359 NIRPS 1569.50 nm "recurrent line" sat on a CO₂ line at
  transmission 0.79 and was present in all eight spectra of two nights).
- **Whole-spectrum contamination needs a spectrum-level veto.** One
  ESPRESSO ross-128 frame carried the Hg I lamp pentad and 31 peaks
  above threshold; the per-feature ladder dispositions each peak, but a
  count of PSF-consistent peaks across the band (> 10 → contaminated)
  is the honest gate. The same pentad also revealed a +84 km/s grid
  error: **ESPRESSO phase-3 `WAVE` is vacuum, HARPS `WAVE` is air** —
  never assume a frame per archive; check a known line (K I 7699,
  Na D) per product type at recon.
- **Thresholds set by the max over 60 × 2×10⁵ pixels are hostage to
  the ugliest pixel.** Cosmic hits (1–3 px), normalisation blow-ups
  next to masked runs, zero-flux orders of an M dwarf normalised by a
  near-zero running median, a frame with SNR 0.95 and a frame with a
  wrong header BERV each put T at 10³–10⁵ on the first pass. The
  remedies, in order of leverage: a PSF-consistent statistic (only
  Gaussian-fit-compliant local maxima count — identically for nulls and
  in-window), the SNR gate at spectrum *and* pixel granularity, each
  spectrum's own photon error as the σ floor (the ensemble-median floor
  over-weights low-SNR spectra: robust z scale 0.6–3.4, r = −0.9 with
  SNR), and a wavelength-solution gate against the template. Four dev
  passes were needed; every one of them was a gate, not a hypothesis.
- **Bookkeeping for "max over n spectra vs max over N nulls" is
  n/(N+n), not 1/(N+1).** Units with 18–40 in-window spectra have null
  exceedance probabilities of 0.23–0.40 per cell by construction; the
  coadd statistic is the sharper test for them, and the expected-count
  line of the report must use the right formula.
- **Unit-conversion arithmetic in a recon note is a freeze input —
  check it with a second route.** The recon's order-of-magnitude
  power floors were 10⁴ too high (erg→W applied as 10⁻³ instead of
  10⁻⁷ J per erg with the cm²→m² factor); the injection chain caught it
  because the numbers disagreed with a hand calculation from the
  measured continuum. Related: **log-log interpolation of Gaia+2MASS
  photometry misses the 1.0–1.1 µm flux peak of an M6 dwarf by 2.3×**;
  NIRPS's absolute `FLUX_CAL` (checked to 0.15 mag against 2MASS J/H)
  is the right F_λ for NIR line cells, and a measured SED should be
  the default whenever a flux-calibrated product exists.
- **Barycentric bookkeeping**: SPIRou APERO and CARMENES caracal store
  observer-frame vacuum wavelengths with `BERV` in the header
  (λ_bary = λ_obs (1 + BERV/c), verified by stellar-line alignment);
  ESO phase 3 is already barycentric (`SPECSYS`). The BERV-sign check
  must be done on stellar-line regions — whole-band residuals are
  telluric-dominated and give the opposite verdict.
- **Product-level traps**: HARPS s1d has no error vector (`FLUXERR
  = -1`); X-shooter IDPs come in four variants per exposure and in nm;
  older SPIRou `t` files lack `OHLine`/`Recon`/`MJDMID`; SOPHIE public
  headers are date-stripped (BJD rounded to the day); CARMENES DR1 zip
  members carry an `_A` suffix and the Karmn `+` must be
  percent-encoded.

## §14 — WISPR survey (2026-09-06): template-based static removal on a fast observer

- **Static content dominates the exceedance budget when the template
  is incomplete.** 9 of 15 blind exceedances were static flux in the
  source aperture that Tycho-2/Hipparcos missed (G 10–12.5 stars,
  bright neighbours 2–3 px off-centre, a high-PM star counted twice
  because Tycho positions were not PM-propagated); Gaia DR3 accounted
  for each to ~20 %. Rule: any aperture-template construction must be
  built from Gaia with proper motions and a measured encircled-fraction
  curve — and the Gaia check belongs in the adjudication ladder.
- **Ridge curvature recurs in every field with a brightness ridge**
  (HI-1 H1 in ecliptic latitude, WISPR in orbit latitude): symmetric
  ±offset controls straddle the ridge and their median under-predicts
  the source's background. Quadratic interpolation through the control
  ladder removes it (van-maanen S1 S_stack 17 → 5).
- **Gate on the quantity that matters**: a per-star ZP-scatter gate
  rejected 87 % of perihelion frames although the ZP itself was
  determined to 0.05–0.14 mag; the uncertainty gate (MAD/√n) kept them.
  Likewise, single-frame pulse statistics on summed-exposure images are
  particle-hit-limited (control thresholds 57–410σ); persistence
  belongs in the statistic.
- **Bright in-patch stars are systematics-limited at a few percent**
  (V 4.7 at S/N 100 → ±10σ per event); self-calibration per exposure
  regime helps but cannot remove monotonic drifts. Declare the class.
- **Process hygiene**: long runs must be launched detached (`setsid
  nohup … &`) — a session restart and a power outage each killed the
  confirmatory run (append-only measurement files made resumption
  lossless); `pkill -f "<pattern>"` from a shell whose command line
  contains the pattern kills that shell (use `"[p]attern"`).

## §15 — SoloHI survey (2026-09-06): a star-fixed baseline on a fast observer, and its edges

- **Measure the field model from the headers, never from the
  literature.** SoloHI's tiles are fixed in the *orbit-plane* frame
  (the spacecraft rolls to keep them there) and sit on the anti-ram
  side — the opposite of WISPR; a solar-north field model would have
  wandered by ±5°. The 0.49° detector seams matter physically: the
  inner seam lies on the orbital plane, exactly where the deepest
  in-plane crossings (ross-128, b 0.13 R☉) put their sources.
- **A fixed source entering a field from larger elongation carries
  its own baseline.** Because the SoloHI source approaches the Sun
  through the field, every arc is preceded by 60–240 h of off-beam
  same-tile frames — the HI-1 star-fixed differential worked on an
  inner-heliosphere substrate and removed the template systematics
  that dominated WISPR's exceedance budget (0 static-content
  exceedances here).
- **Saturation hides as a plateau, not as DSATVAL.** The F-corona
  along the sunward edge saturates at 0.90–0.92 × the header DSATVAL
  in the ≥ 45-s exposure regime and erases stars without any flag;
  mask at 0.85 DSATVAL and let the validity gate remove the epochs.
  The same plateau removes any bright control star's core
  intermittently — control patches need a brighter mask (G ≤ 8) than
  source patches wherever saturation can bite.
- **Gate on measured per-epoch noise, not on geometry.** The inner
  ~120 px (ε ≲ 8°) of every arc has 6–11 units of unresolved corona
  structure per frame (against 0.5 outside), and it biases, not just
  scatters; the annulus noise `err` tracks it, and an absolute gate
  (err ≤ 1.5 units) turned control thresholds of 30–190 into 3–17.
- **Quadratic control interpolation is not robust.** The WISPR v2
  lesson (quadratic in latitude removes ridge curvature) was applied
  from the start here, and on one-sided ladders (125 of 178
  unit-events) it *extrapolates*: one control excursion of −50 units
  became +60 on the source. Half the blind exceedances were this. A
  robust interpolation (median or clipped quadratic) and the
  curvature correction are both needed; declare the ladder geometry
  per unit-event.
- **Movers set the pulse cell.** With 12–48-min cadence, uncatalogued
  objects crossing a fixed patch persist for 2–3 frames and pass the
  pair-min persistence rule; the pixel test (a peak that moves between
  frames) is the discriminator, and belongs in the ladder.
- **Process hygiene, again**: `pkill -f` from a shell whose command
  contains the pattern kills that shell (twice in one session) — kill
  by PID from `pgrep -f '^python …'`.
