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
- Two structural theorems, confirmed independently three+ times: an
  elongation-90° surveyor cannot see channels B / A-0.1 at all
  (elongation gate); and ~annual multi-month windows defeat the
  8-offset pseudo-window temporal control family (ZTF, PS1, WISE,
  TESS channel A alike) — a v2-style null ensemble is the designated
  route if the wide rung is ever searched with calibrated error rates.
- Blended on-star channels are systematics-dominated (PM dipole against
  multi-year references, k factors 10²–10⁴); the parallax-factor
  systematics template + empirical variance rescale recipe is reusable.
  The antipode (downlink pre-lens) channel is the workhorse.
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
- SPHEREx template refit with injections (slow-source absorption
  unmodelled); a six-detector SPHEREx joint cell from the stored
  accumulators.
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
