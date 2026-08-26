# Dev-stage machinery log (started 2026-08-25)

Engineering record of the dev search stage (threshold freeze v1.0 §7
order of work). No survey position has been measured; everything here
is chain construction and validation on the recon frames (2010-06-15,
verified in zero sunward windows — hypotheses §11).

## Chain built (`scripts/lasco_lib.py`)

Load (level 0.5: OFFSET bias subtraction, MID_DATE/MID_TIME midpoint
timing) → azimuthal-median radial-profile detrend (0.1 R☉ bins, the
frozen D3 layer) → synthetic celestial WCS → Hipparcos star fit
(roll resolution + translation + radius-resolved star ZP) → Gaussian-
weighted forced photometry (σ 1.2 px, aperture 2.5 px, local annulus
5–9 px) at ICRS positions. Star catalog: Hipparcos V ≤ 9.5 with PM
and B−V, VizieR ASU-TSV snapshot sha-pinned
(`runs/heliospheric-crossings/starcat/`, 98,887 stars; the TAPVizieR
sync endpoint 403s — use ASU).

**Validations (C3 frame, 19 s, 2010-06-15):**

- Computing the Sun's direction **from SOHO** (Horizons table) instead
  of geocentric collapses the recon's 15.6 px parallax translation to
  ~1–2 px — the star-fit translation is now a residual-pointing term
  only.
- Roll hypothesis resolution works: the frame's CROTA 180 state is
  picked over the +180 alternative by match count (35 vs ~0).
- ζ Tau (V 3.01, B−V −0.15) recovered at V 2.89 through the full
  chain (S/N 108 single frame at r = 9.8 R☉).
- Radius-resolved ZP is mandatory on level 0.5 (not flat-fielded):
  a flat ZP shows ~1.3 mag "scatter" that is really the vignetting
  curve + sub-detection noise. Calibrator gates that fix it:
  S/N ≥ 10, V ≥ 4.5 (brighter saturates/bleeds low — β Tau V 1.65
  measures ~4 mag low), −0.3 ≤ B−V ≤ 1.2, radius inside the adopted
  annulus, one 3×MAD clip.

## Finding L1 — C2 per-frame star gate unachievable (amendment needed)

A representative C2 frame (Orange, 25 s) detects ~500 sources
(cosmic-ray dominated) but only the **2–3 brightest catalog stars
(V ≲ 5.5–6) are measurable** over the inner-corona background — the
frozen per-frame gate (≥ 5 matched stars) can never pass on C2. The
WCS itself is fine: the V 4.9 and V 5.4 stars land 2.4/2.5 px from
detections in the correct roll hypothesis and >35 px in the wrong one.
**Proposed amendment (to be frozen from dev ensemble statistics,
pre-confirmatory):** C2 astrometry = ≥ 2-star translation check
(falling back to the nearest fitted C3 frame's pointing when < 2);
C2 photometric scale = per-window stacked calibrator ZP
(EXPTIME-normalized across the window's frames; the ZP gate applies
to the stack, not the frame).

## Finding L2 — Clear-band color term dominates the ZP scatter

With all calibrator gates the per-frame C3 radial-ZP scatter is
~0.4 mag — the intrinsic star-to-star color spread of a ~400–850 nm
bandpass referenced to V (red stars measure bright: 119 Tau,
B−V 2.06, comes out 1.3 mag bright; B−V > 1.2 calibrators sit ~+0.6).
The frozen 0.2 gate is mis-calibrated for this bandpass.
**Resolution path (dev):** measure the global Clear-vs-V color
coefficient from the full dev calibrator ensemble (the dev driver
records (V, B−V, r, flux) per calibrator per frame), freeze it as a
static correction, re-measure the achievable per-frame scatter, and
re-set the gate from that distribution (with any residual scale
uncertainty declared TESS-style on the depths). Survey targets are
red (M dwarfs) and the S2 excess hypothesis is a laser line — the
color system must be declared either way.

## Status

`dev_fetch.py` running (18,294 (cam,day) entries: 61 full-cadence C2
days / 1,442 subsampled C3 days / 16,791 baseline days; ~36k frames
expected). Next: the dev search driver (per-event position tracks,
series assembly, S_stack/S_event/S_pulse, PA-ring controls), the
L1/L2 ensemble measurements, then the dev run proper on the 5 dev
units.

## Addendum (later 2026-08-25) — dev driver built + smoke-tested

- `scripts/tracks.py`: per-event source tracks. **Validated against the
  frozen list**: min-over-window b_e(t) reproduces `soho_v1` b_min to
  3 decimals on 6 sampled grazing events; ΔPA ring rotation preserves
  radius exactly. Both channels share one construction (apparent
  source at radius b_e/R☉ at the position angle of the observer's
  transverse axis offset).
- Baseline construction decision (documented, not an amendment): the
  unit series lives at solar-frame positions; baseline epochs sample
  the event's median visible radius at the *current* transverse
  azimuth (continuous with the in-window track near the window; the
  27 d solar rotation scrambles PA-stationarity anyway). The S2 star
  leaves the C3 FOV days after t_ca, so a star-fixed baseline cannot
  exist — corona-frame baselines are the only well-defined choice.
- `scripts/dev_series.py`: index → measure → reduce. One load+fit per
  frame serving all mapped positions (source + 8 ring controls);
  resumable jsonl; EXPTIME-normalized fluxes; frozen detrend (30 d
  baseline running median + k rescale) and S_stack/S_event/S_pulse
  with T = max over rings in reduce.
- Frame gates added (dev finding, header-level): synoptic frames only —
  1024², FILTER = Orange (C2) / Clear (C3), POLAR = Clear; 512²
  binned subframes and polarizer-sequence frames are excluded and
  counted (~5 % of C2 in the smoke sample).
- L1 implemented as proposed: C2 astrometric minimum = 2 stars
  (unambiguous against the predicted pattern; wrong-roll rejects at
  >35 px); ZP optional per frame (statistics are unit-free z-scores;
  NaN ZP recorded — depths come from the ensemble/stacked ZP at the
  completeness stage).
- Performance: radial-profile detrend rewritten sort-grouped
  (3 s → 0.12 s/frame); fit ~1–2.5 s/frame → full dev measure ≈ 2–4 h
  at 6 workers.
- Smoke test (724 frames fetched so far): 0 errors, 3,825 photometry
  records, window fluxes null-consistent (median 1.5 ± 7.8
  flux/exptime units). Pipeline validated; full measure + reduce +
  the L2 colour-coefficient ensemble run after the fetch completes.

## 2026-08-26 — full dev measure + statistic iterations (v1.0 → v1.2)

Full measurement pass complete: 32,843 frames (62 GB), 0 worker
errors; astrometric validity **C2 80.3 % (≥2-star, L1 rule) /
C3 99.99 %**. Frozen-gate pass rates confirm L1/L2: C2 0 % (5-star
gate unachievable), C3 12.7 % (median ZP scatter 0.37 = colour term).

**Reduce iterations** (each kept: `results/dev_search_v1.json` =
frozen v1.0 construction, `dev_search_v11.json` = v1.1/v1.2):

- v1.0 (frozen construction, literal): z-scores of 10^29 — epochs
  whose photometry annulus sits in the occulted-core *fill plateau*
  (constant → MAD 0 → infinite weight); the valid-fraction gate
  misses it because the fill is nonzero. Epoch gate added (err ≤ 0
  or ≪ series norm → dropped: physically occulted epochs).
- v1.1: A1 planet mask (frozen §4 mask, previously unimplemented —
  14 masked epochs); A2 ring-differential statistic (source minus
  ring reference on the same frame); A3 empirical night-level σ.
  Systematics collapse (ross-128 C2: z_MAD 0.49) but a +0.6 bias
  remains (window differential not centred on baseline differential
  — the ecliptic-elongated F-corona PA structure).
- v1.2: baseline-median centring; robust ring median; V ≤ 8
  star-proximity mask (3 px, KDTree); night medians. Geometry cache
  vectorized (per-call astropy ephemeris was the reduce bottleneck).

**v1.2 dev picture:** ross-128 S1 2.5 R☉ (C2 wings): quiet null
(z_med +0.04, S_stack −1.4 vs T 0.4) — the C2 grazing machinery
works. gj-908 S1: clean null. ross-154 S1: still heavy-tailed.

**Finding L3 — S2 stellar self-detection (also a positive-control
success).** gj-908 S2 shows a window-locked positive bias
(z_med +0.93, S_stack 6.9 vs T 3.2): **that is the star** —
gj-908 is V 8.98, detectable in night-median stacks; the freeze's
"empty-field" S2 assumption breaks at *stack* depth for the
brightest targets (ross-128 V 11.2 and van-maanen V 12.4 may be
marginal in full recurrence stacks). Amendment: the S2 statistic
must include the PM-propagated stellar flux as a known term
(excess-above-star, not excess-above-zero). Silver lining: the
chain just recovered a V 9 star at the predicted position,
window-locked — an in-situ end-to-end positive control.

**Finding L4 — low-galactic-latitude fields.** ross-154's field is
at b ≈ −4° (Sagittarius): catalogue-star transits through the
aperture dominate (z_MAD 1.8 even after the V ≤ 8 mask + night
medians). Every confirmatory target is high-|b|; options for
ross-154 = deeper star mask (Tycho-2/Gaia to V ~ 12) or a declared
high-background class. Decision at amendment freeze.

**Calibration note:** C2 z-units are under-dispersed (z_MAD 0.09) —
the ±110 d empirical σ over-covers short-window scatter; the
S > max(T, 0) ring comparison is unit-consistent regardless, and
per-unit z-normalisation by the ring ensemble is the candidate fix.

Next: formal amendment freeze v1.1 (L1–L4 + the statistic
construction as implemented in v1.2 + the L2 colour coefficient from
the calibrator ensemble), then the dev re-run under the amended
freeze and the dev verdict.

## 2026-08-26 — amendments v1.1/v1.2 frozen; DEV STAGE CLOSED (0 candidates)

Amendment v1.1 (thresholds.md, config sha 4abc6617…): epoch-validity
rework (astrometry-gated, ZP retired from search inclusion),
ring-differential statistic, L2 colour system measured (c_C3 = 0.429
mag/(B−V), 8,414 records, per-frame scatter 0.53 → 0.38), L3 stellar
template (V ≤ 12.5 targets), L4 high-background class, Tycho-2
VT ≤ 11 mask (sha-pinned, 860,738 stars).

Amendment v1.2 (config sha 32f0eeae…), from dev adjudication:
S2 source-star mask exemption (the mask had deleted ross-154's S2
unit — its source IS a Tycho star; gj-908 escaped only via Tycho-2's
high-PM gap), and the bright-planet in-FOV veto — ross-154 S1's
z = 20–32 events sit on the exact 8-year Venus synodic cycle
(2000/2008/2016/2024 July windows; Venus elongation 5.9–8.0°,
verified vs the SOHO ephemeris): frame-wide stray light, far beyond
the 10 px proximity mask.

**Final dev run (v1.3.2 reduce, `results/dev_search_v11.json`;
iterations v1.0 → v1.3.2 all logged, v1.0 record kept in
`dev_search_v1.json`): 11 trials, 1 exceedance vs 1.2 expected —
0 candidates.**

| unit | S_stack (S/T) | S_event (S/T) | S_pulse (S/T) | verdict |
|---|---|---|---|---|
| ross-128 S1 2.5 R☉ (C2) | −1.35/0.39 | **1.45/1.20 exc** | 16.3/49.6 | S_event = 2015-03-17 St. Patrick's Day CME storm window (2013-03-17 next-ranked, also a documented storm; a ring control shows 1.23 in the same window) → adjudicated control-crossing, retained |
| gj-908 S1 | −4.55/2.23 | 0.06/1.93 | — | null |
| gj-908 S2 | −2.20/2.71 | 1.61/2.85 | — | null; stellar self-detection → template → null = in-situ positive control |
| ross-154 S1 | −3.44/7.57 | 3.71/23.2 | — | null after the Venus veto (was 13/32 — the 8-yr family) |
| ross-154 S2 | 3.55/7.27 | 5.01/24.6 | — | null; high-background class |

Machinery validated on both cameras and channels. Next: the blind
confirmatory run (14 units) — fetch, measure, reduce ONCE under the
frozen v1.0 + v1.1 + v1.2 chain.
