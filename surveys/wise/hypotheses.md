---
title: "WISE/NEOWISE shakedown — baseline hypothesis freeze"
status: "v1.1 — parameters frozen 2026-08-18; endpoint set expanded 2026-08-18 (see addendum)"
date: 2026-08-18
---

> **v1.5 addendum (2026-08-20).** Batch-5 endpoint expansion from the
> universal target list **v2** (`targets/universal_v2`: horizon 10 pc,
> picky-network baskets — Engineering Backbone, Science Interest,
> desirability-filtered Selective Network Neighbor; see
> `notes/picky_network_hypothesis_sgl_seti.md`): 41 endpoints added
> (35 queued overlay_v2 systems incl. GJ 338 A/B, plus GJ 229 A,
> GJ 667 C, LTT 1445 A, GJ 66 A/B unblocked on per-component Gaia
> solutions) — **89 endpoint hypotheses across 77 corridors**, registry
> v1.5 `sha256:794d9f90…`. **Every parameter in sections 2–9 remains
> unchanged from v1.0.** Calibration pipeline v0.2.0 unchanged.
>
> **v1.4 addendum (2026-08-19).** Batch-4 + all previously deferred
> endpoints: 22 added (61 Cyg A/B, Struve 2398 A/B, Groombridge 34
> A/B, GJ 1111, Luyten's Star, Kapteyn's Star, LP 145-141, GJ 1221,
> GJ 9193, GJ 783, ε Ind Ba/Bb photocenter, GJ 11068, WISE 0855,
> EZ Aqr photocenter, Luhman 16 A/B, Procyon A/B) — **48 endpoint
> hypotheses across 38 corridors**, registry v1.4 `sha256:09366624…`.
> Every parameter in sections 2–9 remains unchanged from v1.0. The
> calibration pipeline moves to v0.2.0 (effective-epoch weight cap) —
> an estimator change, recorded in the analysis-run config, not a
> hypothesis change.
>
> **v1.3 addendum (2026-08-19).** Batch-3 endpoint expansion (WISE
> overlay mid-confusion queue): `eps-eri`, `lacaille-9352`, `gj-1061`,
> `gj-12724`, `wolf-1061` — 26 endpoint hypotheses across 23 corridors,
> registry v1.3 `sha256:a476a3e3…`. **Every parameter in sections 2–9
> remains unchanged from v1.0.**
>
> **v1.2 addendum (2026-08-19).** Batch-2 endpoint expansion from the
> universal target list (`targets/universal_v1`) per the WISE overlay
> queue: `ross-128`, `eps-ind-a`, `tau-cet`, `gj-54`, `teegarden`,
> `lacaille-8760`, `van-maanen`, `gj-908`, `gj-784` — 21 endpoint
> hypotheses across 18 corridors, registry v1.2 `sha256:f03b821f…`.
> ε Ind Ba/Bb deferred (photocenter-orbit model needed). **Every
> parameter in sections 2–9 remains unchanged from v1.0.**
>
> **v1.1 addendum (2026-08-18).** Endpoint set expanded after the pilot
> passed (plan §4.2 expansion gate): added `proxima-cen`, `wolf-359`,
> `ross-248` (linear provider) and `gj65-a`, `gj65-b` (component
> endpoints, GRAVITY 2024 orbit) — 12 endpoint hypotheses across 9
> corridors, registry v1.1 `sha256:27895a11…`. **Every parameter in
> sections 2–9 is unchanged from v1.0.** Runs pin the hypothesis
> version they executed under; v1.0-pinned records remain valid (their
> per-target source hashes are unchanged by the registry expansion).

# Baseline hypothesis freeze — v1.0

Plan §3.1 requires every analysis to freeze and identify the items below.
Parameters were decided 2026-08-18. This version is frozen: run configs
record this file's version and content hash; any change creates a new
hypothesis version, never an edit. Per-endpoint solution values and their
provenance are bound separately by the registry version pinned in each run
config.

## 1. Target endpoints and target-state models

Full pilot shortlist, five systems (see
`notes/pilot_registry_shortlist.md`): Barnard's Star, Alpha Cen AB,
Sirius AB, Ross 154, Lalande 21185.

For the multiple systems the endpoint hypotheses are the **stellar
components** (Alpha Cen A, Alpha Cen B, Sirius A, Sirius B), each under
the orbital target-state provider — **not** the barycenter. Barycenter
endpoints are a separate model ID, out of scope for this freeze. That
gives 7 endpoint hypotheses across 5 systems.

Each registry entry uses the best available solution with covariance and
per-value provenance (not Gaia DR3 by default), fixed at registry
curation.

## 2. Role

Rx and Tx evaluated as **separate hypotheses** for every endpoint, per
plan §3.1. Both roles run in the pilot; they share frame discovery, so
the marginal cost is in the precise pass only.

## 3. Relay-distance prior and sampling

- Range: 550–10,000 AU (heliocentric), per plan §3.1.
- Prior: **log-uniform** over the range (governs injection placement and
  completeness statements; discovery geometry is prior-independent).
- Sampling: sglseti tolerance-guaranteed adaptive locus sampling, not a
  fixed grid; covered relay-distance intervals reported exactly.

## 4. Geometry model and observer-state versions

- Geometry: sglseti `tusay2022_eq5_7_v1` role model.
- Observer: Earth center for WISE. **Declared approximation:** WISE is in
  a ~500 km LEO; Earth-center parallax error at 550 AU is bounded by
  ~R_orbit/550 AU ≈ 2.5 milliarcsec — negligible against the search
  padding, but carried in the accuracy budget.
- Pinned sglseti commit, kernels, and IERS per run config.

## 5. Uncertainty confidence level and search padding

- Target-state covariance propagated via sglseti's seeded Monte Carlo;
  loci reported at **99% confidence**.
- Additional fixed search padding: **+10 arcsec**, recorded separately
  from propagated uncertainty per plan §3.1. Revisit (as a new hypothesis
  version) after the pilot measures actual astrometric residuals.

## 6. Source morphology and spectral model

- Morphology: unresolved point source at WISE resolution (6"–12").
- Spectral interpretation is **band-specific** (plan §4): W1/W2 (3.4/4.6
  µm) test reflected sunlight, unusually hot components, or nonthermal
  emission; W3/W4 (12/22 µm, 4-band cryo only) carry the waste-heat
  sensitivity for warm (~300 K) radiators. A 550 AU solar-equilibrium
  blackbody (~12 K, peak ~240 µm) is invisible to all WISE bands — flux
  limits must not be described as generic thermal coverage.

## 7. Persistence / duty-cycle assumption

- Pilot cell: **duty cycle ≥ 0.5** (includes fully persistent sources).
  Injections sample duty cycle over [0.5, 1].
- Lower duty cycles and structured intermittency are separate grid cells
  for later analyses.

## 8. Stationkeeping / residual-motion bounds

- Baseline: station-kept relay on the focal line with bounded residual
  motion, fit jointly with relay distance: **|µ_resid| ≤ 1 arcsec/yr**,
  recorded as an explicit assumption. (Plan §5 keeps the physical
  stationkeeping prior open research; this bound defines the pilot cell.)

## 9. Detection pipeline and decision threshold

Layered per plan §3.4; specific thresholds (catalog-screen radius, forced
photometry S/N, track-fit acceptance, held-out-epoch prediction test) are
frozen in the run config once the pipeline exists. Empirical false-alarm
rates must account for all targets × roles × distances × residual-motion ×
signal-shape trials.

## Out of scope for this freeze

Barycenter endpoints, planetary endpoints, relay swarms, off-axis
infrastructure, inactive relics, duty cycle < 0.5 — separate model IDs
and separate freezes (plan §3.1).
