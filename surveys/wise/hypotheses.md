---
title: "WISE/NEOWISE shakedown — baseline hypothesis freeze"
status: "v1.0 — parameters frozen 2026-08-18"
date: 2026-08-18
---

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
