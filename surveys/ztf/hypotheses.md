---
title: "ZTF pilot — baseline hypothesis freeze"
status: "v1.0 — frozen 2026-08-20 (3-corridor pilot)"
date: 2026-08-20
---

# Baseline hypothesis freeze — ZTF v1.0

Second-adapter pilot (plan §5). Parameters were inherited from the WISE
freeze (`surveys/wise/hypotheses.md` v1.0) wherever the physics is
archive-independent, so that a ZTF null and a WISE null on the same
endpoint close the *same* cell of the §3.1 grid in different bands.
Only the items marked **ZTF-specific** are new. Run configs record this
file's version and content hash; any change creates a new version.

## 1. Target endpoints and target-state models

Pilot corridors (3), reusing registry `registries/pilot_wise_2026.yaml`
entries **unchanged** (the registry is survey-agnostic):

| endpoint | why it is in the pilot |
|---|---|
| `ross-128` | corridor on the ecliptic (β ≈ +0.5°): maximal parallax excursion along one axis; high Galactic latitude (b −60°), clean field; Dec −0.8° |
| `eps-ind-a` | southern star whose corridor is at Dec +57°, β +41°: near-circular parallax ellipse; demonstrates that ZTF reaches southern-hemisphere endpoints; mid-latitude field (b +48°) |
| `proxima-cen` | flagship target; corridor at b +1.9°, Dec +63°: the crowded Galactic-plane stress case for a 1" survey |

Endpoint hypotheses are stellar components under their registry
providers (all three are linear-astrometry endpoints). Component
endpoints of the multiple systems follow in the scale-up, as in WISE.

## 2. Role

Rx and Tx as separate hypotheses, both evaluated (unchanged).

## 3. Relay-distance prior and sampling

550–10,000 AU, log-uniform prior, sglseti adaptive-locus sampling with
exact covered-interval reporting (unchanged).

## 4. Geometry model and observer-state versions — ZTF-specific observer

- Geometry: sglseti `tusay2022_eq5_7_v1` (unchanged).
- Observer: **Palomar P48 terrestrial site** (lon −116.8650°, lat
  +33.3563°, h 1712 m) via `Observer.from_geodetic`. Measured
  topocentric-vs-geocentric locus shift ≤ 0.016" at 550 AU; carried in
  the accuracy budget, negligible against padding.
- Epoch: `obsjd` is exposure start; loci are evaluated at mid-exposure
  (start + 15 s). Corridor drift over a 30 s exposure ≤ 0.003" — the
  midpoint model is adequate (plan §3.3).

## 5. Uncertainty confidence level and search padding

99% propagated locus confidence + **10" fixed padding** (unchanged).
At 1"/pix this is 10 pixels; revisit downward as a new version once
the pilot measures ZTF astrometric residuals along the track.

## 6. Source morphology and spectral model — ZTF-specific interpretation

- Morphology: unresolved point source at ZTF resolution (median seeing
  ~2"). Extended or trailed morphologies are separate cells.
- Bands: ZTF g (0.47 µm), r (0.64 µm), i (0.79 µm); all public frames
  used, constraints reported per filter.
- Physical interpretation of an optical flux limit:
  1. **Reflected sunlight** — the primary cell. For a Lambertian
     sphere of geometric albedo p and diameter D at heliocentric
     distance z (AU), observed near opposition at Δ ≈ z:
     m ≈ m_⊙ − 2.5 log10[ p (D/2)² / (4 z⁴ AU²) ] with m_⊙,r = −26.93.
     A 5σ stack depth m_r converts directly to an upper limit on
     p^{1/2} D at each z. Limits must quote p explicitly.
  2. **Self-luminous optical emission** (beacons, leakage, thermal
     > 2000 K) — reported as a flux limit in Jy without a size
     conversion.
  A 550 AU solar-equilibrium body (12 K) is invisible in g/r/i; these
  limits are never to be described as thermal coverage.

## 7. Persistence / duty-cycle assumption

Duty cycle ≥ 0.5; injections sample [0.5, 1] (unchanged).

## 8. Stationkeeping / residual-motion bounds

|µ_resid| ≤ 1"/yr, fit jointly with relay distance (unchanged). Note
that at ZTF resolution this bound spans ±8 pixels over the 8-year
baseline, so residual motion is a *resolved* nuisance parameter here,
unlike in WISE.

## 9. Detection pipeline and decision threshold — ZTF-specific inputs

Layered per plan §3.4. Frozen data-quality inputs:

- Frame selection: `ipac_gid = 1` (public), `infobits` bit 25 clear
  (ZSDS bad-quality flag), `imgtype = object`. Airmass, seeing, moon
  illumination and `maglimit` are recorded per frame and enter the
  stack as weights, not as hard cuts.
- Usable pixels: mask template 6141 — bits {0,2,3,4,5,6,7,8,9,10,12}
  fatal; bits 1 and 11 (source present) **not** fatal.
- Primary search image: `scimrefdiffimg` (PSF-matched difference),
  which removes the static field; `sciimg` fallback where no
  difference image exists, with a local background model and the
  static-field veto delegated to the parallax-phase test.
- Screening catalogs: per-frame `psfcat`, DR24 objects table. Catalog
  absence is neither detection nor null (plan §3.4).
- Parallax-phase test, offset-trajectory controls (8), effective-epoch
  weight cap (20× median) and 1/z-uniform z-grid carried over from WISE
  calibration v0.2.0.
- Positive control: a numbered main-belt asteroid recovered through the
  shift-and-stack code along its JPL ephemeris — new in this pilot.

Specific thresholds are frozen in the run config when each stage is
first executed.

## Out of scope for this freeze

Barycenter/planetary endpoints, swarms, off-axis infrastructure,
inactive relics, duty cycle < 0.5, trailed sources, proprietary
(gid 2/3) frames.
