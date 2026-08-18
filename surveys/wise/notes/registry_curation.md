---
title: "Pilot registry curation — resolution of the three CURATION flags"
date: 2026-08-18
status: "Resolved; registry v1.0 hash sha256:82743c0976c2…"
---

# Registry curation record (2026-08-18)

How each CURATION flag from the draft registry was resolved. Decision:
revalidate copied data against primary sources; reasonable documented
assumptions for the rest.

## 1. Alpha Cen barycenter astrometry — revalidated (replaced)

The draft carried the sglseti example's "demonstration" barycenter
snapshot. Revalidated against **Kervella et al. 2016, A&A 594, A107**
(ar5iv full text): the demo PM values were their rounded Table 1 values,
but the parallax was not (743 mas from Pourbaix & Boffin vs Kervella's
refined **747.17 ± 0.61 mas** — a 4 mas difference). Adopted the fully
coherent Kervella solution: barycenter position from their Table 2
(14:39:31.813, −60:50:00.01 at epoch 2009.1517 → 219.88255417,
−60.83333611 deg; 30 mas conservative uncertainty at printed precision),
PM (−3619.9 ± 3.9, +693.8 ± 3.9), RV −22.3930 ± 0.0043, and their orbit
(P 79.929, T₀ 1955.604, e 0.5208, a 17.592″, i 79.320°, Ω 205.064°,
ω 232.006°, m_B/m_tot 0.45884 ± 0.00027), replacing the Pourbaix &
Boffin mix. Orbit-element uncertainties entered in provenance so Monte
Carlo can sample them.

## 2. Sirius photocenter-as-barycenter — investigated and validated (kept)

The draft flagged that hip2 HIP 32349 tracks the photocenter (≈ A), so
using it as barycenter astrometry looked like an approximation needing
correction. Investigation showed the opposite:

- A closure fit through sglseti (subtracting A's orbital offset/velocity
  at 1991.25 per the Bond orbit) implied a barycenter PM correction of
  ~(313, 340) mas/yr — yet the century-baseline FK5 proper motion of
  Sirius agrees with the hip2 values at the mas/yr level, which it could
  not if the catalog PM contained A's ~460 mas/yr orbital velocity.
- Discrimination test against the 2MASS position of Sirius A
  (JD 2451541.62 ≈ J2000.0, 8.75 yr from the Hipparcos epoch):

  | Interpretation of hip2 values | Predicted-A miss vs 2MASS |
  | --- | --- |
  | **barycenter (+ Bond orbit → A)** | **0.16″** |
  | instantaneous photocenter (linear) | 1.44″ |
  | closure-corrected barycenter (+ orbit) | 3.62″ |

Conclusion: van Leeuwen's Sirius solution is already orbit-corrected —
the catalog values ARE the barycenter. Adopted as-is, with PM
uncertainties inflated to 20 mas/yr per axis (the systematic allowance
implied by 0.16″ over 8.75 yr), Bond et al. 2017's adopted parallax
378.90 ± 1.27 mas (hip2's 379.21 ± 1.58 agrees at 0.2σ), and Bond's
gravitational-redshift-corrected systemic RV −8.47 km/s (replacing the
Gontcharov −5.5, which is an A-component spectroscopic value). The Bond
Table 4 orbit itself was confirmed against the ar5iv full text — the
paper's single-column new solution matches sglseti's validated fixture
exactly, now with published uncertainties recorded in provenance.

Scripts: scratchpad `sirius_barycenter.py`, `sirius_discriminate.py`.

## 3. Correlation matrices — diagonal by declared assumption

ESA's Gaia archive is unreachable from the dev network and VizieR's
I/355 export omits correlation columns, so full covariance matrices are
not available without manual effort. Adopted assumption: **diagonal
covariance** (per-value 1σ in provenance; sglseti then samples
independent Gaussians). Justification: propagated 99% locus radii are
≤ 0.1″ for every pilot endpoint against a 10″ frozen search padding —
even fully correlated errors could not move any discovery or usability
decision at pilot scale. Revisit only if a constraint ever becomes
padding-limited.

## Downstream effects

- Registry hash changed: runs pin
  `sha256:82743c0976c2f7d8f3a994e56f117380f0454691c903807d8c0ee9c0e610969d`.
- The component-vs-barycenter experiment (`component_vs_barycenter.md`)
  used the draft α Cen orbit; its conclusions are scale-level and
  unaffected by the Kervella refinements (same orbit to ~1%).
