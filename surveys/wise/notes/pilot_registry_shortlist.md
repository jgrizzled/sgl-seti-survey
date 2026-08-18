---
title: "Pilot registry shortlist — WISE shakedown"
date: 2026-08-18
status: "Selection decided 2026-08-18; final registry goes to registries/ with full provenance"
---

# Pilot registry shortlist

Plan §4.2 wants 3–5 deliberately diverse endpoints: an isolated
linear-motion case, a high-proper-motion case, and at least one bright or
multiple system using orbital target-state support — each with the best
available solution (not Gaia DR3 by default), covariance, and per-value
provenance.

The SGL corridor for a target at (α, δ) is centered on the antipode
(α+12h, −δ); antipode coordinates below are exact arithmetic, while
ecliptic/galactic context and WISE depth-of-coverage per corridor should
be computed during curation, not assumed.

| Candidate | Role in the pilot | Antipode (approx) | Notes |
| --- | --- | --- | --- |
| **Barnard's Star (GJ 699)** | High-proper-motion case (~10.4"/yr) | 05h58m, −04°41′ | Isolated M dwarf, superb astrometry; stresses the secular-drift term of the SGL track like nothing else. |
| **Alpha Cen AB** | Bright multiple, orbital target-state model | 02h40m, +60°50′ | Nearest system; orbit solution validated in sglseti; A/B too bright for clean Gaia DR3, so published orbit + barycenter solution is the natural source. Antipode sits near the Galactic plane (b ≈ +0.7°) — a deliberately hard, crowded-background corridor. Prior SGL search art (Tusay et al. 2022) gives a cross-check geometry. |
| **Sirius AB** | Second multiple with validated orbit | 18h45m, +16°43′ | Exercises the orbital provider on a very different mass ratio and period; extremely bright primary is also a latent-artifact stress test if the corridor lands near bright-star ghosts. |
| **Ross 154 (GJ 729)** | Isolated ordinary linear-motion case | 06h50m, +23°50′ | Moderate PM (~0.7"/yr) M dwarf — the "boring" control-like endpoint the grid needs. |
| **Lalande 21185 (GJ 411)** | Alternate/5th: bright isolated M dwarf, high PM (~4.8"/yr) | 23h03m, −35°58′ | Middle ground between Barnard's and Ross 154; excellent solutions available. |

## Selection considerations to check during curation

- **WISE depth-of-coverage** at each antipodal corridor: the scan pattern
  gives many more epochs at high ecliptic latitude; compute actual epoch
  counts from the frame inventory rather than assuming.
- **Corridor background**: galactic latitude, source density, and
  bright-star artifacts at the antipode (matters for the matched control
  corridors in plan §3.4 too).
- **Endpoint hypothesis per system** (plan §5): for Alpha Cen and Sirius,
  component vs. barycenter endpoints are distinct model IDs — decide which
  are in the pilot.

## Decisions (2026-08-18)

1. **Full shortlist selected** — all five systems above.
2. **Component endpoints for the multiples** (Alpha Cen A, Alpha Cen B,
   Sirius A, Sirius B), not barycenters → 7 endpoint hypotheses across 5
   systems. Per-endpoint solution sources are fixed at registry curation
   with per-value provenance.

## Curation mechanics

sglseti's registry format (docs/registry.md in `../sglseti`) is the
target; its validated Alpha Cen AB and Sirius AB orbit fixtures/examples
should seed those entries. Solution values must be entered with per-value
provenance and covariance — fetchable from Gaia DR3 TAP where appropriate
and from published orbit papers otherwise.
