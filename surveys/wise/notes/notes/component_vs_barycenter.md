---
title: "Component vs barycenter endpoints: measured locus separations"
date: 2026-08-18
status: "sglseti experiment, 2026-08-18; supports the component-endpoint decision"
---

# Component vs barycenter: how different are the focal lines?

Question: before committing to component endpoints for the multiples, is
the component-vs-barycenter distinction observationally meaningful for
the WISE search, or would a single barycenter hypothesis cover both?

## Method

Using sglseti (`evaluate_locus`, `tusay2022_eq5_7_v1`, Earth-center
observer), evaluate the Rx locus for each component endpoint from
`registries/pilot_wise_2026.yaml` and for a barycenter variant sharing
the identical barycenter astrometry and orbit block (differing only in
`component`/`endpoint_kind`), at yearly epochs 2010-04 → 2024-04 and
z ∈ {550, 1500, 10000} AU; report on-sky separations. Because both
targets in each pair share the same astrometry block, the comparison is
differential and insensitive to the registry's open CURATION items
(α Cen barycenter PM revalidation, Sirius photocenter approximation).

## Results (arcsec, Rx role)

| Pair | min | max | drift 2010→2024 | shape |
| --- | --- | --- | --- | --- |
| α Cen A vs AB barycenter | 1.86 | 3.81 | 1.94 | minimum in 2016, rising |
| α Cen B vs AB barycenter | 2.17 | 4.44 | 2.27 | minimum in 2016, rising |
| Sirius A vs AB barycenter | 2.93 | 3.74 | 0.82 | slow rise, flattening |
| Sirius B vs AB barycenter | 5.93 | 7.59 | 1.65 | slow rise (near apastron) |

z-dependence is negligible (< 0.1" between 550 and 10,000 AU), as
expected: the relay-direction offset inherits the component's angular
offset from the barycenter almost unchanged.

Context scales:

- Propagated 99% locus uncertainty (`propagate_locus_uncertainty`,
  256 samples): ≤ 0.1" for Sirius A, ≪ 0.01" for α Cen A and Barnard's
  Star. The endpoint offsets are 20–1000× larger — a genuine model
  choice, not something uncertainty propagation absorbs.
- W1/W2 PSF FWHM ≈ 6.1"/6.4" (σ ≈ 2.6"), pixels 2.75"; frozen search
  padding 10"; frames 47'.
- For a Gaussian W1 PSF, forced photometry centered at the wrong track
  retains ~74% of peak flux at 2" offset, ~51% at 3", ~30% at 4", and
  ~1% at 7.6".

## Search-perspective reading

1. **Discovery (coarse pass): no difference.** Offsets of arcseconds are
   invisible to 47' frame selection; both hypotheses select the same
   frames, so component endpoints add no discovery or download cost.
2. **Forced photometry / shift-and-stack: decisive.** A barycenter-track
   search would lose ~25–70% of a component-true source's peak flux at
   α Cen / Sirius A offsets, and ~99% at Sirius B's 6–7.6" — i.e. it
   would simply miss a faint relay on the Sirius B focal line. The
   hypotheses are not interchangeable at W1/W2 resolution. (W4's ~12"
   PSF blurs the distinction, but W4 has only the 2010 cryo epochs.)
3. **Residual-motion fitting cannot rescue the wrong endpoint cleanly.**
   The wrong-endpoint offset drifts at only ~0.06–0.16"/yr — well inside
   the frozen |µ_resid| ≤ 1"/yr bound — so a component-true source found
   in a barycenter search would be partially absorbed as "stationkeeping
   residual," corrupting the model comparison rather than failing
   loudly. Running the correct endpoint hypotheses explicitly keeps the
   residual-motion prior meaningful.
4. **A and B are also distinct from each other** (their loci differ by
   the sum of the offsets, ~4–12") — searching one component does not
   cover the other.

## Conclusion

The component-endpoint decision stands and matters: separations are a
substantial fraction of (to larger than) the W1 PSF, far exceed
propagated astrometric uncertainty, and change the forced-photometry
answer. Barycenter endpoints remain a separate, explicitly-labeled model
ID (plan §3.1) rather than a substitute; per plan §5, component
astrometry must never be used as an unlabeled barycenter stand-in — this
experiment quantifies why.

Script: scratchpad `component_vs_barycenter.py` (session artifact); the
barycenter variants used are recorded inline there and mirror the main
registry values exactly.
