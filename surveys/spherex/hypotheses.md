---
title: "SPHEREx survey — hypothesis freeze v2.0 (decision rule and geometry specification)"
status: "v2.0 — frozen 2026-08-22 before any v2 script ran on the SPHEREx confirmatory set; hash in configs/v2_freeze.json; to be re-run when the next quick release adds a parallax phase"
date: 2026-08-22
---

# SPHEREx hypothesis freeze v2.0

The physics cell of the v1.0 freeze (`hypotheses_v1.md`, + v3
amendments) is unchanged (88 endpoints / 77 corridors; Rx and Tx;
550–10,000 AU log-uniform; unresolved; duty ≥ 0.5; component-wise
|µ| ≤ 1″/yr; Earth-centre observer, ≤ 17 mas). The decision rule, null
model, completeness definition and hold-out are those of the WISE
freeze **v2.1** through `sglsurvey/v2/` with the SPHEREx bindings in
`profile.py`:

1. **Statistic and null.** Cells per detector D1–D6; S(z, µ) over
   96 × 3 × 3 nodes (uniform in 1/z, 3.7″; µ = −1/0/+1″/yr; T0 = MJD
   61000), search image IMAGE − ZODI with the exposure's PSF plane as
   the kernel at 2 × 2 sub-pixel phases, the v3 static-sky template
   subtracted at sampling with its reduced-χ² variance scaling,
   single-epoch clip |S_e| ≤ 5, per-frame cap. 48-offset ring (20/30/40″
   = 3.3–6.6 FWHM), R = S_max/T (8 designated), R̃ = R/q95; phase
   scramble and trajectory randomisation as annotations; heavy-tail /
   radius-dependent cells void. Full (SIP) WCS sampling.
2. **Candidate rule.** R̃ ≥ R̃_FWER (α = 0.05 over the family).
3. **Rejection tests.** None from catalogues: the template is the
   calibrated static-sky treatment (a catalogue flux-consistency test
   on template-subtracted fluxes double-counts). A candidate surviving
   the rule is `retained-ambiguous` with its annotations (phase split,
   p_phase, p_trajectory, nearest 2MASS / CatWISE source, template
   coverage). No epoch hold-out (one quick release).
4. **Quality masks.** Primary = loose = all usable exposures; strict =
   excluding the deep-field collection.
5. **Completeness.** 400 image-level injections per cell with the
   exposure's own PSF plane, flat-Fν AB (the same µJy in every
   detector — the declared SED), magnitude uniform on m90_v1 ± 2,
   continuous z / µ / cross-track, four temporal models (visit-scale ≈
   the 2–3 observing seasons). Caveat: template absorption of slow real
   sources is not modelled.
6. **Geometry.** Monte Carlo envelopes (N = 2,000, seed 20260822,
   2025.5 / 2026.0 / 2026.5); cross-track dimension where
   σ_xt,99 > 0.5 × 6.0″.
7. **Hold-out.** Random corridor split stratified by the WISE
   confusion class, seed 20260822; confirmatory set analysed once.
8. **Scope.** Targeted coverage of the frozen 88-endpoint portfolio in
   QR2; no population inference.
