---
title: "SPHEREx survey — hypothesis freeze v2.0 (decision rule and geometry specification)"
status: "v2.0 — frozen 2026-08-22 before any v2 script ran on the SPHEREx confirmatory set; hash in configs/v2_freeze.json; to be re-run when the next quick release adds a parallax phase"
date: 2026-08-22
---

# SPHEREx hypothesis freeze v2.0

The carried-forward physical-source hypothesis is reproduced directly
in §0. The decision rule, null model, completeness definition, static
template treatment and hold-out are those of the WISE freeze **v2.1**
with the SPHEREx-specific bindings below.

> **Documentation restoration (2026-08-24).** This section restores the
> still-active v1 physics text after the v1 file was retired. It changes
> no search parameter, data product, candidate decision or result. The
> `hypotheses_hash` in `configs/v2_freeze.json` remains the historical
> pre-run hash and is not a hash of this post-run documentation copy.

## 0. Physical-source hypothesis carried forward from v1

1. **Targets and endpoint models.** The frozen registry v1.5 portfolio
   contains 88 endpoint hypotheses in 77 corridors. Each endpoint is the
   registry-defined stellar component or explicitly modelled photocentre
   with fixed target-state provenance and covariance; barycentric and
   planetary endpoints are separate, out-of-scope hypotheses.
2. **Relay and role.** The source is a compact, actively station-kept
   artifact on or near the Sun's focal line for the endpoint. Receive
   (`Rx`) and transmit (`Tx`) geometries are searched separately.
3. **Distance and motion.** Heliocentric distance is 550–10,000 AU with
   a log-uniform physical prior. Residual motion is bounded
   component-wise by |µ*α| ≤ 1″/yr and |µ*δ| ≤ 1″/yr and is fitted with
   distance. The source is unresolved at the roughly 6″ SPHEREx PSF;
   extended or visibly trailed sources are different cells.
4. **Observer.** The analysis uses the Earth centre. SPHEREx is in low
   Earth orbit, so the resulting observer displacement is at most about
   17 milliarcseconds at 550 AU, negligible relative to the PSF but
   retained as a declared geometry approximation.
5. **Bands and physical interpretation.** Results are reported for six
   detector ranges: D1 0.75–1.12, D2 1.10–1.64, D3 1.62–2.42, D4
   2.40–3.82, D5 3.80–4.42 and D6 4.40–5.00 µm. The v2 search and
   injections assume an unresolved flat-Fν source—the same flux density
   in every detector. Reflected sunlight is the primary interpretation
   in D1–D3 and remains possible in D4–D6. D4–D6 also test unusually hot
   components, for which 400, 700 and 1000 K blackbodies are useful
   physical translations, plus self-luminous or nonthermal emission.
   A passive solar-equilibrium body at 550 AU is about 12 K and invisible
   throughout 0.75–5 µm, so this is not cold-thermal coverage. Each
   detector stack spans its detector-wide wavelength range; a
   single-spectral-channel line search is a separate hypothesis.
6. **Visibility in time.** Duty cycle is at least 0.5, instantiated by
   persistent, exposure-flicker, visit and long-block injections. QR2's
   two or three observing seasons limit those tests but do not broaden
   the cell to rare glints or low-duty pulses.
7. **Outside this cell.** Relay swarms, off-axis infrastructure,
   inactive or dark relics, distance or residual motion outside the
   stated bounds, extended/trailed morphologies, single-channel spectral
   lines and duty < 0.5 are not constrained. A detection would establish
   an SGL-consistent moving point source, not its emission mechanism; a
   null makes no population or network-architecture inference.

The executed decision rule was:

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
