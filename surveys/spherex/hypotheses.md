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

# Amendment v2.1 — six-detector joint cell (frozen 2026-09-04)

> Added after the v2.0 confirmatory run and before any joint-cell
> statistic was evaluated on the confirmatory set. The v2.0 rule, its
> freeze (`configs/v2_freeze.json`) and its results are unchanged; this
> amendment adds one dependent family on the same data. Hash in
> `configs/joint6_freeze.json`. Implementation `surveys/spherex/joint6.py`.

1. **Statistic.** One cell per endpoint × role (band label `J6`):
   S_J(z, µ) = Σ_b A_b / √(Σ_b B_b) over D1–D6 from the stored v2
   per-trajectory accumulators (per-frame cap, |S_e| ≤ 5 clip and the
   v3 template subtraction inherited per detector), on the common
   96 × 3 × 3 grid — the declared flat-Fν SED gives equal weight to every
   detector. n_J = Σ_b n_b ≥ 5; a cell exists when ≥ 2 detectors have
   n_epochs_ok ≥ 5. Cross-track variants are combined per trajectory by
   the larger S_max (the v2.0 rule).
2. **Null and candidate rule.** As v2.0: 48-offset ring, R = S_max/T
   (8 designated), R̃ = R/q95, heavy-tail / inner–outer KS void flags;
   phase scramble and trajectory donors as annotations. R̃ ≥ R̃_FWER at
   α = 0.05 over the joint family of the set — a second, dependent
   family on the same data (the PS1 + ZTF joint-stage precedent); the
   survey-wide FWER across the per-detector and joint families is ≤ 0.10.
3. **Rejection tests.** None from catalogues (v2.0). A surviving cell is
   `retained-ambiguous` with annotations: per-detector S and fitted flux
   at the joint node with a flat-Fν χ² (annotation, not a veto), phase
   split, p_phase, p_trajectory, nearest 2MASS / CatWISE source.
4. **Completeness.** The v2 injection chain re-run once with a common
   magnitude window per pair — the v1 six-detector joint m90
   (`calib_v4/m90_curves.npz`, `ALL`, median over z) ± 2 — so the j-th
   injection is one physical source in all six detectors
   (`runs/spherex/v2/injections_joint`); 400 per pair; the six window
   sums combined exactly as the accumulators; bright limit = the
   brightest of the six single-epoch clip magnitudes. The visit / block
   on-patterns are drawn per detector (an inherited property of the
   chain, recorded as a caveat).
5. **Hold-out.** The v2.0 corridor split inherited unchanged; dev set
   first under all three masks, confirmatory once under primary.
6. **Template-absorption control (not a search rule).** Measured
   alongside: for a ladder of persistent sources on the real track (8
   z-interval centres × 3 magnitudes at µ = 0; µ = ±1″/yr at the two
   most distant z), the template is refitted with the source present and
   the absorbed fraction f_abs at the track is recorded; Δm = −2.5
   log10(1 − f_abs) is reported per cell and z interval as a correction
   to the injection-calibrated m90 (report-level; frozen records
   unchanged but flagged). `surveys/spherex/scripts/template_absorb.py`.

**Dev-driven amendment (a), 2026-09-04, before the confirmatory run.**
The dev nulls exposed a degenerate cell: proxima-cen/rx lies in an
over-subtracted Galactic-plane field where every trajectory's joint
S_max is negative, so T = 0.057 and the ring q95 = −12.3; dividing by a
negative normaliser turned the most negative ring values into the
largest R̃ of the family and set R̃_FWER = 5.26 for all 52 usable cells.
Rule added: a cell whose ring normaliser is not a positive scale
(q95 ≤ 0) is void (`null_degenerate`); the heavy-tail test is evaluated
only for q95 > 0. The v2.0 per-detector engine has no such guard (its
degenerate cells happened to be caught by the heavy-tail rule); noted
for the QR3 re-run. Freeze re-hashed with the amendment recorded
(`supersedes_freeze` keeps the pre-amendment hash).
