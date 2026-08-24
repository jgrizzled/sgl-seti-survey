---
title: "Pan-STARRS1 survey — hypothesis freeze v2.0 (decision rule and geometry specification)"
status: "v2.0 — frozen 2026-08-22 before any v2 script ran on the PS1 confirmatory set; hash in configs/v2_freeze.json"
date: 2026-08-22
---

# PS1 hypothesis freeze v2.0

The physics cell of `surveys/panstarrs/hypotheses.md` v1.0 is unchanged
(69 endpoints / 62 corridors with v1 tensors; Rx and Tx; 550–10,000 AU
log-uniform; unresolved point source; duty ≥ 0.5; component-wise
|µ| ≤ 1″/yr; Haleakalā observer). The decision rule, null model,
completeness definition and hold-out are those of the WISE freeze
**v2.1** (`surveys/wise-v2/hypotheses.md` §3–§5, §8) through the
engine `sglsurvey/v2/` with the PS1 bindings in `profile.py`:

1. **Statistic and null.** S(z, µ) over 360 × 5 × 5 nodes (uniform in
   1/z, 1.0″ spacing, T0 = MJD 59800 — the joint-stage common epoch),
   single-epoch clip |S_e| ≤ 5, per-frame weight cap at 20 × the band's
   median frame weight, skycell duplicates collapsed, per-warp zero
   point = the v1 star calibration (DR2 mean-table stars through the
   same matched filter; never the header FPA.ZP). 48-offset ring
   (20/30/40″), R = S_max/T over the 8 designated controls, R̃ = R/q95;
   phase scramble and trajectory randomisation as annotations;
   heavy-tail / radius-dependent cells void.
2. **Candidate rule.** R̃ ≥ R̃_FWER (α = 0.05, 10,000 pseudo-experiments
   over the family), computed from the family's null before any real
   R̃ is examined; BH q-values and rank statements reported.
3. **Calibrated veto.** The flux-consistent catalogued static-source /
   halo test against the PS1 DR2 mean objects (≥ 3 detections) of the
   v1 screening snapshots. The PS1 mission is over (2009–2014): no
   epoch hold-out; the 97:3 parallax-phase split of the 3π cadence
   makes the phase annotation uninformative — the joint stage supplies
   the second phase.
4. **Quality masks.** Primary: `badflag = 0`; strict: seeing ≤ 1.5″;
   loose: all usable warps.
5. **Completeness.** 400 image-level injections per cell, Moffat(β = 3)
   at the warp's CHIP.SEEING, magnitude uniform on the union of this
   archive's and the ZTF partner band's v1 m90 ± 2 (the j-th injection
   of a cell is the same physical source in both archives for the
   joint stage), flat-Fν AB; four temporal models (visit-scale ≈
   persistent for the single-phase cadence).
6. **Geometry.** Monte Carlo envelopes (N = 2,000, seed 20260822,
   2010 / 2012 / 2014); cross-track dimension where σ_xt,99 > 0.5 × 1.0″.
7. **Hold-out.** Random corridor split stratified by the WISE
   confusion class, seed 20260822; confirmatory set analysed once.
8. **Scope.** Targeted coverage of the 69-endpoint PS1 subset; no
   population inference.
