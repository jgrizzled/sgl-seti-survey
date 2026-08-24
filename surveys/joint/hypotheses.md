---
title: "Joint PS1 + ZTF stage — hypothesis freeze v2.0"
status: "v2.0 — frozen 2026-08-22; the joint split is the PS1/ZTF split (identical by construction); hash in configs/v2_freeze.json"
date: 2026-08-22
---

# Joint PS1 + ZTF hypothesis freeze v2.0

The common PS1/ZTF physical-source hypothesis is reproduced directly
below so that this freeze is meaningful without another retired file.

> **Documentation restoration (2026-08-24).** This is a post-run,
> meaning-preserving copy of the physical assumptions already used by
> both archive freezes. It changes no search parameter, data product,
> candidate decision or result; the frozen configuration remains the
> provenance authority for the executed run.

## Physical-source hypothesis

- **Targets and geometry:** the common PS1/ZTF subset contains 69
  registry v1.5 endpoint hypotheses in 62 corridors. Registry-defined
  stellar components or explicitly modelled photocentres are used;
  barycentric and planetary endpoints are separate hypotheses. `Rx` and
  `Tx` are searched separately for a compact, station-kept relay on the
  Sun's focal line at 550–10,000 AU, with a log-uniform distance prior
  and component-wise residual motion |µ*α|, |µ*δ| ≤ 1″/yr.
- **Source and spectrum:** the source is unresolved in both archives.
  The joint g, r and i cells use a flat-Fν AB spectrum and pair PS1
  g/r/i with ZTF zg/zr/zi. Reflected sunlight is the primary physical
  interpretation; self-luminous optical emission, a beacon, or
  sufficiently broad and recurrent leakage are also detectable as flux.
  A reflected-light size constraint must state the assumed albedo. A
  passive roughly 12 K body at 550 AU is invisible in these bands, so
  the joint stage is not thermal coverage.
- **Time behavior:** duty cycle is at least 0.5 under the persistent,
  exposure-flicker, visit and long-block families. Rare glints and
  low-duty pulses are outside the cell.
- **Scope:** swarms, off-axis infrastructure, inactive/dark relics,
  extended or trailed sources, distances or motions outside the stated
  bounds and different endpoint hypotheses require separate searches.
  The result is a targeted flux constraint, not an inference about relay
  prevalence or network architecture.

Operationally, the common reference epoch is MJD 59800, the AB zero
point is 25, and a joint cell is endpoint × role × band pair (g = PS1 g
+ ZTF zg; r; i).

1. **Statistic.** S_joint(z, µ) = (A_ztf + A_ps1) / √(B_ztf + B_ps1)
   from the per-archive v2 accumulated stack sums (per trajectory,
   per mask, capped and clipped per archive), the ZTF sums interpolated
   in 1/z onto the PS1 grid (360 nodes). Trajectory index t is the same
   sky offset in both archives (identical 48-offset ring), so the joint
   ring null is exact; donors are annotations.
2. **Candidate rule.** R = S_max/T over the 8 designated controls,
   R̃ = R/q95(ring), heavy-tail / radius-dependent cells void,
   **one family-wise threshold across archives**: R̃_FWER at α = 0.05
   over the set's joint cells.
3. **Calibrated veto.** Flux-consistent catalogued static source: the
   PS1 DR2 mean and ZTF DR objects pushed through each archive's
   response at the joint node; rejection only if together they account
   for the joint S within a factor 2. Phase split (the joint stage is
   where PS1's single phase meets ZTF's two), p_phase and p_trajectory
   are annotations.
4. **Completeness.** The j-th injection of a cell is the same physical
   source in both archives (per-injection seeded draws on the union
   magnitude window, flat-Fν AB): the per-archive injected window sums
   are combined as in (1) and classified with the joint rule;
   threshold and final-candidate curves per temporal model.
5. **Hold-out.** The PS1/ZTF corridor split (seed 20260822);
   confirmatory set analysed once.
