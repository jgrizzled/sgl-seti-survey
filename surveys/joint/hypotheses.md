---
title: "Joint PS1 + ZTF stage — hypothesis freeze v2.0"
status: "v2.0 — frozen 2026-08-22; the joint split is the PS1/ZTF split (identical by construction); hash in configs/v2_freeze.json"
date: 2026-08-22
---

# Joint PS1 + ZTF hypothesis freeze v2.0

Physics cell as the PS1 and ZTF v2.0 freezes (common T0 = MJD 59800,
AB zero point 25, 69 endpoints / 62 corridors). Joint cell = endpoint
× role × band pair (g = PS1 g + ZTF zg; r; i).

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
