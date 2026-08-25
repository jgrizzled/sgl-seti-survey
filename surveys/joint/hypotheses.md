---
title: "Joint Pipeline A stage — hypothesis freeze v3.0 (PS1 + ZTF stack, WISE colour axis)"
status: "v3.0 — frozen 2026-08-24 before any v3 script touched the confirmatory set; optical candidate statistic unchanged from v2.0; hash in configs/v3_freeze.json (v2.0 kept below as the superseded record, hash in configs/v2_freeze.json)"
date: 2026-08-24
---

# Joint hypothesis freeze v3.0 — PS1 + ZTF + WISE

The three-archive extension of the v2 joint stage (plan §4 open items;
`notes/learnings.md` §10 item 3). The optical candidate-generating
statistic and its family are **unchanged from v2.0**; v3 adds the WISE
W1/W2 colour axis as a calibrated veto and per-cell annotations. This
document is self-contained (learnings §9: the freeze hashes the
complete active hypothesis document).

## 1. Physical-source hypothesis

- **Targets and geometry:** the common PS1/ZTF subset — 69 registry
  v1.5 endpoint hypotheses in 62 corridors (WISE, all-sky, covers every
  one). Registry-defined stellar components or explicitly modelled
  photocentres; barycentric and planetary endpoints are separate
  hypotheses. `Rx` and `Tx` are searched separately for a compact,
  station-kept relay on the Sun's focal line at 550–10,000 AU, with a
  log-uniform distance prior and component-wise residual motion
  |µ*α|, |µ*δ| ≤ 1″/yr about the common reference epoch T0 = MJD 59800.
- **Source and spectrum:** unresolved in all three archives. The frozen
  spectrum is **flat Fν in AB** across PS1 g/r/i, ZTF zg/zr/zi and WISE
  W1/W2: on the common flux scale (AB, ZP 25) a flat-Fν source has
  equal expected tensor flux in every band of every archive. Reflected
  sunlight is the primary physical interpretation; self-luminous
  optical emission, a beacon, or sufficiently broad and recurrent
  leakage are also detectable as flux. A reflected-light size
  constraint must state the assumed albedo. A passive ~12 K body at
  550 AU is invisible in every band used here, so this stage is not
  thermal coverage; sources whose true SED falls steeply from the
  optical to 3–5 µm are in scope for detection but outside the
  colour-veto cell (§4: the veto can only reject candidates whose W
  flux contradicts the frozen flat-Fν cell).
- **Time behavior:** duty cycle ≥ 0.5 under the persistent,
  exposure-flicker, visit and long-block families. Rare glints and
  low-duty pulses are outside the cell.
- **Scope:** swarms, off-axis infrastructure, inactive/dark relics,
  extended or trailed sources, distances or motions outside the stated
  bounds and different endpoint hypotheses require separate searches.
  The result is a targeted flux constraint, not an inference about
  relay prevalence or network architecture.

## 2. Optical candidate statistic (unchanged from v2.0)

Joint cell = endpoint × role × band pair (g = PS1 g + ZTF zg; r; i).

1. **Statistic.** S_joint(z, µ) = (A_ztf + A_ps1) / √(B_ztf + B_ps1)
   from the per-archive v2 accumulated stack sums (per trajectory, per
   mask, capped and clipped per archive), the ZTF sums interpolated in
   1/z onto the PS1 grid (360 nodes). Trajectory index t is the same
   sky offset in both archives (identical 48-offset ring), so the
   joint ring null is exact; donors are annotations.
2. **Candidate rule.** R = S_max/T over the 8 designated controls,
   R̃ = R/q95(ring), heavy-tail / radius-dependent cells void,
   **one family-wise threshold across archives**: R̃_FWER at α = 0.05
   over the set's joint cells. **The W bands never enter S_joint or
   the family** — v3 adds no candidate cells, so the family and its
   threshold calibration are input-identical to v2.
3. **Static veto (unchanged).** Flux-consistent catalogued static
   source: the PS1 DR2 and ZTF DR predictions pushed through each
   archive's response at the joint node, summed against the joint S;
   rejection only if together they account for it within a factor 2.
4. **Annotations (unchanged).** Phase split (p_phase), p_trajectory.

## 3. WISE colour axis — conventions

W1/W2 tensors rebuilt through the engine under the joint conventions
(`surveys/wise/profile.py`, `runs/wise/v3/`; operational freeze
`surveys/wise/configs/v3_freeze.json`):

- T0 = 59800 (v2 WISE tensors at T0 = 57800 do not correspond on the
  grid and are not used);
- AB on the common ZP 25 scale via the frame Vega MAGZP plus the
  Vega→AB offsets W1 +2.699, W2 +3.339;
- µ grid 9 nodes at 0.25″/yr, whose [::2] subgrid is exactly the
  optical 5-node grid (T0 sits ~5 yr off the WISE mid-baseline; the
  0.5″/yr step would quantise tracks by ~0.5×FWHM — the PS1 T0
  lesson); joint µ nodes map exactly, z evaluated at the nearest WISE
  node (64 nodes uniform in 1/z);
- W1/W2 only; W3/W4 excluded (scientific-review rule: threshold
  statements only, in the standalone WISE survey);
- WISE measurement conventions otherwise identical to the frozen v2.1
  survey (quality masks, weight cap, no single-epoch clip, spacecraft
  observer, empirical PRF grid).

## 4. Colour-consistency veto (new; calibrated)

At a candidate's joint peak node (iz, iµ1, iµ2):

- f_pred = A_joint/B_joint, σ_pred = 1/√B_joint (the joint optical
  stack flux estimate, AB ZP 25 units); under the frozen flat-Fν
  spectrum the predicted WISE stack flux equals f_pred.
- Per W band: f_w = A_w/B_w and n_w from the WISE v3 tensor (primary
  mask, real trajectory) at the mapped node;
  **σ_w = max(1/√B_w, 1.4826 × MAD of the 48 ring-trajectory f̂ at the
  same node)** — the ring term is the empirical confusion floor
  (learnings §5: per-pixel uncertainties understate it at ≥6″
  resolution); the ring term requires ≥ 8 ring members with n ≥ 5 and
  B > 0, else 1/√B_w alone.
- Deficit statistic D_b = (f_pred − f_w) / √(σ_w² + σ_pred²).
- A band is **usable** iff n_w ≥ 5 and σ_w finite and positive.
- **The veto fires iff every usable W band has D_b ≥ ν = 5 and at
  least one band is usable.** Direction: the veto can only reject a
  candidate whose WISE flux is *deficient* against the flat-Fν cell —
  excess W flux (e.g. a red blend) never fires it.

Calibration (the three calibrated-veto conditions, learnings §2): the
selection function is measured by pushing the **same injections**
through the WISE chain — per-injection seeded draws are archive- and
band-independent, and the WISE injection runs reproduce each optical
band family's exact PS1/ZTF union magnitude window
(`runs/wise/v3/injections_{g,r,i}`), so the j-th injection of a cell
is one physical source in all three archives. On injections the veto
is evaluated from the injected W1/W2 window sums (ring floor from the
cell's tensor at the node) and **false vetoes are charged to the
final-candidate completeness curve** (the threshold curve stays
veto-free). An in-scope flat-Fν source satisfies E[f_w] = f_pred, so
its false-veto probability is the measured tail rate, expected ≈ 0;
WISE depth (single-band stack ≳ AB 17–19 against optical windows of
AB ~19–25) means the veto has power only against bright candidates —
that asymmetry is measured, not assumed.

Annotations: f_w, σ_w (both terms), n_w and D_b per band are recorded
on **every** cell at its real-trajectory joint node, veto or no veto.

## 5. Hold-out, ordering and provenance

- Split: the PS1/ZTF corridor split (seed 20260822), identical by
  construction across the three per-archive freezes — 20 development /
  42 confirmatory corridors (24 / 45 endpoints).
- Ordering: development set first (validating the measured false-veto
  rate); the confirmatory set analysed once. The v2 confirmatory
  optical statistic is unchanged and its 0-candidate outcome is
  expected to reproduce identically; v3's new confirmatory content is
  the colour axis.
- Rejected-candidate statuses: `static_flux_consistent` (§2.3),
  `colour_inconsistent` (§4); everything else `retained-ambiguous`.
- Pins: ZTF v2.0, PS1 v2.0 and WISE v3.0-joint-conventions freezes by
  content hash; the v2.0 joint freeze is superseded with a recorded
  hash; bulk products under `runs/joint/v3/`.

---

# Superseded: hypothesis freeze v2.0 (2026-08-22)

Kept as the record of the executed v2 run (report
`report/joint_ps1_ztf.md`, products `runs/joint/v2/`, hash
`configs/v2_freeze.json`). The v2 physical-source hypothesis and rule
are reproduced in §§1–2 above with only the v3 colour-axis additions;
no v2 parameter was altered.

> **Documentation restoration (2026-08-24, v2).** The v2 freeze section
> was a post-run, meaning-preserving copy of the physical assumptions
> already used by both archive freezes. It changed no search parameter,
> data product, candidate decision or result; the frozen configuration
> remains the provenance authority for the executed v2 run.

Operationally (v2): common reference epoch MJD 59800, AB zero point
25, joint cell = endpoint × role × band pair; statistic, candidate
rule, static veto and annotations as §2 above; completeness by
same-source per-injection seeded draws on the union magnitude window;
hold-out = the PS1/ZTF corridor split (seed 20260822), confirmatory
set analysed once. Result: 0 candidates (confirmatory 45 endpoints /
230 cells, R̃_FWER 1.778; development 24 / 102).
