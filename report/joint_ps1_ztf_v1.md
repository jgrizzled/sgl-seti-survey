---
title: "Joint Pan-STARRS1 + ZTF stage 2 v1 — both parallax phases, 2009–2026, on 62 corridors"
date: 2026-08-21
status: "v1 (µ = 0) and v2 (full µ grid, common T0) complete; 0 candidates; PS1 retained cells all closed"
---

# Joint PS1 + ZTF stage 2 v1

**Question.** The PS1 survey (`ps1_survey_v1.md`) found that 3π samples
every corridor at one parallax phase, so a static background source at
the track position could not be excluded from PS1 alone, and it left
three marginal cells open. Can the ZTF 2018–2026 stacks on the same
corridors supply the other phase, close those cells, and what does the
combined 17-year, 1″ optical baseline say about station-kept relays?

**Answer.** Yes on all counts. Combining the per-archive sample tensors
on the exact common trajectory (µ_resid = 0) gives 321 of 414
endpoint-role-band cells with both parallax phases populated (PS1
alone: effectively none), and the census of 35 exceedances in 414 cells
(8.5 %, below the 12.5 % chance rate) resolves entirely into
single-phase cells, static sources at one phase position, non-persistent
bumps, and — for the two cells the automatic rules retained — catalogued
stars within 1–2″ of i-band tracks in crowded fields. **No candidate
survives.** The three PS1 marginal cells are vetoed directly: forced
photometry on 900–1,200 ZTF frames along each PS1-fitted (z, µ)
trajectory gives S = 0.2, 1.4, −0.1. Joint µ = 0 depths reach median
m90 = 23.35 (g), 23.26 (r), 21.48 (i) AB, ≈ 2.3 mag beyond PS1 alone.

## 1. Method

- Inputs: PS1 tensors (360-node 1/z grid, per-frame star-calibrated AB
  flux) and ZTF tensors (192-node grid, MAGZP scale) for the 138
  endpoint-roles both archives cover; ZTF interpolated in 1/z and scaled
  by the 0.5 mag asteroid-control throughput so both sit on the PS1
  total-flux scale. Bands paired g/zg, r/zr, i/zi (achromatic within a
  pair assumed; stated on every record).
- Hypothesis: station-kept relay, µ_resid = 0. This is the one
  trajectory identical in both tensor families — their µ reference
  epochs differ (T0 = 56000 vs 59800), so a µ ≠ 0 cell is a different
  track in each archive; |µ| ≤ 1″/yr remains covered by the per-archive
  calibrations.
- Stack: the WISE v0.2.0 estimator (1/σ² weights capped at 20× median,
  epoch floor 5, 5σ single-epoch clip), 8 offset controls shared by both
  archives (thresholds joint by construction), parallax phase
  recomputed for both archives against one reference, µ = 0 analytic
  injections at duty 0.5 for the depths.
- Adjudication: automatic phase veto, split-half persistence, other
  paired bands at the same z; stage-7 catalogued-star test
  (`runs/joint/ps1_ztf_v1/retained_star_check.json`).
- PS1 marginal cells: `marginal_ztf_test.py` — ZTF hybrid diff/sci flux
  maps sampled along the PS1-extrapolated trajectory (z, µ, T0 = 56000)
  and at the 8 offsets, all ZTF bands.

## 2. Results

AnalysisRun `run-fd75b2c982c2`: 3,312 Constraints (3,302 recovery
curves; 2,568 flagged both-phase), 35 Candidates (33 vetoed
automatically, 2 vetoed at stage 7), 0 retained. Details and tables:
`surveys/joint/results/joint_v1_summary.md`.

The two stage-7 vetoes illustrate the failure mode of i-band joint
cells: ZTF i has frames on only 76 of 138 endpoint-roles and those
cluster in weeks, so i cells are effectively PS1-only and crowded
fields put a DR2 star within 2″ of the track (GJ 783 tx at 10,000 AU,
where the track moves only 20″; Ross 154 rx at 1,290 AU, an i = 19.6
star 1–2″ away). In both, the ~1,000-epoch ZTF g/r stacks at the same
z are ≤ 2.7σ.

## 3. Interpretation

For a station-kept relay the optical cell is now closed to m ≈ 23.3 AB
(g, r) over 2009–2026 on 62 corridors with the static-background veto
in force: reflected sunlight (albedo 0.1) D ≲ 3×10⁴ km at 550 AU; self
luminous sources ≲ 2 µJy (r). These remain limits on the reflected /
self-luminous cell only. The |µ| ≤ 1″/yr nuisance family is constrained
per archive (PS1 m90 ≈ 21, ZTF per its own calibration), not jointly.

## 4. v2 — the full (z, µ) family on a common reference epoch

The PS1 tensors were rebuilt with T0 = 59800 (ZTF's) and the joint
stack run over the 5 × 5 µ grid (AnalysisRun `run-ac08543c5b29`,
`surveys/joint/results/joint_v1_summary.md` §v2). 43 of 414 cells
exceed the cube threshold (10.4 %, chance 12.5 %); four pass the
automatic rules — now including the catalogued-static-source test —
and all four are vetoed at stage 7: a field-wide systematic (GJ 1111
rx g, all controls at 17–19σ), faint catalogued stars along the track
at the signal's own brightness (GJ 229 A tx g), and the GJ 783 i-only
pair absent in ZTF g+r along the same track. **Still no candidate.**
On-grid depths are m90 = 23.0 (g), 22.9 (r), 21.1 (i) AB;
off-grid-marginalised 22.6 / 22.5 / 14.2: with T0 = 59800 the PS1
epochs are 7–13 yr from the reference and the 0.5″/yr µ step is
1.7–3.2″ of displacement there, so sources between µ nodes are lost
from PS1. The v1 µ = 0 result remains the cleanest joint statement; a
v3 needs a mid-baseline T0 with ≈ 0.1″/yr µ sampling or analytic
interpolation between nodes.

The catalogued-static-source test now lives in `sglsurvey/vetting.py`
and runs automatically in the PS1 and ZTF calibrations and both joint
scripts (WISE's census tool uses its VizieR CatWISE loader). With a
≥ 3-detection requirement it is conservative; the six PS1 cells it left
were all closed by direct ZTF forced photometry along their tracks.

## 5. Next

1. v3 joint µ family: mid-baseline T0 and finer µ sampling (or
   analytic interpolation) so PS1 epochs contribute off-grid.
2. Extend the joint stage to WISE (thermal cell) on these corridors;
   DECam/NOIRLab or SPHEREx for the 15 southern corridors (plan §7 TODO).

---

## Status note (2026-08-22) — exploratory, pending v2

The WISE scientific review of 2026-08-21
(`surveys/wise/scientific_review.md`) applies to the joint stage: the
per-archive 8-control thresholds are per-search rank statistics and
the joint stage multiplies the family (one family-wise error rate
across archives is needed); injections were tensor-level, so the joint
depths are *threshold sensitivity* (`completeness_kind = threshold`)
without confidence intervals and the same source was never injected
into both archives' images; the phase, split-half and role rules are
heuristic, the catalogued-star and direct forced-photometry tests the
only calibrated-in-principle rejections; loci nominal. No number is
recomputed here. The defensible conclusion is *no compelling candidate
after heuristic review*. The calibrated version is the joint v2 stage
(project plan §10.1 step 6, after ZTF v2 and PS1 v2).
