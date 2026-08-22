# Joint PS1 + ZTF stage 2 v1 — summary (2026-08-21)

AnalysisRun `run-fd75b2c982c2` (`runs/joint/ps1_ztf_v1/`), hypothesis
`joint-ps1-ztf-v1.0` = station-kept relay (µ_resid = 0), 550–10,000 AU,
duty ≥ 0.5, paired bands g/zg, r/zr, i/zi, ZTF on the PS1 star-calibrated
flux scale (+0.5 mag throughput), common parallax-phase reference.

## Coverage

138 endpoint-roles (69 endpoints, 62 corridors) in both archives; 414
endpoint-role-band cells; median epochs per cell 15 (PS1) + 349 (ZTF);
epoch range 2009-02 → 2026-07. **321 of 414 cells have both parallax
phases populated** (≥ 3 epochs each) — versus 4 of 40 in the PS1 pilot
and ≈ 0 across PS1 alone. The 93 single-phase cells are almost all i
band, where ZTF has frames on only 76 of 138 endpoint-roles and those
cluster in short windows.

## Depths (µ = 0, 90 % recovery, duty ≥ 0.5, AB on the PS1 scale)

| band | median m90 | 10–90 % | PS1-alone median |
|---|---|---|---|
| g (+zg) | 23.35 | 21.99–23.80 | 21.06 |
| r (+zr) | 23.26 | 22.04–23.79 | 20.86 |
| i (+zi) | 21.48 | 20.65–21.87 | 20.72 |

3,312 Constraint records (3,302 recovery curves, 10 not-constrainable);
2,568 carry `both_parallax_phases = true`. The ≈ 2.3 mag gain over PS1
alone in g/r is the ~350-epoch ZTF stack at a single trajectory (no µ
trials); these depths apply to the µ = 0 hypothesis only — the
per-archive µ-grid constraints remain the reference for |µ| ≤ 1″/yr.

## Exceedance census

35 of 414 cells exceed the 8-control threshold (8.5 %; chance rate
12.5 %). Automatic vetoes: 14 single-phase (i-band cells), 16 phase-split
disagreement (static source at one phase position: e.g. Barnard rx i
6.0/−0.8, Luhman 16 tx g 9.6/0.9, van Maanen rx r 1.4/7.8), 3 split-half
not persistent (EZ Aqr rx g 3.2/1.4, GJ 832 tx r, GJ 876 tx r), 2 absent
in the other paired bands (GJ 667 C tx r, GJ 674 tx g). Two retained by
the automatic rules, both i-band, vetoed at stage 7
(`retained_star_check.json`, `records/candidate_adjudicated.jsonl`):

- **GJ 783 tx i**, z = 10,000 AU (grid edge), S 12.7 vs T 10.1: DR2
  i ≈ 21 sources 1–2.5″ from the (20″-parallax) track at every epoch;
  control maxima up to 10.1 (crowded); g/r 0.3/2.7σ at the same z.
- **Ross 154 rx i**, z ≈ 1,290 AU, S 6.7 vs T 5.2: DR2 star i = 19.6
  (58–64 detections) 1.0–2.3″ from the track at the median epochs; all
  46 ZTF i epochs in one 2-month window; ZTF g/r (~1,000 epochs, both
  phases) 2.1/1.8σ at the same z.

**No candidate survives.**

## PS1 marginal cells — direct ZTF test (`marginal_tests.json`)

Forced photometry on the ZTF cutouts along the PS1-fitted (z, µ)
trajectories, all ZTF bands, 8 shared offset controls:

| PS1 cell | PS1 S / T | ZTF frames | ZTF S (all) / T | verdict |
|---|---|---|---|---|
| 82 Eri rx y, z 7,231, µ (−0.5, +1.0) | 9.37 / 8.19 | 897 (g 465, r 425) | 0.19 / 2.83 | vetoed |
| Fomalhaut rx y, z 574, µ (+1, 0) | 5.42 / 5.15 | 1,206 | 1.42 / 1.46 | vetoed |
| GJ 526 tx z, z 559, µ (−1, −0.5) | 5.15 / 5.13 | 1,023 | −0.11 / 0.25 | vetoed |

A source at the PS1 amplitude with any plausible colour would reach
S ≫ 10 in these stacks. The PS1 survey closes with **0 candidates**.

## Lessons

- ZTF i is too sparse and too clustered in time to carry the phase
  test; i-band joint cells are effectively PS1-only and should be read
  as such.
- Dense-field i-band stacks are star-contaminated at the ≤ 2″ level;
  the catalogued-star test must be part of the automatic rules, not
  only stage 7.
- The µ reference epoch must be common across archives for any future
  joint µ-grid analysis (rebuild tensors with a shared T0 or store
  offsets).

---

# Joint PS1 + ZTF stage 2 **v2** — full (z, µ) family on a common T0 (2026-08-21)

AnalysisRun `run-ac08543c5b29` (`runs/joint/ps1_ztf_v2/`), hypothesis
`joint-ps1-ztf-v2.0`: relay on the focal line at T0 = 59800 with
|µ| ≤ 1″/yr on the 5 × 5 grid. PS1 tensors rebuilt on T0 = 59800
(`runs/panstarrs/calib_t0_59800/tensors`, `run_common_t0.sh`, 3.3 h,
≈ 180 GB transient); ZTF tensors unchanged. Script
`joint_ps1_ztf_mugrid.py`; stage 7 `adjudicate_v2.py`.

## Depths (90 % recovery, duty ≥ 0.5, AB on the PS1 scale, median over 138 endpoint-roles)

| band | µ on grid | off-grid marginalised | v1 (µ = 0) | PS1 alone |
|---|---|---|---|---|
| g | 23.00 | 22.56 | 23.35 | 21.06 |
| r | 22.90 | 22.48 | 23.26 | 20.86 |
| i | 21.09 | 14.17 | 21.48 | 20.72 |

3,312 Constraints (3,298 recovery curves); `flux_limit.value` is the
on-grid depth, `value_off_grid_marginalised` the uniform-µ one.

**The common-T0 lesson.** With T0 = 59800 the PS1 epochs are 7–13 yr
from the reference, so the 0.5″/yr µ step is 1.7–3.2″ of displacement
at PS1 epochs — larger than the PSF. The search on the grid nodes is
valid, but a source with µ between nodes is lost from PS1 (and
partially from early ZTF), which is why the off-grid i-band depth
(PS1-dominated) collapses and g/r lose ~0.45 mag. A mid-baseline T0
(≈ 58000) halves the lever arm but still needs ≈ 0.1″/yr µ steps —
≈ 25× the tensor volume — or an analytic µ-interpolation between
nodes. Recorded as the v3 item; v1 (µ = 0, exact) remains the
cleanest joint statement.

## Exceedance census

43 of 414 cells exceed the 8-control cube maximum (10.4 %; chance
12.5 %). Automatic vetoes: 22 single-phase, 13 phase-split
disagreement, 1 split-half, 3 absent in the other paired bands. Four
retained, all vetoed at stage 7 (`adjudication.json`,
`direct_ztf_tests.json`, `records/candidate_adjudicated.jsonl`):

| cell | S / T | verdict |
|---|---|---|
| GJ 1111 rx g, z 5,237, µ (−1, +1) | 19.03 / 18.94 | field-wide systematic: all 8 controls at 16.9–18.9σ; 0.5 % margin |
| GJ 229 A tx g, z 739, µ (−0.5, +0.5) | 6.41 / 5.88 | faint catalogued stars along the track through the ZTF era (DR2 objects with 43/13/34 detections at 2.1″/1.35″/2.4″, Gaia G = 21.0 at 2.1″) at the ≈ 23 AB level of the signal; PS1 1.4σ; 9 % margin |
| GJ 783 rx i, z 981 | 8.74 / 8.27 | i-only; absent in ZTF g+r along the track (3.3 vs 3.7); crowded field (controls to 8.3) |
| GJ 783 tx i, z 6,165 | 13.40 / 11.91 | i-only; absent in ZTF g+r along the track (2.8 vs 3.5); controls to 11.9 |

**No candidate survives.**

## PS1 automatic rules (`sglsurvey/vetting.py`, 2026-08-21)

The PS1 calibration was rerun with the catalogued-static-source test in
its automatic rules (AnalysisRun `run-c2806acd56a3` supersedes
`run-eaa6d89a1ec9`; same 5,520 constraints and 78 exceedances). With
the ≥ 3-detection requirement the catalog alone vetoes none of the six
automatically retained cells — the earlier GJ 783 / LHS 1723 "stars"
were 1–2-detection DR2 entries — and all six are vetoed by direct ZTF
forced photometry along their PS1-fitted (z, µ) tracks (S = 0.19,
1.42, −0.11, 1.87, 2.25, 1.77 vs control maxima 2.8, 1.5, 0.25, 6.7,
4.2, 2.8 over 897–1,206 frames; `marginal_tests.json`).
