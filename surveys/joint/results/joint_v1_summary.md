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
