# Joint beam-crossings stage (Pipeline B) — report

Closes the crossings programme's cross-archive stage (plan §11.2
step 6). Executed 2026-08-24 under the pre-registration
`surveys/joint-crossings/plan.md` (v1.0 + documented amendment v1.1).
Two deliverables: the unified covered-window ledger, and the
programme's one designated follow-up — the gj-1276 recurrence test.
**Verdict: the flat-SED persistent-relay interpretation of the PS1
retained-ambiguous exceedance is refuted at 90 %; the crossings
programme closes with 0 candidates.**

## 1. Unified covered-window ledger

`surveys/joint-crossings/results/covered_window_ledger_v1.ecsv` —
593 rows (ZTF 86, PS1 316, WISE 191), one per (archive, channel,
rung, event, band) with coverage, its disposition under that survey's
frozen rules, and its calibrated depth where one exists. 42
target-channels are covered by ≥ 2 archives (the cross-archive
recurrence axis: eras are disjoint, so "shared" means same target and
rung across different years).

| disposition | rows |
|---|---|
| searched_null | 44 |
| exceedance_adjudicated (ZTF, control-crossing statistics) | 3 |
| exceedance_vetoed (PS1 ross-128, stack-catalogued static) | 1 |
| retained_ambiguous (PS1 gj-1276 — resolved below) | 1 |
| single_epoch (reportable class) | 20 |
| constraint_only (frozen gates) | 28 |
| track_masked (PS1, nominal-covered / mask-unusable) | 7 |
| no_usable_data (ZTF) | 4 |
| coverage_only (WISE structural null; PS1 A-1.0 rung) | 485 |

Geometric-null cells (WISE channels B and A-0.1: 0 observable events
by the elongation theorem) carry no rows — they are not coverage.

## 2. The gj-1276 recurrence test

**Input:** PS1 `evt-344c12d32a8d` (retained-ambiguous): a
night-consistent i_AB ≈ 22.8 signal at the z = 550 AU track node,
2014-03-06 — census clean, mover-vetoed, no static stack counterpart
to ~23+, single covered PS1 window. A persistent d = 1 relay at that
amplitude predicts recurrence in *every* covered window.

**Arm R1 (broadband, pre-registered single joint trial).** The frozen
ZTF channel-B stack at z fixed to 550 AU over the 11 ZTF event-bands
with data (7 events × g/r, 2019–2026, retained confirmatory cutouts),
ring controls paired across event-bands. Amendment v1.1 (documented
in the plan before any verdict): the ZTF-recorded gj-1276 r bad-epoch
anomaly (−31/−48 σ instrumental outliers) is removed by the
established ZTF single-epoch clip (|f|/σ > 20), applied symmetrically
to on-track, control, and response samples; the unamended numbers are
retained in the result file.

| quantity | value |
|---|---|
| J (joint statistic) | 0.154 |
| T_J (max paired-ring joint control) | 1.116 |
| joint exceedance | **no** |
| joint recovery at 22.8 AB (flat Fν) | **0.94** |
| joint m90 | 22.93 AB |

**Verdict (frozen rule): flat-SED recurrence refuted at 90 %.** A
persistent z = 550 relay at the PS1 amplitude with a flat SED would
have been recovered with 94 % probability across ZTF's covered
windows; nothing recurred (J consistent with the paired-ring null).

**Arm R2 (752 nm line / red SED): untestable** — ZTF covered exactly
2 in-window i epochs across all gj-1276 wide-rung windows (2-epoch
depth ~21.1 vs the required 22.8). **Static arm: settled** by the PS1
DR2 stack catalog (no counterpart within 6″ to ~23+).

**Final status of the anomaly:** retained in the record, *narrowed*.
Surviving interpretations: (a) a red/line-SED source invisible to
g/r — untestable by any current optical archive at this depth; (b) a
non-persistent transient (which the crossing hypothesis does not
require to recur); (c) the 1/9 control-crossing budget doing what it
was designed to absorb (2 exceedances vs 2.1 expected in PS1). It is
not promotable and is not a candidate.

## 3. Programme-wide crossings tally (steps 1–6)

- **Searched units:** ZTF 36, PS1 11, WISE 0 (structural null) — 47
  searched trials programme-wide.
- **Exceedance budget:** 5 observed vs ~6.1 expected control
  crossings — the 1/9 discipline behaved as designed in every survey.
- **Candidates: 0.** Dispositions: 3 adjudicated control-crossing
  statistics (ZTF), 1 calibrated static veto (PS1), 1
  retained-ambiguous narrowed to near-exhaustion (this stage).
- **Calibrated depths where searched:** ZTF B median m90 ≈ 21.8
  (relay ≳ ~130 W through the 2.5 R☉ cone, 2018–2026 windows); PS1 B
  median 21.6 (≳ ~170–280 W, 2010–2014 windows); joint gj-1276
  ensemble 22.93.
- **Structurally open cells**, stated for the record: photosphere/
  coronal grazing rungs (ZTF: 5 windows on 2 targets searched; PS1:
  the b = 0.28 R☉ event mask-lost); 1064/1550 nm and all line SEDs
  redward of ~900 nm; the mid-IR crossings regime (WISE elongation +
  control-geometry theorems); pulse periods between exposure length
  and window length; transmitters scheduled to avoid Earth crossings;
  the d = 1 wide-beam (1 AU) rung as a *controlled* search anywhere
  (three independent proofs it defeats the frozen temporal-control
  family — the designated future route is a v2-style null ensemble).

## 4. Hand-off

Step 7 (radio-scope decision) is the programme's remaining open step.
The ledger is the coverage product a radio follow-up (or any future
crossings survey) consumes: per target × rung, which windows any
archive has searched, to what depth, and which are open. Yearly
refresh of the crossing lists (2028 window end) and re-running the
SPHEREx decision when 3+ annual windows exist are the standing
maintenance items.

## 5. Products

`surveys/joint-crossings/{plan.md, scripts/*, results/*}`; ledger
`results/covered_window_ledger_v1.ecsv` + `ledger_v1_summary.json`;
recurrence test `results/gj1276_recurrence_v1.json`; discovery
snapshot for the re-derived ZTF epoch map under
`runs/joint-crossings/`.
