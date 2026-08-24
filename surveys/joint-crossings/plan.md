# Joint crossings stage — pre-registration (v1.0)

Frozen 2026-08-24, before any joint statistic is computed. Plan §11.2
step 6. Two deliverables: the unified covered-window ledger and the
pre-registered gj-1276 recurrence test. This stage runs **no new blind
search**: it aggregates the three finished crossings surveys and
executes exactly one designated follow-up test, declared below before
computation.

## 1. Unified covered-window ledger

One row per (archive, channel, rung, event, band) that any survey
covered, from the frozen per-survey products:

| archive | inputs |
|---|---|
| ZTF | `surveys/ztf-crossings/results/{coverage_v1_events.ecsv, confirmatory_v1.json, completeness_v1.json}` |
| PS1 | `surveys/ps1-crossings/results/{coverage_v1_events.ecsv, confirmatory_v1.json, dev_search_v1.json, completeness_v1.json}` |
| WISE | `surveys/wise-crossings/results/coverage_v1_events.ecsv` + `configs/threshold_freeze_v1.json` unit classes |

Status vocabulary (frozen): `searched_null`, `exceedance_vetoed`,
`retained_ambiguous`, `single_epoch`, `single_window`,
`constraint_only`, `track_masked` (nominal-covered / mask-unusable),
`coverage_only` (WISE: no statistic exists by frozen construction).
Geometric-null cells (WISE channels B and A-0.1) are a summary note,
not ledger rows — they are not coverage. Columns: archive, channel,
rung_au, target_id, event_id, t_ca_mjd, window days, band, n_epochs,
status, S, T, margin, m90 (+censoring), notes.

## 2. The gj-1276 recurrence test (the programme's one designated follow-up)

**Input anomaly:** PS1 `evt-344c12d32a8d`, retained-ambiguous —
night-consistent i_AB ≈ 22.7–23.4 (1.1–2.3 σ/epoch × 3 exposures) at
the z = 550 AU track node, 2014-03-06; no stack counterpart to ~23+;
census clean; single covered PS1 window. Hypothesis under test: a
persistent relay at z = 550 AU on the gj-1276 axis with the PS1
amplitude (reference flux **i_AB = 22.8**), in-beam during every
wide-rung window (d = 1). Prediction: it recurs in *every* ZTF-covered
gj-1276 wide-rung window (7 events with g/r data, 2019–2026, cutouts
retained from the ZTF confirmatory run).

**Arm R1 — broadband recurrence (g + r), pre-registered rules:**

- Per (event, band): S₅₅₀ = the frozen ZTF channel-B stack along the
  z = 550 track — z *fixed* at the PS1-fitted value, no nested-z max;
  one trial. Same 8 designated ring controls (20/30/40″), same
  quality mask, same WEIGHT_CAP, from the retained confirmatory
  cutouts.
- Joint statistic J = Σ S₅₅₀ / √N over the event-bands with data;
  joint controls J_c pair ring index c across event-bands;
  T_J = max_c J_c. Exceedance rule J > max(T_J, 0) — a single trial.
- Sensitivity: joint image-level injections through the identical
  chain (the ZTF stamp-response machinery) under a **flat-Fν SED**
  normalized to i_AB = 22.8 (g = r = i = 22.8); joint recovery
  fraction at 22.8 and joint m90 reported with the per-survey
  conventions (≥ 100 draws).
- **Verdicts (frozen):**
  - joint recovery(22.8) ≥ 0.9 and no J exceedance →
    *flat-SED persistent relay at the PS1 amplitude REFUTED at 90 %*
    for z = 550;
  - joint recovery(22.8) < 0.9 and no J exceedance → *inconclusive*
    (depth margin reported);
  - J exceedance → *joint retained anomaly* — reported as
    cross-archive recurrence support, escalating to dedicated
    follow-up outside this stage; never silently absorbed.

**Arm R2 — line/red SED (i), decided by coverage arithmetic (no
pixels):** ZTF covered exactly 2 in-window i epochs across all
gj-1276 wide-rung windows (evt-27f78476fe86, 2026). A 2-epoch stack
reaches ~21.1 (single-epoch ~20.8 + 0.4), 1.7 mag short of the
required 22.8: the 752 nm-line interpretation of the PS1 signal is
**untestable by ZTF** and remains so for every current optical
archive; recorded in the ledger as an open cell.

**Static arm — already settled:** the PS1 DR2 *stack* catalog is
empty within 6″ of the node to ~23+ depth — deeper than any ZTF
reference stack could test, and ZTF difference imaging subtracts
static sources by construction. No new computation; the PS1 result is
cited.

**Trials accounting:** one pre-registered joint trial (Arm R1's J).
No other statistic computed in this stage may be promoted.

**Amendment v1.1 (2026-08-24, after the first R1 computation, before
any verdict was adopted):** the v1.0 ensemble included the two ZTF
r-band units carrying the bad-epoch anomaly *retained in the ZTF
confirmatory record* (gj-1276 r single-epoch outliers at −31/−48 σ);
those instrumental epochs poison both J (−13.5) and ~half the
bootstrap null draws (ring controls reaching −566), flattening the
recovery curve. R1 adopts the established ZTF single-epoch clip
(|f|/σ > 20 → instrumental, excluded; the ZTF pilot lesson and v1
survey rule) applied identically to on-track, control, and response
samples. The exclusion criterion is record-based (the anomaly was
frozen into the ZTF confirmatory record before this stage existed)
and value-symmetric (a +48 σ epoch would be clipped the same way).
Both versions are reported; v1.1 is the primary.

## 3. Report

`report/joint_crossings.md`: the ledger summary, the recurrence-test
verdict, the programme-wide crossings tally (searched units,
exceedance budget vs observed across all three surveys, retained
anomalies, structurally open cells), and the step-7 hand-off.
