# WISE/NEOWISE survey v2

Re-analysis of the v1 WISE survey (`surveys/wise/`) under the frozen
v2.0 decision rule, built per `surveys/wise/v2_plan.md` in response to
`surveys/wise/scientific_review.md`. Same registry (v1.5, 88 endpoints /
77 corridors), same frames, cutouts and physics cell; new statistical
experiment: a 48-offset exchangeable null per cell with a family-wise
error rate, image-level PRF injections with four temporal models,
calibrated vetoes with measured selection functions, covariance
propagation, a spacecraft observer, an asteroid positive control, and
content-hashed products.

- `hypotheses.md` — the v2.0 freeze (decision rule, null model,
  completeness definition, geometry, hold-out).
- `configs/v2_0_freeze.json` — hashed freeze with the stratified
  development / confirmatory corridor split (seed 20260822);
  `configs/cross_track_cells.json` — cells needing a cross-track
  dimension (from `geometry_check.py`).
- `scripts/` — `run_v2.sh` stages A–I; `v2common.py` holds every
  frozen parameter. Shared code lives in `sglsurvey/` (`nulls`,
  `inject`, `manifest`, `vetting`, `geometry`, `corridors`).
- `results/` — stage summaries and the ledger-generated report tables.
- Bulk products under `runs/wise/v2/` (tensors ~25 GB, injections,
  nulls, PRFs, observer table, control, records), all regenerable.

Order of work and gates: `v2_plan.md` §10 (git history —
`git show 8fb226f:surveys/wise/v2_plan.md`; see also
`notes/project_history.md` §5). The confirmatory set is
analysed once (`run_v2.sh G`); `null_ensemble.py` refuses a second run.

## v3 joint-conventions rebuild (2026-08-24, operational)

`profile.py` + `run.py` bind the flattened engine to the same v1 WISE
inputs to rebuild W1/W2 tensors and injections for the joint stage v3
(`surveys/joint/`): T0 = 59800, AB on ZP 25 (Vega→AB W1 +2.699 / W2
+3.339), 9-node µ grid, W1/W2 only — `hypotheses.md` §9. Not a
standalone search: `run.py` refuses the analysis stages; the decision
rule is the joint v3.0 freeze. Operational freeze (joint split):
`configs/v3_freeze.json`; products under `runs/wise/v3/` (reusing the
v2 observer table, PRF grids and covariance envelopes).
