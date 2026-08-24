# joint v2

v2 re-analysis of the v1 joint survey under the frozen rule in
`hypotheses.md` (WISE v2.1 design; v2 programme step 6, `notes/project_history.md` §5), run by
the survey-agnostic engine `sglsurvey/` with this directory's
`profile.py` bindings. `configs/v2_freeze.json` holds the hashed freeze
and the stratified development / confirmatory split;
`configs/cross_track_cells.json` the cells with a cross-track
dimension; `results/` the ledger-generated tables. Bulk products under
`runs/joint/v2/`.

Stages: `uv run python surveys/joint/run.py {freeze,geometry,build,nulls,inject,completeness,adjudicate,report}`
(joint: `surveys/joint/joint.py {freeze,nulls,completeness,adjudicate,report}`).
