# spherex v2

v2 re-analysis of the v1 spherex survey under the frozen rule in
`hypotheses.md` (WISE v2.1 design; v2 programme step 6, `notes/project_history.md` §5), run by
the survey-agnostic engine `sglsurvey/` with this directory's
`profile.py` bindings. `configs/v2_freeze.json` holds the hashed freeze
and the stratified development / confirmatory split;
`configs/cross_track_cells.json` the cells with a cross-track
dimension; `results/` the ledger-generated tables. Bulk products under
`runs/spherex/v2/`.

Stages: `uv run python surveys/spherex/run.py {freeze,geometry,build,nulls,inject,completeness,adjudicate,report}`
(joint: `surveys/joint/joint.py {freeze,nulls,completeness,adjudicate,report}`).

Deferred controls (2026-09-04, `report/spherex_joint6.md`; plan §5.15
O3): `joint6.py {freeze,inject,nulls,completeness,adjudicate,report}`
— the six-detector joint cell (hypotheses v2.1, `configs/
joint6_freeze.json`, products `runs/spherex/joint6/`, injections
`runs/spherex/v2/injections_joint/`, tables `results/joint6/`);
`scripts/template_control.py {fit,verify,validate}` — regenerate the
static templates from the cutouts (the v3 set was lost in the v1
retirement), verify against a rebuilt tensor, image-level validation of
the absorption module; `scripts/template_absorb.py {run,summarise,
correct}` — the template-absorption ladder and the report-level
corrected limits (`results/template_absorption{,_corrected}.md`).
