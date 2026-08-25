# joint v3 — PS1 + ZTF stack, WISE colour axis

v3 extension of the joint stage under the frozen rule in
`hypotheses.md` (v3.0, 2026-08-24): the v2 optical candidate statistic
unchanged (PS1 + ZTF stack sums on the PS1 grid, one FWER across
archives), plus the WISE W1/W2 colour axis — per-cell annotations and
the calibrated colour-consistency veto — from the joint-conventions
WISE tensors (`surveys/wise/profile.py`, `runs/wise/v3/`).
`configs/v3_freeze.json` holds the hashed freeze (split = the PS1/ZTF
corridor split); `configs/v2_freeze.json` is the superseded v2 record
(report `report/joint_ps1_ztf.md`, products `runs/joint/v2/`). Bulk v3
products under `runs/joint/v3/`.

Stages (`surveys/joint/joint.py`):

    uv run python surveys/wise/run.py build            # WISE v3 tensors (once)
    uv run python surveys/joint/joint.py freeze
    uv run python surveys/joint/joint.py wise-inject --set dev|confirmatory
    uv run python surveys/joint/joint.py {nulls,completeness,adjudicate,report} --set dev|confirmatory
