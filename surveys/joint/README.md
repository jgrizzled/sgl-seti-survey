# surveys/joint — cross-archive stage 2 (PS1 + ZTF)

Plan §3.4 motion-model stage across archives. First instance: the
joint Pan-STARRS1 (2009–2014) + ZTF (2018–2026) stack on every corridor
both archives cover, motivated by the PS1 finding that 3π samples each
corridor at one parallax phase — the static-background veto needs the
other phase, which ZTF supplies; the joint baseline is 17 years at 1″.

## Contents

- `scripts/joint_ps1_ztf.py` — joint weighted stack per endpoint × role
  × paired band (g/zg, r/zr, i/zi) for the station-kept relay
  (µ_resid = 0; the one trajectory identical in both tensor families,
  whose µ reference epochs differ), ZTF fluxes put on the PS1
  star-calibrated scale (+0.5 mag ZTF asteroid-control throughput),
  ZTF interpolated onto the PS1 1/z grid, 8 shared offset controls,
  common parallax-phase reference, µ = 0 injection depths, automatic
  phase veto + stage-7 split-half / other-band adjudication. Records
  under `runs/joint/ps1_ztf_v1/`.
- `scripts/marginal_ztf_test.py` — direct forced photometry on the ZTF
  cutouts along the PS1-fitted (z, µ ≠ 0) trajectory of each PS1
  marginal cell (the cross-archive test the PS1 report deferred).
- `results/joint_v1_summary.md` — numbers.

Inputs are the per-archive sample tensors (`runs/panstarrs/calib_v1/tensors`,
`runs/ztf/calib_v1/tensors`) and the ZTF cutouts; nothing is re-downloaded.
