# ptf

PTF/iPTF corridor survey (project plan §4.13): Pipeline A over the
IRSA PTF level-1 epochal images (2009-03 → 2015-01, g + Mould R) on
the 62 Palomar-visible corridors — the PS1-era optical archive from a
different site, cadence and filter set, adopted for the parallax
phases PS1 lacks and an era-independent persistence test. Built
directly on the v2 engine (DECam pattern, no v1 exploratory pass).

- Adapter: `sglsurvey/adapters/irsa_ptf.py` (built for the crossings
  survey; `contains_any` added for the corridor coarse stage).
  Observer: Palomar P48 = `GeometryContext.ztf_default()`.
- Corridors: `scripts/ptf_corridors.py` (ZTF-visible set; searchable
  subset from `targets/overlay_v1.json`, built by
  `targets/build_ptf_overlay.py` from the coarse + precise records).
- Flux scale: `scripts/ptf_calib.py` (PS1 DR2 mean calibrators from
  `runs/panstarrs/screen_v1`, Jordi R + 0.21 AB, dwarf-locus rule) and
  `scripts/calibrate_zeropoints.py` → `runs/ptf/zeropoints.jsonl`.
- Hypotheses: `hypotheses.md` (v0.1 draft → v1.0 at the freeze,
  `configs/v2_freeze.json`).
- Positive control: (60000) Miminko, `scripts/asteroid_control_v2.py`
  (`--stage fetch` → `configs/asteroid_control_v1.json`; `--stage score`).
- Recon checks: `scripts/wcs_linearity_check.py` →
  `results/wcs_linearity_check.json`; `scripts/flux_scale_check.py` →
  `results/flux_scale_check.json`.

Stages (Pipeline A):

1. `scripts/coarse_discovery.py [endpoint ...]` → `runs/ptf/coarse_v1/`
2. `scripts/precise_pass.py [endpoint ...]` → `runs/ptf/precise_v1/`
   (dmask cutouts under `runs/ptf/products/msk/`)
3. `targets/build_ptf_overlay.py` → `targets/overlay_v1.{json,md}`
4. `scripts/fetch_cutouts.py` → `runs/ptf/products/cut/` (+ manifest)
5. `scripts/calibrate_zeropoints.py` → `runs/ptf/zeropoints.jsonl`
6. v2 engine (`profile.py` + `run.py`):
   `uv run python surveys/ptf/run.py
   {freeze,geometry,build,nulls,inject,completeness,adjudicate,report}`
   — development set first, blind confirmatory only after the rule is
   frozen on dev.

Survey state (2026-09-10): **complete** — development (36 cells) and
blind confirmatory (66 cells) both 0 candidates at FWER α = 0.05;
1,632 Constraints; canonical report `report/ptf_survey.md`, tables in
`results/report_tables.md`. Post-freeze driver: `scripts/run_chain.sh
{data,dev,conf}`. The archive is closed (no yearly refresh).
