# surveys/panstarrs — Pan-STARRS1 pilot (fourth archive adapter, first non-IRSA)

Plan §3.6 priority 4: a 3-corridor pilot of Pipeline A against PS1 DR2
single-epoch warp images at MAST, mirroring the ZTF pilot (§5) on the
same corridors so that PS1 (2009–2014) extends each ZTF (2018–) corridor
by a 5–10-year baseline. Started 2026-08-20 while IRSA throughput was
saturated by the ZTF scale-up and the SPHEREx pilot.

## Contents

- `hypotheses.md` — baseline hypothesis freeze v1.0 (inherits WISE/ZTF
  v1.0 physics; PS1-specific observer, bands, mask template, flux scale)
- `notes/mast_recon.md` — verified facts about the warp listing, fitscut,
  mask bits, catalogs API, and the obsTime offset
- `scripts/` — pipeline entry points, in execution order:
  `ps1_corridors.py` (config) → `coarse_discovery.py` → `precise_pass.py`
  → `catalog_screen.py` + `screen_recurrence.py` → `fetch_cutouts.py`
  → `sample_tensor.py` → `injection_calibrate.py` → `make_figures.py`;
  `asteroid_control.py` is the independent positive control;
  `adjudicate_candidates.py` is the stage-7 pass over retained
  candidates; `purge_products.py` + `run_scaleup.sh` drive the batched
  scale-up. Each run
  writes `run_config_*.json` under `../../runs/panstarrs/<run-id>/`
  (gitignored)
- `results/` — compact tracked summaries per stage
- `targets/` — PS1 overlay (`build_ps1_overlay.py` → `overlay_v1.json/.md`)

Adapter: `sglsurvey/adapters/mast_ps1.py`. Photometry:
`sglsurvey.photometry.build_flux_map_ps1`. Geometry:
`GeometryContext.ps1_default()` (Haleakalā observer). Registry: shared
`registries/pilot_wise_2026.yaml` (survey-agnostic).

## Pilot step status

1. Freeze hypotheses — **v1.0 frozen 2026-08-20**
2. Corridors — `ross128`, `epsind`, `proxima` (same as ZTF pilot)
3. Coarse discovery — **run** (`coarse_v1`: 1,945 warps over 11
   skycells; ~50 % geometric hits per filter)
4. Precise pass — **run** (`precise_v1`: full skycell masks; ~75 % of
   hits usable, mean usable fraction 0.5–0.8; 6 masks not served)
5. Catalog screening — **run** (`screen_v1`: DR2 detection table;
   6,040 matches; nothing track-following)
6. Forced photometry / stack — see `results/pilot_v1_summary.md`
7. Injection calibration + asteroid positive control — see summary
8. Report — `../../report/ps1_pilot_v1.md`

Stage details: `results/pilot_v1_summary.md`.

## Scale-up (2026-08-20 → 21, complete)

`targets/build_ps1_overlay.py` → `targets/overlay_v1.{json,md}` (62 of
77 corridors, Dec > −30°; phase split and calibrator density per
corridor). `scripts/run_scaleup.sh` ran the chain over all 62 corridors
in batches (cutouts → tensors → `purge_products.py`), then
`injection_calibrate.py` and `adjudicate_candidates.py`. Results:
`results/scaleup_v1_summary.md`, report `../../report/ps1_survey_v1.md`.
AnalysisRun `run-eaa6d89a1ec9`: 5,520 Constraints, 78 exceedances
(chance rate), 0 detections, 3 marginal cells held for the ZTF
cross-archive test.
