# surveys/spherex — SPHEREx pilot (third archive adapter)

Plan §6: a 3-star pilot of Pipeline A against SPHEREx Level-2 spectral
images at IRSA, mirroring the WISE and ZTF shakedowns, run in parallel
with the ZTF scale-up and the last WISE batch because those two are
IRSA-download-bound and SPHEREx products can be read from the public
S3 mirror.

## Contents

- `hypotheses.md` — baseline hypothesis freeze v1.0 (inherits WISE v1.0
  physics; SPHEREx-specific observer, detectors, FLAGS template, flux
  scale, static-sky template)
- `notes/irsa_recon.md` — verified facts about the TAP tables, the
  CAOM join that yields product URIs + published MD5s, the Level-2 MEF
  layout, S3 byte-range access and the IBE cutout service
- `scripts/` — pipeline entry points, in execution order:
  `spherex_corridors.py` (config) → `coarse_discovery.py` →
  `precise_pass.py` (also assembles the slim cutouts) →
  `star_control.py` (flux-scale positive control) →
  `static_template.py` → `track_screen.py` → `sample_tensor.py` →
  `injection_calibrate.py`; `run_pilot.sh` runs the chain. Each run
  writes `run_config_*.json` under `../../runs/spherex/<run-id>/`
  (gitignored)
- `results/` — compact tracked summaries per stage

Adapter: `sglsurvey/adapters/irsa_spherex.py`. Registry: shared
`registries/pilot_wise_2026.yaml` (survey-agnostic; no SPHEREx entries
needed). Photometry: `sglsurvey.photometry.build_flux_map_spherex`
(exposure-PSF matched filter at 2×2 sub-pixel phases, µJy scale).

## What is SPHEREx-specific in the chain

| stage | SPHEREx specific |
|---|---|
| discovery | IRSA TAP `spherex.obscore` ⋈ `plane` ⋈ `artifact` (URI, MD5); centre-distance superset + client-side ObsCore polygon footprint |
| products | one *slim cutout* per detector-exposure assembled by S3 byte-range reads: IMAGE/FLAGS/VARIANCE/ZODI sub-arrays, the PSF plane of the zone containing the cutout, the 9×9 wavelength table; ~260 kB on disk, ~3 s |
| usable pixels | FLAGS fatal template 708343 (all bits except FULLSAMPLE and SOURCE) |
| search image | IMAGE − ZODI; no difference images exist → per corridor × detector **static-sky template** (per-node linear-in-λ fit across epochs) subtracted at sampling time |
| matched filter | the exposure's own PSF binned to detector sampling, evaluated at 4 sub-pixel phases (2× map) — the undersampled 5.3" PSF on 6.15" pixels otherwise loses up to 0.5 mag at half-pixel offsets |
| flux scale | MJy/sr × OMEGA_MEDIAN → µJy; verified on catalogued 2MASS/CatWISE/Gaia stars (`star_control.py`) |
| bands | D1–D6; per-epoch wavelength and bandwidth carried on every sample; constraints quote the covered λ range; candidate reports carry a 4-bin spectral significance |
| observer | Earth centre (LEO offset ≤ 0.017" carried as a budget term) |

## Pilot step status

1. Freeze hypotheses — **v1.0 frozen 2026-08-20** (+ role-coincidence amendment)
2. Stars — `lalande-21185`, `gj-687`, `sigma-dra` (hypotheses §1)
3. Coarse discovery — **run** (`coarse_v1`: 13,488 detector-exposures, 16,978 hits)
4. Precise pass + cutouts — **run** (`precise_v1`: 16,761 usable/partial evaluations, 0 missing products)
5. Flux-scale control — **run** (`control_v1`: 13,657 star measurements, all detectors within ±0.2 mag)
6. Static template + layer-1 screen — **run** (`screen_v1`: no track-following cluster)
7. Sampling + injection calibration — **run** (`calib_v1`, AnalysisRun `run-55623bb1127a`: 288 Constraints, 8 Candidates all vetoed)
8. Report — **done**: `../../report/spherex_pilot_v1.md`

Stage details: `results/pilot_v1_summary.md`.

## Scale-up (2026-08-20, complete)

All 15 ZTF-inaccessible corridors (19 endpoints) run through the same
chain with `SPHEREX_CALIB_RUN=calib_v2`. Results:
`results/scaleup_v1_summary.md`, `../../report/spherex_survey_v1.md`
(AnalysisRun `run-4f8c9fd41ce9`: 1,824 Constraints, 28 Candidates all
vetoed). Scripts honour `SPHEREX_CALIB_RUN` for the template/tensor/
calibration directory; discovery, precise, control and screen runs
extend their `*_v1` directories in place.

## v2 estimator (2026-08-21, `calib_v3`)

`sample_tensor.py` scales epoch variances by the template's per-node
reduced χ²; `injection_calibrate.py` adds the joint six-detector stack
(band `ALL`) and the template-coverage and bright-static-neighbour
vetoes (hypotheses amendment 2). AnalysisRun `run-4fbed2e52668`:
2,128 Constraints, 37 Candidates all vetoed. `calib_v3/templates` is a
symlink to `calib_v2/templates`.

## All-sky run (2026-08-21, complete)

v3 estimator (16 controls, cross-detector phase test; amendment 3)
then the 62 northern corridors (69 endpoints) through the chain with
`SPHEREX_CALIB_RUN=calib_v4` and `sample_tensor.py --only-missing`.
Results: `results/allsky_v1_summary.md`, `../../report/spherex_survey_v1.md`
(AnalysisRun `run-ddd97a27faad`: 9,808 Constraints, 64 Candidates all
vetoed). 77 corridors / 88 endpoints; 39 GB under `runs/spherex/`.
