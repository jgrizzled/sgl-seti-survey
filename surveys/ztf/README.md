# surveys/ztf — ZTF pilot (second archive adapter)

Plan §5: a 3-corridor pilot of Pipeline A against ZTF public
single-epoch products at IRSA, mirroring the WISE shakedown steps.

## Contents

- `hypotheses.md` — baseline hypothesis freeze v1.0 (inherits WISE v1.0
  physics; ZTF-specific observer, bands, quality inputs)
- `notes/irsa_recon.md` — verified facts about the IBE metadata search,
  product tree, cutout service, mask bits and catalogs
- `scripts/` — pipeline entry points, in execution order:
  `ztf_corridors.py` (config) → `coarse_discovery.py` → `precise_pass.py`
  → `catalog_screen.py` + `screen_recurrence.py` → `fetch_cutouts.py` →
  `sample_tensor.py` → `injection_calibrate.py`; `asteroid_control.py`
  is the independent positive control. Each run writes
  `run_config_*.json` under `../../runs/ztf/<run-id>/` (gitignored)
- `results/` — compact tracked summaries per stage
- `targets/` — ZTF overlay of the universal list (scale-up; not pilot)

Adapter: `sglsurvey/adapters/irsa_ztf.py`. Registry: shared
`registries/pilot_wise_2026.yaml` (survey-agnostic; no ZTF entries
needed for the pilot).

## Pilot step status

1. Freeze hypotheses — **v1.0 frozen 2026-08-20**
2. Corridors — `ross-128`, `eps-ind-a`, `proxima-cen` (see hypotheses §1)
3. Coarse discovery — **run** (`coarse_v1`: 7,930 public quadrant-
   exposures, 1,896 hits; all three antipodes in/beside primary-grid
   CCD gaps)
4. Precise pass — **run** (`precise_v1`: 1,480 usable evaluations;
   ≈10 % of metadata rows have no served products)
5. Catalog screening — **run** (`screen_v1`: 7,786 psfcat matches; one
   static star resolved; nothing track-following)
6. Forced photometry / stack — **run** (797 hybrid diff/sci cutout sets;
   variance model Gaussian to 2 %)
7. Injection calibration + asteroid positive control — **run**
   (`calib_v1`, AnalysisRun `run-1f02c325fd66`: 80 Constraints, 0
   Candidates; control (60000) recovered, flux-scale fix + 0.5 mag
   throughput correction applied)
8. Report — **done**: `../../report/ztf_pilot_v1.md`

Stage details: `results/pilot_v1_summary.md`.
