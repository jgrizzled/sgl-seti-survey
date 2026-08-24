# decam

DECam/NOIRLab southern survey (project plan §4.6, first
archive-family expansion on the v2 design): the 15 southern corridors
(antipode Dec < −30°) that PS1 and ZTF cannot reach, including the
engineering-backbone picks σ Dra, HD 219134 and Lalande 21185.

Recon facts (service probes, coverage sweep, `astro-datalab` client
assessment): `surveys/decam/notes/decam_recon_2026-08-24.md` at repo level.

- Adapter: `sglsurvey/adapters/noirlab_decam.py` (Astro Archive
  adv_search discovery, EXPNUM-joined instcal image/dqmask/wtmap,
  `?hdus=` single-CCD fetches, dqmask-borne exact TPV footprints).
- Focal plane: `configs/decam_focal_plane_v1.json` — static per-CCD
  tangent-plane layout derived from a pinned reference dqmask,
  validated against 2012/2024-era exposures
  (`scripts/build_focal_plane.py`); nominal footprints need no
  per-exposure download.
- Corridors: `scripts/decam_corridors.py` (15 southern; pilot =
  lalande, sigmadra, hd219134).
- Hypotheses: `hypotheses.md` **v1.0 frozen 2026-08-24**
  (`configs/v2_freeze.json`, seed 20260824): dev = lalande + sigmadra
  + hd219134 (forced pilots) + struve2398 (drawn); confirmatory = 11
  corridors / 14 endpoints. First survey with no v1 exploratory
  phase: v2 discipline from day one.
- Positive control: (60000) Miminko
  (`configs/asteroid_control_v1.json`,
  `scripts/asteroid_control_v2.py`).
- Per-frame star ZP: `scripts/calibrate_zeropoints.py` →
  `runs/decam/zeropoints.jsonl` (header MAGZERO never used).

Stages (Pipeline A):

1. `scripts/coarse_discovery.py [endpoint ...]` →
   `runs/decam/coarse_v1/` (Observations + coarse
   IntersectionEvaluations; defaults to the pilot corridors).
2. `scripts/precise_pass.py [endpoint ...]` → `runs/decam/precise_v1/`
   (full dqmask per hit exposure, covered z-intervals, cutout index).
3. `scripts/catalog_screen.py [corridor ...]` → `runs/decam/screen_v1/`
   (NSC DR2 via `astro-datalab`, object-cone → meas-by-objectid;
   NSC is time-partial — see `catalog_stats.json`), then
   `scripts/screen_recurrence.py` (bin occupancy + fixed-z point
   filter triage).
4. `scripts/fetch_cutouts.py` → `runs/decam/products/cut/` (single-CCD
   image + wtmap HDUs per usable exposure; CCDs selected by the
   covered locus — multi-CCD arcs fetch every crossing CCD, and the
   flux-map builder mosaics them per hypotheses §8.7).
5. v2 engine (`profile.py` + `run.py`):
   `uv run python surveys/decam/run.py
   {freeze,geometry,build,nulls,inject,completeness,adjudicate,report}`
   — development set first, blind confirmatory only after the rule is
   frozen on dev.

Survey state (2026-08-24): **complete** — development (30 cells) and
blind confirmatory (86 cells) both 0 candidates at FWER α = 0.05;
1,904 Constraints; canonical report `report/decam_survey.md`, tables
in `results/report_tables.md`. Confirmatory Pipeline A driver:
`scripts/run_confirmatory.sh`. Yearly archive refresh re-runs
discovery with an extended TIME_RANGE stop.
