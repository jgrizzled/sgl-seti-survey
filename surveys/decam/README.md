# decam

DECam/NOIRLab southern survey (project plan §10 step 8, first
archive-family expansion on the v2 design): the 15 southern corridors
(antipode Dec < −30°) that PS1 and ZTF cannot reach, including the
engineering-backbone picks σ Dra, HD 219134 and Lalande 21185.

Recon facts (service probes, coverage sweep, `astro-datalab` client
assessment): `notes/decam_recon_2026-08-24.md` at repo level.

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
- Hypotheses: `hypotheses.md` is a **draft** (v0.1). The geometric
  stages below are hypothesis-light and reusable; the freeze (v1.0,
  with dev/confirmatory split) must happen before any screening
  threshold, stack, injection or candidate rule runs. This survey is
  the first with no v1 exploratory phase: v2 discipline from day one.

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
   image + wtmap HDUs per usable exposure; unit of retrieval is the
   CCD containing the locus centre — multi-CCD loci lose the off-CCD
   arc, to revisit at freeze).
5. Stacking / injection / candidate rules: not yet built — freeze
   `hypotheses.md` v1.0 first (v2 engine profile).

Pilot Pipeline A state (2026-08-24): 306 Observations, 612 coarse /
412 precise evaluations (222 hit exposures, ~all usable), 11,430
ScreenMatch records, recurrence triage clean.
