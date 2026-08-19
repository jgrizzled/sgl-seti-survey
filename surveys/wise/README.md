# surveys/wise — WISE/NEOWISE shakedown

First milestone: a completeness-calibrated end-to-end run of Pipeline A
against the WISE merged L1b products (plan §4).

## Contents

- `hypotheses.md` — the baseline hypothesis freeze (draft until signed off;
  open decisions are marked `[DECISION]`)
- `notes/irsa_recon.md` — reconnaissance of the IRSA tables, services, and
  data-product identifiers this sub-project will query
- `notes/pilot_registry_shortlist.md` — candidate pilot endpoints and what
  each one exercises
- `configs/` — frozen run configurations (target registry ref, sglseti
  pin, table/release identifiers, query parameters)
- `scripts/` — runnable pipeline entry points (e.g.
  `coarse_discovery.py`); every run writes its own `run_config_*.json`
  under the run directory
- `results/` — compact, reviewable summaries of runs (tracked); bulky
  generated products go to `../../runs/wise/<run-id>/` (gitignored)

## Dependency pins

- `sglseti`: adjacent checkout `../../../sglseti`, pin at commit `21f6f3d`
  (v1.1.0 + type fix), clean working tree as of 2026-08-18. Every run
  config must record the commit and dirty state at run time.

## Data policy

Disk budget: up to a few hundred GB locally for downloaded frames.
Everything under `../../runs/` is gitignored; any bulk download must be
reproducible from a tracked script/config plus the snapshotted query, so
the data can always be re-fetched rather than committed.

## Shakedown step status (plan §4)

1. Freeze baseline hypotheses — **frozen v1.0 (2026-08-18)**
2. Curate pilot registry — **curated v1.0 (2026-08-18)**
   (`registries/pilot_wise_2026.yaml`, 5 systems / 7 component
   endpoints, sglseti-validated, hash `sha256:82743c09…`). All CURATION
   flags resolved — Alpha Cen revalidated against Kervella et al. 2016,
   Sirius hip2-as-barycenter interpretation validated against 2MASS to
   0.16", diagonal covariance by declared assumption — see
   `notes/registry_curation.md`. Corridor viability confirmed for all
   five systems (`notes/corridor_coverage.md`: 470–630 W1/W2 epochs over
   ~14 yr each, plus cryo W3/W4).
3. Snapshot frame discovery (IRSA TAP/IBE) — **adapter implemented and
   run** (`sglsurvey/adapters/irsa_wise.py` +
   `scripts/coarse_discovery.py`): sglseti discovery cones → snapshotted
   TAP queries against the merged L1b inventory → Observation records →
   coarse locus-vs-nominal-WCS IntersectionEvaluations (hits and misses
   retained). Fetch path (IBE cutouts + md5-verified full products)
   smoke-tested. Results: `results/coarse_v1_summary.md`.
4. Precise pass (loci × exact WCS + mask usable-pixel tests via
   `covered_z_intervals`) — **implemented**
   (`scripts/precise_pass.py`): -msk products (md5-verified; the mask
   header carries the full frame WCS) define usable pixels with the
   Explanatory Supplement fatal-bit set; emits precise-stage
   IntersectionEvaluations with covered relay-distance intervals
   (disjoint when masks split the locus) and usable-pixel fractions.
   Results: `results/precise_v1_summary.md`.
5. Screen catalog products — **run** (`scripts/catalog_screen.py` +
   `scripts/screen_recurrence.py`): 43,915 ScreenMatch records from 140
   snapshotted queries across all corridors and mission phases; all
   recurrence peaks resolved as static background stars via the
   parallax-phase test → defensible layer-1 null. Results:
   `results/screen_v1_summary.md`.
6. Search images — **run** (`scripts/fetch_cutouts.py`,
   `sglsurvey/photometry.py`, `scripts/forced_stack.py`): matched-filter
   forced photometry on 4,280 cutouts stacked over a joint 64-z × 5×5-µ
   grid per endpoint × role × band, with parallax-phase-split stacks
   and offset-trajectory controls; AnalysisRun `run-2f0f4c6fd24a`.
   No detection; 2 marginal cells flagged for stage-7 adjudication.
   Results: `results/stack_v1_summary.md`.
7. Calibrate via injections and controls — **run**
   (`scripts/sample_tensor.py` + `scripts/injection_calibrate.py`,
   AnalysisRun `run-88887a3ac8d9`): 8 offset-control trajectories per
   endpoint × role define predeclared thresholds (FAR < 1/8 per grid
   search); analytic Gaussian-source injections on a 1/z-uniform grid
   (fixes stack_v1's low-z undersampling) yield 90%-recovery depths and
   **448 ledger-ready Constraint records**; both stack_v1 marginal
   cells vetoed (Candidate records). No surviving candidate. Results:
   `results/calib_v1_summary.md`.
8. Report — not started
