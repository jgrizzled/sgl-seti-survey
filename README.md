# SGL SETI Survey

Archival search for SGL SETI signals with the [sglseti](../sglseti) geometry
engine.

## Layout

| Path                       | Purpose                                                                                                                            | Tracked         |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| `notes/`                   | Cross-project research and planning                                                                                                | yes             |
| `docs/`                    | Methods docs, record schemas, ADRs shared across archives                                                                          | yes             |
| `sglsurvey/`               | Shared Python package: core records, adapter interfaces, sglseti integration, injection framework                                  | yes             |
| `targets/`                 | Universal (survey-agnostic) target portfolio: census snapshots, ranking config + builder, versioned target list with basket labels | yes             |
| `registries/`              | Curated target registries (sglseti registry format), shared across archives                                                        | yes             |
| `surveys/<archive>/`       | One sub-project per archive family: frozen hypotheses, configs, sub-project notes, compact reviewable results                      | yes             |
| `runs/<archive>/<run-id>/` | Generated products: query snapshots, downloaded data, analysis outputs                                                             | no (gitignored) |
| `report/`                  | Aggregate cross-archive results and publication drafts                                                                             | yes             |

Sub-projects are organized **per archive** (`surveys/wise/`, later
`surveys/spherex/`, …), not per pipeline: Pipeline A (corridor discovery +
precise intersection) and Pipeline B (beam-axis proximity) are stages in the
shared `sglsurvey` package that every archive sub-project invokes, and the
archive adapter, downloaded data, and footprint logic are shared between
both pipelines for a given archive.

## Current status

27 archival surveys are complete (7 Pipeline A corridor surveys, 20
Pipeline B crossings surveys) with **zero candidates**; the remaining
plan rows wait on data releases (SPHEREx QR3, Rubin images, Gaia DR4,
WINTER / PGIR). The programme-wide index of every constraint and open
cell — 28,588 rows over 88 endpoints, with the shared Constraint /
AnalysisRun / Candidate records — is `report/programme_ledger.md`
(built 2026-09-10).

See `notes/project_plan.md` for per-survey status and
`notes/project_history.md` for the execution record.
