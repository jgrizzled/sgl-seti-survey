# SGL SETI Survey

Archival search for SGL SETI signals with the [sglseti](../sglseti) geometry
engine. See `notes/project_plan.md` for the full plan.

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

Status (2026-08-24): the v2 programme is complete — WISE, ZTF,
Pan-STARRS1, the joint PS1+ZTF stage and SPHEREx each searched under a
frozen decision rule with a blind confirmatory hold-out (family-wise
α = 0.05); **0 candidates** in every survey. Canonical reports:
`report/{wise,ztf,ps1,spherex}_survey.md`, `report/joint_ps1_ztf.md`.
