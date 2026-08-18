# SGL SETI Survey

Archival search for SGL SETI signals with the [sglseti](../sglseti) geometry
engine. See `notes/project_plan.md` for the full plan.

## Layout

| Path | Purpose | Tracked |
| ---- | ------- | ------- |
| `notes/` | Cross-project research and planning | yes |
| `docs/` | Methods docs, record schemas, ADRs shared across archives | yes |
| `sglsurvey/` | Shared Python package: core records, adapter interfaces, sglseti integration, injection framework | yes |
| `registries/` | Curated target registries (sglseti registry format), shared across archives | yes |
| `surveys/<archive>/` | One sub-project per archive family: frozen hypotheses, configs, sub-project notes, compact reviewable results | yes |
| `runs/<archive>/<run-id>/` | Generated products: query snapshots, downloaded data, analysis outputs | no (gitignored) |
| `report/` | Aggregate cross-archive results and publication drafts | yes |

Sub-projects are organized **per archive** (`surveys/wise/`, later
`surveys/spherex/`, …), not per pipeline: Pipeline A (corridor discovery +
precise intersection) and Pipeline B (beam-axis proximity) are stages in the
shared `sglsurvey` package that every archive sub-project invokes, and the
archive adapter, downloaded data, and footprint logic are shared between
both pipelines for a given archive.

## Current status

WISE/NEOWISE shakedown (plan §4) in setup: see `surveys/wise/README.md`.
