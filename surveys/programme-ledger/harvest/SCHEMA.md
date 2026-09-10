# Programme ledger v2 — harvest schema

One YAML file per survey report (`harvest/<report_stem>.yaml`). The harvest is
a **reviewable index of what each report states**, never a re-analysis: every
number must be quotable from the report or from a results file the report
cites. Unknown → `null`, never a guess.

```yaml
survey_id: ztf_crossings            # report file stem
report: report/ztf_crossings.md
plan_section: "5.1"                 # notes/project_plan.md section
pipeline: A | B                     # A = corridor (persistent relay track), B = crossings
archives: [ztf]                     # lower-case archive ids
hypothesis_version: "v1.0"          # plus amendments, e.g. "v1.0 + A1–A5"
era: "2018-03 → 2026-03"            # as the report states it
run_dir: runs/ztf-crossings         # null if none
records_dir: null                   # e.g. runs/ztf/v2/records if constraint.jsonl exists
candidates: 0                       # final candidate count stated in the report
candidate_notes: null               # one line if any retained/ambiguous items exist
rows:                               # one row per (target, channel, rung, band) cell the report speaks to
  - target_id: teegarden            # registry id (see vocabulary); "programme" for target-agnostic
    channel: B                      # A | B | S1 | S2 | corridor  (see vocabulary)
    rung: "0.1 AU"                  # crossings: "1.2 Rsun" | "2.5 Rsun" | "0.1 AU" | "0.90 AU" | "0.95 AU" | "1.0 AU"
                                    # corridor: z interval string, e.g. "550-10000 AU"
    band: "zr"                      # instrument band / line / energy range as the report names it; "any" if none
    substrate: "ZTF difference images"   # short free text: what data was searched
    status: searched                # vocabulary below
    limit_kind: "m90"               # what the limit is: m90 | flux | eirp | power | rate | fluence | null
    limit_value: 21.6               # number in native unit (or null)
    limit_unit: "AB mag"            # native unit string
    limit_range: null               # "min–max" string if the report gives only a range over cells
    power_mw: null                  # derived transmitter power in MW if the report gives one (number)
    power_mw_range: "0.6–24"        # if the report gives only a range
    n_events: 3                     # windows / events / cells this row aggregates (int or null)
    n_trials: null
    epoch_range: null               # as stated, if the row is time-bound
    source: "report/ztf_crossings.md §5 table 2"   # exact citation: report section, or results file path
    notes: null                     # one line
```

## Status vocabulary (carry the programme's existing states; do not invent new ones)

- `searched` — calibrated search with a measured depth (recovery curve / m90 / threshold-calibrated limit)
- `constraint_only` — searched but no calibrated null distribution / threshold; limit is a threshold statement
- `ledger_only` — coverage exists and is logged, no statistic possible (e.g. no temporal signature)
- `not_constrainable` — the archive cannot observe that cell in principle (Gate I/II, occulter, elongation)
- `structurally_open` — cell is in scope but no usable data (uncovered, epoch-starved, single epoch, track masked, saturated)
- `vetoed_known_source` — exceedance attributed to a catalogued source
- `retained_ambiguous` — exceedance retained, not promotable
- `no_survey` — plan row closed/blocked/deferred without a search

## Channel vocabulary

- `A` — uplink interception at the target star near opposition (Earth channels)
- `B` — downlink pre-lens interception at the star's antipode near the anti-solar point
- `S1` — downlink post-lens, sunward (target side): Earth inside the beam after the solar graze
- `S2` — uplink past the Sun, sunward (anti-target side): star near conjunction
- `corridor` — Pipeline A persistent relay track on the focal line (antipode corridor); use `role`-like detail in `notes` if the report separates tx/rx

Rungs are beam radii at the observer; corridor rows give the z interval instead.

## Target ids

Use the registry ids exactly as the crossing lists spell them (e.g. `wolf-359`,
`gj-1276`, `teegarden`, `van-maanen`, `eps-ind-b`, `gj65-b`, `wise-0855`). A row
whose statement covers several named targets may list them as a YAML list under
`target_id`. A row covering "all N targets" uses `target_id: programme` with
`n_events`/`notes` giving N.

## Granularity

Per-target rows wherever the report (or a results table it cites under
`surveys/<name>/results/`) gives per-target numbers. Aggregate rows otherwise,
with `limit_range` / `power_mw_range` and the cell count. Always include the
rows for the cells the report declares `not_constrainable`, `structurally_open`
or `ledger_only`: the open-cell census is half the point of the ledger.

## Layout

- `surveys/programme-ledger/harvest/*.yaml` — the transcribed files (tracked): all
  Pipeline B crossings reports plus `plan_no_survey.yaml`.
- `runs/programme-ledger/v2/harvest_generated/*.yaml` — the seven Pipeline A files
  (one row per Constraint record, ~24 MB), regenerated by
  `scripts/harvest/pipeline_a_*.py`; not tracked.
- Generators for the 22 generated files live in `scripts/harvest/`; each takes an
  optional output directory argument and reproduces its files byte for byte.
- `scripts/build_ledger_v2.py` reads both directories; `scripts/render_tables.py`
  renders the report tables.
