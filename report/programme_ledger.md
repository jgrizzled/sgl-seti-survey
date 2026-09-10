# Programme ledger v2 — every survey's constraints and open cells in one index

*Built 2026-09-10 from the 27 completed survey reports plus the plan's
closed / blocked / deferred rows. Ledger version 2.0.0. Products under
`surveys/programme-ledger/` (harvest, scripts, results) and
`runs/programme-ledger/v2/` (records, manifest). Plan goal 5; supersedes
nothing — the joint covered-window ledger v1 (`report/joint_crossings.md`,
593 rows over ZTF/PS1/WISE) stays as pre-registered and is folded in here
as one of 28 inputs.*

## 1. Question and answer

**Question.** After 27 surveys over two pipelines and 88 endpoints, which
cells of the search grid — target × channel × beam-radius rung (or
corridor z-interval) × band — carry a calibrated constraint, which carry
only coverage, which are open, and which no archive can observe? And are
those statements available as the stable Constraint / AnalysisRun /
Candidate records the plan promised, rather than as 27 differently shaped
results directories?

**Answer.** The ledger holds **28,588 rows from 28 sources** (27 reports +
the plan's no-survey rows); **0 candidates** programme-wide. By coverage
state: 20,541 `searched` (a calibrated depth), 191 `constraint_only`, 344
`ledger_only`, 7,301 `structurally_open`, 87 `not_constrainable`, 97
`no_survey`, plus 16 `retained_ambiguous` and 11 `vetoed_known_source`
exceedance rows. The Pipeline A corridor half is dense and uniform: all 88
endpoints have searched corridor cells, 61 of them in six surveys, and the
per-endpoint depths are in Table 5. The Pipeline B crossings half is
sparse by construction: **102 distinct (target, channel, rung) cells are
searched across 37 targets**, 35 of those targets carry at least one
published transmitter-power limit, and the per-target matrix (Table 3) is
the first place every sunward, antipode and uplink cell of a given star
can be read together. Records: 28 AnalysisRun, 21,163 Constraint and 27
Candidate records in the shared schema, content-addressed, each citing
its report section or results file.

Three things the merge shows that no single report could:

1. **Fourteen registry targets have no per-target crossings row in any
   report** — 82-eri, alpha-cen-a/b, gj-1061, gj-367, gj-66-a/b, gj-687,
   kapteyn-star, lp-145-141, proxima-cen, sigma-dra, struve-2398-a/b.
   They are covered only by programme-wide statements (WISE Gate I,
   SPHEREx deferred, the DASCH/Rubin/ATLAS no-pull rows): far-southern or
   ecliptic-pole systems the northern optical crossings surveys never
   reached and the sunward surveys' unit cuts never selected. A further 37
   named targets have coverage rows but no searched crossings cell. The
   88-endpoint list is therefore fully covered in Pipeline A and 37/88
   covered in Pipeline B.
2. **The A 1.0 AU rung is ledger-only almost everywhere** (249 ledger-only
   rows; the one searched unit is ZTF's wolf-1069 g-band window, S 0.47
   against T 2.75, at m90 17.6). Every other crossings survey declared it
   coverage-without-statistic, as the plan's standing deferral says; the
   ledger now shows the size of that cell.
3. **The S1 grazing rungs are `not_constrainable` in every substrate but
   one**: LASCO C2 at 2.5 R☉ (four targets searched, MW-class); the 1.2 R☉
   S1 rung is behind an occulter or mission floor in LASCO, STEREO, WISPR
   and SoloHI alike. The B grazing rungs, by contrast, are searched on
   4–5 targets at W-class optical depths (ZTF, PS1, PTF, DASCH, TESS) and
   at 0.1–0.3 MW in gamma rays (LAT).

The twelve retained-ambiguous cells (Table 7) are the complete list of
exceedances the programme carries forward; none is promotable and each
report's disposition is quoted on its row.

## 2. Design

**Harvest, not re-analysis.** One YAML file per report (schema in
`surveys/programme-ledger/harvest/SCHEMA.md`) transcribes what the report
states — per-target where the report or a results file it cites gives
per-target numbers, aggregate otherwise. The 21 crossings / plan files
(0.9 MB) are tracked under `harvest/`; the seven Pipeline A files carry
one row per Constraint record (27,383 rows, 24 MB) and are regenerated
into `runs/programme-ledger/v2/harvest_generated/` by
`scripts/harvest/pipeline_a_*.py`. Every generated file — 22 of the 28 —
has a generator under `scripts/harvest/` that reproduces it byte for
byte; the six others (TESS, GALEX, DASCH, Rubin, ATLAS/ASAS-SN, plan
rows) are hand transcriptions. Every row carries `source` (report section, or the results /
records path plus record id), the native limit (`limit_value`,
`limit_unit`, `limit_kind`), and the derived transmitter power in MW only
where the report itself gives one. Nothing is recomputed; unit conversions
are numeric only (W/kW/GW → MW, noted per row). Six harvest passes ran in
parallel over report groups; their mapping decisions are listed in §5.

**States** (the programme's existing vocabulary, unchanged): `searched`
(calibrated depth), `constraint_only` (threshold statement, no calibrated
null or no depth), `ledger_only` (coverage without statistic),
`not_constrainable` (unobservable in principle — Gate I/II, occulter,
elongation, mission floor), `structurally_open` (in scope, no usable data:
uncovered, epoch-starved, single epoch, track-masked, saturated,
void null), `vetoed_known_source`, `retained_ambiguous`, `no_survey`.

**Channels and rungs.** A (uplink interception at the star near
opposition), B (downlink pre-lens at the antipode), S1 (downlink
post-lens, sunward), S2 (uplink past the Sun, conjunction), `corridor`
(Pipeline A focal-line track; rung = z interval). Rungs 1.2 R☉ / 2.5 R☉ /
0.1 AU / 0.90 AU / 0.95 AU / 1.0 AU. The LAT survey is the one place
channel A carries grazing rungs.

**Merge** (`scripts/build_ledger_v2.py`): list-valued target / channel /
rung fields explode into one row per combination (the harvest row id and
group size are kept); vocabulary and registry ids are validated (build
fails on a violation; unusual channel–rung pairs warn); the flat table is
written as ECSV with a `rung_class` column; the per-target matrix takes,
per (target, channel, rung), the best state in the order searched >
retained-ambiguous > vetoed > constraint-only > ledger-only > open >
no-survey > not-constrainable, and the lowest published MW limit among
searched / constraint-only rows. Records: one AnalysisRun per harvest
file (pipeline `programme-ledger-harvest`, observation-set hash = the
report's sha256, environment = both repo commits); a Constraint per row in
a constraining state (kind `recovery_curve` for searched rows with a
limit, `qualified_null` for constraint-only / ledger-only / searched
without a limit, `not_constrainable`), with `z_interval_au` = the corridor
interval or the 550–10,000 AU prior and the full harvest row in `extra`;
a Candidate per retained-ambiguous / vetoed row. Structurally-open and
no-survey rows are coverage states, not constraints, and live only in the
table. Manifest with sha256 of every harvest file, record file and the
ledger.

**Render** (`scripts/render_tables.py`) writes the seven tables below to
`results/programme_ledger_v2_tables.md`; Tables 1–3 and 5–7 are inlined
here, Table 4 (the 236 programme-wide statements) is in the results file.

**Refresh.** Re-run the generators (or edit a hand-transcribed YAML),
then `build_ledger_v2.py` and `render_tables.py`; a re-run survey needs
only its generator's source paths pointed at the new records / results.
The full 28,588-row table lives with the records under `runs/`; the
tracked results directory keeps the 1,182-row crossings + programme-wide
subset, the matrix, the summary and the tables.

## 3. Contents

### Table 1 — survey index

| survey | plan § | pipeline | rows | candidates | states |
|---|---|---|---|---|---|
| `wise_survey` | 4.1 | A | 5640 | 0 | searched 4831, structurally_open 809 |
| `ztf_survey` | 4.2 | A | 2683 | 0 | searched 2376, structurally_open 288, no_survey 19 |
| `spherex_joint6` | 4.3 | A | 1412 | 0 | searched 1332, structurally_open 80 |
| `spherex_survey` | 4.3 | A | 8462 | 0 | searched 8150, structurally_open 312 |
| `ps1_survey` | 4.4 | A | 5552 | 0 | searched 971, structurally_open 4562, no_survey 19 |
| `joint_ps1_ztf_wise` | 4.5 | A | 2667 | 0 | searched 2030, ledger_only 1, structurally_open 615, no_survey 21 |
| `decam_survey` | 4.6 | A | 967 | 0 | searched 512, structurally_open 455 |
| `plan_no_survey` | 4.7–4.12, 5.15, 5.20–5.22, 5.24 | – | 42 | 0 | ledger_only 5, structurally_open 4, no_survey 23, not_constrainable 10 |
| `ztf_crossings` | 5.1 | B | 76 | 0 | searched 17, constraint_only 27, structurally_open 27, no_survey 2, not_constrainable 3 |
| `ptf_crossings` | 5.10 | B | 26 | 0 | searched 7, constraint_only 1, ledger_only 1, structurally_open 13, no_survey 2, not_constrainable 2 |
| `lasco_crossings` | 5.11 | B | 26 | 0 | searched 14, retained_ambiguous 1, constraint_only 6, no_survey 1, not_constrainable 4 |
| `galex_crossings` | 5.12 | B | 28 | 0 | searched 11, retained_ambiguous 1, constraint_only 4, ledger_only 1, structurally_open 11 |
| `rubin_crossings` | 5.13 | B | 14 | 0 | constraint_only 1, structurally_open 11, no_survey 1, not_constrainable 1 |
| `dasch_crossings` | 5.14 | B | 21 | 0 | searched 18, structurally_open 1, no_survey 1, not_constrainable 1 |
| `stereo_hi_crossings` | 5.16 | B | 18 | 0 | searched 9, retained_ambiguous 1, vetoed_known_source 1, constraint_only 4, ledger_only 1, not_constrainable 2 |
| `spectral_archives` | 5.17 | B | 83 | 0 | searched 38, constraint_only 28, structurally_open 17 |
| `wispr_crossings` | 5.18–5.19 | B | 46 | 0 | searched 25, retained_ambiguous 1, vetoed_known_source 5, constraint_only 4, ledger_only 2, structurally_open 1, no_survey 1, not_constrainable 7 |
| `ps1_crossings` | 5.2 | B | 39 | 0 | searched 11, retained_ambiguous 1, vetoed_known_source 1, constraint_only 6, ledger_only 1, structurally_open 15, no_survey 2, not_constrainable 2 |
| `solohi_crossings` | 5.23 | B | 37 | 0 | searched 18, constraint_only 13, ledger_only 3, structurally_open 1, not_constrainable 2 |
| `skirt_crossings` | 5.25 | B | 59 | 0 | searched 26, retained_ambiguous 3, ledger_only 24, structurally_open 2, no_survey 4 |
| `highenergy_crossings` | 5.26 | B | 218 | 0 | searched 81, vetoed_known_source 3, constraint_only 57, ledger_only 53, structurally_open 12, not_constrainable 12 |
| `wise_crossings` | 5.3 | B | 104 | 0 | ledger_only 66, structurally_open 1, not_constrainable 37 |
| `joint_crossings` | 5.4 | B | 238 | 0 | searched 26, retained_ambiguous 1, vetoed_known_source 1, constraint_only 29, ledger_only 151, structurally_open 26, not_constrainable 4 |
| `radio_crossings` | 5.5 | B | 14 | 0 | ledger_only 6, structurally_open 8 |
| `tess_crossings` | 5.6 | B | 21 | 0 | searched 11, retained_ambiguous 1, constraint_only 2, ledger_only 1, structurally_open 6 |
| `radio_crossings_ext` | 5.7 | B | 26 | 0 | ledger_only 14, structurally_open 12 |
| `radio_quicklook` | 5.8 | B | 17 | 0 | ledger_only 14, structurally_open 3 |
| `atlas_asassn_crossings` | 5.9 | B | 52 | 0 | searched 27, retained_ambiguous 6, constraint_only 9, structurally_open 9, no_survey 1 |

### Table 2 — row census by coverage state and channel

| state | corridor | A | B | S1 | S2 | total |
|---|---|---|---|---|---|---|
| searched | 20202 | 100 | 147 | 35 | 57 | 20541 |
| retained_ambiguous | 0 | 1 | 9 | 1 | 5 | 16 |
| vetoed_known_source | 0 | 0 | 5 | 3 | 3 | 11 |
| constraint_only | 36 | 107 | 21 | 16 | 11 | 191 |
| ledger_only | 53 | 249 | 12 | 2 | 28 | 344 |
| structurally_open | 7121 | 60 | 116 | 0 | 4 | 7301 |
| no_survey | 65 | 8 | 12 | 3 | 9 | 97 |
| not_constrainable | 0 | 46 | 20 | 18 | 3 | 87 |

### Table 3 — per-target crossings matrix (best state; best transmitter-power limit where one is published)

Glyphs: S searched (calibrated depth), C constraint-only, L ledger-only, open = structurally open, – no survey, × not constrainable; a trailing * marks a retained-ambiguous exceedance on the cell and v an exceedance vetoed as a known source (S* / Sv alone: the cell has only that disposition). The power is the lowest published transmitter-power limit among the cell's searched or constraint-only rows, whatever the band or statistic. Blank = no report speaks to the cell for that target (programme-wide rows are in Table 4).

| target | A 0.1 AU | A 1.0 AU | B 1.2 Rsun | B 2.5 Rsun | B 0.1 AU | S1 1.2 Rsun | S1 2.5 Rsun | S1 0.1 AU | S2 0.1 AU | S2 0.90 AU | S2 0.95 AU | S2 1.0 AU |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 61-cyg-a |  | L |  |  |  |  |  |  |  |  |  |  |
| 61-cyg-b |  | L |  |  |  |  |  |  |  |  |  |  |
| 61-vir |  | L |  |  |  |  |  | C | S* |  |  |  |
| barnard-star |  | L |  |  |  |  |  |  |  |  |  |  |
| eps-eri |  | L |  |  |  |  |  |  |  |  |  |  |
| eps-ind-a |  | L |  |  |  |  |  |  |  |  |  |  |
| eps-ind-b |  | L |  |  |  |  |  |  |  | S 850 kW | S 600 kW | L |
| ez-aqr |  | L |  |  |  |  |  | S 48 MW | S 70 MW |  |  |  |
| fomalhaut |  | L |  |  |  |  |  |  | open |  |  |  |
| gj-1002 |  | C |  |  |  |  |  | S 83 MW | S 1.4 GW |  |  |  |
| gj-1087 |  | L |  |  |  |  |  | × | Sv | S 8.5 MW | S 6 MW | L |
| gj-11068 |  | L |  |  |  |  |  |  |  | S | S | L |
| gj-1111 |  | L |  |  |  |  |  | S 3.9 GW | S | open | open | L |
| gj-11547 |  |  |  |  |  |  |  |  |  | S 14 MW | S 15 MW | L |
| gj-1221 |  | C |  |  |  |  |  |  |  |  |  |  |
| gj-12724 |  |  |  |  |  |  |  |  |  | S 8.5 MW | S | L |
| gj-1276 | S 32 kW | C | S 49 kW | S 210 kW | S* 6 kW | × | S 6.7 MW | S 12 MW | S* 61 MW | S 1.4 MW | S 1.5 MW | L |
| gj-13157 |  | C |  |  |  |  |  |  |  | L | S 24 MW | L |
| gj-2012 |  | C |  |  |  |  |  |  |  | S | S 2.4 MW | L |
| gj-2066 |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-229-a |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-251 |  | L |  |  |  |  |  | S 520 MW | S 350 MW |  |  |  |
| gj-293 |  |  |  |  |  |  |  |  |  | L | L | L |
| gj-3112 |  |  |  |  |  |  |  |  |  | L | S 15 MW | L |
| gj-318 |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-3306 |  | C |  |  |  |  |  |  |  | S 3.4 MW | S | L |
| gj-338-a |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-338-b |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-3512 |  | C |  |  |  |  |  |  |  | S 14 MW | S* 6 MW | L |
| gj-514 |  | L |  |  |  |  |  | S 1.7 GW | S 1.2 GW |  |  |  |
| gj-518 |  | L |  |  |  |  |  | S 500 MW | S 1.4 GW | S | S 3.8 MW | L |
| gj-526 |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-54 |  | L |  |  |  |  |  | × | S 84 MW |  |  |  |
| gj-581 |  | L |  |  |  |  |  | Sv 710 MW | S 1.1 GW |  |  |  |
| gj-588 |  |  |  |  |  |  |  | S 140 MW |  |  |  |  |
| gj-625 |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-667-c |  | L |  |  |  |  |  | Sv 340 MW | S 3.2 GW |  |  |  |
| gj-674 |  | L |  |  |  |  |  | S 23 MW |  |  |  |  |
| gj-682 |  | L |  |  |  |  |  | S 470 MW |  |  |  |  |
| gj-783 |  | L |  |  |  |  |  | S 14 MW | Sv |  |  |  |
| gj-784 |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-832 |  | L |  |  |  |  |  |  |  |  |  |  |
| gj-876 |  | L |  |  |  |  |  | S | S 440 MW |  |  |  |
| gj-908 | S 8.5 kW | L |  |  | S 26 kW |  |  | C | C |  |  |  |
| gj-915 |  | L |  |  |  |  |  |  |  | S 8.5 MW | S 6 MW | L |
| gj-9193 |  | L |  |  |  |  |  |  |  | S 2.1 MW | S 3.8 MW | L |
| gj65-a |  | L |  |  |  |  |  |  |  |  |  |  |
| gj65-b |  | C |  |  |  |  |  |  |  |  |  |  |
| groombridge-34-a |  | L |  |  |  |  |  |  |  |  |  |  |
| groombridge-34-b |  | L |  |  |  |  |  |  |  |  |  |  |
| hd-219134 |  | L |  |  |  |  |  |  |  |  |  |  |
| lacaille-8760 |  | L |  |  |  |  |  |  | S 950 MW |  |  |  |
| lacaille-9352 |  | L |  |  |  |  |  |  |  |  |  |  |
| lalande-21185 |  | L |  |  |  |  |  |  |  |  |  |  |
| lhs-1723 |  | L |  |  |  |  |  |  |  |  |  |  |
| ltt-1445-a |  | L |  |  |  |  |  |  |  |  |  |  |
| luhman16-a |  |  |  |  |  |  |  |  |  | – | – |  |
| luhman16-b |  |  |  |  |  |  |  |  |  | – | – |  |
| luyten-star |  | L |  |  |  |  |  | × |  |  |  |  |
| procyon-a |  | L |  |  |  |  |  | × |  |  |  |  |
| procyon-b |  | L |  |  |  |  |  | × |  |  |  |  |
| ross-128 | S 870 W | L | × | S 110 W | Sv 1.3 kW |  | C | Sv | S 9.5 MW |  |  |  |
| ross-154 | S 4.4 kW | L |  |  | S 34 kW |  |  | C | Sv |  |  |  |
| ross-248 |  | L |  |  |  |  |  |  |  |  |  |  |
| sirius-a |  | L |  |  |  |  |  |  |  |  |  |  |
| sirius-b |  | L |  |  |  |  |  |  |  |  |  |  |
| tau-cet |  | L |  |  |  |  |  |  |  |  |  |  |
| teegarden | S 130 W | C | S 840 W | S 420 W | S* 35 kW | × | S 290 MW | S 12 MW | S 7.4 MW | S* | S* | L |
| van-maanen | S 3 MW | L | S*v 160 W | S*v 190 W | Sv 8 kW | × | S* 530 MW | S 29 MW | S 16 MW |  |  |  |
| wise-0855 |  | C 200 W |  |  |  |  |  |  |  |  |  |  |
| wolf-1061 |  | L |  |  |  |  |  | S 120 MW | S 110 MW |  |  |  |
| wolf-1069 |  | S |  |  |  |  |  |  |  | L | L | L |
| wolf-359 | S* 3.9 kW | C | S 90 W | S* 2.4 kW | S 41 kW | × | S 41 MW | S 17 MW | S 2.9 MW |  |  |  |
| wolf-437 |  | L |  |  |  |  |  | S 1.2 GW | S 790 MW |  |  |  |

Table 4 (programme-wide statements, 232 rows) is in `surveys/programme-ledger/results/programme_ledger_v2_tables.md`.

### Table 5 — corridor (Pipeline A) depth per endpoint: deepest published persistent-source m90 per survey (band)

WISE depths are Vega magnitudes; all others AB. A blank means the endpoint has no searched corridor cell in that survey.

| endpoint | `wise_survey` | `ztf_survey` | `ps1_survey` | `joint_ps1_ztf_wise` | `spherex_survey` | `spherex_joint6` | `decam_survey` | n surveys |
|---|---|---|---|---|---|---|---|---|
| 61-cyg-a | 12.1 (W1) |  |  |  | 19.6 (D4) | 20.2 (J6) | 22.2 (z) | 4 |
| 61-cyg-b | 12.9 (W1) |  |  |  | 19.7 (D4) | 20.0 (J6) | 22.3 (i) | 4 |
| 61-vir | 13.4 (W1) | 22.9 (zg) | 20.7 (z) | 22.8 (g) | 20.2 (D2) | 20.6 (J6) |  | 6 |
| 82-eri | 12.6 (W1) | 22.4 (zr) | 21.7 (g) | 22.3 (r) | 20.4 (D2) | 21.1 (J6) |  | 6 |
| alpha-cen-a | 11.6 (W1) | 22.8 (zr) | 21.1 (r) | 22.8 (r) | 18.8 (D6) | 18.6 (J6) |  | 6 |
| alpha-cen-b | 12.0 (W1) | 22.7 (zg) | 20.9 (r) | 22.8 (r) | 18.0 (D2) | 18.3 (J6) |  | 6 |
| barnard-star | 12.3 (W1) | 22.3 (zg) | 22.6 (g) | 22.3 (g) | 19.6 (D2) | 19.9 (J6) |  | 6 |
| eps-eri | 13.0 (W1) | 21.9 (zr) | 21.7 (g) | 22.0 (g) | 20.1 (D4) | 20.6 (J6) |  | 6 |
| eps-ind-a | 13.1 (W1) | 20.7 (zr) | 21.9 (g) | 21.7 (r) | 20.5 (D4) | 21.0 (J6) |  | 6 |
| eps-ind-b | 10.9 (W3) | 20.6 (zr) | 21.5 (r) |  | 19.8 (D4) | 19.9 (J6) |  | 5 |
| ez-aqr | 14.0 (W1) | 22.0 (zg) | 21.8 (g) | 21.9 (g) | 19.9 (D4) | 20.2 (J6) |  | 6 |
| fomalhaut | 13.6 (W1) | 22.4 (zg) | 21.5 (g) | 22.3 (g) | 20.2 (D4) | 21.0 (J6) |  | 6 |
| gj-1002 | 13.1 (W1) | 22.4 (zg) | 20.5 (z) | 22.1 (g) | 20.3 (D4) | 20.8 (J6) |  | 6 |
| gj-1061 | 14.2 (W1) | 21.7 (zr) | 20.6 (z) | 22.4 (g) | 20.9 (D1) | 21.2 (J6) |  | 6 |
| gj-1087 | 11.4 (W1) | 22.3 (zg) |  | 22.0 (r) | 19.1 (D1) | 19.5 (J6) |  | 5 |
| gj-11068 | 10.0 (W1) | 20.6 (zr) |  |  | 17.7 (D1) | 18.0 (J6) |  | 4 |
| gj-1111 | 12.3 (W1) | 20.4 (zr) | 21.7 (g) |  | 19.9 (D4) | 20.2 (J6) |  | 5 |
| gj-11547 | 12.6 (W1) | 23.4 (zg) | 21.0 (i) | 23.3 (r) | 20.3 (D3) | 21.2 (J6) |  | 6 |
| gj-1221 | 12.4 (W1) |  |  |  | 19.9 (D4) | 20.5 (J6) | 21.8 (g) | 4 |
| gj-12724 | 12.6 (W1) | 21.6 (zr) | 21.5 (g) | 21.9 (g) | 20.3 (D4) | 20.6 (J6) |  | 6 |
| gj-1276 | 12.9 (W1) | 22.2 (zr) | 20.7 (z) | 22.4 (r) | 20.2 (D1) | 20.7 (J6) |  | 6 |
| gj-13157 | 11.4 (W1) |  |  |  | 19.0 (D5) | 19.1 (J6) | 22.4 (i) | 4 |
| gj-2012 | 13.0 (W1) | 22.4 (zg) | 21.4 (i) | 22.2 (g) | 20.4 (D2) | 21.3 (J6) |  | 6 |
| gj-2066 | 12.6 (W1) | 22.6 (zr) | 21.8 (g) | 22.3 (g) | 19.0 (D1) | 19.4 (J6) |  | 6 |
| gj-229-a | 12.2 (W1) | 22.7 (zg) | 21.7 (i) | 22.4 (g) | 20.1 (D4) | 20.7 (J6) |  | 6 |
| gj-251 | 11.3 (W1) |  |  |  | 19.3 (D3) | 19.6 (J6) |  | 3 |
| gj-293 | 12.0 (W1) | 22.8 (zr) | 21.7 (g) | 22.8 (r) | 20.1 (D5) | 20.5 (J6) |  | 6 |
| gj-3112 | 12.5 (W1) | 22.7 (zg) | 21.7 (g) | 22.8 (g) | 20.8 (D2) | 21.3 (J6) |  | 6 |
| gj-318 | 11.6 (W1) | 22.1 (zg) |  | 22.3 (g) | 19.4 (D5) | 20.0 (J6) |  | 5 |
| gj-3306 | 13.0 (W1) | 22.5 (zr) | 22.0 (g) | 22.5 (r) | 20.2 (D3) | 20.8 (J6) |  | 6 |
| gj-338-a | 13.4 (W1) |  |  |  | 20.3 (D4) | 21.4 (J6) | 23.6 (g) | 4 |
| gj-338-b | 13.3 (W1) |  |  |  | 20.5 (D4) | 21.6 (J6) | 23.8 (g) | 4 |
| gj-3512 | 12.4 (W1) |  |  |  | 19.9 (D3) | 19.8 (J6) | 23.6 (g) | 4 |
| gj-367 | 11.9 (W1) |  | 21.8 (g) |  | 19.0 (D3) | 19.9 (J6) |  | 4 |
| gj-514 | 14.1 (W1) | 22.2 (zg) | 21.3 (r) | 22.2 (g) | 20.2 (D1) | 20.8 (J6) |  | 6 |
| gj-518 | 13.3 (W1) | 22.4 (zr) | 21.7 (g) | 22.1 (g) | 20.1 (D2) | 20.8 (J6) |  | 6 |
| gj-526 | 13.4 (W1) | 22.5 (zg) | 21.6 (g) | 22.1 (g) | 20.3 (D4) | 21.4 (J6) |  | 6 |
| gj-54 | 13.8 (W1) | 22.1 (zg) | 22.2 (g) | 22.0 (r) | 20.5 (D4) | 21.1 (J6) |  | 6 |
| gj-581 | 12.2 (W1) | 22.5 (zg) | 19.6 (y) | 22.4 (r) | 19.9 (D1) | 20.4 (J6) |  | 6 |
| gj-588 | 12.6 (W1) | 22.7 (zg) | 21.2 (i) | 22.4 (r) | 19.8 (D3) | 19.9 (J6) |  | 6 |
| gj-625 | 13.9 (W1) |  |  |  | 21.1 (D3) | 21.4 (J6) | 24.0 (g) | 4 |
| gj-66-a | 14.4 (W1) | 22.8 (zg) | 21.7 (g) | 22.5 (g) | 20.8 (D4) | 21.3 (J6) |  | 6 |
| gj-66-b | 14.4 (W1) | 22.7 (zr) | 21.3 (i) | 22.6 (r) | 20.8 (D4) | 21.5 (J6) |  | 6 |
| gj-667-c | 11.4 (W1) | 22.4 (zg) | 21.9 (g) | 22.1 (r) | 18.7 (D3) | 19.1 (J6) |  | 6 |
| gj-674 | 12.1 (W1) | 22.2 (zr) | 20.7 (i) | 22.0 (r) | 19.4 (D4) | 20.0 (J6) |  | 6 |
| gj-682 | 12.7 (W1) | 22.1 (zg) | 21.3 (g) | 22.4 (r) | 19.6 (D4) | 20.1 (J6) |  | 6 |
| gj-687 | 12.1 (W1) |  |  |  | 20.7 (D3) | 21.2 (J6) | 22.5 (r) | 4 |
| gj-783 | 13.3 (W1) | 22.5 (zg) | 21.4 (g) | 22.1 (r) | 19.7 (D4) | 19.9 (J6) |  | 6 |
| gj-784 | 12.6 (W1) | 22.2 (zr) | 21.4 (g) | 22.1 (r) | 20.0 (D4) | 20.8 (J6) |  | 6 |
| gj-832 | 13.2 (W1) | 22.8 (zr) | 23.0 (g) | 22.6 (r) | 20.6 (D4) | 20.9 (J6) |  | 6 |
| gj-876 | 12.7 (W1) | 22.3 (zr) | 21.5 (r) | 22.2 (r) | 20.0 (D4) | 21.2 (J6) |  | 6 |
| gj-908 | 13.8 (W1) | 22.1 (zg) | 21.9 (g) | 22.0 (r) | 20.2 (D2) | 20.8 (J6) |  | 6 |
| gj-915 | 12.1 (W1) | 21.8 (zg) | 21.2 (i) | 21.9 (g) | 20.2 (D3) | 20.8 (J6) |  | 6 |
| gj-9193 | 11.4 (W2) | 22.5 (zg) | 21.5 (r) | 22.2 (g) | 19.9 (D3) | 20.2 (J6) |  | 6 |
| gj65-a | 13.0 (W1) | 22.1 (zg) | 23.8 (i) | 22.1 (g) | 20.2 (D4) | 20.7 (J6) |  | 6 |
| gj65-b | 13.2 (W1) | 21.9 (zg) | 21.5 (i) | 22.1 (g) | 20.5 (D3) | 21.1 (J6) |  | 6 |
| groombridge-34-a | 12.2 (W1) |  |  |  | 20.2 (D4) | 20.6 (J6) | 23.0 (g) | 4 |
| groombridge-34-b | 12.2 (W1) |  |  |  | 20.1 (D3) | 21.2 (J6) | 23.2 (g) | 4 |
| hd-219134 | 11.3 (W1) |  |  |  | 19.3 (D4) | 19.3 (J6) | 24.0 (g) | 4 |
| kapteyn-star | 14.1 (W1) | 22.9 (zg) | 21.1 (i) | 23.0 (g) | 20.3 (D3) | 21.4 (J6) |  | 6 |
| lacaille-8760 | 14.1 (W1) | 22.6 (zg) | 21.9 (g) | 22.3 (r) | 20.2 (D3) | 20.7 (J6) |  | 6 |
| lacaille-9352 | 13.9 (W1) | 22.6 (zg) | 21.0 (i) | 22.6 (g) | 20.4 (D2) | 20.8 (J6) |  | 6 |
| lalande-21185 | 14.5 (W1) |  |  |  | 20.4 (D4) | 21.3 (J6) | 23.9 (g) | 4 |
| lhs-1723 | 7.9 (W3) | 22.7 (zg) | 21.9 (r) | 22.4 (g) | 20.3 (D2) | 20.3 (J6) |  | 6 |
| lp-145-141 | 12.0 (W1) | 22.7 (zg) |  | 22.4 (r) | 19.2 (D1) | 20.2 (J6) |  | 5 |
| ltt-1445-a | 12.9 (W1) | 22.1 (zg) | 21.6 (r) | 22.5 (r) | 20.4 (D2) | 21.0 (J6) |  | 6 |
| luhman16-a | 11.7 (W1) | 21.8 (zg) | 21.3 (r) | 22.1 (g) | 19.2 (D4) | 19.4 (J6) |  | 6 |
| luhman16-b | 10.3 (W2) | 21.9 (zg) | 21.4 (r) | 22.4 (g) | 19.2 (D4) | 19.6 (J6) |  | 6 |
| luyten-star | 12.1 (W1) | 22.4 (zr) | 20.4 (z) | 22.2 (g) | 19.1 (D4) | 20.0 (J6) |  | 6 |
| procyon-a | 11.8 (W1) | 21.7 (zr) | 20.0 (z) | 22.1 (g) | 19.4 (D4) | 19.4 (J6) |  | 6 |
| procyon-b | 12.0 (W1) | 22.3 (zg) | 20.1 (z) | 22.2 (g) | 19.3 (D4) | 19.5 (J6) |  | 6 |
| proxima-cen | 11.9 (W1) | 23.0 (zg) | 22.0 (g) | 22.0 (r) | 18.4 (D4) |  |  | 5 |
| ross-128 | 12.6 (W1) | 21.0 (zg) | 21.5 (g) | 21.6 (r) | 19.9 (D3) | 20.3 (J6) |  | 6 |
| ross-154 | 12.2 (W1) | 22.6 (zg) | 21.3 (g) | 22.3 (g) | 19.5 (D5) | 19.4 (J6) |  | 6 |
| ross-248 | 12.6 (W1) |  |  |  | 20.3 (D4) | 20.7 (J6) | 23.3 (g) | 4 |
| sigma-dra | 13.1 (W1) |  |  |  | 20.8 (D3) | 21.1 (J6) | 23.7 (g) | 4 |
| sirius-a | 11.2 (W1) | 22.1 (zg) | 21.9 (g) | 22.0 (g) | 19.4 (D4) | 20.0 (J6) |  | 6 |
| sirius-b | 11.7 (W1) | 22.3 (zg) | 21.3 (r) | 21.9 (g) | 19.3 (D4) | 20.1 (J6) |  | 6 |
| struve-2398-a | 12.4 (W1) |  |  |  | 20.9 (D1) | 21.8 (J6) | 24.1 (r) | 4 |
| struve-2398-b | 12.2 (W1) |  |  |  | 20.7 (D1) | 21.6 (J6) | 24.1 (r) | 4 |
| tau-cet | 13.2 (W1) | 22.2 (zr) | 21.5 (i) | 22.2 (r) | 19.9 (D2) | 20.3 (J6) |  | 6 |
| teegarden | 12.5 (W1) | 21.5 (zg) | 21.6 (g) | 22.5 (r) | 20.2 (D3) | 20.6 (J6) |  | 6 |
| van-maanen | 12.7 (W1) | 21.2 (zr) | 22.0 (g) | 21.9 (g) | 20.2 (D3) | 20.9 (J6) |  | 6 |
| wise-0855 | 13.3 (W1) | 22.5 (zr) | 21.6 (r) | 22.7 (g) | 20.5 (D4) | 21.0 (J6) |  | 6 |
| wolf-1061 | 13.2 (W1) | 22.2 (zr) | 21.4 (g) | 22.1 (r) | 19.9 (D4) | 20.8 (J6) |  | 6 |
| wolf-1069 | 12.1 (W1) |  |  |  | 20.6 (D2) | 21.1 (J6) | 23.8 (g) | 4 |
| wolf-359 | 13.6 (W1) | 21.6 (zg) | 21.6 (r) | 21.8 (r) | 20.3 (D3) | 20.5 (J6) |  | 6 |
| wolf-437 | 12.8 (W1) | 21.4 (zg) | 20.5 (z) | 21.8 (g) | 20.2 (D3) | 20.6 (J6) |  | 6 |

### Table 6 — open and unconstrainable cells (crossings channels), by rung

| channel | rung | structurally_open rows | targets | not_constrainable rows | no_survey rows |
|---|---|---|---|---|---|
| A | 0.1 AU | 41 | 7 | 4 | 3 |
| A | 1.0 AU | 3 | 0 | 32 | 5 |
| B | 1.2 Rsun | 24 | 4 | 5 | 3 |
| B | 2.5 Rsun | 25 | 4 | 4 | 3 |
| B | 0.1 AU | 58 | 7 | 5 | 4 |
| S1 | 1.2 Rsun | 0 | 0 | 8 | 0 |
| S1 | 2.5 Rsun | 0 | 0 | 4 | 0 |
| S1 | 0.1 AU | 0 | 0 | 6 | 0 |
| S2 | 0.1 AU | 2 | 2 | 1 | 0 |
| S2 | 0.90 AU | 1 | 1 | 0 | 2 |
| S2 | 0.95 AU | 1 | 1 | 0 | 2 |
| S2 | 1.0 AU | 0 | 0 | 0 | 2 |

### Table 7 — retained-ambiguous and vetoed exceedances

| survey | target | channel | rung | band | state | notes |
|---|---|---|---|---|---|---|
| `atlas_asassn_crossings` | gj-1276 | B | 0.1 AU | o | retained_ambiguous | statistic: S_event; X2: S 5.09 > T 3.03 at z = 10000 AU (all 28 in-window epochs enter); in-window nightly means consistent with zero (−11…+10 µJy ± 5–9) — the  |
| `atlas_asassn_crossings` | gj-1276 | B | 0.1 AU | o | retained_ambiguous | statistic: S_pulse; X2: S 3.04 > T 2.93 — a +31 ± 15 µJy exposure whose quad reads −8, +22, +31, −3 (intra-quad rung 2 fail); same adjudication as the S_event r |
| `atlas_asassn_crossings` | van-maanen | B | 1.2 Rsun | o | retained_ambiguous | statistic: S_event (per-window chord); X1: S 2.78 > T 1.97 (4 valid controls) — three Haleakala exposures at dt +0.27…+0.30 d read 11, 66, 43 ± 30 µJy (quad-con |
| `atlas_asassn_crossings` | van-maanen | B | 2.5 Rsun | o | retained_ambiguous | statistic: S_event; X1 nested: the same 2017-04-03 epoch set as the 1.2 Rsun row, S 2.78 > T 2.64 (7 valid controls) — one epoch set counted in two nested rungs |
| `atlas_asassn_crossings` | wolf-359 | B | 2.5 Rsun | o | retained_ambiguous | dev unit (not blind); statistic: S_event; dev exceedance S 3.55 > T 1.45 on one epoch — a single 3.5σ Sutherland exposure in a degraded-sky quad, quad-inconsist |
| `atlas_asassn_crossings` | wolf-359 | B | 2.5 Rsun | o | retained_ambiguous | dev unit; statistic: S_pulse; the same single-exposure epoch: S 3.18 > T 2.14 (second of the 2 dev exceedances on one epoch), retained, non-promotable; pulse ce |
| `galex_crossings` | wolf-359 | A | 0.1 AU | NUV | retained_ambiguous | S_burst 4.68 > T 4.13 (5 photons in 0.5 s, 0.80 expected); S_rate 1.61 < 3.84, S_period 27.9 < 32.6; veto ladder cannot fire (FRED unresolved at 5 photons; 0 FU |
| `highenergy_crossings` | van-maanen | B | 0.1 AU | 100 MeV–1 GeV (lane L) + 0.1–300 GeV burst | vetoed_known_source | Lane L and burst trials constraint_only at the threshold freeze (amendment v1.2: 3C 279 = 4FGL J1256.1−0547, variability index 33,299, 1.8° from the antipode in |
| `highenergy_crossings` | van-maanen | B | 1.2 Rsun | 100 MeV–1 GeV (lane L) + 0.1–300 GeV burst | vetoed_known_source | Lane L and burst trials constraint_only at the threshold freeze (amendment v1.2: 3C 279 = 4FGL J1256.1−0547, variability index 33,299, 1.8° from the antipode in |
| `highenergy_crossings` | van-maanen | B | 2.5 Rsun | 100 MeV–1 GeV (lane L) + 0.1–300 GeV burst | vetoed_known_source | Lane L and burst trials constraint_only at the threshold freeze (amendment v1.2: 3C 279 = 4FGL J1256.1−0547, variability index 33,299, 1.8° from the antipode in |
| `joint_crossings` | gj-1276 | B | 0.1 AU | i | retained_ambiguous | ledger status retained_ambiguous; 1 ledger row(s), 1 event(s); m90 median 22.0, best 22.0, n with m90 1, 1 grid-censored; S 3.31; T 2.3; n_epochs [11]; PS1 gj-1 |
| `joint_crossings` | ross-128 | B | 0.1 AU | r | vetoed_known_source | ledger status exceedance_vetoed; 1 ledger row(s), 1 event(s); m90 median 21.08, best 21.08, n with m90 1, 0 grid-censored; S 5.05; T 3.61; n_epochs [2]; PS1 ros |
| `lasco_crossings` | van-maanen | S1 | 2.5 Rsun | C2 Orange | retained_ambiguous | S_event 4.51 vs T 1.19: one window, ~3 h all-position-angle annulus disturbance whose onset follows a CDAW-catalogued C2 CME (20:48) by 48 min; both ±25° rings  |
| `ps1_crossings` | gj-1276 | B | 0.1 AU | i | retained_ambiguous | evt-344c12d32a8d: S 3.315 vs T 2.303; night-consistent i ~22.7–23.4 signal (1.1–2.3 sigma per epoch, 3 usable exposures over ~35 min) at the z = 550 node; censu |
| `ps1_crossings` | ross-128 | B | 0.1 AU | r | vetoed_known_source | evt-3f7ad50e70e1: S 5.055 vs T 3.613; vetoed — flux-consistent catalogued static (DR2 stack sources 1.46 arcsec r 22.5–22.8 and 1.99 arcsec r 23.5 from the trac |
| `skirt_crossings` | gj-3512 | S2 | 0.95 AU | o | retained_ambiguous | S_skirt retained_ambiguous_colour_unmatched, S_year retained_ambiguous_colour_unmatched; S_skirt 0.85/−0.12, S_year 2.77/2.14, Δ skirt/opp +0.3 / +0.3 %; colour |
| `skirt_crossings` | teegarden | S2 | 0.90 AU | o | retained_ambiguous | dev unit, o median 12.55; S_sym 2.00/-0.29 exc, S_skirt 1.22/-0.55 exc, S_year 3.98/2.97 exc (S/T); dispositions: S_sym retained_ambiguous_colour_unmatched, S_s |
| `skirt_crossings` | teegarden | S2 | 0.95 AU | o | retained_ambiguous | dev unit, o median 12.55; S_sym 1.96/-0.41 exc, S_skirt 0.47/-0.87 exc, S_year 4.49/3.49 exc (S/T); dispositions: S_sym retained_ambiguous_colour_unmatched, S_s |
| `stereo_hi_crossings` | gj-1276 | S2 | 0.1 AU | HI-1 630–730 nm | retained_ambiguous | S_event 8.22 vs T 7.75: a +0.39-unit plateau over 54 epochs (V_eq 12.1) with no point source in the pixels at any epoch (gj-1276 is V 16); patch's own baseline  |
| `stereo_hi_crossings` | ross-154 | S2 | 0.1 AU | HI-1 630–730 nm | vetoed_known_source | Dev-stage unit (no injection depth): S_pulse exceedance 2010-10-06 single-frame +6.4-unit flare-class event on V1216 Sgr, vetoed by the persistence rule; crowde |
| `tess_crossings` | teegarden | B | 0.1 AU | T | retained_ambiguous | statistic: chord S_c; S_c 79.7 > T_c 52.1 retained-ambiguous, non-promotable (shape test weakened by sector-end truncation — no egress; shares final-day cadence |
| `wispr_crossings` | 61-vir | S2 | 0.1 AU | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | retained_ambiguous | Dev-stage unit: S_event exceedance = monotonic 0.95 → 0.88 calibration drift of the V 4.7 star; retained, non-promotable; dev S_stack 18.70 / 20.05, S_event 31. |
| `wispr_crossings` | gj-1087 | S2 | 0.1 AU | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | vetoed_known_source | S_stack 33.34 / 17.87 exc, S_event 72.13 / 8.76 exc, S_pulse 3.84 / 2.51 exc (S/T); static content (template incompleteness): z > 0 in all 10 events, Gaia DR3 p |
| `wispr_crossings` | gj-581 | S1 | 0.1 AU | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | vetoed_known_source | S_stack 24.99 / 21.17 exc, S_event 20.88 / 51.75, S_pulse 54.28 / 68.50 (S/T); latitude curvature + static content: z +5.0 median; curvature-corrected S_stack 1 |
| `wispr_crossings` | gj-667-c | S1 | 0.1 AU | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | vetoed_known_source | S_stack 44.20 / 10.25 exc, S_event 21.32 / 5.22 exc, S_pulse 1.57 / 3.58 (S/T); static content (bright-neighbour template error): z +14.5 in every event; Gaia 4 |
| `wispr_crossings` | gj-783 | S2 | 0.1 AU | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | vetoed_known_source | S_stack 19.32 / 10.77 exc, S_event 34.00 / 29.73 exc, S_pulse 17.00 / 8.90 exc (S/T); bright-star class (declared v1.3 limit): V 5.3 star; template counted the  |
| `wispr_crossings` | ross-128 | S1 | 0.1 AU | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | vetoed_known_source | S_stack 24.30 / 16.69 exc, S_event 13.35 / 25.75, S_pulse 22.29 / 102.05 (S/T); static content: z +4.9 median, 96 % of events positive; curvature leaves +4.0; G |

## 4. Reading the matrix

- **Best-constrained systems.** gj-1276, teegarden, van-maanen, wolf-359
  and ross-128 carry searched cells in every channel: W-class limits on
  the B grazing and 0.1 AU rungs (optical pre-lens downlink; ZTF, PS1,
  PTF, DASCH, TESS, ATLAS), kW-class on A 0.1 AU (uplink interception at
  opposition), MW-class on the sunward S1 / S2 0.1 AU rungs (LASCO C3,
  STEREO HI-1, WISPR, SoloHI) and on the S1 2.5 R☉ rung (LASCO C2), and
  MW-class on the S2 sub-1 AU skirt (ATLAS). The LAT adds 0.1–0.3 MW
  recurrence-stacked limits on the grazing rungs of its 29 units — the
  only substrate where channel A has grazing rungs.
- **The sunward 0.1 AU rungs are the widest-covered crossings cells** (S1
  19 targets, S2 18 targets searched), because the heliospheric imagers
  see every target's axis on every orbit; their depths span 10 MW to a few
  GW depending on the imager and the star's brightness.
- **The B 0.1 AU rung** is searched on 7 targets with 58 structurally-open
  rows over the same 7 (single-epoch, track-masked, uncovered windows):
  the cell is epoch-starved, not unobservable, and every later optical
  release re-opens it.
- **Corridor depths** (Table 5) are per-endpoint deepest persistent-source
  m90 per survey. For SPHEREx they are the template-absorption-corrected
  values (the v2 frozen depths are overstated by ~0.7 mag for slow
  sources); for PS1 only 1,336 of 7,200 confirmatory records are recovery
  curves and the rest are structurally open (no 90 %-complete regime in
  the stack); for WISE W3/W4 the rows are searched at a calibrated
  threshold with no exclusion claimed.

## 5. Harvest judgement calls (review list)

Each pass recorded where the report's wording did not map one-to-one onto
the state vocabulary; the mapping is stated on the row's `notes`.

- **Record kind `not_constrainable` → `structurally_open`** for Pipeline A
  void nulls (`null_unstable`) and unsupported fits
  (`insufficient_recovery_fit`): the records use the word for "no usable
  statistic", the schema reserves it for physical impossibility.
- **Instrumental / defect vetoes have no state**: TESS teegarden B 2.5 R☉
  (chord-shape veto), DASCH wolf-359 A 0.1 AU (plate defects), ATLAS X3
  teegarden A (proper-motion dipole) stay `searched` with the veto in
  notes; `vetoed_known_source` is used only for catalogued sources.
- **"Retained, non-promotable" → `retained_ambiguous`** (ATLAS X2 gj-1276
  B 0.1 AU, the ATLAS dev wolf-359 B 2.5 R☉ exposure, the WISPR dev 61-vir
  S2 drift, the STEREO gj-1276 S2 plateau). Only ATLAS X1 (van-maanen),
  TESS teegarden 0.1 AU, GALEX wolf-359 NUV, LASCO van-maanen S1 2.5 R☉
  and the skirt colour-unmatched rows are worded "retained-ambiguous" in
  their reports.
- **Dev-stage units** (statistics, no injection depth) are
  `constraint_only`; spectral-archive dev units with measured depths are
  `searched` and flagged "dev unit (not blind)".
- **Coverage-without-statistic wording** ("constraint-only at freeze,
  zero trials", quick-look rms floors, TESS's programme-wide A 1.0 AU) is
  `ledger_only`; the quick-look rms sits in `limit_value` with a unit
  string that says it is a noise floor.
- **WISE Gate I/II units** the report calls constraint-only are
  `not_constrainable` (no pixel was searched); its saturation-excluded
  coverage rows are `ledger_only`.
- **Spectral cells "unconstrained at v1"** (A90 beyond the 200 % cap) are
  `constraint_only` with a null limit. The high-energy corridor screen
  (`catalogue_screen_only` in the report) is `constraint_only` where an
  eROSITA upper limit exists and `ledger_only` for the 52 eastern
  corridors.
- **Aggregate rows spanning rungs** use `rung: any` (36 rows).
- **Absent endpoints** (19 of 88 not in the ZTF / PS1 / joint corridor
  records) are listed as `no_survey` rows derived from the record set
  difference; the reports only say "69-endpoint subset".
- **Report-internal inconsistencies transcribed as stated, not resolved**:
  ZTF corridor 228 vs 230 cells (report vs records) and the joint the
  other way round; SPHEREx dev 323 vs 324; DECam confirmatory 86 vs 89;
  PS1 crossings 33 vs 5 uncovered wide-rung events (§1 vs §3); the
  spectral report's teegarden ESPRESSO 107 W (text) vs 129 W (table); the
  high-energy per-unit uncovered-window arithmetic (65 vs 62).
- **Not representable**: the skirt survey's 0.85 AU rung (no rung
  string); per-unit LAT best-window powers (notes only); BAT at 0.1 AU
  (the report does not speak to it).

## 6. Files

| product | path |
|---|---|
| schema | `surveys/programme-ledger/harvest/SCHEMA.md` |
| harvest, transcribed (21 YAML) | `surveys/programme-ledger/harvest/*.yaml` |
| harvest, generated Pipeline A (7 YAML) | `runs/programme-ledger/v2/harvest_generated/*.yaml` |
| harvest generators (5) | `surveys/programme-ledger/scripts/harvest/*.py` |
| build / render | `surveys/programme-ledger/scripts/{build_ledger_v2,render_tables}.py` |
| full ledger table (28,588 rows) | `runs/programme-ledger/v2/programme_ledger_v2.ecsv` |
| crossings + programme-wide subset | `surveys/programme-ledger/results/programme_ledger_v2_crossings.ecsv` |
| per-target matrix | `surveys/programme-ledger/results/programme_ledger_v2_matrix.json` |
| summary | `surveys/programme-ledger/results/programme_ledger_v2_summary.json` |
| all seven tables | `surveys/programme-ledger/results/programme_ledger_v2_tables.md` |
| records | `runs/programme-ledger/v2/records/{analysis_run,constraint,candidate}.jsonl` |
| manifest | `runs/programme-ledger/v2/manifest.json` |
