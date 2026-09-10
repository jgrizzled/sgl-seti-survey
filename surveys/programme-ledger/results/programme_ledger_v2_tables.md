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

### Table 4 — programme-wide statements (rows that speak to the whole target list)

| survey | channel | rung | band | state | limit | power | n | source |
|---|---|---|---|---|---|---|---|---|
| `wise_survey` | corridor | 550-10000 AU | W1 | searched | 13.44 Vega mag |  | 784 | report/wise_survey.md §4 table; surveys/wise/results/report_tables.md |
| `wise_survey` | corridor | 550-10000 AU | W2 | searched | 12.28 Vega mag |  | 888 | report/wise_survey.md §4 table; surveys/wise/results/report_tables.md |
| `wise_survey` | corridor | 550-10000 AU | W3 | searched | 9.7 Vega mag |  | 751 | report/wise_survey.md §4 table; surveys/wise/results/report_tables.md |
| `wise_survey` | corridor | 550-10000 AU | W4 | searched | 7.35 Vega mag |  | 880 | report/wise_survey.md §4 table; surveys/wise/results/report_tables.md |
| `wise_survey` | corridor | 550-10000 AU | W3 | structurally_open |  |  |  | report/wise_survey.md §4, §10 (open) |
| `wise_survey` | corridor | 550-10000 AU | W4 | structurally_open |  |  |  | report/wise_survey.md §4, §10 (open) |
| `wise_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 70 | report/wise_survey.md §1, §2; surveys/wise/results/report_tables.md |
| `wise_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 27 | report/wise_survey.md §1, §3; surveys/wise/results/report_tables.md |
| `ztf_survey` | corridor | 550-10000 AU | zg | searched | 23.21 AB mag |  | 674 | report/ztf_survey.md §Completeness table; surveys/ztf/results/report_tables.md |
| `ztf_survey` | corridor | 550-10000 AU | zr | searched | 23.08 AB mag |  | 690 | report/ztf_survey.md §Completeness table; surveys/ztf/results/report_tables.md |
| `ztf_survey` | corridor | 550-10000 AU | zi | searched | 20.96 AB mag |  | 329 | report/ztf_survey.md §Completeness table; surveys/ztf/results/report_tables.md |
| `ztf_survey` | corridor | 550-10000 AU | zg | structurally_open | 20.37–21.32 AB mag |  | 135 | report/ztf_survey.md (layered-search rule); runs/ztf/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ztf_survey` | corridor | 550-10000 AU | zi | structurally_open | 20.07–20.98 AB mag |  | 60 | report/ztf_survey.md (layered-search rule); runs/ztf/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ztf_survey` | corridor | 550-10000 AU | zr | structurally_open | 20.41–21.36 AB mag |  | 137 | report/ztf_survey.md (layered-search rule); runs/ztf/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ztf_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 4 | report/ztf_survey.md §Results; surveys/ztf/results/report_tables.md |
| `ztf_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 1 | report/ztf_survey.md §Results; surveys/ztf/results/report_tables.md |
| `spherex_joint6` | corridor | 550-10000 AU | J6 | searched | 20.02 AB mag |  | 944 | report/spherex_joint6.md §Summary + 'Corrected limits' table; surveys/spherex/results/joint6/report_tables.md |
| `spherex_joint6` | corridor | 550-10000 AU | J6 | searched | 20.69 AB mag |  | 386 | report/spherex_joint6.md §Joint cell — results table; surveys/spherex/results/template_absorption_corrected.md §joint6 |
| `spherex_joint6` | corridor | 550-10000 AU | J6 | structurally_open |  |  | 3 | report/spherex_joint6.md §Joint cell — results table |
| `spherex_joint6` | corridor | 550-10000 AU | J6 | structurally_open |  |  | 4 | report/spherex_joint6.md §Joint cell — results table + 'Dev-driven amendment (a)' |
| `spherex_survey` | corridor | 550-10000 AU | D1 | searched | 19.3 AB mag |  | 952 | report/spherex_survey.md §Completeness table + §Results 'Deferred controls' paragraph; surveys/spherex/results/report_tables.md |
| `spherex_survey` | corridor | 550-10000 AU | D2 | searched | 19.2 AB mag |  | 920 | report/spherex_survey.md §Completeness table + §Results 'Deferred controls' paragraph; surveys/spherex/results/report_tables.md |
| `spherex_survey` | corridor | 550-10000 AU | D3 | searched | 19.43 AB mag |  | 952 | report/spherex_survey.md §Completeness table + §Results 'Deferred controls' paragraph; surveys/spherex/results/report_tables.md |
| `spherex_survey` | corridor | 550-10000 AU | D4 | searched | 19.44 AB mag |  | 960 | report/spherex_survey.md §Completeness table + §Results 'Deferred controls' paragraph; surveys/spherex/results/report_tables.md |
| `spherex_survey` | corridor | 550-10000 AU | D5 | searched | 18.75 AB mag |  | 952 | report/spherex_survey.md §Completeness table + §Results 'Deferred controls' paragraph; surveys/spherex/results/report_tables.md |
| `spherex_survey` | corridor | 550-10000 AU | D6 | searched | 18.24 AB mag |  | 928 | report/spherex_survey.md §Completeness table + §Results 'Deferred controls' paragraph; surveys/spherex/results/report_tables.md |
| `spherex_survey` | corridor | 550-10000 AU | D1 | searched | 19.34 AB mag |  | 416 | surveys/spherex/results/report_tables.md §Development set; surveys/spherex/results/template_absorption_corrected.json families.per-detector.summary dev/D1 |
| `spherex_survey` | corridor | 550-10000 AU | D2 | searched | 19.29 AB mag |  | 408 | surveys/spherex/results/report_tables.md §Development set; surveys/spherex/results/template_absorption_corrected.json families.per-detector.summary dev/D2 |
| `spherex_survey` | corridor | 550-10000 AU | D3 | searched | 19.6 AB mag |  | 384 | surveys/spherex/results/report_tables.md §Development set; surveys/spherex/results/template_absorption_corrected.json families.per-detector.summary dev/D3 |
| `spherex_survey` | corridor | 550-10000 AU | D4 | searched | 19.68 AB mag |  | 392 | surveys/spherex/results/report_tables.md §Development set; surveys/spherex/results/template_absorption_corrected.json families.per-detector.summary dev/D4 |
| `spherex_survey` | corridor | 550-10000 AU | D5 | searched | 18.93 AB mag |  | 384 | surveys/spherex/results/report_tables.md §Development set; surveys/spherex/results/template_absorption_corrected.json families.per-detector.summary dev/D5 |
| `spherex_survey` | corridor | 550-10000 AU | D6 | searched | 18.44 AB mag |  | 392 | surveys/spherex/results/report_tables.md §Development set; surveys/spherex/results/template_absorption_corrected.json families.per-detector.summary dev/D6 |
| `spherex_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 18 | report/spherex_survey.md §Results table; surveys/spherex/results/report_tables.md §Confirmatory set |
| `spherex_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 9 | report/spherex_survey.md §Results table; surveys/spherex/results/report_tables.md §Development set |
| `ps1_survey` | corridor | 550-10000 AU | g | searched | 21.71 AB mag |  | 107 | report/ps1_survey.md §Completeness table; surveys/panstarrs/results/report_tables.md |
| `ps1_survey` | corridor | 550-10000 AU | r | searched | 21.56 AB mag |  | 64 | report/ps1_survey.md §Completeness table; surveys/panstarrs/results/report_tables.md |
| `ps1_survey` | corridor | 550-10000 AU | i | searched | 20.99 AB mag |  | 71 | report/ps1_survey.md §Completeness table; surveys/panstarrs/results/report_tables.md |
| `ps1_survey` | corridor | 550-10000 AU | z | searched | 20.46 AB mag |  | 60 | report/ps1_survey.md §Completeness table; surveys/panstarrs/results/report_tables.md |
| `ps1_survey` | corridor | 550-10000 AU | y | searched | 19.38 AB mag |  | 104 | report/ps1_survey.md §Completeness table; surveys/panstarrs/results/report_tables.md |
| `ps1_survey` | corridor | 550-10000 AU | g | structurally_open | 18.15–22.42 AB mag |  | 138 | report/ps1_survey.md (layered-search rule); runs/panstarrs/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ps1_survey` | corridor | 550-10000 AU | i | structurally_open | 20.62–21.91 AB mag |  | 138 | report/ps1_survey.md (layered-search rule); runs/panstarrs/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ps1_survey` | corridor | 550-10000 AU | r | structurally_open | 21.04–22.17 AB mag |  | 138 | report/ps1_survey.md (layered-search rule); runs/panstarrs/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ps1_survey` | corridor | 550-10000 AU | y | structurally_open | 19.1–20.14 AB mag |  | 138 | report/ps1_survey.md (layered-search rule); runs/panstarrs/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ps1_survey` | corridor | 550-10000 AU | z | structurally_open | 20.39–21.18 AB mag |  | 138 | report/ps1_survey.md (layered-search rule); runs/panstarrs/v2/records/constraint.jsonl extra.bright_limit_mag |
| `ps1_survey` | corridor | 550-10000 AU | any | structurally_open |  |  |  | report/ps1_survey.md §Completeness (structural result); surveys/panstarrs/results/report_tables.md §Ledger |
| `ps1_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 8 | report/ps1_survey.md §Results; surveys/panstarrs/results/report_tables.md |
| `ps1_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 4 | report/ps1_survey.md §Results; surveys/panstarrs/results/report_tables.md |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | g | searched | 22.83 AB mag |  | 607 | report/joint_ps1_ztf_wise.md §Completeness table; surveys/joint/results/report_tables.md |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | r | searched | 22.86 AB mag |  | 639 | report/joint_ps1_ztf_wise.md §Completeness table; surveys/joint/results/report_tables.md |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | i | searched | 21.09 AB mag |  | 78 | report/joint_ps1_ztf_wise.md §Completeness table; surveys/joint/results/report_tables.md |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | W1/W2 (colour axis) | ledger_only |  |  | 230 | report/joint_ps1_ztf_wise.md §Summary, §Colour veto calibration, §Results |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | W3 | no_survey |  |  |  | report/joint_ps1_ztf_wise.md §Scope and provenance notes |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | W4 | no_survey |  |  |  | report/joint_ps1_ztf_wise.md §Scope and provenance notes |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | any | structurally_open |  |  | 5 | report/joint_ps1_ztf_wise.md §Results; surveys/joint/results/report_tables.md |
| `joint_ps1_ztf_wise` | corridor | 550-10000 AU | any | structurally_open |  |  | 1 | report/joint_ps1_ztf_wise.md §Results; surveys/joint/results/report_tables.md |
| `decam_survey` | corridor | 550-10000 AU | g | searched | 23.18 AB mag |  | 78 | report/decam_survey.md §Completeness table; surveys/decam/results/report_tables.md §Confirmatory set |
| `decam_survey` | corridor | 550-10000 AU | r | searched | 22.13 AB mag |  | 39 | report/decam_survey.md §Completeness table; surveys/decam/results/report_tables.md §Confirmatory set |
| `decam_survey` | corridor | 550-10000 AU | i | searched | 22.73 AB mag |  | 64 | report/decam_survey.md §Completeness table; surveys/decam/results/report_tables.md §Confirmatory set |
| `decam_survey` | corridor | 550-10000 AU | z | searched | 22.08 AB mag |  | 44 | report/decam_survey.md §Completeness table; surveys/decam/results/report_tables.md §Confirmatory set |
| `decam_survey` | corridor | 550-10000 AU | Y | searched | 21 AB mag |  | 27 | report/decam_survey.md §Completeness table; surveys/decam/results/report_tables.md §Confirmatory set |
| `decam_survey` | corridor | 550-10000 AU | g | searched | 23.54 AB mag |  | 57 | surveys/decam/results/report_tables.md §Development set |
| `decam_survey` | corridor | 550-10000 AU | r | searched | 23.42 AB mag |  | 48 | surveys/decam/results/report_tables.md §Development set |
| `decam_survey` | corridor | 550-10000 AU | i | searched | 22.89 AB mag |  | 41 | surveys/decam/results/report_tables.md §Development set |
| `decam_survey` | corridor | 550-10000 AU | z | searched | 22.56 AB mag |  | 5 | surveys/decam/results/report_tables.md §Development set |
| `decam_survey` | corridor | 550-10000 AU | any | structurally_open |  |  | 378 | report/decam_survey.md §Summary + §Coverage and data notes |
| `plan_no_survey` | A | 0.1 AU | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | A | 0.1 AU | J | structurally_open |  |  | 14 | notes/project_plan.md §5.20; surveys/gattini-crossings/notes/gattini_winter_recon_2026-09-04.md |
| `plan_no_survey` | A | 0.1 AU | Y/J/Hs (1064 nm needs WINTER Y) | no_survey |  |  |  | notes/project_plan.md §5.20 |
| `plan_no_survey` | A | 0.1 AU | G/BP/RP epoch photometry | no_survey |  |  |  | notes/project_plan.md §5.22 |
| `plan_no_survey` | A | 0.1 AU | 0.75–5 µm | no_survey |  |  | 1130 | notes/project_plan.md §5.24 |
| `plan_no_survey` | A | 1.0 AU | G/BP/RP epoch photometry | no_survey |  |  |  | notes/project_plan.md §5.22 |
| `plan_no_survey` | A | 1.0 AU | 0.75–5 µm | no_survey |  |  | 1130 | notes/project_plan.md §5.24 |
| `plan_no_survey` | B | 0.1 AU | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | B | 0.1 AU | J | structurally_open |  |  | 14 | notes/project_plan.md §5.20; surveys/gattini-crossings/notes/gattini_winter_recon_2026-09-04.md |
| `plan_no_survey` | B | 0.1 AU | Y/J/Hs (1064 nm needs WINTER Y) | no_survey |  |  |  | notes/project_plan.md §5.20 |
| `plan_no_survey` | B | 0.1 AU | G/BP/RP epoch photometry | no_survey |  |  |  | notes/project_plan.md §5.22 |
| `plan_no_survey` | B | 0.1 AU | 0.75–5 µm | no_survey |  |  | 1130 | notes/project_plan.md §5.24 |
| `plan_no_survey` | B | 1.0 AU | G/BP/RP epoch photometry | no_survey |  |  |  | notes/project_plan.md §5.22 |
| `plan_no_survey` | B | 1.0 AU | 0.75–5 µm | no_survey |  |  | 1130 | notes/project_plan.md §5.24 |
| `plan_no_survey` | B | 1.2 Rsun | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | B | 1.2 Rsun | J | structurally_open |  |  | 14 | notes/project_plan.md §5.20; surveys/gattini-crossings/notes/gattini_winter_recon_2026-09-04.md |
| `plan_no_survey` | B | 1.2 Rsun | Y/J/Hs (1064 nm needs WINTER Y) | no_survey |  |  |  | notes/project_plan.md §5.20 |
| `plan_no_survey` | B | 1.2 Rsun | G/BP/RP epoch photometry | no_survey |  |  |  | notes/project_plan.md §5.22 |
| `plan_no_survey` | B | 1.2 Rsun | 0.75–5 µm | no_survey |  |  | 1130 | notes/project_plan.md §5.24 |
| `plan_no_survey` | B | 2.5 Rsun | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | B | 2.5 Rsun | J | structurally_open |  |  | 14 | notes/project_plan.md §5.20; surveys/gattini-crossings/notes/gattini_winter_recon_2026-09-04.md |
| `plan_no_survey` | B | 2.5 Rsun | Y/J/Hs (1064 nm needs WINTER Y) | no_survey |  |  |  | notes/project_plan.md §5.20 |
| `plan_no_survey` | B | 2.5 Rsun | G/BP/RP epoch photometry | no_survey |  |  |  | notes/project_plan.md §5.22 |
| `plan_no_survey` | B | 2.5 Rsun | 0.75–5 µm | no_survey |  |  | 1130 | notes/project_plan.md §5.24 |
| `plan_no_survey` | S1 | 0.1 AU | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | S1 | 1.2 Rsun | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | S1 | 2.5 Rsun | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | S2 | 0.1 AU | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | S2 | 1.2 Rsun | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | S2 | 2.5 Rsun | Kepler | not_constrainable |  |  | 1435 | notes/project_plan.md §5.15 |
| `plan_no_survey` | corridor | 550-10000 AU | any | no_survey |  |  |  | notes/project_plan.md §4.7 |
| `plan_no_survey` | corridor | 550-10000 AU | any | no_survey |  |  |  | notes/project_plan.md §4.8 |
| `plan_no_survey` | corridor | 550-10000 AU | G/BP/RP astrometry | no_survey |  |  |  | notes/project_plan.md §4.9 |
| `plan_no_survey` | corridor | 550-10000 AU | thermal IR / sub-mm | no_survey |  |  |  | notes/project_plan.md §4.10 |
| `plan_no_survey` | corridor | 550-10000 AU | JHK / photographic | no_survey |  |  |  | notes/project_plan.md §4.11 |
| `plan_no_survey` | corridor | 550-10000 AU | any | no_survey |  |  |  | notes/project_plan.md §4.12 |
| `ztf_crossings` | A | 0.1 AU | any | structurally_open |  |  | 15 | surveys/ztf-crossings/results/coverage_v1_summary.json; report/ztf_crossings.md §1 item 3 |
| `ztf_crossings` | A | 1.0 AU | any | not_constrainable |  |  | 318 | surveys/ztf-crossings/results/coverage_v1_summary.json; report/ztf_crossings.md §1 item 3 |
| `ztf_crossings` | A | 1.0 AU | g/r | constraint_only |  |  |  | report/ztf_crossings.md §1 item 6; §4 lesson 1 |
| `ztf_crossings` | A | any | g/r | constraint_only | 14.4–18.4 AB mag |  | 18 | report/ztf_crossings.md §2 depth table |
| `ztf_crossings` | A | any | g/r/i | structurally_open |  |  |  | report/ztf_crossings.md §1 item 4 |
| `ztf_crossings` | A | any | 1064/1550 nm | not_constrainable |  |  |  | report/ztf_crossings.md §3 |
| `ztf_crossings` | B | 0.1 AU | g/r | searched | 21.8 AB mag | 130 W | 31 | report/ztf_crossings.md §2 depth table; §3 |
| `ztf_crossings` | B | 0.1 AU | g/r | structurally_open | 21.2 AB mag |  | 8 | report/ztf_crossings.md §2 depth table |
| `ztf_crossings` | B | 0.1 AU | any | structurally_open |  |  | 27 | surveys/ztf-crossings/results/coverage_v1_summary.json; report/ztf_crossings.md §1 item 3 |
| `ztf_crossings` | B | 1.2 Rsun | any | structurally_open |  |  | 34 | surveys/ztf-crossings/results/coverage_v1_summary.json; report/ztf_crossings.md §1 item 3 |
| `ztf_crossings` | B | 2.5 Rsun | any | structurally_open |  |  | 38 | surveys/ztf-crossings/results/coverage_v1_summary.json; report/ztf_crossings.md §1 item 3 |
| `ztf_crossings` | B | any | 1064/1550 nm | not_constrainable |  |  |  | report/ztf_crossings.md §3 |
| `ztf_crossings` | B | any | any | structurally_open |  |  |  | report/ztf_crossings.md §3 |
| `ztf_crossings` | S1 | any | any | no_survey |  |  |  | report/ztf_crossings.md §1 item 2; §3 |
| `ztf_crossings` | S2 | any | any | no_survey |  |  |  | report/ztf_crossings.md §1 item 2; §3 |
| `ptf_crossings` | A | 0.1 AU | g/R | structurally_open |  |  | 38 | surveys/ptf-crossings/results/coverage_v1_summary.json; report/ptf_crossings.md §1 item 3 |
| `ptf_crossings` | A | 1.0 AU | g/R | ledger_only |  |  | 76 | surveys/ptf-crossings/results/coverage_v1_summary.json; report/ptf_crossings.md §1 item 3 |
| `ptf_crossings` | A | 1.0 AU | g/R | structurally_open |  |  | 437 | surveys/ptf-crossings/results/coverage_v1_summary.json; report/ptf_crossings.md §1 item 3 |
| `ptf_crossings` | A | any | 1064/1550 nm | not_constrainable |  |  |  | report/ptf_crossings.md §1 item 2; §3 |
| `ptf_crossings` | B | 0.1 AU | g/R | structurally_open |  |  | 35 | surveys/ptf-crossings/results/coverage_v1_summary.json; report/ptf_crossings.md §1 item 3 |
| `ptf_crossings` | B | 1.2 Rsun | g/R | structurally_open |  |  | 24 | surveys/ptf-crossings/results/coverage_v1_summary.json; report/ptf_crossings.md §1 item 3 |
| `ptf_crossings` | B | 2.5 Rsun | g/R | structurally_open |  |  | 27 | surveys/ptf-crossings/results/coverage_v1_summary.json; report/ptf_crossings.md §1 item 3 |
| `ptf_crossings` | B | any | g/R | searched | 20.58–22 AB mag | 110 W–30 kW |  | report/ptf_crossings.md §2 |
| `ptf_crossings` | B | any | 1064/1550 nm | not_constrainable |  |  |  | report/ptf_crossings.md §1 item 2; §3 |
| `ptf_crossings` | B | any | g/R | structurally_open |  |  |  | report/ptf_crossings.md §3 |
| `ptf_crossings` | S1 | any | any | no_survey |  |  |  | surveys/ptf-crossings/hypotheses.md §1 |
| `ptf_crossings` | S2 | any | any | no_survey |  |  |  | surveys/ptf-crossings/hypotheses.md §1 |
| `lasco_crossings` | S2 | 1.0 AU | any | no_survey |  |  |  | report/lasco_crossings.md §1 |
| `galex_crossings` | B | 1.2 Rsun | any | structurally_open |  |  | 39 | report/galex_crossings.md §1, §3 caveats; surveys/galex-crossings/results/coverage_v1.md; surveys/galex-crossings/results/coverage_v1_summary.json |
| `galex_crossings` | B | 2.5 Rsun | any | structurally_open |  |  | 49 | report/galex_crossings.md §1, §3 caveats; surveys/galex-crossings/results/coverage_v1.md; surveys/galex-crossings/results/coverage_v1_summary.json |
| `rubin_crossings` | A | 1.0 AU | any | no_survey |  |  | 85 | report/rubin_crossings.md §3; surveys/rubin-crossings/results/coverage_v1_summary.json |
| `dasch_crossings` | A | 1.0 AU | any | no_survey |  |  |  | report/dasch_crossings.md §1.3 (ladder: + A 1.0 AU constraint-only), §5 |
| `stereo_hi_crossings` | S1 | 1.2 Rsun | HI-1 630–730 nm | not_constrainable |  |  |  | report/stereo_hi_crossings.md §1, §5 |
| `stereo_hi_crossings` | S1 | 2.5 Rsun | HI-1 630–730 nm | not_constrainable |  |  |  | report/stereo_hi_crossings.md §1, §5 |
| `stereo_hi_crossings` | S2 | 1.0 AU | HI-1 630–730 nm | ledger_only |  |  | 289 | report/stereo_hi_crossings.md §1 |
| `wispr_crossings` | S1 | 0.1 AU | WISPR-O (ε 50°–108°) | ledger_only |  |  |  | report/wispr_crossings.md §6 item 4 |
| `wispr_crossings` | S1 | 1.2 Rsun | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | not_constrainable |  |  |  | report/wispr_crossings.md §1 |
| `wispr_crossings` | S1 | 2.5 Rsun | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | not_constrainable |  |  |  | report/wispr_crossings.md §1 |
| `wispr_crossings` | S2 | 0.1 AU | WISPR-O (ε 50°–108°) | ledger_only |  |  |  | report/wispr_crossings.md §6 item 4 |
| `wispr_crossings` | S2 | 1.0 AU | WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band) | no_survey |  |  |  | report/wispr_crossings.md §1 |
| `ps1_crossings` | A | 0.1 AU | any | structurally_open |  |  | 28 | surveys/ps1-crossings/results/coverage_v1_summary.json; report/ps1_crossings.md §1 item 3 |
| `ps1_crossings` | A | 0.1 AU | grizy | structurally_open |  |  |  | report/ps1_crossings.md §1 item 4 |
| `ps1_crossings` | A | 1.0 AU | any | ledger_only |  |  | 294 | surveys/ps1-crossings/results/coverage_v1_summary.json; report/ps1_crossings.md §1 item 3 |
| `ps1_crossings` | A | 1.0 AU | any | structurally_open |  |  | 269 | surveys/ps1-crossings/results/coverage_v1_summary.json; report/ps1_crossings.md §1 item 3 |
| `ps1_crossings` | A | any | 1064/1550 nm | not_constrainable |  |  |  | report/ps1_crossings.md §1 item 2; §3 |
| `ps1_crossings` | B | 0.1 AU | g/r/i | searched | 21.6 AB mag | 170 W | 11 | report/ps1_crossings.md §2 depth table; §3 |
| `ps1_crossings` | B | 0.1 AU | g | searched | 21.1 AB mag | 280 W | 4 | report/ps1_crossings.md §2 depth table; §3 |
| `ps1_crossings` | B | 0.1 AU | r | searched | 21.1 AB mag |  | 3 | report/ps1_crossings.md §2 depth table |
| `ps1_crossings` | B | 0.1 AU | i | searched | 22 AB mag | 37 W | 4 | report/ps1_crossings.md §2 depth table; §3 |
| `ps1_crossings` | B | 0.1 AU | any | structurally_open |  |  | 33 | surveys/ps1-crossings/results/coverage_v1_summary.json; report/ps1_crossings.md §1 item 3 |
| `ps1_crossings` | B | 1.2 Rsun | any | structurally_open |  |  | 26 | surveys/ps1-crossings/results/coverage_v1_summary.json; report/ps1_crossings.md §1 item 3 |
| `ps1_crossings` | B | 2.5 Rsun | any | structurally_open |  |  | 33 | surveys/ps1-crossings/results/coverage_v1_summary.json; report/ps1_crossings.md §1 item 3 |
| `ps1_crossings` | B | any | any | structurally_open |  |  | 7 | report/ps1_crossings.md §2 |
| `ps1_crossings` | B | any | 1064/1550 nm | not_constrainable |  |  |  | report/ps1_crossings.md §1 item 2; §3 |
| `ps1_crossings` | B | any | any | structurally_open |  |  |  | report/ps1_crossings.md §3 |
| `ps1_crossings` | S1 | any | any | no_survey |  |  |  | surveys/ps1-crossings/hypotheses.md §1 |
| `ps1_crossings` | S2 | any | any | no_survey |  |  |  | surveys/ps1-crossings/hypotheses.md §1 |
| `solohi_crossings` | S1 | 0.1 AU | SoloHI λ_c 610 nm, W_eff 250 nm (532 nm in band) | ledger_only |  |  | 49388 | report/solohi_crossings.md §2, §5, §6 item 2 |
| `solohi_crossings` | S1 | 1.2 Rsun | SoloHI λ_c 610 nm, W_eff 250 nm (532 nm in band) | not_constrainable |  |  |  | report/solohi_crossings.md §1 |
| `solohi_crossings` | S1 | 2.5 Rsun | SoloHI λ_c 610 nm, W_eff 250 nm (532 nm in band) | not_constrainable |  |  |  | report/solohi_crossings.md §1 |
| `solohi_crossings` | S2 | 0.1 AU | SoloHI λ_c 610 nm, W_eff 250 nm (532 nm in band) | ledger_only |  |  | 49388 | report/solohi_crossings.md §2, §5, §6 item 2 |
| `solohi_crossings` | S2 | 1.0 AU | SoloHI λ_c 610 nm, W_eff 250 nm (532 nm in band) | ledger_only |  |  | 214 | report/solohi_crossings.md §1, §6 item 5 |
| `highenergy_crossings` | A | 0.1 AU | 0.1–300 GeV | structurally_open |  |  | 1 | surveys/highenergy-crossings/results/coverage_v1.md (coverage-without-statistic list, counted per channel and rung); report/highenergy_crossings.md §2 coverage |
| `highenergy_crossings` | A | 0.1 AU | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | A | 0.1 AU | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | A | 0.1 AU | any | structurally_open |  |  |  | report/highenergy_crossings.md §2 (geometry) |
| `highenergy_crossings` | A | 1.2 Rsun | 0.1–300 GeV | structurally_open |  |  | 14 | surveys/highenergy-crossings/results/coverage_v1.md (coverage-without-statistic list, counted per channel and rung); report/highenergy_crossings.md §2 coverage |
| `highenergy_crossings` | A | 1.2 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | A | 1.2 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | A | 1.2 Rsun | any | structurally_open |  |  |  | report/highenergy_crossings.md §2 (geometry) |
| `highenergy_crossings` | A | 2.5 Rsun | 0.1–300 GeV | structurally_open |  |  | 13 | surveys/highenergy-crossings/results/coverage_v1.md (coverage-without-statistic list, counted per channel and rung); report/highenergy_crossings.md §2 coverage |
| `highenergy_crossings` | A | 2.5 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | A | 2.5 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | A | 2.5 Rsun | any | structurally_open |  |  |  | report/highenergy_crossings.md §2 (geometry) |
| `highenergy_crossings` | B | 0.1 AU | 0.1–300 GeV | structurally_open |  |  | 7 | surveys/highenergy-crossings/results/coverage_v1.md (coverage-without-statistic list, counted per channel and rung); report/highenergy_crossings.md §2 coverage |
| `highenergy_crossings` | B | 0.1 AU | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | B | 0.1 AU | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | B | 0.1 AU | any | structurally_open |  |  |  | report/highenergy_crossings.md §2 (geometry) |
| `highenergy_crossings` | B | 1.2 Rsun | 0.1–300 GeV | structurally_open |  |  | 13 | surveys/highenergy-crossings/results/coverage_v1.md (coverage-without-statistic list, counted per channel and rung); report/highenergy_crossings.md §2 coverage |
| `highenergy_crossings` | B | 1.2 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | B | 1.2 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | B | 1.2 Rsun | any | structurally_open |  |  |  | report/highenergy_crossings.md §2 (geometry) |
| `highenergy_crossings` | B | 2.5 Rsun | 0.1–300 GeV | structurally_open |  |  | 14 | surveys/highenergy-crossings/results/coverage_v1.md (coverage-without-statistic list, counted per channel and rung); report/highenergy_crossings.md §2 coverage |
| `highenergy_crossings` | B | 2.5 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | B | 2.5 Rsun | any | not_constrainable |  |  |  | report/highenergy_crossings.md §2 (geometry); §4 |
| `highenergy_crossings` | B | 2.5 Rsun | any | structurally_open |  |  |  | report/highenergy_crossings.md §2 (geometry) |
| `wise_crossings` | A | 0.1 AU | W1/W2/W3/W4 | not_constrainable |  |  | 101 | surveys/wise-crossings/configs/elongation_gate_v1.json; report/wise_crossings.md §2 table |
| `wise_crossings` | A | 1.0 AU | W1/W2 | not_constrainable |  |  | 1272 | surveys/wise-crossings/configs/elongation_gate_v1.json; report/wise_crossings.md §2 table; §3 |
| `wise_crossings` | A | 1.0 AU | W1/W2 | ledger_only |  |  | 91 | report/wise_crossings.md §4 |
| `wise_crossings` | A | 1.0 AU | W1/W2 | structurally_open |  |  |  | report/wise_crossings.md §1 item 5; §5 lesson 2 |
| `wise_crossings` | A | any | 3–5 um | not_constrainable |  |  |  | report/wise_crossings.md §4 |
| `wise_crossings` | B | 0.1 AU | W1/W2/W3/W4 | not_constrainable |  |  | 102 | surveys/wise-crossings/configs/elongation_gate_v1.json; report/wise_crossings.md §2 table |
| `wise_crossings` | B | 1.2 Rsun | W1/W2/W3/W4 | not_constrainable |  |  | 50 | surveys/wise-crossings/configs/elongation_gate_v1.json; report/wise_crossings.md §2 table |
| `wise_crossings` | B | 2.5 Rsun | W1/W2/W3/W4 | not_constrainable |  |  | 67 | surveys/wise-crossings/configs/elongation_gate_v1.json; report/wise_crossings.md §2 table |
| `wise_crossings` | B | any | 3–5 um | not_constrainable |  |  |  | report/wise_crossings.md §4 |
| `joint_crossings` | A | 1.0 AU | any | constraint_only |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | A | any | any | searched |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | A | any | 1064/1550 nm and line SEDs redward of ~900 nm | not_constrainable |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | A | any | mid-IR | not_constrainable |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | A | any | any | structurally_open |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | B | 0.1 AU | g/r | searched | 21.8 AB mag | 130 W |  | report/joint_crossings.md §3 |
| `joint_crossings` | B | 0.1 AU | g/r/i | searched | 21.6 AB mag | 170 W–280 W |  | report/joint_crossings.md §3 |
| `joint_crossings` | B | any | any | searched |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | B | any | any | structurally_open |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | B | any | 1064/1550 nm and line SEDs redward of ~900 nm | not_constrainable |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | B | any | mid-IR | not_constrainable |  |  |  | report/joint_crossings.md §3 |
| `joint_crossings` | B | any | any | structurally_open |  |  |  | report/joint_crossings.md §3 |
| `radio_crossings` | A | 0.1 AU | any | structurally_open |  |  |  | report/radio_crossings.md §2 table; §3 item 1 |
| `radio_crossings` | A | 1.0 AU | any | ledger_only |  |  | 22 | report/radio_crossings.md §2 table; §3 item 2; surveys/radio-crossings/results/coverage_v1_summary.json |
| `radio_crossings` | A | 1.0 AU | any | ledger_only |  |  | 69 | surveys/radio-crossings/results/coverage_v1_summary.json (vlass A_1.0) |
| `radio_crossings` | B | 0.1 AU | any | structurally_open |  |  |  | report/radio_crossings.md §2 table; §3 item 3; §5 |
| `radio_crossings` | B | 1.2 Rsun | any | structurally_open |  |  |  | report/radio_crossings.md §2 table; §3 item 3; §5 |
| `radio_crossings` | B | 1.2 Rsun | any | structurally_open |  |  |  | report/radio_crossings.md §2 table |
| `radio_crossings` | B | 2.5 Rsun | any | structurally_open |  |  |  | report/radio_crossings.md §2 table; §3 item 3; §5 |
| `radio_crossings` | B | 2.5 Rsun | any | structurally_open |  |  |  | report/radio_crossings.md §2 table |
| `tess_crossings` | A | 1.0 AU | any | ledger_only |  |  | 292 | report/tess_crossings.md §4, §6; surveys/tess-crossings/results/coverage_gate_v1_summary.json |
| `radio_crossings_ext` | A | 0.1 AU | any | ledger_only |  |  | 50 | report/radio_crossings_ext.md §2 table; §2.2 |
| `radio_crossings_ext` | A | 0.1 AU | 144 MHz | ledger_only |  |  | 50 | report/radio_crossings_ext.md §2 table |
| `radio_crossings_ext` | A | 1.0 AU | any | ledger_only |  |  | 606 | report/radio_crossings_ext.md §2 table; §2.3 |
| `radio_crossings_ext` | A | 1.0 AU | 144 MHz | ledger_only |  |  | 300 | report/radio_crossings_ext.md §2 table; §2.3 |
| `radio_crossings_ext` | B | 0.1 AU | any | ledger_only |  |  | 50 | report/radio_crossings_ext.md §2 table; §2.1; §3 item 2 |
| `radio_crossings_ext` | B | 0.1 AU | 144 MHz | structurally_open |  |  | 30 | report/radio_crossings_ext.md §2 table; §2.1 |
| `radio_crossings_ext` | B | 1.2 Rsun | any | structurally_open |  |  | 29 | report/radio_crossings_ext.md §2 table; §3 item 3 |
| `radio_crossings_ext` | B | 1.2 Rsun | 144 MHz | structurally_open |  |  | 10 | report/radio_crossings_ext.md §2 table |
| `radio_crossings_ext` | B | 2.5 Rsun | any | structurally_open |  |  | 36 | report/radio_crossings_ext.md §2 table; §3 item 3 |
| `radio_crossings_ext` | B | 2.5 Rsun | 144 MHz | structurally_open |  |  | 20 | report/radio_crossings_ext.md §2 table |
| `atlas_asassn_crossings` | A | 0.1 AU | o | structurally_open |  |  | 9 | report/atlas_asassn_crossings.md §1.4; surveys/atlas-asassn-crossings/results/coverage_v1_summary.json |
| `atlas_asassn_crossings` | A | 1.0 AU | any | no_survey |  |  | 957 | report/atlas_asassn_crossings.md §2 (constraint-only ledger); surveys/atlas-asassn-crossings/results/era_scope_v1.json |
| `atlas_asassn_crossings` | B | 0.1 AU | o | structurally_open |  |  | 3 | report/atlas_asassn_crossings.md §1.4; surveys/atlas-asassn-crossings/results/coverage_v1_summary.json |
| `atlas_asassn_crossings` | B | 1.2 Rsun | o | structurally_open |  |  | 36 | report/atlas_asassn_crossings.md §1.4, §3 (declared unconstrained), §4; surveys/atlas-asassn-crossings/results/coverage_v1_summary.json |
| `atlas_asassn_crossings` | B | 2.5 Rsun | o | structurally_open |  |  | 30 | report/atlas_asassn_crossings.md §1.4, §3; surveys/atlas-asassn-crossings/results/coverage_v1_summary.json |

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
