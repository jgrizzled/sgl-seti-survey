# confirmatory search v1 — unit table

Construction: threshold_freeze_v1.2; frames usable 27973, unusable {'astrometry_or_zp': 83}; median matched stars 252.0, astrometric rms 0.561 px, ZP0 -19.929, ZP MAD 0.106 mag. Epoch vetoes: {'recon_exclusion': 35, 'near_jupiter': 689, 'near_earth': 356, 'near_mercury': 9448, 'near_venus': 565, 'prox_venus': 3}. Stamp response 0.953 ± 0.067 (n 12247); colour coefficient 0.371 (n 685523).

**Trials 54, exceedances 8, expected control crossings 6.0.**

| unit | status | n_ev | z med / MAD | S_stack S / T | S_event S / T | S_pulse S / T |
|---|---|---|---|---|---|---|
| ez-aqr:S1 | searched | 5 | +0.85 / 2.97 | 2.52 / 4.13 | 4.92 / 12.48 | 6.99 / 14.80 |
| ez-aqr:S2 | constraint_only | 1 | included, insufficient_after_veto, insufficient_epochs, no_baseline | | | |
| gj-1002:S1 | searched | 3 | +2.56 / 2.12 | -0.63 / 6.56 | 3.99 / 15.37 | 2.78 / 10.11 |
| gj-1002:S2 | constraint_only | 2 | included, insufficient_after_veto, insufficient_epochs | | | |
| gj-1111:S1 | constraint_only | 0 | insufficient_epochs | | | |
| gj-1111:S2 | searched | 6 | +8.68 / 9.60 | 19.16 / 13.33 **EXC** | 16.93 / 12.81 **EXC** | 38.43 / 30.16 **EXC** |
| gj-1276:S1 | searched | 3 | -3.34 / 0.16 | -7.59 / 14.21 | -3.23 / 19.59 | 2.17 / 7.22 |
| gj-1276:S2 | constraint_only | 2 | included, insufficient_epochs | | | |
| gj-251:S2 | searched | 4 | -0.74 / 0.56 | -1.02 / 9.79 | 0.84 / 8.56 | 5.48 / 5.42 **EXC** |
| gj-581:S1 | searched | 4 | +2.43 / 0.89 | 3.44 / 13.91 | 3.48 / 12.14 | 1.47 / 8.43 |
| gj-588:S1 | searched | 3 | +2.13 / 3.21 | 2.75 / 5.74 | 4.30 / 6.41 | 2.63 / 10.70 |
| gj-667-c:S1 | searched | 4 | -0.87 / 1.31 | -1.43 / 10.67 | 1.45 / 13.46 | 3.56 / 4.39 |
| gj-674:S1 | searched | 3 | +2.60 / 0.22 | 4.71 / 5.31 | 3.10 / 5.00 | 3.14 / 3.74 |
| gj-682:S1 | searched | 5 | -0.87 / 1.04 | -0.83 / 7.85 | 3.56 / 6.49 | 2.50 / 7.64 |
| gj-783:S1 | searched | 7 | +0.60 / 1.94 | 2.89 / 3.47 | 6.67 / 8.80 | 6.79 / 10.13 |
| gj-876:S1 | searched | 5 | +1.20 / 6.97 | 4.60 / 3.71 **EXC** | 7.87 / 4.85 **EXC** | 2.10 / 17.96 |
| gj-876:S2 | constraint_only | 2 | included, insufficient_after_veto, insufficient_epochs | | | |
| ross-128:S1 | constraint_only | 2 | included, no_baseline | | | |
| ross-128:S2 | searched | 3 | +1.10 / 4.48 | 0.82 / 8.36 | 4.12 / 15.33 | 2.47 / 17.75 |
| ross-154:S2 | constraint_only | 0 | insufficient_epochs | | | |
| teegarden:S1 | constraint_only | 2 | included, insufficient_epochs | | | |
| teegarden:S2 | searched | 8 | -0.35 / 2.45 | 1.78 / 4.37 | 7.80 / 6.77 **EXC** | 15.41 / 6.37 **EXC** |
| van-maanen:S1 | searched | 3 | -1.40 / 6.50 | -1.92 / 7.76 | 3.86 / 12.07 | 4.58 / 5.90 |
| van-maanen:S2 | searched | 3 | +0.59 / 1.08 | 0.75 / 8.89 | 1.32 / 9.20 | 4.26 / 33.89 |
| wolf-1061:S1 | searched | 3 | -2.62 / 0.28 | -2.66 / 6.80 | 0.81 / 12.73 | 1.74 / 8.26 |
| wolf-359:S1 | constraint_only | 1 | included, insufficient_epochs, no_baseline | | | |
| wolf-359:S2 | searched | 3 | -1.07 / 0.86 | -4.34 / 1.12 | -0.48 / 3.74 | 2.44 / 8.22 |

## Completeness (m90, V-equivalent) and power

| unit | S_stack m90 / P | S_event m90 / P | S_pulse m90 / P |
|---|---|---|---|
| ez-aqr:S1 | 12.57 / 4.80e+07 W | 9.66 / 7.06e+08 W | 6.18 / 1.73e+10 W |
| gj-1002:S1 | 10.47 / 3.34e+08 W | 9.32 / 9.59e+08 W | 6.68 / 1.09e+10 W |
| gj-1276:S1 | 9.07 / 1.21e+09 W | 8.44 / 2.16e+09 W | 6.75 / 1.02e+10 W |
| gj-251:S2 | 10.09 / 3.48e+08 W | 9.66 / 5.19e+08 W | 7.19 / 5.04e+09 W |
| gj-581:S1 | 9.65 / 7.07e+08 W | 9.53 / 7.92e+08 W | 6.79 / 9.91e+09 W |
| gj-588:S1 | 11.39 / 1.43e+08 W | 10.85 / 2.35e+08 W | 6.79 / 9.85e+09 W |
| gj-667-c:S1 | 10.46 / 3.38e+08 W | 10.09 / 4.73e+08 W | 7.54 / 4.97e+09 W |
| gj-674:S1 | 13.36 / 2.33e+07 W | 11.40 / 1.42e+08 W | 7.88 / 3.61e+09 W |
| gj-682:S1 | 10.09 / 4.74e+08 W | 10.47 / 3.33e+08 W | 6.78 / 9.97e+09 W |
| gj-783:S1 | 13.89 / 1.43e+07 W | 10.84 / 2.37e+08 W | 6.86 / 9.31e+09 W |
| ross-128:S2 | 10.47 / 8.92e+07 W | 9.66 / 1.89e+08 W | 6.37 / 3.91e+09 W |
| teegarden:S2 | 11.83 / 3.30e+07 W | — | — |
| van-maanen:S1 | 9.68 / 6.93e+08 W | 8.98 / 1.31e+09 W | 7.19 / 6.83e+09 W |
| van-maanen:S2 | 9.62 / 3.20e+08 W | 8.91 / 6.15e+08 W | 4.75 / 2.84e+10 W |
| wolf-1061:S1 | 9.95 / 5.39e+08 W | 9.17 / 1.10e+09 W | 7.08 / 7.59e+09 W |
| wolf-359:S2 | 10.48 / 4.53e+07 W | 10.46 / 4.58e+07 W | 6.58 / 1.63e+09 W |

## Adjudication (frozen §6 ladder; blind statistics stand)

**8 exceedances vs 6.0 expected, 0 candidates.** Evidence:
`adjudicate.py`, the pixel stamps at the flagged epochs (re-fetched
frames), the known-object census (Horizons `@-144`, H < 7
asteroids within 1°), aligned mean stacks of the arcs, and one
post-blind *diagnostic* reduce with the robust median of the controls
in place of the frozen quadratic interpolation
(`confirmatory_search_v1_diag_median.json`; `SOLOHI_INTERP=median`;
not a re-search — the frozen S/T stand).

| unit / statistic | S / T | class | evidence |
|---|---|---|---|
| gj-1111 S2 — S_stack 19.2 / 13.3, S_event 16.9 / 12.8, S_pulse 38.4 / 30.2 | | **control-interpolation artifact** | The +1.5° control (one-sided ladder, source at orbit latitude +4.4…+5.0°) holds a bright star whose core is intermittently removed by the v1.1 saturation mask: its E′ drops to −51, −54, −47 units on single frames (2022-03-31 08:38–09:14 and the like in P07, P11); the frozen quadratic, extrapolated to the source offset, turns that into D_src = +59, +60, +53. The source pixels at those epochs are flat (re-measured flux −3.6 … +0.4 units, 9 × 9 stamps ≤ 0.5 units/px). Diagnostic median: S_stack −3.5 / 5.6, S_event 1.0 / 8.3, S_pulse 2.4 / 5.4. DX Cnc's own baseline 0.09 units (V 14.9). |
| gj-251 S2 — S_pulse 5.48 / 5.42 | | **control-interpolation artifact** (marginal) | Same mechanism (one-sided ladder at +6.5°); diagnostic median 2.38 / 3.12; S_stack, S_event null (−1.0 / 9.8, 0.8 / 8.6). |
| teegarden S2 — S_event 7.80 / 6.77 (P12, 16 epochs), S_pulse 15.4 / 6.4 | | **uncatalogued moving objects through the fixed patch** | P12 2026-02-23 20:02 and 20:26: a compact 2.6–2.9 units/px peak displaced by 4 px (5′) between consecutive frames, absent at 19:38 and 20:50 (flux 12.3 → 14.6 → −2.1); P10 2025-03-24 05:13–05:37 (pulse 11.9σ): a diffuse 0.5–0.8 units/px blob over ~5 px, gone by 06:00. No H < 7 asteroid within 1° at either epoch; Uranus 10°–14° away. The P12 z is carried by those 3 of 16 epochs (halves −0.05 / +11.1); the star's baseline is 0.14–0.19 units (V 15.1). Persists under the median diagnostic (movers are real flux) — rule 2 (a fixed-position source) fails. |
| gj-876 S1 — S_stack 4.60 / 3.71, S_event 7.87 / 4.85 | | **extended level offset, alternating sign** | Antipode patch, starless (Gaia template −0.07, baseline −0.2 … +0.1): per-orbit z +1.2 (P04), +0.7 (P08), **+5.9 (P09), −5.4 (P10), +7.9 (P11)**; split halves consistent within each orbit. Aligned mean stacks (169 / 79 / 86 frames): uniform +0.02–0.04, −0.01–0.03, +0.00–0.02 units/px over the whole 15 × 15 stamp, centre pixel 0.036 / −0.026 / 0.024 against the ≥ 0.1 a PSF of the measured 0.4-unit excess would give — no point source; a residual background gradient (bottom half +0.033 vs top −0.002 in P09). Marginal (1.2–1.6 × T, controls spanning −4.2 … +3.7; median diagnostic 3.9 / 3.7, 5.6 / 5.3); no consistent-sign recurrence (rule 6). |

The LASCO / HI-1 / WISPR cross-substrate family: van-maanen S1
(S_stack −1.9 / 7.8) and S2 (0.8 / 8.9), wolf-359 S2 (−4.3 / 1.1),
teegarden S2 (1.8 / 4.4; the S_event/S_pulse movers above), ross-128
S2 (0.8 / 8.4), gj-1276 S1 (−7.6 / 14.2) — null; wolf-359 S1,
ross-128 S1, teegarden S1, gj-1276 S2, ross-154 S2 fell to
`constraint_only` (< 3 events survive the edge gates); gj-908 S1/S2
and ross-154 S1 null in dev.

**Not searched (9 confirmatory units, `constraint_only`)**: ez-aqr S2,
gj-1002 S2, gj-1111 S1, gj-1276 S2, gj-876 S2, ross-128 S1, ross-154
S2, teegarden S1, wolf-359 S1 — their arcs lie mostly inside the
sunward-edge band that v1.1/v1.2 exclude (saturation plateau in the
long-exposure regime, unresolved corona structure), or lack the
same-tile baseline; 0–2 events survive.

**Post-blind lesson (v2, not applied)**: the frozen quadratic-in-
latitude control interpolation is not robust — a single control
excursion on a one-sided ladder (125 of 178 unit-events) extrapolates
into the source differential; the diagnostic median reduce gives 4
exceedances vs 6.0 expected (the two movers and gj-876 S1), the same
0 candidates, and lower thresholds (S_stack T median 6.8 vs 8.1).
Control patches should also be masked for Gaia G ≤ 8 stars (not 6)
wherever the saturation mask can bite.
