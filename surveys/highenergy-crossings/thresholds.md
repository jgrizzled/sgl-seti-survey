# Threshold freeze v1.0 — high-energy (Fermi-LAT) crossings (2026-09-07)

Analytic Poisson thresholds (hypotheses D6): T = 3.451 in −log10 p for 145 confirmatory trials (29 units × 5 statistics; amendment v1.1) at FWER 0.05 (Šidák per-trial α = 3.54e-04). Seed 20260907. sha256 of `results/thresholds_v1.json` in `results/threshold_freeze_v1_sha.txt`.

Validation on the pseudo-window ensembles (off-window photons only): dispersion gate (≤ 10 % of members with p < 0.05) and the per-trial exceedance gate (amendment v1.2: ≤ max(3, 3 × expected) members above T).

| Unit | dev | n_searchable | L: n_pw / disp / exc T / stack-ens max | H: n_pw / disp / exc T / stack-ens max | constraint-only trials |
|---|---|---|---|---|---|
| A gj-1276 1.2Rsun | dev | 15 | 1246 / 0.030 / 0 / 3.0 | 1246 / 0.012 / 0 / 1.0 | — |
| A gj-1276 2.5Rsun |  | 16 | 1113 / 0.034 / 0 / 2.6 | 1113 / 0.017 / 1 / 2.7 | — |
| A gj-1276 0.1AU |  | 18 | 291 / 0.055 / 0 / 0.9 | 291 / 0.024 / 0 / 0.8 | — |
| A gj-908 0.1AU | dev | 18 | 365 / 0.052 / 0 / 2.4 | 365 / 0.027 / 1 / 0.8 | — |
| A ross-128 2.5Rsun |  | 17 | 1685 / 0.037 / 0 / 2.1 | 1685 / 0.027 / 0 / 2.2 | — |
| A ross-128 0.1AU |  | 18 | 307 / 0.065 / 0 / 1.3 | 307 / 0.039 / 0 / 1.0 | — |
| A ross-154 0.1AU |  | 18 | 279 / 0.072 / 2 / 2.3 | 279 / 0.050 / 0 / 1.3 | — |
| A teegarden 1.2Rsun |  | 14 | 1318 / 0.036 / 0 / 1.5 | 1318 / 0.021 / 0 / 1.7 | — |
| A teegarden 2.5Rsun |  | 15 | 1269 / 0.047 / 0 / 1.7 | 1269 / 0.009 / 0 / 1.7 | — |
| A teegarden 0.1AU |  | 17 | 295 / 0.088 / 0 / 0.9 | 295 / 0.024 / 0 / 1.2 | — |
| A van-maanen 1.2Rsun |  | 15 | 1546 / 0.027 / 0 / 1.3 | 1546 / 0.016 / 0 / 1.7 | — |
| A van-maanen 2.5Rsun |  | 15 | 1131 / 0.041 / 1 / 2.4 | 1131 / 0.026 / 0 / 2.6 | — |
| A van-maanen 0.1AU |  | 18 | 310 / 0.042 / 1 / 1.7 | 310 / 0.023 / 0 / 1.1 | — |
| A wolf-359 1.2Rsun |  | 15 | 1506 / 0.024 / 0 / 2.7 | 1506 / 0.024 / 0 / 2.4 | — |
| A wolf-359 2.5Rsun |  | 15 | 1192 / 0.031 / 0 / 1.7 | 1192 / 0.017 / 0 / 1.9 | — |
| A wolf-359 0.1AU |  | 18 | 318 / 0.053 / 1 / 2.0 | 318 / 0.028 / 0 / 1.0 | — |
| B gj-1276 1.2Rsun |  | 15 | 1522 / 0.027 / 0 / 2.3 | 1522 / 0.022 / 0 / 2.0 | — |
| B gj-1276 2.5Rsun |  | 15 | 1166 / 0.036 / 1 / 1.8 | 1166 / 0.022 / 0 / 1.3 | — |
| B gj-1276 0.1AU |  | 18 | 318 / 0.057 / 1 / 2.4 | 318 / 0.022 / 0 / 0.8 | — |
| B gj-908 0.1AU | dev | 18 | 373 / 0.059 / 0 / 2.6 | 373 / 0.038 / 0 / 1.0 | — |
| B ross-128 2.5Rsun |  | 15 | 1528 / 0.035 / 0 / 2.5 | 1528 / 0.022 / 0 / 1.6 | — |
| B ross-128 0.1AU |  | 18 | 304 / 0.066 / 0 / 1.1 | 304 / 0.016 / 0 / 1.0 | — |
| B ross-154 0.1AU |  | 15 | 263 / 0.084 / 0 / 1.2 | 263 / 0.015 / 0 / 1.7 | — |
| B teegarden 1.2Rsun |  | 12 | 1118 / 0.033 / 0 / 1.6 | 1118 / 0.028 / 1 / 2.0 | — |
| B teegarden 2.5Rsun |  | 13 | 985 / 0.035 / 0 / 1.4 | 985 / 0.012 / 0 / 1.6 | — |
| B teegarden 0.1AU |  | 14 | 234 / 0.060 / 0 / 1.6 | 234 / 0.043 / 0 / 1.1 | — |
| B van-maanen 1.2Rsun |  | 16 | 1444 / 0.073 / 17 / 26.5 | 1444 / 0.021 / 0 / 2.8 | S_stack_L, S_event_L, S_burst |
| B van-maanen 2.5Rsun |  | 16 | 1084 / 0.101 / 21 / 27.0 | 1084 / 0.017 / 0 / 2.2 | S_stack_L, S_event_L, S_burst |
| B van-maanen 0.1AU |  | 16 | 260 / 0.204 / 33 / 76.5 | 260 / 0.054 / 1 / 1.0 | S_stack_L, S_event_L, S_burst |
| B wolf-359 1.2Rsun |  | 16 | 1482 / 0.025 / 1 / 1.9 | 1482 / 0.013 / 0 / 1.6 | — |
| B wolf-359 2.5Rsun |  | 17 | 1192 / 0.035 / 3 / 2.1 | 1192 / 0.021 / 0 / 1.6 | — |
| B wolf-359 0.1AU |  | 18 | 293 / 0.061 / 0 / 2.0 | 293 / 0.038 / 0 / 0.7 | — |

Ensemble: 55474 members, 86 above T vs 19.6 expected; 71 of the 86 are the B van-maanen lane-L members (3C 279 at 1.8° from the antipode), the rest 15 vs 18.6.

Constraint-only trials: 9 of 145.