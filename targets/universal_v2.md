---
title: "Universal SGL target list — v2"
date: 2026-08-20
status: "76 systems / 95 component tracks; 39 systems not yet searched in WISE"
---

# Universal target list v2

Survey-agnostic target portfolio built by `scripts/build_universal_list.py` from `config.yaml` and the curated `science_interest.yaml`. Census: CNS5 within 10 pc (`census/cns5_10pc_2026-08-20.csv`), engineering joins in `census/engineering_2026-08-20.csv`. The network prior, the engineering/desirability columns, and archive searchability are separate columns, never merged.

Baskets: **Distance Core** (nearest 25 systems, protected), **Robust Geometric Neighbor** (≥3 of 7 sparse-graph families among all systems), **Selective Network Neighbor** (same vote on the *desirable* subgraph — level 1 drops close multiples, substellar objects and evolved stars; level 2 also drops intermediate multiples, detected accelerations, active and rapidly rotating stars), **Engineering Backbone** (top 12 rank-sum of compactness M/R², quietness, dynamical cleanliness, lifetime and distance among level-2 dwarfs), **Science Interest** (curated tier 1 / tier 2), **Historical Neighbor** (closest approach within ±1 Myr), **Compact Lens Wildcard** (isolated white dwarfs).

Desirability is the network's presumed taste and is distinct from the solution gate (our ability to predict the track).

| System | d (pc) | Tracks | Mult. | SpT | M/R² | Activity | Accel. | Votes all / L1 / L2 | Eng. score | Gate | WISE | Baskets |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Proxima Cen | 1.30 | 3 | close | M5.5Ve | 5.3 | active | clean | 6 / 0 / 0 | — | clean | done | Distance Core, Historical Neighbor, Robust Geometric Neighbor |
| Barnard's Star | 1.83 | 1 | single | M4V | 4.3 | quiet | clean | 6 / 6 / 7 | 0.73 | clean | done | Distance Core, Engineering Backbone, Robust Geometric Neighbor |
| Luhman 16 AB | 2.00 | 2 | close | — | — | unknown | unknown | 2 / 0 / 0 | — | clean | done | Distance Core |
| WISE 0855-0714 | 2.28 | 1 | single | — | — | unknown | unknown | 2 / 0 / 0 | — | hard | done | Distance Core |
| Wolf 359 | 2.41 | 1 | single | — | 6.0 | active | clean ruwe only | 2 / 3 / 0 | — | clean | done | Distance Core, Selective Network Neighbor |
| Lalande 21185 | 2.55 | 1 | single | M2+V | 2.6 | quiet | clean | 2 / 2 / 5 | 0.68 | clean | done | Distance Core, Engineering Backbone, Selective Network Neighbor |
| Sirius AB | 2.67 | 2 | close | — | 3840.0 | active | accel detected | 2 / 0 / 0 | — | clean | done | Distance Core |
| GJ 65 AB | 2.70 | 2 | close | — | 4.3 | active | accel detected | 3 / 0 / 0 | — | clean | done | Distance Core, Robust Geometric Neighbor |
| Ross 154 | 2.98 | 1 | single | M3.5Ve | 4.0 | active | clean | 1 / 1 / 0 | — | clean | done | Distance Core |
| Ross 248 | 3.16 | 1 | single | — | 4.7 | active | clean ruwe only | 3 / 3 / 0 | — | clean | done | Distance Core, Historical Neighbor, Robust Geometric Neighbor |
| eps Eri | 3.22 | 1 | single | K2V | 1.4 | active | clean | 1 / 3 / 0 | — | clean | done | Distance Core, Science Interest, Selective Network Neighbor |
| Lacaille 9352 | 3.29 | 1 | single | M2V | 2.1 | quiet | clean | 2 / 2 / 3 | 0.69 | clean | done | Distance Core, Engineering Backbone, Science Interest, Selective Network Neighbor |
| Ross 128 | 3.37 | 1 | single | M4V | 4.1 | quiet | clean | 0 / 0 / 1 | 0.60 | clean | done | Distance Core, Science Interest |
| EZ Aqr | 3.41 | 3 | close | — | — | unknown | unknown | 2 / 0 / 0 | — | hard | done | Distance Core |
| 61 Cyg AB | 3.50 | 2 | intermediate | K5V | 1.4 | active | accel detected | 1 / 1 / 0 | — | clean | done | Distance Core |
| Procyon AB | 3.51 | 1 | single | — | — | active | unknown | 0 / 0 / 0 | — | clean | done | Distance Core |
| Struve 2398 AB | 3.52 | 2 | close | M3V | 2.7 | quiet | accel detected | 1 / 0 / 0 | — | orbit_needed | done | Distance Core, Science Interest (tier 2) |
| Groombridge 34 AB | 3.56 | 2 | intermediate | M2V | 2.4 | quiet | accel detected | 0 / 0 / 0 | — | clean | done | Distance Core, Science Interest (tier 2) |
| GJ 1111 | 3.58 | 1 | single | — | 6.5 | unknown | clean ruwe only | 1 / 1 / 1 | — | clean | done | Distance Core, Historical Neighbor |
| eps Ind | 3.64 | 2 | wide | K5V | 1.3 | active | accel detected | 0 / 1 / 0 | — | clean | done | Distance Core, Science Interest |
| tau Cet | 3.65 | 1 | single | G8V | 1.1 | quiet | clean | 0 / 1 / 2 | 0.71 | clean | done | Distance Core, Engineering Backbone, Science Interest |
| GJ 1061 | 3.67 | 1 | single | — | 5.3 | quiet | clean ruwe only | 0 / 2 / 2 | 0.52 | clean | done | Distance Core, Science Interest |
| YZ Cet | 3.72 | 1 | single | M4.0Ve | 4.9 | active | clean | 0 / 1 / 0 | — | clean | done | Distance Core, Science Interest |
| Luyten's Star | 3.79 | 1 | single | M3.5V | 3.0 | quiet | clean | 0 / 2 / 2 | 0.69 | clean | done | Distance Core, Engineering Backbone, Science Interest |
| Teegarden's Star | 3.83 | 1 | single | — | 6.7 | unknown | clean ruwe only | 1 / 1 / 2 | — | clean | done | Distance Core |
| Kapteyn's Star | 3.93 | 1 | single | M1VIp | 3.1 | quiet | clean | 0 / 1 / 1 | 0.65 | clean | done | Engineering Backbone, Historical Neighbor, Science Interest (tier 2) |
| Lacaille 8760 | 3.97 | 1 | single | M1V | 1.7 | active | clean | 0 / 0 / 0 | — | clean | done | Historical Neighbor |
| GJ 12724 | 4.00 | 1 | single | — | 6.8 | unknown | clean ruwe only | 0 / 1 / 2 | — | clean | done | Historical Neighbor |
| GJ 11547 | 4.04 | 1 | single | — | 6.9 | unknown | clean ruwe only | 0 / 1 / 2 | — | clean | queue | Historical Neighbor |
| Wolf 1061 | 4.31 | 1 | single | M3V | 2.9 | quiet | clean | 0 / 0 / 0 | 0.63 | clean | done | Engineering Backbone, Historical Neighbor, Science Interest (tier 2) |
| van Maanen's Star | 4.31 | 1 | single | DZ7.5 | 3840.0 | unknown | clean | 0 / 0 / 1 | — | clean | done | Compact Lens Wildcard |
| GJ 687 | 4.55 | 1 | single | M3.0V | 2.4 | quiet | clean | 0 / 1 / 1 | 0.64 | clean | queue | Engineering Backbone, Science Interest (tier 2) |
| GJ 674 | 4.55 | 1 | single | M3V | 2.6 | active | clean | 0 / 0 / 0 | — | clean | queue | Science Interest (tier 2) |
| LP 145-141 (WD) | 4.64 | 1 | single | DQ | 3840.0 | unknown | clean | 0 / 1 / 1 | — | clean | done | Compact Lens Wildcard |
| GJ 876 | 4.67 | 1 | single | M3.5V | 2.7 | quiet | clean | 0 / 0 / 0 | 0.63 | clean | queue | Engineering Backbone, Science Interest |
| GJ 1002 | 4.85 | 1 | single | — | 5.7 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Science Interest |
| GJ 832 | 4.97 | 1 | single | M2/3V | 2.2 | quiet | accel detected | 0 / 0 / 0 | — | clean | queue | Science Interest (tier 2) |
| GJ 682 | 5.01 | 1 | single | M3.5 | 3.1 | quiet | clean | 0 / 0 / 0 | 0.58 | clean | queue | Science Interest (tier 2) |
| LHS 1723 | 5.37 | 1 | single | — | 4.2 | active | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Science Interest |
| GJ 526 | 5.43 | 1 | single | M2V | 2.0 | quiet | clean | 0 / 0 / 0 | 0.64 | clean | queue | Engineering Backbone |
| GJ 251 | 5.58 | 1 | single | M3V | 2.6 | quiet | clean | 0 / 0 / 0 | 0.61 | clean | queue | Engineering Backbone, Science Interest |
| GJ 229 AB | 5.76 | 2 | close | M1V | 1.8 | active | accel detected | 0 / 0 / 0 | — | orbit_needed | queue | Science Interest (tier 2) |
| sigma Dra | 5.76 | 1 | single | K0V | 1.2 | quiet | clean | 0 / 0 / 0 | 0.64 | clean | queue | Engineering Backbone |
| GJ 908 | 5.91 | 1 | single | M1VFe-1 | 2.4 | quiet | clean | 0 / 0 / 0 | 0.62 | clean | done | Engineering Backbone, Historical Neighbor |
| GJ 588 | 5.92 | 1 | single | M2.5V | 2.2 | quiet | clean | 0 / 0 / 0 | 0.61 | clean | queue | Engineering Backbone |
| GJ 783 | 6.01 | 2 | close | K2.5V | 1.5 | quiet | clean | 0 / 0 / 0 | — | orbit_needed | done | Historical Neighbor |
| 82 Eri | 6.04 | 1 | single | G6V | 1.1 | quiet | clean | 0 / 0 / 0 | 0.64 | clean | queue | Engineering Backbone, Historical Neighbor, Science Interest |
| GJ 784 | 6.16 | 1 | single | M0V | 1.8 | quiet | clean | 0 / 0 / 0 | 0.57 | clean | done | Historical Neighbor |
| GJ 1221 | 6.21 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | done | Compact Lens Wildcard |
| GJ 581 | 6.30 | 1 | single | M3V | 2.9 | quiet | clean | 0 / 0 / 0 | 0.54 | clean | queue | Science Interest |
| GJ 338 AB | 6.33 | 2 | intermediate | K7V | 1.6 | active | accel detected | 0 / 0 / 0 | — | clean | queue | Science Interest (tier 2) |
| GJ 9193 | 6.44 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | done | Compact Lens Wildcard |
| GJ 625 | 6.48 | 1 | single | M1.5V | 2.8 | quiet | clean | 0 / 0 / 0 | 0.60 | clean | queue | Engineering Backbone, Science Interest (tier 2) |
| HD 219134 | 6.54 | 1 | single | K3V | 1.4 | quiet | clean | 0 / 0 / 0 | 0.61 | clean | queue | Engineering Backbone, Science Interest |
| GJ 11068 | 6.80 | 1 | single | — | — | unknown | unknown | 0 / 0 / 0 | — | hard | done | Historical Neighbor |
| LTT 1445 ABC | 6.86 | 3 | close | M3.5+M3.0 | 3.3 | active | accel detected | 0 / 0 / 0 | — | orbit_needed | queue | Science Interest |
| GJ 667 ABC | 7.24 | 3 | close | — | 2.8 | active | clean ruwe only | 0 / 0 / 0 | — | orbit_needed | queue | Science Interest |
| GJ 514 | 7.62 | 1 | single | M1.0Ve | 2.0 | quiet | clean | 0 / 0 / 0 | 0.52 | clean | queue | Science Interest (tier 2) |
| Fomalhaut | 7.70 | 1 | single | — | — | unknown | unknown | 0 / 0 / 0 | — | clean | queue | Science Interest |
| Wolf 437 | 8.08 | 1 | single | M3.5Ve | 2.9 | active | clean | 0 / 0 / 0 | — | clean | queue | Science Interest |
| GJ 1087 (WD) | 8.11 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 293 (WD) | 8.17 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 66 | 8.19 | 2 | close | — | 1.4 | quiet | clean ruwe only | 0 / 0 / 0 | — | orbit_needed | queue | Historical Neighbor |
| GJ 915 | 8.33 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 518 (WD) | 8.35 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 13157 | 8.46 | 1 | single | — | 3840.0 | active | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 318 | 8.51 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 1276 | 8.53 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| 61 Vir | 8.54 | 1 | single | G6.5V | 1.0 | quiet | clean | 0 / 0 / 0 | 0.53 | clean | queue | Science Interest |
| GJ 2066 | 8.94 | 1 | single | M2.0V | 2.2 | quiet | clean | 0 / 0 / 0 | — | clean | queue | Historical Neighbor |
| GJ 2012 | 9.09 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 3306 | 9.40 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |
| GJ 367 | 9.42 | 1 | single | M1.0 | 2.2 | quiet | clean | 0 / 0 / 0 | — | clean | queue | Historical Neighbor, Science Interest |
| GJ 3512 | 9.49 | 1 | single | — | 5.2 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Science Interest |
| Wolf 1069 | 9.57 | 1 | single | — | 4.3 | quiet | accel detected | 0 / 0 / 0 | — | clean | queue | Science Interest |
| GJ 3112 | 9.71 | 1 | single | — | 3840.0 | unknown | clean ruwe only | 0 / 0 / 0 | — | clean | queue | Compact Lens Wildcard |

## Science-interest rationale

- **eps Eri** (tier 1): giant planet + multi-belt debris system (resource + architecture) — Mawet 2019; Su 2017
- **Lacaille 9352** (tier 1): 2-3 super-Earths, quiet host, HZ candidate — Jeffers 2020
- **Ross 128** (tier 1): temperate Earth-mass planet b — Bonfils 2018
- **Struve 2398 AB** (tier 2): B hosts 2 planets — Feng 2020
- **Groombridge 34 AB** (tier 2): GJ 15A b super-Earth + c cold planet — Howard 2014; Pinamonti 2018
- **eps Ind** (tier 1): cold Jupiter imaged (JWST) + T-dwarf binary — Matthews 2024
- **tau Cet** (tier 1): 4 candidate planets + massive debris disk — Feng 2017; MacGregor 2016
- **GJ 1061** (tier 1): 3 low-mass planets, d in HZ — Dreizler 2020
- **YZ Cet** (tier 1): 3 Earth-mass planets, compact system; innermost terrestrial system known — Astudillo-Defru 2017; Stock 2020
- **Luyten's Star** (tier 1): GJ 273b low-mass planet in/near HZ — Astudillo-Defru 2017
- **Kapteyn's Star** (tier 2): halo-star planet claims (disputed); oldest nearby planetary candidate — Anglada-Escude 2014; Bortle 2021
- **Wolf 1061** (tier 2): 3 planets, c near HZ — Wright 2016
- **GJ 687** (tier 2): 2 Neptune-mass planets — Burt 2014; Feng 2020
- **GJ 674** (tier 2): hot Neptune — Bonfils 2007
- **GJ 876** (tier 1): Laplace-resonant 4-planet system; science overrides lens cleanliness — Rivera 2010; Millholland 2018
- **GJ 1002** (tier 1): two Earth-mass temperate planets around a quiet M5.5 — Suarez Mascareno 2023
- **GJ 832** (tier 2): cold Jupiter analogue — Bailey 2009
- **GJ 682** (tier 2): 2 candidate planets incl. HZ (disputed) — Tuomi 2014
- **LHS 1723** (tier 1): 2 planets, b in HZ — Astudillo-Defru 2017
- **GJ 251** (tier 1): super-Earth b + HZ candidate c; quiet M3 — Stock 2020; Burt 2025
- **GJ 229 AB** (tier 2): first T dwarf (now a T+T binary) + planet candidates — Nakajima 1995; Xuan 2024
- **82 Eri** (tier 1): 3 planets incl. HZ super-Earth d; quiet G6V — Pepe 2011; Nari 2025
- **GJ 581** (tier 1): 3 confirmed planets incl. HZ-edge e/c; archetype RV system — Mayor 2009; Robertson 2014
- **GJ 338 AB** (tier 2): wide binary, super-Earth around B — Gonzalez-Alvarez 2020
- **GJ 625** (tier 2): super-Earth at HZ inner edge — Suarez Mascareno 2017
- **HD 219134** (tier 1): 6 planets incl. 2 transiting rocky; quiet K3V; nearest transiting system — Motalebi 2015; Gillon 2017
- **LTT 1445 ABC** (tier 1): A hosts 2 transiting rocky planets; nearest transiting M-dwarf triple — Winters 2019, 2022
- **GJ 667 ABC** (tier 1): C hosts multi-planet system with HZ super-Earth c; triple (lens = C) — Anglada-Escude 2013
- **GJ 514** (tier 2): eccentric super-Earth crossing the HZ — Damasso 2022
- **Fomalhaut** (tier 1): eccentric debris ring + wide triple; exceptional resource/architecture system (A-type, short-lived) — Kalas 2008; Gaspar 2023
- **Wolf 437** (tier 1): transiting rocky planet with best-characterised M-dwarf-planet atmosphere constraints — Trifonov 2021; Moran 2023
- **61 Vir** (tier 1): 3 low-mass planets + debris disk; solar twin-ish G5V — Vogt 2010; Wyatt 2012
- **GJ 367** (tier 1): ultra-short-period iron sub-Earth (transiting) + 2 outer planets — Lam 2021; Goffo 2023
- **GJ 3512** (tier 1): Jupiter-mass planet around 0.12 Msun star — anomalous architecture — Morales 2019
- **Wolf 1069** (tier 1): Earth-mass planet in HZ of quiet M5 — Kossakowski 2023

## Deferred by track budget

- Vega: Science Interest (tier 2)
- GJ 1151: Science Interest (tier 2)
- GJ 686: Science Interest (tier 2)
- GJ 849: Science Interest (tier 2)
- GJ 357: Science Interest (tier 2)
- GJ 176: Science Interest (tier 2)
- GJ 436: Science Interest (tier 2)

## Changes from v1

Added (39):

- GJ 11547 (4.04 pc): Historical Neighbor
- GJ 687 (4.55 pc): Engineering Backbone, Science Interest (tier 2)
- GJ 674 (4.55 pc): Science Interest (tier 2)
- GJ 876 (4.67 pc): Engineering Backbone, Science Interest
- GJ 1002 (4.85 pc): Science Interest
- GJ 832 (4.97 pc): Science Interest (tier 2)
- GJ 682 (5.01 pc): Science Interest (tier 2)
- LHS 1723 (5.37 pc): Science Interest
- GJ 526 (5.43 pc): Engineering Backbone
- GJ 251 (5.58 pc): Engineering Backbone, Science Interest
- GJ 229 AB (5.76 pc): Science Interest (tier 2)
- sigma Dra (5.76 pc): Engineering Backbone
- GJ 588 (5.92 pc): Engineering Backbone
- 82 Eri (6.04 pc): Engineering Backbone, Historical Neighbor, Science Interest
- GJ 581 (6.30 pc): Science Interest
- GJ 338 AB (6.33 pc): Science Interest (tier 2)
- GJ 625 (6.48 pc): Engineering Backbone, Science Interest (tier 2)
- HD 219134 (6.54 pc): Engineering Backbone, Science Interest
- LTT 1445 ABC (6.86 pc): Science Interest
- GJ 667 ABC (7.24 pc): Science Interest
- GJ 514 (7.62 pc): Science Interest (tier 2)
- Fomalhaut (7.70 pc): Science Interest
- Wolf 437 (8.08 pc): Science Interest
- GJ 1087 (WD) (8.11 pc): Compact Lens Wildcard
- GJ 293 (WD) (8.17 pc): Compact Lens Wildcard
- GJ 66 (8.19 pc): Historical Neighbor
- GJ 915 (8.33 pc): Compact Lens Wildcard
- GJ 518 (WD) (8.35 pc): Compact Lens Wildcard
- GJ 13157 (8.46 pc): Compact Lens Wildcard
- GJ 318 (8.51 pc): Compact Lens Wildcard
- GJ 1276 (8.53 pc): Compact Lens Wildcard
- 61 Vir (8.54 pc): Science Interest
- GJ 2066 (8.94 pc): Historical Neighbor
- GJ 2012 (9.09 pc): Compact Lens Wildcard
- GJ 3306 (9.40 pc): Compact Lens Wildcard
- GJ 367 (9.42 pc): Historical Neighbor, Science Interest
- GJ 3512 (9.49 pc): Science Interest
- Wolf 1069 (9.57 pc): Science Interest
- GJ 3112 (9.71 pc): Compact Lens Wildcard

Dropped (0): none

Renamed: GJ 54 → YZ Cet (GJ 54.1; v1 truncated decimal GJ ids).

## Notes

- Mass/radius provenance per system is in the JSON (`mass_radius_source`: mann_mks / tic / flame / wd_nominal / unknown).
- Activity evidence strings record which indicator decided the class; `unknown` is not `quiet`. ROSAT non-detections are recorded as Lx/Lbol upper limits and count as quiet only when the limit is below the active threshold.
- Adjacency votes are computed against all systems in the census (all / level-1 / level-2 subgraphs), so only ~15 systems can hold Sun-Delaunay edges at once; a vote of ≥3 is meaningful.
