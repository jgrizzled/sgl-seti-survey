---
title: "Universal SGL target list — v1"
date: 2026-08-19
status: "37 systems / 49 component tracks; survey-agnostic"
---

# Universal target list v1

Survey-agnostic target portfolio built by
`scripts/build_universal_list.py` from `config.yaml` (method:
`notes/sgl_seti_star_ranking_methods.md`, simplified staged design).
Census: CNS5 within 8 pc (200 objects → 176 systems, snapshot in
`census/`). The network prior and archive searchability are separate
columns, never merged; per-survey overlays (e.g.
`surveys/wise/targets/`) consume this list and order their own work
queues without changing membership.

Baskets: **Distance Core** (nearest 25 systems, protected, no
exclusions — includes substellar systems like Luhman 16 and
WISE 0855), **Robust Geometric Neighbor** (≥3 of 7 sparse-graph
families: Delaunay, Gabriel, RNG, MST, mutual-kNN k=1–3), **Selective
Network Neighbor** (new Sun-edges after dropping flagged hosts),
**Historical Neighbor** (linear-motion closest approach within
±1 Myr), **Compact Lens Wildcard** (isolated white dwarfs).

Solution gates: *clean* (Gaia solution or curated orbit), *orbit-needed*
(close multiple without a usable published orbit — deferred with a
named unblocking condition), *hard* (no usable Gaia astrometry, e.g.
ultracool dwarfs needing special solutions).

| System | d (pc) | Tracks | Multiplicity | Adjacency votes | d_min ±1 Myr (pc) | Gate | Baskets |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Proxima Cen | 1.30 | 3 | close | 6 | 0.95 | clean | Distance Core, Historical Neighbor, Robust Geometric Neighbor |
| Barnard's Star | 1.83 | 1 | single | 6 | — | clean | Distance Core, Robust Geometric Neighbor |
| Luhman 16 AB | 2.00 | 2 | close | 2 | — | clean | Distance Core |
| WISE 0855-0714 | 2.28 | 1 | single | 2 | — | hard | Distance Core |
| Wolf 359 | 2.41 | 1 | single | 2 | 2.27 | clean | Distance Core, Selective Network Neighbor |
| Lalande 21185 | 2.55 | 1 | single | 2 | — | clean | Distance Core |
| Sirius AB | 2.67 | 2 | close | 2 | — | clean | Distance Core, Selective Network Neighbor |
| GJ 65 AB | 2.70 | 2 | close | 3 | 2.26 | clean | Distance Core, Robust Geometric Neighbor |
| Ross 154 | 2.98 | 1 | single | 1 | — | clean | Distance Core |
| Ross 248 | 3.16 | 1 | single | 3 | 0.93 | clean | Distance Core, Historical Neighbor, Robust Geometric Neighbor |
| eps Eri | 3.22 | 1 | single | 1 | — | clean | Distance Core |
| Lacaille 9352 | 3.29 | 1 | single | 2 | 3.28 | clean | Distance Core |
| Ross 128 | 3.37 | 1 | single | 0 | — | clean | Distance Core |
| EZ Aqr | 3.41 | 3 | close | 2 | — | hard | Distance Core |
| 61 Cyg AB | 3.50 | 2 | intermediate | 1 | 2.8 | clean | Distance Core |
| Procyon AB | 3.51 | 1 | single | 0 | — | clean | Distance Core |
| Struve 2398 AB | 3.52 | 2 | close | 1 | — | orbit-needed | Distance Core |
| Groombridge 34 AB | 3.56 | 2 | intermediate | 0 | 3.47 | clean | Distance Core |
| GJ 1111 | 3.58 | 1 | single | 1 | 1.33 | clean | Distance Core, Historical Neighbor |
| eps Ind | 3.64 | 2 | wide | 0 | 3.26 | clean | Distance Core |
| tau Cet | 3.65 | 1 | single | 0 | — | clean | Distance Core |
| GJ 1061 | 3.67 | 1 | single | 0 | — | clean | Distance Core |
| GJ 54 | 3.72 | 1 | single | 0 | — | clean | Distance Core |
| Luyten's Star | 3.79 | 1 | single | 0 | — | clean | Distance Core |
| Teegarden's Star | 3.83 | 1 | single | 1 | 3.17 | clean | Distance Core |
| Kapteyn's Star | 3.93 | 1 | single | 0 | 2.16 | clean | Historical Neighbor |
| Lacaille 8760 | 3.97 | 1 | single | 0 | 3.78 | clean | Historical Neighbor |
| GJ 12724 | 4.00 | 1 | single | 0 | 3.77 | clean | Historical Neighbor |
| Wolf 1061 | 4.31 | 1 | single | 0 | 3.22 | clean | Historical Neighbor |
| van Maanen's Star | 4.31 | 1 | single | 0 | — | clean | Compact Lens Wildcard |
| LP 145-141 (WD) | 4.64 | 1 | single | 0 | — | clean | Compact Lens Wildcard |
| GJ 908 | 5.91 | 1 | single | 0 | 2.81 | clean | Historical Neighbor |
| GJ 783 | 6.01 | 2 | close | 0 | 2.05 | orbit-needed | Historical Neighbor |
| GJ 784 | 6.16 | 1 | single | 0 | 3.51 | clean | Historical Neighbor |
| GJ 1221 | 6.21 | 1 | single | 0 | — | clean | Compact Lens Wildcard |
| GJ 9193 | 6.44 | 1 | single | 0 | — | clean | Compact Lens Wildcard |
| GJ 11068 | 6.80 | 1 | single | 0 | 0.33 | hard | Historical Neighbor |

Notes:

- All 9 systems searched in the WISE shakedown re-selected themselves
  into the distance core — grandfathering required no overrides.
- GJ 11068 (d_min 0.33 pc within ±1 Myr) is the standout
  historical-approach addition; Kapteyn's Star and GJ 783 also enter on
  past/future proximity.
- Component counts follow CNS5; close pairs merged into single CNS5
  rows (α Cen AB, Procyon AB) count their letters as tracks. Known
  limitation: companions lacking independent CNS5 parallaxes
  (e.g. Procyon B) are folded into their primary's entry.
- Adjacency votes are computed against all 176 systems, so only ~15
  systems can hold Sun-Delaunay edges; a vote of ≥3 is meaningful.
