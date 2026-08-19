---
title: "WISE overlay v1 — work queue for the universal target list"
date: 2026-08-19
status: "24 queued / 8 searched / 5 deferred; no exclusions (all-sky)"
---

# WISE overlay v1

Applies WISE-specific grades to `targets/universal_v1` without changing
membership. Built by `build_wise_overlay.py`; per-corridor measurements
(CatWISE2020 source count in r=0.2°, brightest 2MASS Ks in r=0.3°) are
snapshotted in `corridor_measurements.json`. WISE is all-sky, so no
corridor is excluded; coverage uses |ecliptic β| as proxy (exact epoch
counts arrive with each corridor's coarse-discovery run). Queue order:
solution gate, then measured confusion (bright-star corridors last —
the shakedown showed a Ks < 4 star costs ~4 mag of depth), then
distance.

**Grandfathered/searched (8 system nodes = 9 corridors, incl. both
Proxima and α Cen):** all previously searched systems re-selected
themselves into the universal list; none needed a grandfather override.
They are marked `searched` rather than re-queued.

| System | d (pc) | Tracks | CatWISE ρ (deg⁻²) | Confusion | Coverage | Status |
| --- | --- | --- | --- | --- | --- | --- |
| Proxima Cen | 1.30 | 3 | 74,336 | high+bright-star | mid | searched |
| Barnard's Star | 1.83 | 1 | 46,054 | mid | mid | searched |
| Luhman 16 AB | 2.00 | 2 | 76,038 | high | mid | queued |
| WISE 0855-0714 | 2.28 | 1 | 45,386 | mid | baseline | deferred |
| Wolf 359 | 2.41 | 1 | 27,836 | low | baseline | searched |
| Lalande 21185 | 2.55 | 1 | 30,382 | low | mid | searched |
| Sirius AB | 2.67 | 2 | 83,174 | high | mid | searched |
| GJ 65 AB | 2.70 | 2 | 31,448 | low | mid | searched |
| Ross 154 | 2.98 | 1 | 48,266 | mid | baseline | searched |
| Ross 248 | 3.16 | 1 | 52,418 | mid | mid | searched |
| eps Eri | 3.22 | 1 | 37,422 | mid | mid | queued |
| Lacaille 9352 | 3.29 | 1 | 32,411 | mid | mid | queued |
| Ross 128 | 3.37 | 1 | 28,528 | low | baseline | queued |
| EZ Aqr | 3.41 | 3 | 35,115 | mid | baseline | deferred |
| 61 Cyg AB | 3.50 | 2 | 69,292 | high | mid | queued |
| Procyon AB | 3.51 | 1 | 76,826 | high | baseline | queued |
| Struve 2398 AB | 3.52 | 2 | 52,904 | high | high | deferred |
| Groombridge 34 AB | 3.56 | 2 | 52,944 | high | mid | queued |
| GJ 1111 | 3.58 | 1 | 39,125 | mid+bright-star | baseline | queued |
| eps Ind | 3.64 | 2 | 32,323 | low | mid | queued |
| tau Cet | 3.65 | 1 | 28,719 | low | baseline | queued |
| GJ 1061 | 3.67 | 1 | 37,733 | mid | high | queued |
| GJ 54 | 3.72 | 1 | 29,602 | low | baseline | queued |
| Luyten's Star | 3.79 | 1 | 82,904 | high | baseline | queued |
| Teegarden's Star | 3.83 | 1 | 30,517 | low | baseline | queued |
| Kapteyn's Star | 3.93 | 1 | 46,460 | mid+bright-star | high | queued |
| Lacaille 8760 | 3.97 | 1 | 31,806 | low | baseline | queued |
| GJ 12724 | 4.00 | 1 | 36,006 | mid | mid | queued |
| Wolf 1061 | 4.31 | 1 | 33,962 | mid | baseline | queued |
| van Maanen's Star | 4.31 | 1 | 30,398 | low | baseline | queued |
| LP 145-141 (WD) | 4.64 | 1 | 80,159 | high | mid | queued |
| GJ 908 | 5.91 | 1 | 29,356 | low | baseline | queued |
| GJ 783 | 6.01 | 2 | 32,331 | low | baseline | deferred |
| GJ 784 | 6.16 | 1 | 32,339 | low | baseline | queued |
| GJ 1221 | 6.21 | 1 | 98,616 | high | high | queued |
| GJ 9193 | 6.44 | 1 | 73,604 | high | mid | queued |
| GJ 11068 | 6.80 | 1 | 78,202 | high | mid | deferred |

**Deferred (named unblocking conditions):** WISE 0855-0714 and the
GJ 11068 approacher need non-Gaia astrometric solutions ("hard");
EZ Aqr needs a usable triple-system solution; Struve 2398 AB and GJ 783
need published orbits curated into the registry (the GJ 65 playbook).

Next batch recommendation: the nine `low`-confusion clean corridors
(Ross 128 → GJ 784) as WISE batch 2 — 10 tracks, all single-Gaia-row
curation, i.e. the cheap kind.
