# Universal historical beam-crossing list (Pipeline B)

Produced by `python -m sglsurvey.crossings` (project plan §3.5). One row
per local minimum of the observer–axis impact parameter for every
registry endpoint (`registries/pilot_wise_2026.yaml`, 88 endpoints) and
both link directions (`inbound` = star→relay uplink, apparent axis at
`t_o`; `outbound` = relay→star downlink, axis at `t_o + 2d/c`), model
`tusay2022_eq5_7_v1` / axis `sun_star_axis_v1`, astropy built-in
ephemeris.

| product | observer | window | events | notes |
| --- | --- | --- | --- | --- |
| `universal_v1/` (`xng-a09e2db7681d`) | Earth center | 1980-01-01 → 2028-01-01 | 16,586 | canonical; ground surveys and LEO |
| `wise_v1/` (`xng-ffc922f7bbca`) | WISE spacecraft (L1b-header table, `runs/wise/v2/observer/`, sha-pinned) | 2010-01-01 → 2024-08-01 | 5,170 | exact at frame epochs; Earth-center fallback in gaps > 2 d (incl. the 2011–2013 hibernation); bounded by the ~7,000 km LEO radius = 0.010 R☉ |
| `tess_v1/` (`xng-a943f0f3dbe4`) | TESS spacecraft (Horizons −95 SSB vectors, 6 h, `observers/`, sha-pinned) | 2018-07-04 → 2026-08-24 | 3,390 | HEO apogee 0.54 R☉ — Earth-center **invalid** at grazing b: matched events shift by \|Δb\| up to 0.33 R☉ (median 0.18) and \|Δt_ca\| up to 3.2 h (van-maanen validation); the 13.7-d wobble adds shallow extra minima for the ecliptic-pole b ≈ 1 AU family (sigma-dra 16 → 98 in-era events) — physical, kept by the every-minimum rule |
| `spherex_v1/` (`xng-a2773287b29d`) | Earth center (SPHEREx budget) | 2025-03-15 → 2028-01-01 | 1,130 | LEO ~650 km carried as a declared budget: \|Δb\| ≤ 0.010 R☉, \|Δt_ca\| ≤ ~4 min — the pipeline-A observer convention; refresh with the era as quick releases extend |
| `universal_1885_v1/` (`xng-e2d1063af9d0`) | Earth center | 1885-01-01 → 1993-01-01 | 37,096 | the pre-1980 backward extension (plan §5.8 item 5 prerequisite, DASCH era); overlaps `universal_v1` 1980–92 **by design** — consistency check: all 4,104 shared-era events match 1:1, max \|Δt_ca\| 23 s, max \|Δb\| 3.4e-7 R☉. 13,053 events degraded, dominated by `long_propagation_span` (the declared 75-yr linear-motion bound from the 2016.0 catalog epoch — conservative bookkeeping; the measured 1885 budget is ≤ 0.0002 R☉ from RV/µ sensitivity + ≤ 1e-5 R☉ ephemeris vs DE440S, see `surveys/dasch-crossings/notes/dasch_recon_2026-08-26.md`) |
| `soho_v1/` (`xng-298d55d0ce6b`) | SOHO spacecraft (Horizons −21 SSB vectors, 6 h, `observers/`, sha-pinned) | 1996-01-01 → 2026-10-01 (stop = Horizons SPK end 2026-10-05; extend at refresh) | 10,714 | L1 halo cross-track ~0.9 R☉ → vs Earth-center: grazing-family \|Δb\| median 0.12 / max 0.21 R☉ (large-b events shift up to 2.7 R☉ via the projected 0.01 AU radial offset); sunward grazing-rung membership unchanged (117 ev ≤ 1.2 R☉ / 153 ≤ 2.5 R☉, same 4/5 targets); the van-maanen annual October S1 family deepens to b 0.13–0.24 R☉. Census: `surveys/heliospheric-crossings/results/soho_census_v1.json` |

All four share the product schema (events/windows ecsv, result.json,
manifest.json, summary.json) and the survey-independent construction: no
beam radius, no `report_max_b_au` cut, ephemeris-coverage failures kept
as `invalid` rows (0 in all four), window-boundary minima flagged
`degraded`. Each `surveys/<name>-crossings/` sub-project intersects the
matching `events.ecsv` with that archive's exposure coverage and applies
beam-radius / wavelength / duty-cycle hypotheses there. Spacecraft
tables and raw Horizons responses live under `observers/` with content
hashes recorded in each run's manifest inputs.

Impact-parameter uncertainty is **not** propagated (`uncertainty_not_propagated`
warning); `sglseti.crossing_uncertainty` exists for per-event follow-up.
