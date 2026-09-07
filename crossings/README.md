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
| `kepler_v1/` (`xng-07012bdf124d`) | Kepler spacecraft (Horizons −227 SSB vectors, 6 h, `observers/`, sha-pinned) | 2009-05-01 → 2018-11-01 | 3,352 | Earth-trailing heliocentric orbit (period 372.5 d): 0.04 AU (2009) → 1.14 AU (2018) from Earth, so Earth-center is **invalid at every rung** — matched events shift by \|Δt_ca\| median 33 d / max 117 d and \|Δb\| up to 1 AU vs `universal_v1`. Built for the plan §5.8 item 8 footprint intersect (`surveys/kepler-crossings/`); channel positions at the searched rungs sit ≤ 6° from the spacecraft's anti-sun point, which no Kepler/K2 field ever contained |
| `stereoa_v1/` (`xng-3dd5b764776e`) | STEREO-A spacecraft (Horizons −234 SSB vectors, 6 h, `observers/`, sha-pinned) | 2007-01-01 → 2026-12-01 (stop = Horizons predicted-trajectory end 2026-12-18; extend at refresh) | 7,516 | Heliocentric 0.96 AU orbit leading Earth by ~22°/yr (full circuit 2006–2023): Earth-center is **invalid at every rung** — matched events shift by \|Δt_ca\| median 51 d (1,241 events have no Earth counterpart within 120 d). Built for the plan §5.15 O2 HI-1 sunward survey (`surveys/stereo-hi-crossings/`): S1/S2 0.1 AU rungs 147/148 events over the same 7 systems as `soho_v1` (van-maanen S1 min b 0.09 R☉); grazing rungs 84/105 events, never visible to HI-1 (4° inner edge). Census: `surveys/stereo-hi-crossings/results/stereoa_census_v1.json` |
| `psp_v1/` (`xng-c339df63a608`) | Parker Solar Probe (Horizons −96 SSB vectors, **10 min** in yearly chunks, `observers/`, sha-pinned) | 2018-08-15 → 2026-12-01 (stop = Horizons SPK end; extend at refresh) | 7,448 | Heliocentric 0.046–0.73 AU orbit (88–150 d): built with `--coarse-step-days 0.5` because the star-side and antipode-side minima of one axis are ~2 d apart around perihelion, and 10-min observer sampling because PSP sweeps ~34°/6 h at perihelion. Earth-center is a different population, not a correction (3,724 sunward events; nearest Earth pairs |Δt_ca| median 62 d, |Δb| median 0.35 AU, 1,319 unpaired). PSP crosses every axis twice per orbit at a per-target repeating radius: grazing family from PSP = wolf-359 / ross-128 / gj-1111 / gj-251 / gj-518 / wolf-437 (perihelion side, b 0.5–2.1 R☉ at r 0.047–0.056 AU) + teegarden / van-maanen / ross-154. Built for the plan §5.15 O5 PSP/WISPR geometry pass (`surveys/wispr-crossings/`): 0.1 AU sunward rungs 1,513 / 1,701 events (S1 / S2), 278 degraded (window-boundary minima of a fast observer) |
| `solo_v1/` (`xng-dea7100a725c`) | Solar Orbiter (Horizons −144 SSB vectors, **10 min** in yearly chunks, `observers/`, sha-pinned) | 2020-05-01 → 2028-01-01 (SoloHI first light → universal window end; SPK to 2030-11-20) | 4,624 | Heliocentric 0.28–1.0 AU Venus-resonant orbit (150–230 d): built with the PSP settings (`--coarse-step-days 0.5`, 10-min sampling; 39 min). Earth-center is a different population (2,313 sunward events; nearest Earth pairs |Δt_ca| median 48 d, |Δb| median 0.20 AU, 910 unpaired). Crosses every axis twice per orbit at a per-target repeating radius; impact parameters step at each Venus flyby (ross-128 b 0.13 R☉ in P05–P09 → 11 R☉ from P10). Built for the plan §5.23 SoloHI geometry pass (`surveys/solohi-crossings/`): 0.1 AU sunward rungs 268 / 182 events (S1 / S2), 98 degraded |

All nine share the product schema (events/windows ecsv, result.json,
manifest.json, summary.json) and the survey-independent construction: no
beam radius, no `report_max_b_au` cut, ephemeris-coverage failures kept
as `invalid` rows (0 in all nine), window-boundary minima flagged
`degraded`. Each `surveys/<name>-crossings/` sub-project intersects the
matching `events.ecsv` with that archive's exposure coverage and applies
beam-radius / wavelength / duty-cycle hypotheses there. Spacecraft
tables and raw Horizons responses live under `observers/` with content
hashes recorded in each run's manifest inputs.

The lists themselves carry the `uncertainty_not_propagated` warning. For
`universal_v1` the target-state uncertainty is propagated as a companion
product by `python -m sglsurvey.crossings_uncertainty --product universal_v1`
(2026-09-06; plan §5.27): `uncertainty.ecsv` (one row per event, joinable
on `event_id`: 95 % bounds, sigmas and median shift on `b_min`, `t_ca` and
`v_perp`, side consistency, window-edge count, the list's nominal
validity), `uncertainty_samples.npz` (the 128 `b_min` draws per event) and
`uncertainty_summary.json` (parameters, input hashes, per-target maxima,
rung-membership census). 128 Monte Carlo draws of the schema-v2 registry
uncertainties (positions, proper motions, parallax, RV, orbital elements)
per event, re-minimised inside ±45 d of the nominal `t_ca` at 1-s
tolerance, seed 20260906 (+ per-event `event_id` hash); observer state,
ephemeris and the geometry-model floor stay `not_propagated`, as the
library labels them. **Result: no rung membership changes** — over the
16,381 nominally valid events the 95 % half-width on `b_min` is
≤ 0.012 R☉ everywhere (eps-ind-b, the 200 mas/yr allocation; median
σ 2 × 10⁻⁶ R☉) and ≤ 1.0 × 10⁻⁴ R☉ on the 960 grazing-family events
(≤ 2.5 R☉); `t_ca` σ median 0.2 s, p99 61 s, max 264 s; side-of-axis
100 % consistent; 0 events ambiguous at 1.2 R☉ / 2.5 R☉ / 0.1 AU / 1 AU.
The target-state term is therefore at or below the declared LEO observer
floor (0.010 R☉) — under it for every target but eps-ind-b, whose worst
event reaches 0.012 R☉ at b = 0.65 AU, far from any rung edge. The 205 interval-boundary rows the list
flags `degraded` get their refined minimum from the same run (121 resolve
inside ±45 d of the boundary, 84 lie further out and stay `degraded`) —
a bookkeeping aid for the yearly refresh, not an uncertainty. The
spacecraft lists share the registry, so the same bound applies to their
target-state term; a per-list run is `--product <name>` (~3 h on 30
workers for 16.6k events).
