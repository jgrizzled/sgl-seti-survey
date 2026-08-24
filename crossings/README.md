# Universal historical beam-crossing list (Pipeline B)

Produced by `python -m sglsurvey.crossings` (project plan §3.5). One row
per local minimum of the observer–axis impact parameter for every
registry endpoint (`registries/pilot_wise_2026.yaml`, 88 endpoints) and
both link directions (`inbound` = star→relay uplink, apparent axis at
`t_o`; `outbound` = relay→star downlink, axis at `t_o + 2d/c`), model
`tusay2022_eq5_7_v1` / axis `sun_star_axis_v1`, astropy built-in
ephemeris.

| product | observer | window | notes |
| --- | --- | --- | --- |
| `universal_v1/` | Earth center | 1980-01-01 → 2028-01-01 UTC | canonical; ground surveys and LEO |

Survey-independent by construction: no beam radius, no `report_max_b_au`
cut, ephemeris-coverage failures kept as `invalid` rows, window-boundary
minima flagged `degraded`. Each `surveys/<name>-crossings/` sub-project
intersects `events.ecsv` with that archive's exposure coverage and
applies beam-radius / wavelength / duty-cycle hypotheses there. Spacecraft
observers (WISE, SPHEREx, TESS…) get derivative runs keyed by
`observer_id`; the impact parameter at the interesting (≲ R☉) events is
sensitive to ~0.01 AU observer offsets, so do not reuse the Earth-center
list for them without checking.

Impact-parameter uncertainty is **not** propagated (`uncertainty_not_propagated`
warning); `sglseti.crossing_uncertainty` exists for per-event follow-up.
