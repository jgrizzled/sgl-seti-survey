# SoloHI dev machinery log (2026-09-06)

Dev units (D8): gj-908 S1+S2, 61-vir S1+S2, ross-154 S1 — 39 unit-events,
9,313 inner-tile frames (5,842 arc + 3,471 baseline) over 173 days;
index 54 s (Gaia DR3 cone queries for 351 patch positions, cached);
run 1.3 h at ~7,000 frames/h (8 NRL + 8 SOAR fetch threads, 12
measure workers; fetch-bound). Tile 2 frames are 960 × 1024
(`RECTROTA` 6) and tile 1 1024 × 960 — the library was made
shape-agnostic before any frame was measured (a v1.0 gate rejected
the first 359 tile-2 frames as `shape_960x1024`; that partial run was
discarded and dev restarted).

## v1.0 reduce (`results/dev_search_vtest.json`, evidence only)

Chain sane per frame (9,306/9,313 usable; Hipparcos match 248 stars,
rms 0.71 px; ZP −19.92 with MAD 0.106; stamp response 0.96 ± 0.06,
n 4,178; colour coefficient 0.366 on 2.7e5 calibrators) but the
star-fixed differential was broken on the sunward end of every arc:
61-vir S2 read 14–18 units in the arc against 101–106 in the
baseline, and control thresholds were S_stack T 20–190.

**Cause 1 — the saturation plateau.** Along the sunward edge of tile 1
(x ≳ 935, ε ≲ 9°) in the ≥ 45-s exposure regime the F-corona reaches a
flat plateau at 0.90–0.92 × the header `DSATVAL` (column medians
3.0 → 6.7 → 6.8–6.9e-11 on 2024-10-26 frames with DSATVAL 7.5e-11;
2021-12-27 145 s: 2.57e-11 vs 2.84e-11; 2026-02-01 65 s: 5.7e-11 vs
6.4e-11); stars vanish there (a ±80-px search around the predicted
61-vir position found no peak). The 16-s perihelion frames show no
plateau (2022-10-10: the column medians rise to 12.7e-11 and fall
again from vignetting). → **v1.1**: pixels ≥ 0.85 DSATVAL are NaN
before the high-pass (masked fraction 0.6 % at 29 s, 7 % at 65 s,
10 % at 145 s; recorded per frame as `sat_masked`); the ≥ 90 %-finite
aperture gate drops the saturated epochs; saturated cores leave the
calibrator set. Dev re-measured from scratch under v1.1.

**Cause 2 — the unresolved inner-corona structure.** Under v1.1 the D
series are flat over the outer three-quarters of every arc and swing
by ±10–100 units in the last quarter. 42,057 dev arc epochs binned by
distance to the sunward edge: local flux scatter 0.46–0.53 units
(V_eq 9.1 per frame at 5σ) beyond 220 px, 1.8–2.0 at 120–220 px,
6.4–10.8 (p90 29–32) inside 120 px (2.5°, ε ≲ 8°) — the recon's
depth-vs-ε table in aperture form. The per-epoch annulus noise `err`
tracks it (0.41–0.46 → 0.92–1.00 → 2.6–6.0 units; scatter/err
1.1–2.4). → **v1.2** (reduce-stage): a patch epoch with err > 1.5
units is dropped (10,258 of ~60k dev arc records, 382 baseline);
colour calibrators V ≥ 4.5 (V 3–4 stars read 0.30 mag faint, V 4–5
0.06–0.08: bright-star nonlinearity). The v1.1 measurements stand.

## v1.2 dev result (`results/dev_search_v1.json`, `dev_v1.md`)

12 trials, **0 exceedances** (1.33 expected). Control thresholds
S_stack T 3.1–16.8, S_event 8.9–18.8, S_pulse 14.5–28.2 (WISPR dev:
12–43 / 7–33 / 33–70).

| unit | events | z median / MAD | S_stack S / T | note |
|---|---|---|---|---|
| gj-908 S1 | 5 | 0.00 / 2.37 | −0.11 / 3.15 | antipode null |
| gj-908 S2 | 4 | −1.20 / 3.51 | −2.38 / 11.77 | **baseline null on the V 9.0 star** (baseline 3.9–4.5 units = V 8.5–8.7 colour-corrected; z −4.1 … +1.8) |
| 61-vir S1 | 4 | −2.14 / 5.57 | −4.63 / 6.59 | antipode; 2 events `no_baseline`, 3 lost to the edge gates |
| 61-vir S2 | 1 | — | constraint_only | **bright-star class, saturated**: the V 4.7 star's aperture fails the finite/valid gates in 5 of 6 events (`bright_star_saturated`); the one surviving event reads z −12 (baseline 85 units = V 5.2, 0.5 mag faint: nonlinearity) |
| ross-154 S1 | 6 | +4.23 / 4.78 | 10.59 / 16.75 | positive in 5 of 6 events (z +1 … +9.6); controls span −13 … +17: the star-content × field-response systematic that T calibrates (the antipode field is in Gemini near the Galactic anticentre) |

The Gaia-template route (`z_template`) is poorly calibrated on this
substrate (response b 0.37, i.e. the template over-predicts faint-star
aperture flux ~2.7×; z_template −49 … +27 where the star-fixed z is
within ±10) — it stays a secondary annotation, never S.

## Positive controls

- **Vesta 2022-08-20 → 09-14** (tile 1, ε 6°–26°, V 7.5–7.9, r
  0.84 → 0.61 AU, 65–130 s): 118 frames (every 10th), **Δmag −0.05 ±
  0.14** (V_eq − Horizons V with B−V 0.75); one +0.69 outlier on a
  roll-campaign frame (2022-09-03, field rotated). Passes the ≤ 0.2 mag
  rule.
- **Ceres 2024-06-30** (tile 1, ε 5°–15°, V 8.8–9.0): 0 usable frames —
  the whole passage lies in the edge band that v1.1/v1.2 exclude in the
  long-exposure regime (a consistency check of the gates, not a
  recovery).
- Stamp-response injections: ratio 0.960 ± 0.058 (n 4,178).

## Lessons

- `pkill -f <pattern>` matches the calling shell: kill by PID from
  `pgrep -f '^python …'` (the WISPR lesson, re-learned twice).
- The SOAR TAP `v_sc_data_item` inventory and the NRL day tree are the
  same file set; alternating fetches doubles throughput (~140/min).
- Header `DSATVAL` is the saturation reference, but the plateau sits at
  0.90–0.92 of it: mask at 0.85.
