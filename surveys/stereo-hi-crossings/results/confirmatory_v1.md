# confirmatory search v1 — unit table

Construction: threshold_freeze_v1.0; frames usable 38830, unusable {'n_images_25': 136, 'nmissing_1.0': 89, 'nmissing_2.0': 8, 'nmissing_4.0': 8, 'nmissing_7.0': 3, 'nmissing_9.0': 3, 'nmissing_3.0': 7, 'nmissing_5.0': 11, 'nmissing_14.0': 1, 'astrometry_or_zp': 17, 'nmissing_10.0': 1, 'nmissing_8.0': 1, 'nmissing_6.0': 2, 'exptime_1499.94': 2, 'exptime_1499.93': 1, 'exptime_1119.94': 1, 'n_images_15': 1329}; median matched stars 700.0, astrometric rms 0.436 px, ZP0 11.068, ZP MAD 0.054 mag. Epoch vetoes: {'near2deg_venus': 1265, 'prox_mercury': 83, 'recon_exclusion': 6, 'near2deg_earth': 854, 'near2deg_moon': 3, 'prox_earth': 21, 'near2deg_jupiter': 265, 'prox_venus': 19}. Stamp response 0.997 ± 0.002 (n 5390); differential flat 0.9990.

**Trials 27, exceedances 3, expected control crossings 3.0.**

| unit | status | n_ev | z med / MAD | S_stack S / T | S_event S / T | S_pulse S / T |
|---|---|---|---|---|---|---|
| gj-1276|S1_0.1AU | searched | 20 | -0.02 / 1.60 | 0.11 / 9.66 | 4.06 / 7.75 | 5.73 / 14.65 |
| gj-1276|S2_0.1AU | searched | 18 | -0.25 / 2.36 | -0.21 / 17.64 | 8.22 / 7.75 **EXC** | 15.71 / 11.16 **EXC** |
| ross-128|S1_0.1AU | constraint_only | 0 | | | | |
| ross-128|S2_0.1AU | searched | 19 | -1.87 / 2.31 | -9.16 / 12.41 | 3.43 / 37.20 | 8.43 / 73.08 |
| teegarden|S1_0.1AU | searched | 17 | -1.80 / 2.12 | -6.66 / 3.72 | 3.58 / 5.78 | 19.57 / 7.71 **EXC** |
| teegarden|S2_0.1AU | searched | 18 | +0.38 / 1.73 | 1.54 / 10.94 | 6.45 / 10.91 | 10.59 / 29.22 |
| van-maanen|S1_0.1AU | searched | 18 | +0.32 / 2.12 | 1.51 / 21.08 | 4.81 / 15.73 | 9.68 / 106.35 |
| van-maanen|S2_0.1AU | searched | 18 | +0.39 / 3.38 | -1.44 / 9.83 | 3.75 / 12.54 | 3.94 / 18.28 |
| wolf-359|S1_0.1AU | searched | 18 | -0.19 / 1.67 | -3.27 / 6.85 | 5.62 / 6.51 | 5.53 / 13.31 |
| wolf-359|S2_0.1AU | searched | 20 | +0.13 / 2.25 | 2.88 / 13.94 | 6.82 / 14.04 | 5.68 / 14.70 |

## Completeness (m90, V-equivalent) and power

| unit | S_stack m90 / P | S_event m90 / P | S_pulse m90 / P |
|---|---|---|---|
| gj-1276|S1_0.1AU | 12.88 / 1.17e+07 W | 12.45 / 1.74e+07 W | 7.95 / 1.09e+09 W |
| gj-1276|S2_0.1AU | 11.91 / 6.09e+07 W | — | 7.84 / 2.57e+09 W |
| ross-128|S2_0.1AU | 11.91 / 9.52e+06 W | 10.70 / 2.89e+07 W | 6.04 / 2.13e+09 W |
| teegarden|S1_0.1AU | 12.88 / 1.17e+07 W | 13.20 / 8.66e+06 W | 7.67 / 1.42e+09 W |
| teegarden|S2_0.1AU | 12.45 / 7.44e+06 W | 11.52 / 1.75e+07 W | 7.42 / 7.69e+08 W |
| van-maanen|S1_0.1AU | 11.91 / 2.86e+07 W | 11.52 / 4.07e+07 W | 5.59 / 9.62e+09 W |
| van-maanen|S2_0.1AU | 11.91 / 1.56e+07 W | 11.52 / 2.22e+07 W | 7.34 / 1.04e+09 W |
| wolf-359|S1_0.1AU | 12.45 / 1.73e+07 W | 12.88 / 1.17e+07 W | 7.85 / 1.20e+09 W |
| wolf-359|S2_0.1AU | 12.45 / 2.94e+06 W | 11.52 / 6.91e+06 W | 7.71 / 2.33e+08 W |

## Adjudication (frozen §6 ladder; 2026-09-04)

27 searched trials over 9 units (ross-128 S1 constraint-only: its
antipode patch is Tycho-masked — a VT 7.2 star 1.7 px away, the same
permanently blended fixed sky point LASCO found); **3 exceedances vs
3.0 expected control crossings; every S_stack null** (max S_stack 2.9
vs T ≥ 3.7). No planet or bright asteroid (Horizons `@-234`, H < 7
list) within 0.5° of any exceedance patch at its epoch.

- **gj-1276 S2 S_pulse 15.71 vs T 11.16** — 2009-01-14 21:29–23:29
  (`evt-9dbcb8c29217`, 4 consecutive frames: D +11.7, +17.6, +8.7,
  +3.0). Pixels: a broad diagonal band of +1–2 DN/s, ~10 px wide,
  sweeping across the 13 × 13 stamp at ~2–3 px per frame over 2.5 h
  — an **extended moving front (CME/streamer class)**, not a point
  source; the same-frame control median stays at +0.26. Persistence
  passed, point-source morphology failed. Vetoed, coronal-transient
  class.
- **teegarden S1 S_pulse 19.57 vs T 7.71** — 2010-08-14 15:29–16:09
  (`evt-7ae71bb63446`, the last two arc frames: D +6.4, +10.0, back to
  −2.6 at 16:49). Pixels: the entire stamp lifts uniformly by ~0.8
  DN/s for two frames around a static field star at (+4, +1) px and
  returns — an **extended transient front**, not a point source.
  Vetoed, coronal-transient class.
- **gj-1276 S2 S_event 8.22 vs T 7.75** — 2012-10-22/24 arc
  (`evt-333d548b27c7`, 54 epochs, mean D +0.39 units ≈ V_eq 12.1,
  split halves 5.85 / 5.78). Pixels at four epochs across the arc:
  **no point source** (gj-1276 itself, V 16, is invisible; the stamp
  residuals are ±0.1–0.3 DN/s). The plateau is the star-fixed
  differential's baseline offset: this patch's own baseline median is
  −0.27 units while its arc level is ~0, i.e. a regional L2 residual-
  background difference between the ε 6°–12° and 4.4°–6° regions,
  not shared by the flank controls. Non-recurrent: the other 17
  windows of the same patch give S_stack −0.21. **Adjudicated
  background-offset systematic, retained, non-promotable** (the same
  class as the dev gj-908 S1 2022-10 case).

**Result: 0 candidates.** Completeness for gj-1276 S2 S_event is
degenerate (the null exceeds T) and is reported as not measurable.
