# DASCH crossings — dev stage record (2026-08-26)

Chain: hypotheses v1.0 (FROZEN) + amendments v1.1/v1.2/v1.3 (all
adopted at dev, before any confirmatory contact) +
`threshold_freeze_v1.json`. Numeric results: `dev_v1.json`;
controls: `control_rycnc_v1.json`, `control_iris_v1.json`;
background rate: `teegarden_offwindow_rate_v1.json`.

## Verdict

**0 exceedances on the 6 dev trials** (expected control crossings
0.67) under the final chain; the one intermediate exceedance
(teegarden S_stack, v1.2 scoring) was adjudicated to two vetoed
hits and closed structurally by amendment v1.3. Both positive
controls passed. Machinery validated for all three unit types.

## Units

| Unit | Statistic | S_event / T_event | S_stack / T_stack | Verdict |
|---|---|---|---|---|
| wolf-359 B 2.5 R☉ | track hits (15 win / 23 exp) | 0 / 0 | 0 / 0 | clean null; ring clean |
| teegarden A 0.1 AU | limits-only hits (61 win / 244 exp) | 1 / 1 | 1 / 1 | no exceedance; the surviving S hit is the defect-class 1941 row (adjudication below) |
| van-maanen A 0.1 AU (forced_dev) | lightcurve z (ATLAS routing, 1,365 rows, 97 in-window) | 1.93 / 9.65 | −1.09 / 1.20 | clean null; one control offset carries a defect-tail z = 9.65 — the exchangeable-null construction absorbing photographic tails as designed |

## Positive controls

- **A-chain (RY Cnc, P = 1.092943 d):** recovered through the
  identical querycat → lightcurve → usable-row chain — 1,008 usable
  rows; 85 faint excursions of which **75 % phase-lock at the
  eclipse phase** (uniform null 30 %); depth 0.79 mag; out-of-
  eclipse scatter 0.18 mag (the documented ~0.15 class).
- **B-chain ((7) Iris, V 9.9, plate me00943, 1911, 30 min):**
  recovered as an uncatalogued mag-10.3 detection **0.2″ from the
  timing-extended Horizons trail**, elongated exactly as a trailed
  mover should be (fwhm 6.3 px vs field 3.3, ellip 0.38 — the
  morphology veto rung fires on it correctly) — after amendment
  v1.2, which its first pass forced (below).

## Dev findings → amendments (each before any confirmatory contact)

1. **v1.1 — refcat routing.** The DR7 APASS refcat carries no PM
   for van-maanen or wolf-359 (dummy entries): the APASS
   "lightcurve" of van Maanen holds 9 spurious detections vs the
   true 1,552-detection ATLAS lightcurve. Routing rule adopted
   (PM-match test, APASS-preferred), era-local ±5 yr baselines,
   Stouffer S_stack. Routing: van-maanen→ATLAS, ross-128→APASS,
   gj-1276/wolf-359/teegarden→limits-only.
2. **v1.2 — SUSPECTED_DEFECT demoted.** DASCH's defect classifier
   stamps real single-plate moving transients (the Iris control:
   1/1 killed) — exactly the survey's signal class. The bit is now
   an annotation feeding the cutout-inspection rung, not a silent
   fatal gate.
3. **v1.3 — quiescent-star gate (limits-only).** The teegarden
   exceedance's 1945 hit (mag 16.83 at limit 16.91, 2.2″) is the
   quiescent star marginally extracted: the off-window same-locus
   measurement found the identical signature (mag 16.89/17.32 on
   2 of 12 deep off-window plates; median 17.11 ≈ the expected
   B ≈ 17.3). Spatial rings cannot see this background; the gate
   (hit ≥ 1.0 mag brighter than the measured quiescent level)
   closes it, cost charged.

## Exceedance adjudication (teegarden A, v1.2 scoring)

- **1945 / mc34305** (mag 16.83, 2.2″, lim 16.9): vetoed
  `quiescent_star` — flux-consistent with the measured off-window
  quiescent level; the calibrated-veto evidence is the 12-plate
  off-window rate measurement.
- **1941 / ir05302** (mag 14.64, 5.4″, lim 14.7, 96-min exposure):
  vetoed **defect-class by cutout inspection** — the plate pixels
  at the position show only a ≤ 3.4σ diffuse background bump (the
  21σ field star on the same cutout shows what mag 14.6 looks
  like); non-PSF morphology (fwhm 7.3 px, ellip 0.29), at-limit
  photometry, SUSPECTED_DEFECT annotation concurring, and no
  recurrence in the unit's 60 other covered windows.

## Substrate facts fixed at dev

- platephot has **no row cap** (counts 0–122 observed) and serves a
  ~13′ subregion; ring offsets ±60–120″ ride the same pull.
- The v1.3 quiescent measurement doubles as a **DASCH-B photometry
  check**: teegarden's measured median 17.11 vs the expected
  B ≈ 17.3 from V = 15.1, B−V ≈ 2.2.
- One lightcurve pseudo-window control offset reached z = 9.65
  (defect-tail row inside a control window): photographic A-chain
  nulls are heavy-tailed; the max-of-8-controls threshold absorbs
  them at depth cost, exactly the v2 heavy-tail behaviour.
