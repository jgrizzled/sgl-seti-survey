# ATLAS + ASAS-SN beam-crossings survey (Pipeline B) — report

Tenth Pipeline-B (§3.5) survey: a search of ATLAS difference-image
forced photometry (four 0.5-m units, o/c bands, nightly quads,
2015-07 → 2026-09) at the Earth-crossing windows of hypothesized
Sun–star relay beam axes, with ASAS-SN Sky Patrol v2 as a window
coverage-fraction ledger. The survey was adopted for the **window
coverage-fraction / recurrence cell** — nightly all-sky cadence over
every in-era window of the narrow-rung targets, including the
van-maanen b = 0.24 R☉ photosphere-grazing family no previous survey
could search. Executed 2026-09-03/04 on the new server (first survey
run there). **No candidates: 27 blind confirmatory trials, 5
exceedances (3.2 expected control crossings), all adjudicated.**

## 1. Construction and provenance chain

1. **Universal crossing list** `crossings/universal_v1`
   (`xng-a09e2db7681d`), Earth-center, 1980→2028, both link
   directions; site-vs-geocenter ≤ R⊕ = 0.8 % of the tightest rung
   (declared budget term).
2. **Reachability recon** (2026-08-25,
   `surveys/atlas-asassn-crossings/notes/atlas_asassn_recon_2026-08-25.md`):
   ATLAS forced photometry is a token REST API, serial per account,
   arbitrary positions, all-sky (southern units from 2022); ASAS-SN
   Sky Patrol v2 serves catalogued sources only (no channel-B point
   photometry) over plain HTTP; v1 is reCaptcha-gated (manual only).
3. **Hypothesis freeze v1.0** (`hypotheses.md`, D1–D8 approved
   2026-08-25; coverage deferred to the server). **Pre-data amendments
   at coverage start (§12, 2026-09-03):** A1 — the freeze's "22
   semiannual windows per target" counted the sunward axis crossing
   of each year (channel position at solar elongation 0.1–0.5°),
   which the universal declaration excludes; the searchable
   population is **11 annual windows per target per channel** (B
   1.2 R☉ 44 / 2.5 R☉ 55 / 0.1 AU 77; A 0.1 AU 77), verified by
   elongation at every task position; A2 — ATLAS era end MJD 61286,
   ASAS-SN v2 ceiling re-measured unchanged at 2025-06-16
   (reprocessing lag); A3 — no radec-list batching exists (schema
   snapshot), one windowed task ≈ 2.2 min.
4. **Task list and coverage** (`results/task_list_v1.ecsv`,
   `coverage_v1_*`): 308 windowed (±110 d) difference-flux tasks —
   231 B (t_ca position + D3 mini-track ingress/egress at z = 550,
   37.5″ apart) + 77 A (the star has no z-track; the mini-track
   collapses by identity); drained in 7.4 h, 0 failures, every result
   sidecar-snapshotted (`runs/atlas-asassn-crossings/lc/`). FAQ mask
   keeps 93 %. **The grazing rungs are attrition-dominated**: 0.3–1.3 d
   windows against nightly-at-best sampling with weather and sun-gap
   losses give B 1.2 R☉ 8/44 windows covered (van-maanen 4/11), 2.5 R☉
   25/55, against 0.1 AU 74/77 (B) and 68/77 (A); grazing stacks hold
   14–30 exposures, not the pre-declared 100–180. ASAS-SN v2
   field-level epoch ledger (`asassn_ledger_v1_*`, union of quality-G
   image ids of catalogued sources within 3′; coverage records only)
   shows the same shape at ~2 d cadence.
5. **Threshold freeze v1.0** (`configs/threshold_freeze_v1.json`,
   hash `3fe21df0…`, hypotheses §13; seed 20260904; user-approved
   recommendations 2026-09-04): unit = (target, channel, rung, o
   band), c annotation only; statistics S_stack (pooled recurrence
   chord, primary cell) / S_event (max window chord) / S_pulse (max
   in-window exposure), all on D5-detrended (30-d off-window running
   median), k-rescaled, response-weighted residuals with best-R
   mini-track position per epoch and the nested-z max; controls = 8
   temporal pseudo-windows (±23/47/71/97 d) at the same fixed
   position, **mirror-gated validity** (an offset counts only where it
   carries data like the real window; pseudo-stacks need ≥ 3 events),
   T = max valid control, per-trial crossing probability
   1/(n_valid + 1) from measured validity; gates: S_pulse/S_event at
   ≥ 1 covered window, S_stack at ≥ 3, **≥ 4 valid controls for a
   trial**. D2 saturation on 20-d reduced-mode series: gj-1276 ok,
   teegarden marginal (o 12.59), van-maanen (12.30), wolf-359,
   ross-128, ross-154, gj-908 excluded. D8 dev = wolf-359 (D1 forced)
   + gj-908. Population: **14 searched units / 37 trials** (dev 10,
   confirmatory 27), 4.55 expected control crossings; **S_stack is
   constraint-only at every grazing rung** (≤ 3 valid pseudo-stacks);
   grazing S_pulse/S_event searched for van-maanen (1.2 + 2.5 R☉),
   teegarden and wolf-359 (2.5 R☉).
6. **Dev stage** (`results/dev_search_v1.json`,
   `dev_adjudication_v1.json`): 10 trials, 2 exceedances on one epoch
   (wolf-359 2.5 R☉ — a single 3.5σ Sutherland exposure in a
   degraded-sky quad, quad-inconsistent, SkyBoT clean, no recurrence:
   retained, non-promotable); a chi/N-45 single-frame artefact in a
   pseudo-window sets one S_pulse threshold at 15.8 (frozen max rule;
   the trial stays valid but insensitive). KS of standardized
   baselines p 0.24–0.97, k 1.0–1.5, engine validity counts match the
   freeze; no amendments.
7. **Blind confirmatory** (`results/confirmatory_search_v1.json`,
   `confirmatory_adjudication_v1.json`; §2).
8. **Completeness + positive control** (`results/completeness_v1.json`,
   `power_limits_v1.json`, `mpc_control_v1.json`; §2–3): D6
   response-model injections into the real baseline series against
   the frozen thresholds; (15000) CCD through the identical chain vs
   JPL Horizons.

## 2. Results

**0 candidates.** Blind confirmatory run over 27 trials (10 units):
**5 exceedances vs 3.23 expected control crossings** (Poisson
P(≥ 5) ≈ 0.23), on three distinct events, none surviving the frozen
veto ladder:

- **X1 — van-maanen B 1.2 R☉ + 2.5 R☉, S_event 2.78** (T 1.97 / 2.64;
  evt-fc686a, 2017-04-03, b = 0.24 R☉). Three Haleakala exposures at
  dt +0.27…+0.30 d read 11, 66, 43 ± 30 µJy — a quad-consistent
  +40 µJy (2.8σ) chord; the mini-track positions 37.5″ away are
  nominal noise; SkyBoT clean (Jupiter's system 5.4° away); pre-dates
  both template steps; **no recurrence** in the other three (1.2 R☉)
  / four (2.5 R☉) covered windows (S_w −1.33, −0.84, −0.48). Disposition:
  **retained-ambiguous, non-promotable** — the one 2.8σ chord in four
  covered photosphere-grazing windows; promotion requires recurrence
  at a second covered window (the frozen rule). One epoch set counted
  in two nested rungs.
- **X2 — gj-1276 B 0.1 AU, S_event 5.09 / S_pulse 3.04** (T 3.03 /
  2.93; evt-9c6a08, 2025-03-31, at z = 10000 where all 28 in-window
  epochs enter). In-window nightly means are consistent with zero
  (−11…+10 µJy ± 5–9); the excess is manufactured by a −10 to
  −25 µJy off-window baseline depression at dt −10…−25 d (and −55 at
  +15 d) that the 30-d running median cannot remove, so each in-window
  epoch carries a +0.8σ residual. The unit's other ten windows span
  S_w −3.4…+1.5 — the same baseline-structure class, wider than the
  8 controls resolve. S_pulse: a +31 ± 15 exposure whose quad reads
  −8, +22, +31, −3 (rung 2 fail). SkyBoT clean at peak and t_ca.
  Disposition: **adjudicated difference-image baseline systematic**,
  retained, non-promotable.
- **X3 — teegarden A 0.1 AU, S_stack 11.8** (T 6.9). The
  "difference flux" at the star is the star's full o-band flux
  (34–59 mJy; m 12.59): Teegarden's 5.1″/yr proper motion has carried
  it > 1 FWHM off its template position, so the difference image
  contains the star (the PM dipole of learnings §8); chi/N 60–5600,
  k 79 → 2.9 × 10⁷ across events, controls −14…+7. Disposition:
  **adjudicated channel-A blended systematic** (marginal saturation +
  PM dipole; the freeze carried the variance rescale but not the
  ZTF v1.1 parallax/PM systematics template); the unit is reported
  systematics-limited with its measured null.

Dev + confirmatory together: 37 trials, 7 exceedances vs 4.55
expected, 0 promotable.

**Constraint-only ledger.** S_stack at all grazing rungs; B 1.2 R☉
teegarden, gj-1276 and B 2.5 R☉ gj-1276, ross-128 (1–3 valid
controls); A 1.0 AU (no ATLAS pull — window ≈ season; deferred to a
ledger addendum); the five saturation-excluded A units.

### 90 %-recovery depths (`results/completeness_v1.json`; o-band AB, D6 response-model injections into the real baseline series, 100 draws per magnitude, seed 20260904, z = 550 track geometry, frozen thresholds)

| unit | S_event m90 (T) | S_pulse m90 (T) | S_stack m90 (T) |
| --- | --- | --- | --- |
| B 1.2 R☉ van-maanen (conf) | 20.55 (1.97) | 19.78 (3.09) | constraint-only |
| B 1.2 R☉ wolf-359 (dev) | 20.39 (−0.66) | 19.94 (0.88) | constraint-only |
| B 2.5 R☉ teegarden (conf) | 20.31 (1.94) | 19.72 (2.33) | constraint-only |
| B 2.5 R☉ van-maanen (conf) | 19.95 (2.64) | 19.50 (3.20) | constraint-only |
| B 2.5 R☉ wolf-359 (dev) | 20.75 (1.45) | 20.17 (2.14) | constraint-only |
| B 0.1 AU gj-1276 (conf) | 20.23 (3.03) | 19.90 (2.93) | 20.35 (2.80) |
| B 0.1 AU gj-908 (dev) | 20.07 (3.03) | 19.73 (2.73) | 20.50 (1.49) |
| B 0.1 AU ross-128 (conf) | — (412, artefact-set) | — (188, artefact-set) | 19.79 (4.74) |
| B 0.1 AU ross-154 (conf) | 20.05 (3.28) | 19.30 (3.32) | 20.23 (2.94) |
| B 0.1 AU teegarden (conf) | 19.71 (4.22) | 20.15 (2.36) | 20.19 (3.12) |
| B 0.1 AU van-maanen (conf) | 20.12 (3.37) | 17.12 (22.4, artefact-set) | 20.50 (2.35) |
| B 0.1 AU wolf-359 (dev) | 19.69 (4.14) | 17.57 (15.8, artefact-set) | 20.02 (4.75) |
| A 0.1 AU gj-1276 (conf) | — (99.9) | — (74.6) | 16.55 (52.4) |
| A 0.1 AU teegarden (conf) | — (18.5) | — (65.1) | — (6.9; recovery ≤ 0.5 at 16.0) |

Reading: the B units recover to o ≈ 19.7–20.75 per window and 19.8–20.5
recurrence-stacked at 0.1 AU — a nightly quad (≈ 19.9 at 5σ) read at a
2–4σ calibrated threshold — i.e. the pre-declared per-window depth,
while the recurrence stack the survey was adopted for adds only
~0.3 mag at 0.1 AU (10–11 windows of 14–17 epochs each, against
systematics-set T of 1.5–4.7) and is not calibratable at the grazing
rungs. Five trials are **artefact-set**: their max-rule threshold is a
single-frame outlier in one pseudo-window (T 15.8–412), so they remain
valid but insensitive trials (m90 17.1–17.6 or none) — the price of
the frozen max rule without a chi/N mask term. Both channel-A units
are systematics-limited (T 7–100 in k-rescaled units; injections
recover at ≤ 50 % even at o = 16): the on-star cell delivers no
constraint here. D6 positive control (15000) CCD: scale gate PASS
(apparition medians −0.05 / −0.10 mag after one colour term, slope
−0.03 mag/mag; nightly scatter 0.32 mag is photon + rotation at
S/N 5–8, not calibration); ±0.1 mag declared systematic on every
depth.

## 3. Physical interpretation (line through the o band, 560–820 nm; monochromatic flux via the 260 nm effective width)

**Downlink pre-lens — the grazing cones (the survey's new cell).** A
relay at z = 550 AU beaming through the Sun-grazing cone of radius r
toward the target (P = F_line · π r² · (d/z)²; the 690 nm o-band
centre stands in for any line in 560–820 nm, 752 nm included) would
have been detected at 90 % in a covered window above:

| target (b) | 1.2 R☉ cone | 2.5 R☉ cone | pulse cell (≥ 1 exposure, 30 s) |
| --- | --- | --- | --- |
| van-maanen (0.24 R☉) | **78 W** (4 windows) | 590 W (5 windows) | 160 W / 890 W |
| wolf-359 (0.69 R☉, dev) | 90 W (1 window) | 280 W (4) | 140 W / 480 W |
| teegarden (0.96 R☉) | constraint-only | 420 W (4) | — / 730 W |

These are **per-window, single-transmission** limits (S_event), not
recurrence-stacked: the van-maanen 0.24 R☉ photosphere-grazing cone —
unconstrained by every previous survey — is now excluded above ~80 W
during four of its eleven 2015–2026 windows (2.8σ retained-ambiguous
chord in one of them, §2), and the 2.5 R☉ cones above 280–590 W in 4–5
windows each. For comparison ZTF reached ~130 W through the 2.5 R☉
cone recurrence-stacked (g, 21.8); ATLAS's nightly cadence buys the
first *per-window* grazing coverage of the deepest family, at
~4× ZTF's per-window flux. The 1.2 R☉ numbers are lower than the
2.5 R☉ ones purely by cone area (r²).

**0.1 AU cone (B).** 26–55 kW per window / 26–41 kW recurrence-stacked
(S_stack m90 19.8–20.5) across all seven targets, 9–11 covered windows
each; ross-128's S_event/S_pulse and the van-maanen/wolf-359 S_pulse
trials are artefact-set (≥ 0.4–0.6 MW). Pulse cell: 36–79 kW per
≥ 30 s exposure.

**Uplink interception (A, 10-m diffraction-limited beam at the
target).** gj-1276 S_stack 2.2 MW (o 16.55) — systematics-limited;
teegarden none. The on-star cell contributes no useful constraint in
ATLAS difference imaging for these proper motions (§4).

**Declared unconstrained.** 1064/1550 nm (outside o/c); sub-exposure
pulses beyond ×d scaling; schedules avoiding Earth-crossing windows;
recurrence at the grazing rungs (S_stack constraint-only there);
single-window transmissions below per-window depth at the 82 % / 55 %
of grazing windows ATLAS did not cover.

## 4. Lessons for the next crossings adapters

- **Count windows with the axis-side filter.** Link direction alone
  double-counts: every target has two axis crossings per year and one
  is sunward. Caught pre-data by the elongation check in the task
  builder; `era_scope` now carries both counts.
- **A window shorter than the revisit interval is covered by luck.**
  ATLAS's nightly cadence delivers 0.1 AU windows (10–12 d) at ~95 %
  but 0.3–1.3 d grazing windows at 18–45 %; the recurrence stack the
  survey was adopted for is not calibratable at the grazing rungs
  (≤ 3 valid pseudo-stacks). The coverage-fraction number is the
  deliverable; the stack needs a higher-cadence or wider-window
  substrate (or many more years).
- **Mirror-gated control validity + 1/(n+1) accounting** kept the
  budget honest where the nominal 8 controls do not exist (grazing
  rungs: 4–7 valid).
- **Single-exposure and single-frame classes** dominate the
  exceedances: the FAQ mask has no chi/N term (a chi/N-45 cosmic-ray
  class artefact passed), and S_event on a one-exposure window is
  just that exposure's S/N. The intra-quad rung of the ladder handles
  both; a v2 strict mask should add chi/N < 5, and S_event could
  require ≥ 2 in-window epochs (PTF's rule).
- **A window-scale baseline test belongs in the ladder**: the D5
  30-d median cannot remove ~20-d baseline depressions, which
  manufacture positive chords over 10-d windows (X2). Test the
  in-window mean against zero, not only against the local median.
- **Channel A on a high-PM star in difference imaging is the PM
  dipole**, not a blended-star residual: at 5″/yr against a
  multi-year template the "difference flux" is the whole star. The
  ZTF v1.1 parallax/PM systematics template (omitted from this
  freeze) is mandatory for any future ATLAS channel-A unit; better,
  request reduced-mode photometry with a per-epoch PM-propagated
  position and fit the systematics there.
- **ASAS-SN v2 catalog matches must be PM-matched**: nearest-source
  matches sit 15–36″ off the high-PM stars (epoch offset).
- Ops: the `comment`-keyed resumable drain + sidecar done-list ran 308
  tasks unattended (median 54 s server time, max 8 min); the
  per-account serial queue is the only throughput limit.

## 5. Products

`surveys/atlas-asassn-crossings/`: `hypotheses.md` (v1.0 + §12
amendments + §13 freeze record), `configs/threshold_freeze_v1.json`,
`results/` (task_list_v1, coverage_v1, asassn_ledger_v1,
asassn_currency_2026-09-03, era_scope_v1, saturation_cut_v1,
dev_search_v1 + dev_adjudication_v1, confirmatory_search_v1 +
confirmatory_adjudication_v1, completeness_v1, power_limits_v1,
mpc_control_v1), `scripts/` (atlas_api, build_task_list, atlas_drain,
coverage_ledger, asassn_currency, asassn_ledger, aux_pulls,
saturation_cut, freeze_thresholds, search_core, run_search,
completeness_v1, mpc_control). `runs/atlas-asassn-crossings/`: 308
light-curve results + sidecars (`lc/`), auxiliary pulls (`aux/`),
ASAS-SN snapshots (`asassn/`), API docs snapshot, SkyBoT census
files, drain log.
