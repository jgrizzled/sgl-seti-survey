# ZTF beam-crossings survey (Pipeline B) — report

First Pipeline-B (§3.5) survey: a search of ZTF public single-epoch
data for optical emission during Earth's crossings of hypothesized
Sun–star relay beam axes. Executed 2026-08-23 → 24. **No candidates.**

## 1. Construction and provenance chain

1. **Universal crossing list** `crossings/universal_v1`
   (`xng-a09e2db7681d`): every b(t) minimum for 88 endpoints × 2 link
   directions, Earth-center, 1980→2028 (16,586 events; `sglseti`
   `sun_star_axis_v1` / `tusay2022_eq5_7_v1`).
2. **Hypothesis freeze v1.0** `surveys/ztf-crossings/hypotheses.md`:
   two ZTF-observable channels — A: uplink interception at the star
   near opposition (radii 0.1 / 1.0 AU); B: downlink pre-lens
   interception at the star's antipode near the anti-solar point
   (radii 1.2 R☉ / 2.5 R☉ / 0.1 AU; relay track parameterized by
   z ∈ 550–10,000 AU). Sunward-arriving combinations declared out of
   scope. 400–900 nm; duty d = 1 primary (channel B is window-locked
   by geometry alone).
3. **Coverage intersection** (`results/coverage_v1.*`): era MJD
   58178–61275. A: 402/720 events covered (all 29 uncovered targets
   δ < −30°). B: grazing radii epoch-starved as the cadence predicts
   (5/43 covered, gj-1276 + Teegarden only); 0.1 AU rung 33/60.
4. **Channel-A saturation cut** (`results/saturation_cut_v1.*`):
   excluded < 13.0 mag; searchable populations g 27 / r 20 / i 14
   targets — all M dwarfs, white dwarfs, and optically-dark late-T/Y.
5. **Threshold freeze v1.0 + amendments**
   (`configs/threshold_freeze_v1*.json`, `thresholds.md`): diff-image
   forced photometry; A ≥2-window stack, B max over the nested z
   family; 8 controls (A temporal pseudo-windows, B spatial ring
   20/30/40″); exceedance S > max(T, 0); dev/confirmatory split, seed
   20260822, grazing rung dev-protected. Amendments from the dev run:
   v1.1 (parallax-factor systematics template, same-rung pseudo-window
   exclusion, ratio retired), v1.2 (empirical variance rescale
   k = median(r²/v)/0.4549, the v1 confusion-not-noise precedent).
6. **Dev search** (`results/dev_search_v1.md` + v11/v12 json): channel
   B validated clean; channel A's two frozen-construction failures
   found, amended, and re-verified; the 1.0-AU rung shown to be
   constraint-only by geometry (windows ≈ observing season).
7. **Confirmatory search** (`results/confirmatory_v1.*`): 8,818/9,590
   epochs, 27 A rows + 59 B rows, per-epoch propagated star positions.
8. **Completeness** (`results/completeness_v1.json`): exact
   stamp-response injections (Moffat β=3 at SEEING), ≥100 (200)
   bootstrap draws per event-window cell, chord profiles, magnitude
   grid 14–22 (values at 22.0 are grid-censored lower limits).

## 2. Results

**36 searchable units (3 A + 33 B): 3 exceedances vs 4.0 expected
control crossings — the frozen null rate exactly. All three
adjudicated as control-crossing statistics** (amplitudes 0.6–2.2σ, no
amplitude-consistent recurrence; gj-908 g ×1, gj-1276 g ×2; two of six
gj-1276 events exceeding has null probability ≈ 0.11). Every
constraint-only and single-epoch row is reported with its statistic;
772 archive-404 epochs and all anomalies (one 5/8-control row, two
gj-1276 r bad-epoch outliers) are counted in the record.

### 90 %-recovery depths (AB mag; band-monochromatic line hypothesis)

| channel / rung | population | median m90 |
|---|---|---|
| B 0.1 AU (searchable units) | 31 | 21.8 (14 grid-censored ≥22) |
| B 0.1 AU (single-epoch) | 8 | 21.2 |
| B 2.5 R☉ grazing (single-epoch) | 9 | 21.6 |
| B 1.2 R☉ grazing (single-epoch) | 2 | 21.9 |
| A searchable (contrast-limited) | 3 | 16.6 (g) |
| A constraint-only (with template) | 18 | 14.4–18.4 |
| A wise-0855 (empty field) | 2 | ≥22 (grid-censored) |

Channel-B depths reach the ~1–2σ exceedance floor because false-alarm
control lives in the 1/9 control-crossing budget and veto ladder, not
in a 5σ cut. Channel-A blended depths are contrast-limited 3–7 mag
above the sky limit, as the k factors (10²–10⁴) dictate.

## 3. Physical interpretation (532 nm line in g unless noted)

**Downlink pre-lens (the survey's distinctive constraint).** A relay at
550 AU beaming through the 2.5 R☉ solar-grazing cone (gain ≈ 9×10⁹)
toward the target star would have been detected at 90 % confidence in
the covered windows at transmitter power ≳ **~130 W** (wide-rung
median m90 = 21.8 → EIRP 1.1×10¹² W at 549 AU; grazing single-epoch
median 21.6 → ~150 W). The constraint applies only during the covered
windows: 33 wide-rung events on 5 targets and 5 grazing windows on 2
targets (gj-1276, Teegarden) — a per-event statement, not a continuous
monitor.

**Uplink interception.** A transmitter in the target system beaming at
our Sun through a 10-m diffraction-limited aperture is excluded at
90 % during covered windows above ~**56 kW** (Teegarden, m90 = 17.0,
constraint-only), ~66 kW (Wolf 359) to ~400 kW (gj-1276) for
searchable units; 1-m-class transmitters scale ×100. The deepest
uplink limit is the optically-dark wise-0855 field (m90 ≥ 22):
10-m-class ≳ 200 W — an empty-field limit at 2.28 pc.

**Not constrained:** 1064/1550 nm links (outside ZTF response);
transmitters scheduled to avoid Earth crossings; pulse periods between
30 s and the window length; δ < −30° endpoints; the anti-target
downlink outside covered windows; post-lens downlink light (arrives
sunward).

## 4. Lessons for the next crossings adapters

1. The crossing geometry phase-locks to the sidereal year; for
   wide-beam rungs the window ≈ the observing season, so temporal
   controls fail by construction and those rungs are constraint-only —
   design narrow-rung-first for optical surveys.
2. Blended-channel (on-star) searches are systematics-dominated:
   PM-dipole vs multi-year references, nonlinear seeing coupling, and
   variance floors 10²–10⁴ above background. The parallax-factor
   template + empirical k rescale (v1.1/v1.2) is the reusable recipe.
3. The antipode channel is the workhorse: empty fields, geometry-forced
   transience, spatial ring controls immune to the phase-lock, and
   ~watt-to-hectowatt relay power limits from modest optical depth.
4. Exceedance budgeting worked as designed: 3 observed vs 4.0 expected,
   every one adjudicated, none promoted.

## 5. Products

`surveys/ztf-crossings/{hypotheses.md, thresholds.md, configs/*,
results/*}`, universal list under `crossings/universal_v1/`, cutouts +
snapshots under `runs/ztf-crossings/` (dev 202 MB + conf 770 MB),
scripts under `surveys/ztf-crossings/scripts/`.
