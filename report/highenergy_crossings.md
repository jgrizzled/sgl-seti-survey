# High-energy (Fermi-LAT + Swift/BAT) beam-crossings survey (Pipeline B) — report

**Result: 0 candidates.** Blind confirmatory run over 29 units / 145
frozen trials (136 calibrated): **9 exceedances, all on the three B
van-maanen units whose lane-L aperture contains the blazar 3C 279, all
vetoed `known_source` at step 1 of the frozen ladder and all already
`constraint_only` at the threshold freeze**; the 136 calibrated trials
show 0 exceedances vs 0.05 expected, with 10,778 lane-L photons
observed against 10,517 expected. The survey opens the programme's
high-energy pulse/burst cell — 0.1–300 GeV photon-event coincidence
with the crossing windows of the seven deep-family targets, both
channels, three rungs, 2008–2026 — with recurrence-stacked 90 %
limits of ~1.5 × 10⁻⁷ ph cm⁻² s⁻¹ (100 MeV–1 GeV) through the
1.2 R☉ cone, i.e. **~0.15 MW-class through the grazing cones, ~6 MW
through 0.1 AU**, and a Swift/BAT 15–150 keV constraint-only arm at
~30 mCrab per grazing window (1.6–6 MW). The eleventh crossings
survey; the first to use a γ-ray archive.

## 1. Construction and provenance chain

Plan §5.26. Substrate: Fermi-LAT P8R3 SOURCE-class photons from
era-long data-server pulls at the 14 channel positions (3° radius,
0.1–300 GeV; 16 queries incl. two controls, ~250 MB, sha-pinned
manifests under `runs/highenergy-crossings/v1/photons/`), gated by
the 944 weekly spacecraft files (θ ≤ 60°, DATA_QUAL > 0,
LAT_CONFIG = 1, Moon/Sun > 8°) which also define the livetime and the
CALDB (P8R3_SOURCE_V2) effective-area exposure of every window and
pseudo-window. Chain, all 2026-09-07 on the server: reachability recon
(`surveys/highenergy-crossings/notes/highenergy_recon_2026-09-07.md`)
→ user approval of seven recommendations → hypothesis freeze v1.0
(`hypotheses.md` D1–D11; pre-freeze contact declared, remedy D1: the
van-maanen B 2020-04-02 event `forced_dev`) → **amendment v1.1**
(Jeffreys floor on the pseudo-window rate; A gj-1276 1.2 R☉ demoted
after its statistics were printed in a machinery test → 145 trials,
T = 3.451) → coverage (`results/coverage_v1.*`: 32 units, 580
windows, 518 searchable at ≥ 1 ks gated livetime) → threshold freeze
(`thresholds.md`, `results/thresholds_v1.json`, sha in
`results/threshold_freeze_v1_sha.txt`; 55,474 pseudo-windows) →
**amendment v1.2** (per-trial ensemble exceedance gate; 9 B van-maanen
trials `constraint_only`: 3C 279 at 1.8° from the antipode) → dev
(`results/dev_v1.*`: gj-908 A/B, the demoted unit, the D1 window, two
positive controls) → blind confirmatory (`results/confirmatory_v1.*`,
per-unit `results/units_v1/`) → completeness
(`results/completeness_v1.*`, seed 20260907) → the Swift/BAT
constraint-only arm (`results/bat_v1.*`; HEASoft 6.24 `batsurvey` in
the `chbrandt/heasoft` container with a local Swift BAT CALDB, 238
pointings, 8 workers) → the pre-registered Swift/XRT look
(`results/xrt_look_v1.*`) → the Pipeline-A corridor catalogue screen
(`corridor_screen.md`, `results/corridor_screen_v1.*`).

Statistics (per unit = target × channel × rung): `S_stack` and
`S_event` in lane L (100 MeV–1 GeV, r = 2.0°) and lane H (1–300 GeV,
r = 0.7°) — Poisson upper tails of the in-window count against the
exposure-matched pseudo-window rate (±60 d grazing, ±120 d 0.1 AU) —
and `S_burst` (max over 1-ks boxcars in the union band, r = 1.0°,
look-elsewhere inside). FWER α = 0.05 over 145 trials → T = 3.451.

## 2. Results

**Geometry that shaped the survey.** Both channels sit at solar
elongation 180° in every window. XMM-Newton (solar aspect 70–110°) and
eROSITA (scan at 90° elongation) therefore never observe a crossing —
Gate I of §5.3 now has four members — and Chandra has no in-window
observation in 27 years; Swift XRT has one (below). The monitors carry
the cell: LAT and BAT.

**Coverage** (`results/coverage_v1.md`). 580 windows in the LAT era:
1.2 R☉ 119 searchable / 27 uncovered, 2.5 R☉ 155 / 27, 0.1 AU 244 / 8
(uncovered = < 1 ks gated livetime, all post-2018-03, the
solar-array-anomaly survey profile). Confirmatory units: 479
searchable windows, 29.8 Ms gated livetime, 9.6 × 10⁵ m² s lane-L
exposure. The Moon/Sun 8° exclusion removed 8.0 Ms of 41.6 Ms (137 of
580 windows touched): the full Moon sits at the anti-Sun point — the
channel positions — once a month, a hazard specific to this geometry.

**Threshold freeze** (`thresholds.md`). Pseudo-window p-values uniform
to the 10 % dispersion gate on 30 of 32 units (2–9 % below 0.05); the
B van-maanen units are over-dispersed in lane L (7–20 %; 17/21/33
members above T at the three rungs, pseudo-stack values to 76)
because 4FGL J1256.1−0547 = 3C 279 (variability index 33,299) sits
1.8° from the antipode, inside the 2° aperture, with two more variable
4FGL sources at 0.4° and 1.9°. Amendment v1.2 makes those 9 trials
constraint-only before any confirmatory statistic was formed. The
other 52,686 pseudo-windows: 15 above T vs 18.6 expected — the
analytic Poisson thresholds are calibrated.

**Dev** (`results/dev_v1.md`). gj-908 A/B 0.1 AU and the demoted
A gj-1276 1.2 R☉: all 15 statistics below T (max S_event_L 3.16). D1
window (B van-maanen 2020-04-02): n_L 4 / 12 / 105 vs λ 4.1 / 12.4 /
110 at the three rungs. Positive controls through the identical
chain: GRB 130427A (240 lane-L photons vs 4.7 expected in a
1.2 R☉-length window; S_event_L, S_event_H, S_burst all at the 300
cap) and the 3C 454.3 November-2010 flare (7,353 vs 816; S_event_L at
the cap) — **PASS**. Machinery: 70/70 gated photons in the recon
data-server pull re-identified in the era-long pull; gate census raw
→ SOURCE+zenith → interval-gated recorded per position.

**Blind confirmatory** (`results/confirmatory_v1.md`).

| Trial family | Calibrated trials | Exceedances | Expected |
|---|---|---|---|
| S_stack_L / S_event_L (26 units) | 52 | 0 | 0.018 |
| S_stack_H / S_event_H (29 units) | 58 | 0 | 0.021 |
| S_burst (26 units) | 26 | 0 | 0.009 |
| constraint-only (B van-maanen L / burst) | 9 | 9 | — |

Highest calibrated values: A ross-154 0.1 AU S_event_L 3.42 (T 3.451),
A van-maanen 1.2 R☉ S_event_L 3.04 / S_event_H 2.91, A ross-154 0.1 AU
S_stack_L 2.71. Totals over the calibrated units: lane L 10,778 photons
vs 10,516.5 expected, lane H 212 vs 198.2. The nine B van-maanen
exceedances (S_stack_L 21–136, S_event_L 203–257, S_burst 5.0–5.6) all
point at the 2014-04-03 window — the 3C 279 April-2014 flare, one of
the brightest LAT blazar flares on record, 1.8° from the antipode:
the excess-photon centroids sit 0.78–0.84° from 3C 279 and 0.4–0.9°
from 4FGL J1249.3−0545 (also variable); disposition `vetoed_known_source`
(three variable 4FGL sources inside the aperture; ladder step 1). The
strict re-run (θ ≤ 50°, zenith ≤ 90°, Moon/Sun > 16°) leaves them in
place, as a real source should. No retained-ambiguous exceedance.

### 90 %-recovery depths (`results/completeness_v1.md`; lane L ⟨E⟩ = 256 MeV, lane H 5.7 GeV, Γ = 2; containment 0.61 / 0.92 / 0.40 in L / H / U)

| Rung | F90 stack L (ph cm⁻² s⁻¹) | F90 event L best window | F90 stack H | Φ50 / Φ90 burst (1 ks, ph cm⁻²) | P90 stack L (W) | P90 event L (W) | P90 stack H (W) |
|---|---|---|---|---|---|---|---|
| 1.2 R☉ (6 units) | 1.4–3.6 × 10⁻⁷ (med 1.6) | 6–13 × 10⁻⁷ | 0.9–2.7 × 10⁻⁸ | 2.7 × 10⁻³ / 1.2 × 10⁻² | 1.2–3.3 × 10⁵ (med 1.4) | 5.5–12 × 10⁵ | 1.7–5.5 × 10⁵ |
| 2.5 R☉ (9) | 0.5–1.3 × 10⁻⁷ (med 0.86) | 2–5 × 10⁻⁷ | 4–10 × 10⁻⁹ | 3.2 × 10⁻³ / 1.0 × 10⁻² | 1.8–5.1 × 10⁵ (med 3.4) | 0.8–1.9 × 10⁶ | 3.5–8.9 × 10⁵ |
| 0.1 AU (11) | 0.7–3.1 × 10⁻⁸ (med 2.0) | 5–16 × 10⁻⁸ | 1.0–2.4 × 10⁻⁹ | 3.7 × 10⁻³ / 1.5 × 10⁻² | 2.0–8.9 × 10⁶ (med 5.8) | 1.4–4.6 × 10⁷ | 6.5–15 × 10⁶ |

Injections are Poisson counts on the real per-window exposures
(containment from the CALDB PSF at cosθ = 0.8); the burst injections
add PSF-scale photons to the real photon list at a random gated time
of the best-exposure window and re-run the boxcar statistic — Φ90 is
reached only at 1–2 × 10⁻² ph cm⁻² because a 1-ks pulse starting near
the end of a good-time run is mostly lost (Φ50 ≈ 3 × 10⁻³). B
van-maanen lane L and burst: `not_constrainable` (the real statistic
already exceeds T through 3C 279). Systematics on the power
conversion: ⟨E⟩ × 1.13 / × 0.89 for Γ = 1.7 / 2.3 in lane L (× 1.88 /
× 0.62 in H), exposure ± 10 % (φ-dependence and livetime-efficiency
ignored).

**Swift/BAT arm** (`results/bat_v1.md`; constraint-only, 14–195 keV).
391 grazing-window rows (217 events × rungs), 175 with ≥ 1 usable
pointing (offset ≤ 20°, partial coding ≥ 0.1; 61 at 1.2 R☉, 114 at
2.5 R☉; median exposure 1.3 / 2.0 ks). 3σ upper limits median
32 mCrab (1.2 R☉; 10–100) and 27 mCrab (2.5 R☉; 7–84), i.e.
~7 × 10⁻¹⁰ erg cm⁻² s⁻¹ and **1.6 MW / 5.9 MW median through the
cones**. The Crab, carried in the same input catalogue, was measured
in 24 pointings (0.0408 units) — the mCrab scale is empirical, not
assumed. Two windows at combined SNR 3.1–3.2 (A gj-1276 2.5 R☉
2021-09-04 and 2025-09-04; 2 and 4 pointings) are recorded, not
adjudicated — no threshold exists in this arm; the BAT survey's
few-σ systematics floor is the reason it was frozen constraint-only.

**Swift/XRT look** (`results/xrt_look_v1.md`; wolf-359 A 0.1 AU,
obsid 00010119025, 2018-03-06, +2.9 d): 471 net counts in 4,160 s,
mean 0.112 ct s⁻¹ (0.3–10 keV, 47″ aperture), flat across the visit
(41 bins, median 0.108, max 0.219, no excursion above 3× the median).
The star's quiescent level is 0.03 ct s⁻¹ (2017 campaign, 23 visits,
median 0.031; 2021-12 campaign, 80 visits, median 0.029), so the
window visit sits at an **elevated but smooth level, ~3.6× quiescent**
— inside the range of the 2021 monitoring (visits up to 0.67 ct s⁻¹
during a flare episode) and of the 2017 campaign's brightest visit
(0.145). A high-activity state or a decaying flare, nothing
pulse-like; a look, no constraint. The UVOT UV-grism spectrum of the
same visit is handed to §5.17.

**Corridor catalogue screen** (`results/corridor_screen_v1.md`; Pipeline
A step 1, `catalogue_screen_only`). 88 anti-star corridors: 46 have no
catalogued X-ray/γ-ray source within 7′; 42 have 281 sources in all,
**212 `fixed_sky`** by the SGL-track test (positions agreeing to
< 10″ across epochs whose predicted relay displacement at
z = 10,000 AU is ≥ 15″ — 5XMM merged detections, LSXPS/2SXPS spans,
the eRASS1 → eRASS:3 pair — or compact multi-scan stacks) and 69
`single_epoch_unresolved` (gj-687 18 and gj-1221 9 in the LMC
direction, luhman16 16 CSC-only sources, van-maanen 6, sigma-dra 4, a
few eRASS1-only detections). The teegarden-antipode source at 1.0′
(LSXPS J145306.5−165245 = 3eRASS J145306.7−165249) is `fixed_sky`
(LSXPS span 2005–2025 and the compact eRASS:3 stack). eROSITA-DE
upper limits at all 36 western-hemisphere corridors: DR1 eRASS1
0.2–2.3 keV median 7.1 × 10⁻¹⁴ erg cm⁻² s⁻¹ (exposure median 132 s),
DR2 eRASS:3 3.4 × 10⁻¹⁴ (347 s); the 52 eastern corridors
`no_public_data`.

## 3. Physical interpretation

Power convention: P = F_E · π b_rung² (flat-top footprint of radius
b_rung, the §5.17 convention), F_E = F × ⟨E⟩ for Γ = 2. Through the
1.2 R☉ cone the recurrence-stacked lane-L limits are 1.2–3.3 × 10⁵ W
(single best window 0.5–1.2 MW); through 2.5 R☉ 0.2–0.5 MW (stack) /
0.8–1.9 MW (event); through 0.1 AU 2–9 MW / 14–46 MW. Lane H (1–300
GeV) gives comparable powers (0.2–0.9 MW grazing, 6–15 MW at 0.1 AU)
from 10–20× fewer photons because its mean photon energy is 22×
higher. The BAT 14–195 keV arm sits at 1.6 / 5.9 MW per grazing
window. Compared with the optical/UV crossings cells (ZTF/PS1/PTF
~100 W–10 kW, ATLAS 26–55 kW at 0.1 AU, GALEX kW-class UV pulse
trains), the high-energy cell is 10²–10⁴× shallower in power but is
the only constraint at all on relays emitting above 15 keV; per
photon, a γ-ray cell is ~10⁶× more energetic, so the *energy per
crossing* excluded (F90 × exposure × ⟨E⟩) is comparable. Persistence:
`S_stack` tests a relay active at every crossing 2008–2026 (15–18
windows per unit); `S_event` a relay active at one crossing; `S_burst`
a single 1-ks pulse of ≥ 3 × 10⁻³ ph cm⁻² (50 %) — ~10⁷ W-s-class at
1.2 R☉. Pulses shorter than 1 ks and coherent trains are declared
unconstrained (hypotheses §4).

## 4. Lessons for the next adapters

- **The elongation gate is a mission-constraint statement, not an
  archive-scan result**: XMM-Newton's 70–110° solar aspect and
  eROSITA's 90° scan remove every crossing window before any query.
  Check the pointing law first; the recon's XMM slew and eRASS1
  scan-date measurements were confirmations, not discoveries.
- **The anti-Sun point is a crowded place**: the full Moon (a LAT
  source) passes through the channel positions monthly (8.0 Ms of
  livetime excluded here), the opposition asteroid stream crosses
  them (optical surveys), and any bright variable source within a
  PSF radius owns the aperture — 3C 279 at 1.8° from the van-maanen
  antipode cost 9 of 145 trials. A PSF-scale aperture on a
  fixed point needs a per-position variable-source census at the
  freeze, not at adjudication.
- **Fermi's survey profile changed in 2018-03**: window-by-window
  livetime from the spacecraft files is mandatory; 62 of 580 windows
  are empty post-anomaly. "Continuous all-sky" is a pre-2018 statement.
- **Analytic Poisson thresholds work when validated per trial**: the
  Jeffreys floor (amendment v1.1) and the per-trial ensemble exceedance
  gate (v1.2) were both needed; 52,686 pseudo-windows then show 15 vs
  18.6 exceedances at T.
- **HEASoft in a container with a local CALDB is a one-hour setup**
  (image `chbrandt/heasoft` 6.24, Swift BAT CALDB 23 MB; `batsurvey`
  with the channel positions as `incatalog` gives per-pointing rates
  at arbitrary positions in ~35 s each; run with `--user`, `HOME=/tmp`,
  `HEADASNOQUERY=1`, `:z` mounts). The Crab in the same catalogue
  calibrates the scale for free.
- The LAT data server answers era-long (18-yr) 3° photon pulls in
  ~90 s; the weekly photon files are never needed.

## 5. Products

`surveys/highenergy-crossings/`: `hypotheses.md` (v1.0 + amendments
v1.1, v1.2), `thresholds.md`, `corridor_screen.md`, `notes/` (recon),
`scripts/` (recon_scan, recon_summary, recon_erodat_rows,
fetch_photons, hlib, coverage, freeze_thresholds, dev_stage,
confirmatory, completeness, bat_stage, xrt_look, corridor_screen),
`results/` (recon_*_v0, coverage_v1, thresholds_v1, dev_v1,
confirmatory_v1 + units_v1/, completeness_v1, bat_v1, xrt_look_v1,
corridor_screen_v1). Runs under `runs/highenergy-crossings/` (recon
1.2 GB; v1: photons 250 MB + sha manifests, spacecraft cache 1.1 GB,
CALDB, BAT obs 2 GB + batsurvey out 86 GB — the per-pointing images are
regenerable from `obs/` and can be deleted after the collect stage —
XRT event lists, catalogue snapshots).
Nothing committed to git.
