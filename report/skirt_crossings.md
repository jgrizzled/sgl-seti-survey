---
title: "S2 1 AU outer-skirt blended search — report"
date: 2026-09-10
status: "COMPLETE — 0 candidates; the sub-1 AU uplink step cell opened at MW class on 16 targets"
---

# S2 1 AU outer-skirt blended search (plan §5.25)

Companion files: `surveys/skirt-crossings/hypotheses.md` (freeze v1.0,
D1–D8 approved 2026-09-07, amendments A1–A4), `notes/skirt_recon_2026-09-07.md`,
`configs/threshold_freeze_v1.json`, `results/*.json`. Execution log:
`notes/project_history.md` §27.

## 1. Question and answer

**Question.** The 2026-08-25 sunward geometry study left one sunward
cell outside the coronagraph substrate: the S2 uplink (star → relay
beam intercepted by Earth on the anti-target side of the Sun) at the
1 AU rung, whose outer skirt at solar elongation ε ≳ 30–40° is
night-sky visible. Can a night-sky archive test it?

**Answer.** Not as the 1 AU rung — Earth's distance from the Sun–star
axis is exactly r_E sin ε, so a 1 AU beam contains Earth all year and
a persistent transmitter is a constant blended with the star (the
rung is a coverage ledger, §5). What *is* testable is the
**elongation-locked step of sub-1 AU uplink beams**: a beam of radius
r_b switches off at ε_b = arcsin(r_b/r_E) and on again past
180° − ε_b at opposition, so the star should be brighter by a constant
in both bands than in the quadrature band between them, every year.
For r_b = 0.90 / 0.95 AU the S2-side band lies between the measured
ATLAS solar-elongation floor (~50°) and the edge (64° / 72°). The
survey searched that step, on both sides, on the ATLAS reduced-mode
substrate with a null ensemble of field stars in the same exposures.

**Result: 0 candidates.** 81 calibrated trials on 28 units (16
targets × 2 rungs, less 4 ledger-only), 24 exceedances vs 9.0
expected — over budget by the two red targets and one crowded one
(§3) — all adjudicated non-promotable or retained-ambiguous under a
pre-declared rule (none promotable). 90 % depths on the 57
non-exceedance trials: **0.6–24 MW through the 0.90–0.95 AU beams
(S_sym median 6–8 MW; 0.6–1.5 MW on the three faintest targets)**,
i.e. a 1.5–6.5 mag contrast below each star.

## 2. Design (frozen 2026-09-07; §11 of the hypotheses)

| item | frozen choice |
| --- | --- |
| substrate | ATLAS reduced-mode (target-image tphot) full-history forced photometry, o searched / c annotation, **server-side proper motion** (`propermotion_ra` = μα cos δ, verified to 0.03″ over 55″ of motion on teegarden — gate G1) |
| cell | symmetric step S_sym (primary), S_skirt (S2 side only — the plan's cell), S_year (single-cycle max); S_opp annotation |
| rungs | 0.90 / 0.95 AU (edges 64° / 72°); 0.85 AU not searched (2–5 skirt nights/yr); 1.0 AU ledger |
| epochs | nightly-unit medians (≥ 2 exposures); response gate R ≥ 0.9; airmass layer f = a + b(X − 1) per star |
| controls | null ensemble: 8 Gaia DR3 field stars within 1° (same exposures), magnitude/colour-matched where the field allows, drawn by seed 20260907 from a pool of 16; T = max valid control; exceedance S > max(T, 0); p = 1/(n_valid + 1) |
| gate | ≥ 3 cycles with ≥ 5 IN and ≥ 5 OUT nights (skirt ≥ 3); ≥ 4 valid controls |
| population | 18 unsaturated-by-estimate targets → 16 after the measured cut (gj-1111 12.44, wolf-1069 12.38 excluded); luhman16 excluded at recon (no verified Gaia counterpart) |
| split | dev = teegarden + gj-2012; 14 confirmatory targets; blind run once |

**Pre-confirmatory amendments** (all recorded before any
confirmatory-target statistic was formed): A1 G2 wording; A2 the
airmass layer is neutral on the ensemble floor and retained for the
target-specific colour term; **A3 parallax response correction**
(§4.1); **A4 colour-unmatched rule** — an exceedance on a target with
no control within 1.0 mag of its BP−RP is `retained_ambiguous_colour_unmatched`,
not promotable.

## 3. Results

**Dev (12 trials).** 7 exceedances vs 1.3 expected, all teegarden
(six trials; S_sym 2.0 at both rungs, Δ̄ +0.5 %, symmetric; S_year
4.0–4.5 from the 2020-21 cycle at +4 %) plus one gj-2012 S_sym at
0.001 vs T −1.02 (vetoed asymmetric). teegarden is an M7 flare star
2.7 mag redder than its nearest control: `retained_ambiguous_colour_unmatched` (A4),
one S_skirt trial `non_promotable_single_cycle`.

**Blind confirmatory (69 trials, 14 targets).** **17 exceedances vs
7.7 expected**, 0 promotable:

| target (o) | rung | trial | S / T | Δ skirt / opp | disposition |
| --- | --- | --- | --- | --- | --- |
| gj-915 (12.98) | 0.90, 0.95 | S_skirt | 0.60/0.30, 1.60/1.11 | +0.6 / −0.1 %, +1.1 / −0.1 % | asymmetric (S2 side only; ensemble-shared) |
| gj-3512 (13.17) | 0.95 | S_skirt, S_year | 0.85/−0.12, 2.77/2.14 | +0.3 / +0.3 % | **retained_ambiguous_colour_unmatched** (gap 1.8 mag) |
| gj-3306 (13.59) | 0.95 | S_sym, S_year | 1.13/0.75, 3.23/2.54 | +0.1 / +0.05 % | single-cycle; ensemble-shared |
| gj-13157 (14.69) | 0.95 | S_skirt | 1.30/0.55 | +9.8 / +1.9 % | asymmetric; **crowded** (§4.2) |
| gj-518 (14.22) | 0.90 | S_sym, S_year | 1.02/1.02, 2.57/2.36 | −0.2 / +0.2 % | asymmetric (opposite signs) |
| gj-12724 (14.48) | 0.90, 0.95 | S_year, S_sym | 0.96/0.45, 0.63/−0.73 | −1.1 / −0.03 %, +0.02 / +1.2 % | asymmetric |
| gj-11547 (14.40) | 0.90, 0.95 | S_skirt | 3.00/1.44, 4.02/3.27 | +4.5 / −0.6 %, +3.8 / −0.8 % | asymmetric (S2 side only, red, gap 2.5 mag) |
| gj-11068 (15.79) | 0.90, 0.95 | S_sym, S_skirt | 0.13/−0.29 … 1.82/1.40 | +0.5 / −0.25 %, +1.1 / −0.2 % | asymmetric (S2 side only) |

Reading. Thirteen of the 17 fail the symmetry veto — the excess sits
on the S2 (twilight, high-airmass) side alone or with the opposite
sign at opposition, which a beam cannot do. The four S2-side-only
excesses on red targets (gj-11547 +4 %, gj-11068 +1 %, gj-915 +1 %,
gj-3512) are the colour-dependent twilight/airmass systematic the
freeze anticipated: the five reddest targets have no colour-matched
controls (recon §4), and the two with in-range exceedances that
survive the ladder (teegarden, gj-3512) end under A4. The over-budget
count (24 vs 9.0) is therefore not a calibration failure of the max
rule on ordinary stars — the 9 targets with colour-matched controls
produced 6 exceedances on 45 trials vs 5.0 expected — but the
declared limit of a null ensemble that cannot match a BP−RP > 3.5
star in a 1° field.

**Ledger.** 6 ledger-only units (wolf-1069, gj-293 at both rungs;
gj-13157, gj-3112 at 0.90 — |β| ≥ ε_b); 2 saturation-excluded
targets; no constraint-only unit — every searchable unit had 8 valid
controls.

## 4. Two findings about the substrate

### 4.1 Parallax is an elongation-locked photometric step

The dev search put every target above all eight of its controls,
and the one property every target has and no control shares is a
parallax. The parallactic displacement peaks at quadrature — exactly
the out-of-beam band — so a fixed-position PSF fit loses flux there
and nowhere else: a step of the signal's sign and shape, ~1 % for
teegarden (π = 0.26″, forced-vs-apparent offset 0.26″ at ε 80–100°
vs 0.07″ at opposition). The response of a fixed-position fit to an
offset d is not the peak factor exp(−d²/2σ²) but was **measured**
(six offset tasks on two quiet controls, 9,219 matched exposures):
R = exp(−d²/(k σ²)) with k = 2.42 (`results/response_calib_v1.json`).
Amendment A3 corrects every flux by 1/R_i against the apparent
position (PM + parallax). Residual after correction (star-to-star
spread in k): ±0.2 % on Δ for the π ≈ 0.26″ targets. This is
general: any fixed-position forced-photometry series of a nearby
star carries an annual parallax signature at the ~1 % level.

### 4.2 Crowding defeats the single-PSF fit beyond the 10″ cone

gj-13157 sits 14″ from a G 9.4 star and 11″ from a G 12.1 star
(both outside the recon's 10″ isolation cone). Its tphot chi/N is
35,800 vs 13 for its controls and its nightly fractional MAD is
0.44 vs 0.017 for its magnitude bin: the series is wing
contamination, not photometry of the star. Its exceedance was vetoed
by the ladder regardless (asymmetric) and its depths (24–150 MW) are
reported systematics-limited. Post-blind annotation only — no
disposition was changed. Lesson: the isolation test needs a
brightness-weighted radius (≳ 30″ for stars ≥ 5 mag brighter).

## 5. Depths and physical interpretation

**Injection completeness** (`results/completeness_v1.json`): a
top-hat step of flux f added to every in-beam exposure of the real
series, 100 draws per magnitude with 20 % nightly dropout, frozen
chain, recovery S > max(T, 0). For the 24 exceedance trials the real
statistic already exceeds the threshold, so injection recovery is
the real excess — no depth is claimed for them (a reporting rule
recorded at completeness, hypotheses §12 A5). The 57 constraining
trials:

| statistic | 0.90 AU beam | 0.95 AU beam |
| --- | --- | --- |
| S_sym (recurrent, both sides) | 0.85–13 MW, median 8.5 (9 trials) | 0.6–24 MW, median 6.0 (11) |
| S_skirt (S2 side only) | 3.4–34 MW, median 6.9 (6) | 1.5–95 MW, median 9.5 (9) |
| S_year (single cycle) | 2.1–34 MW, median 11 (10) | 0.95–95 MW, median 12 (12) |

Per-target S_sym depths (o-band m90 → line power through the beam,
690 nm / 260 nm effective width, ±0.1 mag scale, ±0.2 % parallax
residual for π ≥ 0.25″): gj-1276 21.0 → 1.4 MW; eps-ind-b 21.5–22.0
→ 0.6–0.85 MW; gj-9193 20.0–20.5 → 2.1–3.8 MW; gj-2012 20.5 → 2.4 MW;
gj-518 20.0 → 3.8 MW; gj-3306 20.0 → 3.4 MW; gj-1087 19.0–19.5 →
6–8.5 MW; gj-915, gj-3512, gj-3112, gj-12724, gj-11547 18.5–19.5 →
6–15 MW; gj-13157 18.0 → 24 MW (systematics-limited). Contrast:
1.5–6.5 mag below the star (median 5.7).

**Physical reading.** A transmitter at the target system feeding the
Sun's lens with a beam 0.90–0.95 AU wide at 1 AU (a sub-metre
aperture at 3–5 pc, or a deliberately broad beam) radiating a
persistent line in 560–820 nm would have been detected at 90 % in the
symmetric annual step above **0.6–24 MW** (0.6–1.5 MW at the three
faintest targets). Compared with the antipode channels this is 10³×
shallower in flux but through a beam 80–90× the 0.1 AU rung's area,
and it is the first constraint of any kind on the sub-1 AU uplink
rungs, which the antipode channels cannot see and the heliospheric
imagers only see at ε ≲ 6°.

**1 AU ledger** (`results/ledger_1au_v1.json`): every conjunction
cycle of every target in the era has S2-side nights at ε ≥ 50°
(7–31 per cycle, 77–206 per target over 5–12 cycles) and 200–650
A-side nights; the minimum observed elongation is 39–49° for the
low-|β| targets. This is the coverage record for the 1 AU rung —
not a constraint (no temporal signature exists).

**Declared unconstrained.** The 1 AU rung; Gaussian-profile beams
(a soft edge lowers the step); 1064/1550 nm (outside o); sub-day
pulses; schedules avoiding Earth's in-beam months; the five
colour-unmatched targets at their exceedance trials
(systematics-limited); gj-13157 (crowded).

## 6. Lessons for the next adapters

- **A fixed-position forced-photometry series of a nearby star has a
  ~1 % annual parallax step** (§4.1) peaking at quadrature. Correct
  it with a measured response exponent — the two textbook forms
  differ by the full size of the effect.
- **A null ensemble cannot be colour-matched for BP−RP > 3.5 targets
  in a 1° field**; a colour-dependent twilight/airmass systematic then
  survives every symmetric test on the S2 side. Either widen the
  control field (several degrees) and accept a different exposure
  set, or declare the reddest targets systematics-limited up front
  (A4 did the latter).
- **Injection into the real in-beam set with an S > max(T, 0) rule
  returns the real excess for exceedance units**; a step-statistic
  completeness must either inject into pseudo-sets or test
  S_injected − S_real against T.
- **Isolation cones must scale with neighbour brightness** (§4.2).
- The ATLAS queue applies proper motion server-side (μα cos δ);
  one task per star replaces per-interval schemes at any PM.
- ATLAS Obs-prefix sites: 03 = Sutherland, 04 = El Sauce (the
  2026-08-25 recon note had them swapped).
