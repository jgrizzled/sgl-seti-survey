---
title: "S2 1 AU outer-skirt blended search — recon"
date: 2026-09-07
status: "recon complete; hypotheses drafted; freeze decisions pending user approval"
---

# Outer-skirt recon (plan §5.25)

Goal: turn the geometry study's caveat 3 — the S2 uplink's 1 AU rung
has an outer skirt at ε ≳ 30–40° that is night-sky visible — into a
searchable, pre-registrable cell on the ATLAS substrate, or close it
for the record. Every number below is measured (scripts in
`scripts/`, outputs in `results/`); the ATLAS-survey conventions are
inherited where they apply.

**Headline.** The cell is searchable, but not as the 1 AU rung: that
rung has no temporal signature (Earth is inside a 1 AU beam all
year). What is searchable is the **elongation-locked step of sub-1 AU
rungs** (r_b = 0.90 / 0.95 AU, edges at ε 64° / 72°) whose in-beam
band on the S2 side sits between the measured ATLAS floor (~50°) and
the edge — recurring annually, mirrored at opposition, and tested
against a null ensemble of field stars in the same exposures rather
than the temporal pseudo-windows that the season-scale bands defeat.
Population 18 targets (one blended pair excluded); 162 full-history
reduced-mode tasks at 12–15 min each ≈ 1.5–2 days of serial queue
drain.

## 1. Geometry (`skirt_geometry.py` → `skirt_geometry_v1.json`)

- For a star at infinity the Sun–star axis distance of Earth is
  **b_e = r_E sin ε** exactly, S2 side ⇔ ε < 90°, channel-A side ⇔
  ε > 90°. Daily curves 2015-07 → 2026-09 for all 85 crossing-list
  targets; the yearly conjunction minima reproduce the frozen list's
  942 in-era S2 `b_min` values to a median 0.0003 AU (max 0.025 AU —
  daily-grid quantisation for β ≈ 0 stars). Conjunction elongation =
  |ecliptic latitude|: from 0.1° (van-maanen) to 81.6° (Struve 2398).
- A top-hat beam of radius r_b is on for r_E sin ε < r_b: **both**
  ε < ε_b (S2 skirt) and ε > 180° − ε_b (opposition). ε_b = 30° (0.5
  AU), 44° (0.7), 53° (0.8), 58° (0.85), 64° (0.90), 72° (0.95),
  89° (1.0). Targets with |β| ≥ ε_b never enter the beam (ledger-
  only); within ~2° of ε_b the membership flips with Earth's
  perihelion/aphelion r_E (edge-riding — treated as ledger-only).
- The 1 AU rung is on for every ε: **no temporal signature**; it
  becomes a coverage ledger (nights per S2 window at ε ≥ floor),
  exactly the status the ATLAS survey gave channel A's 1.0 AU rung.

## 2. Substrate facts (`recon_probe.py`, `probe_analysis.py` → `probe_elongation_v1.json`)

Two full-history **reduced-mode** tasks at non-target Gaia DR3 field
stars on the ecliptic (G 14.1, high Galactic latitude, PM < 5 mas/yr):
RA 180.05 / Dec +0.00 and RA 0.23 / Dec +0.10. Declared pre-freeze
data contact: none at any survey position.

| quantity | ecl180 | ecl000 |
|---|---|---|
| rows / FAQ-mask pass | 4,043 / 3,659 (91 %) | 4,193 / 3,773 (90 %) |
| MJD span | 57230 → 61252 | same era |
| execution wall time (full history, reduced) | **14.7 min** | **11.7 min** |
| min solar elongation (good rows) | 40.6° | 40.6° |
| 1st / 5th percentile ε | 51.6° / 66.2° | 53.8° / 65.5° |
| o precision per exposure (fractional MAD, G 14) | 1.6 % (photon 0.4 %) | 1.6 % |
| … at ε < 70° / ε > 120° | 1.8 % / 1.6 % | 1.8 % / 1.6 % |

- **Night efficiency vs elongation** (good nights per available day,
  pooled): 40–45° 2 %, 45–50° 5 %, 50–55° 8 %, 55–60° 11 %, 60–65°
  18 %, 65–70° 21 %, 70–80° 21–28 %, 80–95° 25–30 %, > 100° 33–45 %.
  ATLAS goes to nautical twilight (Sun altitude −12° at the lowest-ε
  epochs, verified per site) but the ε < 55° sky is covered on ~1
  night in 15; the **practical floor is ~50°** (declared), ~55–60° for
  the southern units (El Sauce min 55°, Sutherland 58° on ecl000).
- **Obs-prefix → site mapping corrected**: from star/Sun altitudes at
  each epoch, **03 = Sutherland, 04 = El Sauce** (the 2026-08-25 ATLAS
  recon note had them swapped; harmless there — only unit-era dates
  were used — but the airmass layer here depends on it).
- Median airmass of the ε < 70° epochs 1.5–2.0 by site; the
  opposition band sits near 1.0–1.3: the S2 skirt is the high-airmass,
  twilight regime the geometry study anticipated.
- **The queue applies proper motion**: `propermotion_ra`,
  `propermotion_dec` (mas/yr) and `radec_epoch_year` are task fields
  (OpenAPI snapshot; present in the task JSON). One full-history task
  per star instead of the per-interval scheme the 2-parameter forced
  fit would otherwise demand at 1–5″/yr; the output RA/Dec columns
  verify the application (dev gate G1). Whether `propermotion_ra`
  is μα or μα cos δ is undocumented — G1 tests both readings.

## 3. Expected in-beam coverage (`expected_coverage.py` → `expected_coverage_v1.json`)

Geometry × measured efficiency, nights per year (skirt side /
opposition side / out-of-beam quadrature band):

| target | β | o_est | r_b 0.90 | r_b 0.95 |
|---|---|---|---|---|
| gj-1276, teegarden, gj-1111, gj-518, gj-1087 | 0–18° | 12.5–15.4 | 4 / 50 / 31 | 7–8 / 55 / 21 |
| gj-2012, gj-9193, gj-3306, gj-11068 | 25–31° | 13.5–15.7 | 5 / 47 / 34 | 9 / 54 / 24 |
| gj-915, gj-3512, gj-12724, eps-ind-b, gj-11547 | 39–43° | 12.7–18.5 | 7–8 / 45 / 41 | 11–13 / 54 / 28 |
| gj-13157, gj-3112 | 64–66° | 13.6–14.3 | edge / — | 18 / 30–36 / 52–60 |
| wolf-1069, gj-293 | 72°, 79° | 12.6–13.7 | — | edge / — |

The skirt side alone is thin (4–13 nights/yr, 45–140 over the era);
the symmetric statistic adds the opposition band's ~50 nights/yr.
The 0.85 AU rung (edge 58°) would have 2–5 skirt nights/yr — not
proposed.

## 4. Target population (`target_scope.py` → `target_scope_v1.json`; `control_pool.py`)

- Gaia DR3 photometry for all 88 registry targets (by source id; 12
  without a Gaia id — 7 bright stars and 5 CNS5 entries — by 36″ cone).
  Screening transform o ≈ G − 0.22 + 0.13 (BP−RP), anchored on the four
  substrate-measured unsaturated ATLAS-survey targets (rms 0.13 mag;
  the three saturated anchors were dropped — their measured o is
  compressed by 1–2.5 mag). **The measured reduced-mode value decides
  at the cut stage (ATLAS D2 rule)**.
- Result: 59 excluded (o_est < 12.5), 7 bright non-Gaia stars
  excluded, WISE 0855 (no optical counterpart), van-maanen excluded on
  its measured 12.30; **20 candidates** (marginal: gj-1111 12.5,
  wolf-1069 12.6, teegarden 12.6 measured, gj-915 12.7; ok: gj-3512 …
  eps-ind-b 18.5). Blends: gj-13157 has a G 17.2 neighbour at 4.7″
  (7 % constant dilution — flagged, kept); **luhman16-a/b excluded**:
  the registry position's Gaia cone match is a G 16.0 background star
  25″ away, a second G 16.8 star sits 9″ away, and neither carries the
  pair's 2.8″/yr motion — no verified Gaia counterpart, and the pair
  sweeps ~30″ through that field over the era (tphot single-PSF
  blending). Registry note: the luhman16 epoch-2016 position needs
  checking against Gaia DR3 (flagged for the registry maintainers).
- Sensitivity direction: absolute contrast improves for fainter
  stars (1 % of an o 13.5 star is 18.5 AB; of gj-1276 20.2; of
  eps-ind-b ~23) — faint targets are kept.
- **Null-ensemble pools** (Gaia DR3, 1.0° cone = same ATLAS exposures,
  isolated, non-variable, RUWE < 1.4, PM < 0.1″/yr): 16 per target.
  Blue/moderate targets match within |ΔG| ≤ 0.3, |Δ(BP−RP)| ≤ 0.4.
  The five reddest targets (BP−RP 3.7–5.2: teegarden, gj-1111,
  gj-3512, gj-12724, gj-11547) have **no colour-matched field stars**
  within 1° (pools reach only Δ(BP−RP) 1.8–3.5) — the colour-
  dependent airmass residual is therefore not shared by their
  controls, which is why the airmass-regression layer is proposed
  as a frozen construction rather than an option.

## 5. Cost

Full-history reduced task 12–15 min; 18 targets + 8 × 18 controls =
**162 tasks ≈ 32–40 h of serial drain** (resumable, `atlas_drain.py`
pattern). Ledger-only targets (wolf-1069, gj-293 at both rungs) can
be pulled without controls: −16 tasks.

## 6. What this means for the design (→ `hypotheses.md`)

1. The plan's cell is delivered as the **sub-1 AU elongation-locked
   step** (rungs 0.90 / 0.95 AU), with the S2 skirt as one side of the
   step and per-side statistics reported; the 1 AU rung is a ledger.
2. Controls are a **null ensemble** (8 field stars in the same
   exposures) — the learnings' designated route for wide rungs.
3. Reduced mode + server-side PM + airmass regression + nightly means
   are the systematics stack; injections run through the same stack.
4. Expected depth: m90 ≈ o_star + 4 to 5.5 (0.6–1 % excess), i.e.
   3–15 MW through a 0.9–0.95 AU beam for the o 13–15 targets.

## 7. Freeze decisions presented for approval

Listed in the session summary (most consequential first): D1 cell
definition (step search incl. the A-side mirror), D2 rung ladder,
D3 population and exclusions, D4 statistics/trials, D5 airmass layer,
D6 null-ensemble controls and validity gate, D7 task budget scope,
D8 dev/confirmatory split and seed.
