# PTF beam-crossings survey (Pipeline B) — report

Sixth Pipeline-B (§3.5) survey: a search of Palomar Transient Factory
level-1 epochal images for optical emission during Earth's crossings
of hypothesized Sun–star relay beam axes, 2009–2015 — the era
concurrent with PS1's but from a different site, cadence and filter
set, giving the covered-window record its first **independent
re-observations** of PS1-era windows, including the van-maanen
deep-graze family PS1 lost to chip-gap masks. Executed 2026-08-26.
**No candidates: 10 blind confirmatory trials, 0 exceedances (1.11
expected control crossings).**

## 1. Construction and provenance chain

1. **Universal crossing list** `crossings/universal_v1`
   (`xng-a09e2db7681d`): every b(t) minimum for 88 endpoints × 2 link
   directions, Earth-center, 1980→2028.
2. **Hypothesis freeze v1.0** (`surveys/ptf-crossings/hypotheses.md`,
   D1–D8 approved 2026-08-26): the PS1 crossings construction with
   PTF substitutions — era MJD 54891–57051 (measured global span of
   `ptf.ptf_procimg`; late iPTF never publicly released);
   **scie-direct substrate** (no public difference images; level-2
   reference coadds are annotation inputs only); g + Mould R primary
   (532 nm line in g, 658 nm generic leakage in R; 1064/1550 nm
   outside both bands, unconstrained); unit = (target, radius, band)
   with **S_event primary and S_stack only where ≥ 2 covered events**
   (D4 — coverage here is single-event dominated, the reverse of
   ZTF's stack-primary situation); A 1.0 AU constraint-only by the
   standing window≈season theorem; per-frame field-star calibration
   primary for all frames (measured `photcalflag ≈ 0`); dmask fatal
   template 65533 (Laher et al. 2014 Table 15, only the
   object-detected bit non-fatal); Palomar-vs-Earth-center offset the
   standard ≤ 0.8 %-of-rung budget term.
3. **Coverage intersection** (`results/coverage_v1.*`; 93
   snapshot-disciplined IBE discovery boxes): B 1.2 R☉ **0/24**
   (structurally uncovered — the 0.3–0.6 d windows fall between
   epochs; closest 0.36 d vs a ±0.3 d edge); B 2.5 R☉ **3/30** —
   van-maanen 2011 at **b = 0.28 R☉** with a same-night g pair
   in-window, ross-128 ×2; B 0.1 AU **7/42** — the complete
   van-maanen April recurrence family (2009/2010/2011/2013,
   including `evt-4c4ea2c39734`, the 2010 event PS1 covered and lost
   entirely to masks, here holding 2 independent R epochs) plus
   ross-128 ×2 and wolf-359 ×1; A 0.1 AU 4/42; A 1.0 AU 76/513
   (constraint-only). Campaign cadence: teegarden, gj-908 and
   ross-154 antipodes have zero PTF epochs.
4. **Channel-A saturation cut** (`results/saturation_cut_v1.*`;
   E_R 14.5 / E_g 15.0, Mould-R via Jordi 2006): of the four covered
   A units only gj-1276 R survives — teegarden R (13.7), van-maanen g
   (12.4) and ross-128 R (9.9) are excluded-class. Verified at dev:
   R saturation bracketed 12.98 (bits fire) / 14.06 (clean) against
   the frozen 14.0 — no amendment; the ross-128 star core fires
   bits 8+6 (exclusion confirmed).
5. **Threshold freeze v1.0** (`thresholds.md`,
   `configs/threshold_freeze_v1.json`, sha256:26e235ee…, bound to
   hypotheses + coverage + saturation hashes): 7 searched units, 11
   trials; 8 ring-trajectory controls (B) / 8 temporal pseudo-windows
   (A); S > max(T, 0); PS1-DR2 calibration gate ≥ 5 stars, scatter
   ≤ 0.2 mag; dev/confirmatory split frozen verbatim at D8 (both
   headline B families blind).
6. **Dev stage** (`results/dev_search_v1.md`): machinery validated
   end-to-end; realizations fixed (dwarf-locus calibrator
   restriction; enlarged cutouts); wolf-359 B single-epoch class
   clean (S 0.10 vs T 2.10); **gj-1276 A resolved constraint-only at
   the frozen offset-validity gate** (1/12 valid — campaign cadence)
   → 0 searched dev trials; positive control (8971) Leucocephala
   recovered position-locked, internal RMS 0.025 mag, absolute
   offset +0.24 within Horizons-prediction systematics.
7. **Amendment v1.1** (frozen before re-reduction, documented in
   `scripts/confirmatory_search.py`): confirmatory cutouts
   256 → 384 px — the 256-px support was inoperable at the ross-128
   antipode (142/144 frames < 5 in-frame calibrators; a
   surface-density failure). No statistic, threshold, gate or
   selection element changed; all units re-fetched and re-measured
   uniformly. The superseded pass (kept:
   `confirmatory_v1_cut256_superseded.json`) had measured only the
   van-maanen units — all null, margins −0.7 to −15.7 — and nothing
   in the amendment was conditioned on those values (the amended
   van-maanen statistics moved by ≤ 0.13).

## 2. Results

**0 candidates.** Blind confirmatory run
(`results/confirmatory_v1.json`; 600/600 cutout pairs fetched, zero
404s): **10 trials, 0 exceedances** vs 1.11 expected control
crossings (P(0) ≈ 0.31 — unremarkable). Gate attrition at 384 px was
tiny (0–6 frames/unit); every coverage-stage event retained its
epochs at measure time; k rescales 1.00–1.10 (the robust-background
+ Poisson variance model is honest at these field densities).

| unit | events (epochs) | S_event vs T | S_stack vs T |
|---|---|---|---|
| van-maanen B 2.5 R☉ g — **b = 0.28 R☉, 2011** | 1 (2, same-night pair) | 1.45 / 2.16 | — |
| van-maanen B 0.1 AU g | 2 (4 + 12) | 1.19 / 3.16 | 0.95 / 2.98 |
| van-maanen B 0.1 AU R — incl. the PS1-lost 2010 event | 3 (2 + 4 + 3) | 2.02 / 17.81 | 1.51 / 8.14 |
| ross-128 B 2.5 R☉ R | 2 (7 + 2) | 1.58 / 2.61 | 1.25 / 2.24 |
| ross-128 B 0.1 AU R | 2 (18 + 2) | 2.76 / 5.48 | 2.56 / 5.26 |
| ross-128 B 0.1 AU g | 1 (44) | 3.86 / 5.12 | — |

The van-maanen 0.1 AU R threshold (17.8) is confusion-dominated — a
bright static source rides one ring trajectory; the threshold absorbs
it by construction (the v1 lesson), and the track itself passes
3.0″ from the nearest catalogued static source (> 2″, no
annotation). Constraint-only and ledger classes: gj-1276 A R (dev,
offset gate), wolf-359 B single-epoch (dev, S 0.10 vs T 2.10), 3
saturation-excluded A units, the B 1.2 R☉ rung
(coverage-without-statistic), and the zero-coverage positions.

### 90 %-recovery depths (`results/completeness_v1.json`; 532 nm line in g, 658 nm in R)

Exact Moffat-β3 stamp-response injections through the identical
star-calibrated chain, d = 1 persistent-recurrent chord (present at
every covered recurrence), 200 draws × 15 mags, recovery against
each unit's measured threshold; all five z grid points recoverable
in every unit (no dead z).

| unit | m90 (AB) |
|---|---|
| van-maanen B 2.5 R☉ g | 21.37 |
| van-maanen B 0.1 AU g | ≥ 22.0 (grid-censored) |
| van-maanen B 0.1 AU R | 20.58 (confusion threshold) |
| ross-128 B 2.5 R☉ R | ≥ 22.0 (grid-censored) |
| ross-128 B 0.1 AU R | 21.61 |
| ross-128 B 0.1 AU g | ≥ 22.0 (grid-censored) |

## 3. Physical interpretation (532 nm line in g unless noted)

Rescaling the ZTF anchor (m90 21.8 → 130 W transmitter power through
the 2.5 R☉ solar-grazing cone at 550 AU, gain ≈ 9×10⁹):

- **The photosphere-grazing cell, pre-2015 era, first constraint:**
  a relay beaming through the 2.5 R☉ cone during the van-maanen 2011
  b = 0.28 R☉ window would have been seen at 90 % confidence above
  **≈ 190 W** (g). The ross-128 2.5 R☉ windows constrain **≲ 110 W**
  (R, censored depth). These are per-window statements for the
  covered 2010–2012 windows — an era no other archive in the
  programme could search at these b values (PS1's masks took its
  only instance).
- **Wide-rung (0.1 AU beam) downlink:** ≈ 8 kW (van-maanen g,
  censored), ≈ 30 kW (van-maanen R, confusion-limited), ≈ 11 kW
  (ross-128 R), ≈ 8 kW (ross-128 g, censored) through the 0.1 AU
  cone — kW-class limits on pointing-margin downlinks across the
  2009–2013 April/September recurrences.
- The recurrence-stack cell (S_stack) closes clean on 4 units — the
  first multi-instance recurrence test of the van-maanen April
  family (2009–2013 in g and R) and of ross-128 (2010–2012).
- Declared unconstrained: 1064/1550 nm (outside both bands),
  sub-exposure pulse schedules beyond ×d scaling, transmitters
  avoiding Earth-crossing windows, the B 1.2 R☉ rung (structurally
  uncovered), and everything at the zero-coverage positions.

## 4. Lessons for the next crossings adapters

1. **Calibrator support is a field property, not a constant**: the
   same frozen ≥ 5-calibrator gate that passed at 256 px on one
   antipode was inoperable at another of similar galactic latitude —
   size cutouts for the sparsest field, or gate on measured
   calibrator density at the plan stage (amendment v1.1's origin).
2. The dwarf-locus color restriction on transform-based calibrators
   is load-bearing: off-locus red stars inflate per-frame ZP scatter
   ~5× (0.25 → 0.05).
3. PTF serves every listed product (710/710 cutout fetches across
   dev + confirmatory, zero 404s) with published MD5s — the
   highest-integrity archive interface in the programme so far;
   `achecksum` truncation (prefix-verify) and blank in-header MAGZPT
   on non-photometric frames are the only quirks.
4. Campaign cadence defeats temporal pseudo-window controls at
   on-star positions (gj-1276 A: 1/12 valid offsets) — for
   campaign-era archives, channel-A units need the offset-validity
   gate assessed at the threshold freeze, not assumed (the D7
   phrasing "assessed per unit, not assumed" earned its place).
5. The universal-list side-of-axis condition matters at scoping:
   the recon-stage proxy scan omitted it and doubled the event
   denominators (covered sets unchanged).

## 5. Products

`surveys/ptf-crossings/`: `hypotheses.md` (freeze v1.0),
`thresholds.md` + `configs/threshold_freeze_v1.json` (+ amendment
v1.1 in `scripts/confirmatory_search.py`), `results/era_scope_v0.json`,
`coverage_v1*`, `saturation_cut_v1*`, `dev_search_v1*`,
`saturation_verify_v1.json`, `asteroid_control_v1.json`,
`confirmatory_v1.json` (+ superseded 256-px pass),
`completeness_v1.json`; adapter `sglsurvey/adapters/irsa_ptf.py`;
`build_flux_map_ptf` in `sglsurvey/photometry.py`; snapshots and
query records under `runs/ptf-crossings/`. Recon:
`notes/ptf_recon_2026-08-26.md`.
