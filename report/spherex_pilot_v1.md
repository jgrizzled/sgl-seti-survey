---
title: "SPHEREx pilot v1.0 — third-adapter shakedown of the SGL relay search"
date: 2026-08-20
status: "pilot complete; 3 stars, 288 constraints, 0 candidates"
---

# SPHEREx pilot v1.0

**Question.** Does the archive-facing SGL pipeline survive a third,
spectral-image-shaped archive read from a space-based observer, and
what do 14 months of SPHEREx 0.75–5 µm imaging say about relays on the
Sun's focal lines toward three nearby stars that ZTF cannot reach?

**Answer.** The adapter interface (`sglsurvey/adapters/base.py`)
absorbed the multi-extension Level-2 product, the per-sample
wavelength axis and the S3 byte-range transport without change. Two
SPHEREx-specific pieces were required in the shared estimator — a
sub-pixel-phase matched filter for the under-sampled PSF and a static-
sky template in place of difference images — and each was worth
0.4 mag and ~3 mag of depth respectively. No relay candidate survives
on any of the three corridors. 90 %-recovery depths are m ≈ 19.3–20.8
AB in all six detectors for duty ≥ 0.5 and |µ| ≤ 1″/yr, i.e. reflectors
≳ 1.3–2.4 × 10⁵ km (albedo 0.1) or 700 K emitters ≳ 11–17 km at 550–
720 AU. The flux scale was verified on 13,657 catalogued-star
measurements to within ±0.2 mag in every detector before any depth was
quoted.

## 1. Setup

- Hypothesis freeze `surveys/spherex/hypotheses.md` v1.0: WISE v1.0
  physics (550–10,000 AU log-uniform, Rx/Tx, 99 % + 10″, duty ≥ 0.5,
  |µ| ≤ 1″/yr) with the Earth-centre observer (LEO offset ≤ 0.017″
  budgeted), detectors D1–D6 as reporting bands with the per-epoch
  wavelength carried on every sample, FLAGS fatal template 708343, and
  explicit reflected-light / hot-component interpretations.
- Stars (all ZTF-inaccessible, Dec < −28°): Lalande 21185 (nearest
  such star, clean high-latitude field, 349 exposures), GJ 687
  (corridor 3° from the south ecliptic pole inside the SPHEREx deep
  field, 7,308 exposures), σ Dra (engineering-backbone top pick, 643
  exposures). Registry entries reused unchanged (v1.5).
- Adapter `sglsurvey/adapters/irsa_spherex.py`: IRSA TAP discovery
  through the CAOM plane/artifact join (product URI + published MD5),
  ObsCore-polygon nominal footprint, slim cutouts assembled by HTTPS
  byte-range reads from the public S3 mirror (IMAGE/FLAGS/VARIANCE/
  ZODI sub-arrays, the PSF plane of the containing zone, the 9 × 9
  wavelength table; 260 kB kept per exposure), SIP + FLAGS exact
  footprint. Recon: `surveys/spherex/notes/irsa_recon.md`. Running
  from S3 left IRSA's bandwidth to the concurrent WISE and ZTF runs.

## 2. Coverage

13,488 public detector-exposures discovered (2025-05-24 → 2026-08-11),
16,978 coarse hits, 16,761 usable precise evaluations, 0 missing
products. Every corridor covers the full 550–10,000 AU range at every
visit (6 / 4 / 7 visit clusters for Lalande / GJ 687 / σ Dra) — the
3.5° field has no chip-gap geometry at the 6′ scale of a locus. In the
deep field the detectors observe in different seasons: at the GJ 687
corridor D1/D4 epochs fall almost entirely in one parallax phase and
D3/D6 in the other. Tables: `surveys/spherex/results/pilot_v1_summary.md`.

## 3. Search and calibration

Layer 1: 1.26 million single-epoch 5σ peaks, 50,846 within 10″ of a
track, 1,290 position clusters, 1,032 of them static (≥ 3 epochs or a
CatWISE2020 counterpart); the 258 unresolved clusters are 1–3-
detection singletons at one or two epochs and none follows the
parallax track. Layer 2: exposure-PSF matched-filter photometry at
2 × 2 sub-pixel phases, static-sky template subtracted (per corridor ×
detector, linear in λ per 3″ sky node), 96-node 1/z grid × 3 × 3 µ
grid, 8 offset controls, phase-split stacks. Thresholds T = 1.6–5.9.
Eight real-track maxima exceed T (≈ 4 expected from the empirical
FAR); seven fail the parallax-phase split and the eighth (σ Dra rx D4,
S = 4.3) is absent from the coincident Tx stack (S = 0.9) and from
every other detector, vetoed by the role-coincidence rule added at
this stage. **0 Candidate records survive.**

Positive control: 129 isolated 2MASS/CatWISE/Gaia stars measured
through the same estimator at the exposure wavelength; median Δm per
detector −0.20 … +0.17 mag (scatter 0.07–0.12), 0.05 mag all-band
offset applied conservatively. The same test read −0.35 … −0.50 mag
before the sub-pixel-phase estimator — the pilot's first lesson.

## 4. Constraints

288 Constraint records (3 endpoints × 2 roles × 6 detectors × 8 z
intervals), all `recovery_curve`. Median m90 (duty ≥ 0.5):

| endpoint | D1 | D2 | D3 | D4 | D5 | D6 |
|---|---|---|---|---|---|---|
| Lalande 21185 | 19.5 | 19.3–19.8 | 20.0–20.3 | 20.1–20.5 | 19.5 | 19.0–19.3 |
| GJ 687 | 19.4–20.2 | 19.1–19.7 | 19.7–19.9 | 20.3–20.8 | 19.5 | 20.2–20.7 |
| σ Dra | 19.8 | 19.5 | 20.1–20.3 | 19.9 | 20.0 | 18.9–19.3 |

Each record carries the wavelength range covered, the reflected-light
diameter (albedo 0.1, solar colour) and blackbody-equivalent diameters
at 400/700/1000 K at the interval's geometric-mean distance. At 550–
720 AU: D ≳ 1.3–2.4 × 10⁵ km reflected (D1–D4); ≳ 11–17 km at 700 K
and ≳ 4–6 km at 1000 K (D4–D6); ≳ 40–90 km at 400 K (D5–D6). The
deep-field corridor is not materially deeper than the all-sky ones
because its stack is confusion-limited after the template (per-cutout
noise 3–4× the pipeline model, static, mostly removed).

## 5. Lessons

1. Under-sampled PSFs need sub-pixel-phase evaluation of the matched
   filter; interpolating a pixel-centred map loses up to 0.5 mag.
2. Without difference images the stack is confusion-limited at
   16–17 AB; the static-sky template recovers ~3 mag. v2 should weight
   epochs by the template residual variance rather than the per-cutout
   MAD (11–18× the pipeline model in the deep field).
3. Deep-field seasons are detector-dependent, so add a cross-detector
   phase test (D1/D4 vs D3/D6 at the same cell).
4. Rx/Tx coincidence is a free veto for linear endpoints (added).
5. Raise the control count (8 → 16) or use the empirical cross-cell
   null: 8 exceedances against ≈ 4 expected.
6. S3 byte-range reads (≈ 20 MB/s, 0.15 s latency) make SPHEREx the
   one archive whose scale-up is not IRSA-download-bound; the whole
   3-star pilot fetched 8,505 cutouts in under two hours.

## 6. Next

Remaining 12 ZTF-inaccessible corridors (`spherex_corridors.QUEUE`,
including the component endpoints of 61 Cyg, Struve 2398, Groombridge
34 and GJ 338), then the full universal list; re-run as each quick
release extends the baseline (a third parallax phase arrives with the
2026-11 data). Records: `runs/spherex/{coarse_v1,precise_v1,screen_v1,
control_v1,calib_v1}`; AnalysisRun `run-55623bb1127a`.

---

## Status note (2026-08-22) — exploratory, pending v2

The WISE scientific review of 2026-08-21
(`surveys/wise/scientific_review.md`) applies to this report: the
16-offset control maximum is a per-search rank statistic (a noise-only
cell exceeds it with probability ≈ 1/17, so "FAR < 1/16" is a rank
statement, not a survey-wide false-alarm rate); injections were
analytic and tensor-level, so the quoted "90 %-recovery" depths are
*threshold sensitivity* (`completeness_kind = threshold`) without
confidence intervals; the phase / cross-detector season, template-
coverage, bright-static-neighbour and role-coincidence vetoes are
heuristic review rules with unmeasured selection functions; the locus
was evaluated at its nominal position. No number is recomputed here.
The defensible conclusion is *no compelling candidate after heuristic
review*. The calibrated version is the SPHEREx v2 survey (project plan
§10.1 step 6, last in the order because the next quick release adds a
parallax phase).
