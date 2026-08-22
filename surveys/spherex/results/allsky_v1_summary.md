# SPHEREx all-sky run v1 — 77 corridors / 88 endpoints (2026-08-21)

Hypotheses `spherex-hypotheses-v1.0` + amendments 1–3 (v3 estimator:
sub-pixel-phase exposure-PSF filter, static template, template-
calibrated variances, joint six-detector stack, 16 controls, five
adjudication vetoes). Registry v1.5. Adds the 62 universal-list
corridors (69 endpoints) at Dec > −28° to the 15 southern ones.
Calibration `runs/spherex/calib_v4`, AnalysisRun `run-ddd97a27faad`.
Data under `runs/spherex/` (39 GB; 33,879 slim cutouts, 0 missing).

## Coverage
60,426 detector-exposures discovered; 76,065 precise evaluations:
70,422 usable, 4,340 partial, 1,303 unusable. Every corridor covers
550–10,000 AU at every visit.

## Controls
Flux scale on 208,541 catalogued-star measurements: median Δm per
detector −0.19 / 0.00 / +0.04 / +0.10 / +0.07 / +0.14 (D1–D6);
0.053 mag applied. Layer-1 screen over 77 corridors: no unresolved
cluster with > 2 epochs. 1,226 searches × 16 controls: **64
exceedances (5.2 %), all vetoed** — 56 phase split, 4 role
coincidence, 3 cross-detector season complement, 1 template coverage.
**0 Candidates.**

## Constraints
**9,808 records** (88 endpoints × 2 roles × 7 stacks × 8 z intervals;
all `recovery_curve`). Joint-stack m90 (duty ≥ 0.5, |µ| ≤ 1″/yr):
median 20.75 AB over the 88 rx stacks, range 16.9–21.8.

Deepest: GJ 625 21.8, 82 Eri 21.6, GJ 1061 21.6, GJ 66 A/B 21.5,
GJ 338 A 21.5. Shallowest: α Cen A/B 16.9–17.0 and GJ 11068 17.9
(Galactic-plane corridors at |b| ≲ 2°), Proxima 19.0, GJ 667 C 19.4,
GJ 13157 19.4. The per-endpoint table is in
`runs/spherex/calib_v4/threshold_report.json` / `m90_curves.npz`.

Physical limits at 550–720 AU on a median corridor (m90 ≈ 20.7):
grey reflectors ≳ 1.3 × 10⁵ km (albedo 0.1), 700 K emitters ≳ 17 km,
1000 K ≳ 5 km; the α Cen corridor is ×5–6 worse.

## Caveats
Single-season or season-weighted corridors (Struve 2398 and others
with phase-1 counts < 10) rely on the cross-detector complement or
remain single-season until the next quick release. Binary components
share corridors. Depths are confusion-set: Galactic-plane corridors
are 1.5–4 mag shallower than high-latitude ones regardless of epoch
count.
