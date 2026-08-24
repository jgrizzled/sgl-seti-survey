# Confirmatory run v1 — results

Run 2026-08-23/24 on the confirmatory split: channel B under freeze
v1.0, channel A under v1.0+v1.1+v1.2. Data: 8,818/9,590 planned epochs
complete (772 archive 404s and 4 transient 503s, counted; 770 MB under
`runs/ztf-crossings/products/conf/`). Statistics:
`confirmatory_v1.json`. One implementation correction over the dev
scripts, declared: channel-A photometry and cutouts use per-epoch
PROPAGATED star positions (PM + parallax displacement from registry
astrometry); dev's fixed-position shortcut is noted as a caveat on the
dev diagnostics (its qualitative conclusions — offsets, brightness
systematics — are position-independent).

## Headline

**No candidates.** 36 searchable units (3 A + 33 B) produced **3
exceedances against 4.0 expected control crossings** (36/9) — exactly
the frozen null rate. All three adjudicated as control-crossing
statistics (below); every other searchable unit is null. All
constraint-only and single-epoch rows are reported with their S values;
none dropped.

## Channel A

- Searchable (gate + 8 offsets + ≥2 windows): gj-1276 g 0.1 AU
  (S = −2.39), wolf-359 g 0.1 AU (S = −0.67), wolf-1069 g 1.0 AU
  (S = +0.47). All null.
- Constraint-only: 18 units + 6 singles, as the dev-validated gates
  dictated (2–4 offsets on the 1.0-AU rung; |median control| > 1 for
  the brighter/marginal targets — teegarden shows the familiar
  bright-star negative systematic, S ≈ −18…−41 with medC −1.5…−10).
  k factors span 1.0 (wise-0855, genuinely empty field) to 3.2×10⁴.

## Channel B

33 searchable units across gj-908, gj-1276, ross-128, teegarden,
van-maanen; typical S ≈ 0.0–1.4, T ≈ 1.1–94. Exceedances and
adjudication:

| unit | S | T | n_ep | adjudication |
|---|---|---|---|---|
| gj-908 g evt-1ff54c0eec27 | 1.33 | 1.09 | 2 | amplitude 1.3σ; no recurrence (4 other events, max R = 0.48) → control crossing |
| gj-1276 g evt-57e0c494a587 | 0.59 | 0.28 | 3 | sub-1σ amplitude; exceeds only via freakishly low T → control underdispersion |
| gj-1276 g evt-9c6a08c6df77 | 2.23 | 2.05 | 2 | 2.2σ; no amplitude-consistent recurrence (other 4 events R ≤ 0.77) → control crossing |

Two gj-1276 g events exceeding out of 6 at ~1/9 each has null
probability ≈ 0.11 — unremarkable. Frozen promotion requires
recurrence on the recomputed track, which none shows. MPC census
lookups for the in-window epochs remain open as a formality (the
predicted ~″/day rates cannot match catalogued movers per plan §3.6);
disposition: **vetoed, not candidates, retained in the record.**

Grazing rung (fully confirmatory, gj-1276 + teegarden): all
single-epoch class as expected from coverage; |S| ≤ 2.1 except one
r-band artifact (below). Constraint-only statements stand.

## Anomalies retained (not exceedances)

- teegarden r evt-d6a207b21d89: R = 1.14 but only 5/8 finite ring
  controls → invalid unit under the freeze; reported as
  `controls_incomplete`, not searched.
- gj-1276 r shows two large negative outliers (S = −30.7 at
  evt-57e0…, −47.7 in the grazing singles) — a bad difference-image
  epoch (subtraction artifact); flagged for the data-quality appendix.
- gj-1276/gj-2012/wise-0855 i-band and no-template singles report raw
  (unstandardized) S — constraint-only by construction.

## Remaining programme steps

Injection-calibrated completeness (frozen parameters: 100/event-window
cell, chord profiles, 532/650 nm lines) and the survey report with
depth statements per hypothesis rung.
