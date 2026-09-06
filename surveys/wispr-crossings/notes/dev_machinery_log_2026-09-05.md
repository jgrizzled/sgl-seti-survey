# Dev machinery log — PSP/WISPR sunward survey (2026-09-05)

Chronological record of the dev stage (freeze v1.0 → amendment v1.1),
the positive controls, and the machinery facts the report draws on.
Scripts: `scripts/series.py` (index / run / reduce), `wispr_lib.py`
(chain), `psp_geometry.py`, `moving_body_control.py`, `adjudicate.py`,
`completeness.py`, `summarize.py`.

## 1. Pre-dev amendment v1.1 (before any survey frame was measured)

Recorded in `thresholds.md`: affine WCS refinement (a 0.3° residual
rotation of the header solution on the 2024-12-27 recon frame), and
calibrator S/N ≥ 15 / ≥ 30 calibrators (linear ZP surface below 60) /
MAD gate 0.35 — the ZP shifts by 0.35 mag between S/N cuts of 15 and 8
on a perihelion frame (noise-boosted selection of V 6–8 stars at
15-s exposures).

## 2. Machinery facts found at the dev launch (no threshold content)

- **E2 (orbit 02) L3 frames are unbinned** (1920 × 2048, 15,750,720
  bytes; every other orbit is 2 × 2 binned, 3,954,240 bytes). The
  chain 2 × 2-averages them on load (MSB is a per-pixel surface
  brightness, so the mean reproduces the on-board "divided by 4"
  binned product) and slices the WCS; astrometry on a binned E2 frame:
  829 matches, 0.63 px rms, identity affine. The first launch treated
  the 15.7-MB payloads as bad and retried each three times, which
  looked like a server stall.
- **Fetcher**: `requests` with (15 s, 45 s) connect/read timeouts and a
  90-s streamed-body cap; the first launch's `urllib` fetch hung in a
  TLS handshake with no timeout path. Failed fetches are written as
  `<frame>_missing` and retried by a rerun. Worker pool uses the
  `forkserver` start method (fetch threads and forked workers do not
  mix safely).
- **ZP0 by processing version**: −21.79 (E1, V3, NSUMEXP 8, 20.5 s),
  −21.76 (E2, V3), −21.35 (E22, V1) — a 0.4-mag calibration-factor
  difference between the V3 and V1 products, absorbed by the per-frame
  ZP; the colour system is re-measured from the calibrator dump at
  reduce.
- **Vesta is the moving-body control**: it transits WISPR-I in 24 of
  29 encounters at V 6.7–7.6 from PSP (`results/moving_body_scan_v1.json`,
  Horizons `@-96`; Ceres at V 7.8–8.0 in ten); E22 (58 h, V 7.05–7.15,
  2024-12-21 → 23) and E10 (17 h, V 6.76) are the dev passages.
- A process-management note for the record: `pkill -f "<pattern>"`
  issued from a shell whose own command line contains the pattern
  kills that shell — two edits and two debug runs were lost that way
  before the pattern was anchored (`"[s]eries.py"`).

## 3. Dev findings → amendments v1.2–v1.5 (`thresholds.md`)

Partial dev reduce (6k frames): Tycho-2 extract lacks bright and
high-PM stars (61 Vir, gj-908) → template = Tycho ∪ Hipparcos + S2
star table (v1.2); stored series rounded to zero at the 1e-13 scale →
flux unit 1e-12 (v1.2); frame-wide bright-body veto flagged 36 % of
arc epochs with no distance-independent effect except Mercury/Venus
within 10° → proximity veto (v1.3); 61-vir S2 reads 0.88–0.95 of its
template → bright-star class with self-calibration (v1.3); control-
patch ensemble reads a + b·T with b ≈ 0.64–0.7 → template response
(v1.4). Full dev reduce (24.6k frames): the ZP MAD gate rejected 87 %
of frames inside 0.06 AU → uncertainty gate (v1.5); single-frame
particle hits set S_pulse control thresholds at 57–410σ → pair-min
persistence in the statistic (v1.5). 5,310 gate-rejected frames were
re-measured under v1.5 (the reduce keeps the last record per frame).

## 4. Dev result (v1.5; `results/dev_v1.md`)

5 units, 15 trials, **1 exceedance vs 1.67 expected** — 61-vir S2
S_event, adjudicated bright-star calibration drift (retained,
non-promotable). gj-908 S2 template null passes (z median +0.3).
Vesta E22 recovered at −0.19 ± 0.32 mag. Declared for completeness:
colour coefficient 0.286, flux-scale systematic ±0.3 mag (magnitude-
dependent ZP residual), pulse cell systematics-limited (pair-min T
32–70 on dev units).
