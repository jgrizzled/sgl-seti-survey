---
title: "Pilot corridor depth-of-coverage in the WISE merged L1b inventory"
date: 2026-08-18
status: "Proxy measurement from live IRSA TAP queries, 2026-08-18"
---

# Pilot corridor depth-of-coverage (proxy)

Frame-band counts from `wise.neowiser_merge_p1bm_frm` (union of all four
mission phases, verified in `irsa_recon.md`) with frame center within
0.55° of each system's fixed J2000 antipode point — 0.55° ≈ the L1b
frame half-diagonal, so this approximates "frames plausibly covering the
antipode." Query:

```sql
SELECT band, COUNT(*) AS n, MIN(mjd_obs) AS mjd_min, MAX(mjd_obs) AS mjd_max
FROM wise.neowiser_merge_p1bm_frm
WHERE CONTAINS(POINT('ICRS',crval1,crval2), CIRCLE('ICRS',<ra>,<dec>,0.55))=1
GROUP BY band ORDER BY band
```

| Corridor (antipode) | W1 | W2 | W3 | W4 | W1/W2 span (MJD) |
| --- | --- | --- | --- | --- | --- |
| Barnard's Star (89.45, −4.69) | 555 | 555 | 40 | 21 | 55269–60368 (2010-03 → 2024-02) |
| Alpha Cen AB (39.90, +60.83) | 625 | 625 | 50 | 24 | 55239–60339 |
| Sirius AB (281.29, +16.72) | 581 | 582 | 22 | 22 | 55289–60406 |
| Ross 154 (102.45, +23.84) | 472 | 472 | 22 | 22 | 55279–60379 |
| Lalande 21185 (345.83, −35.97) | 561 | 564 | 23 | 22 | 55339–60453 (→ 2024-05) |

## Reading

- **Every pilot corridor is viable**: 470–630 W1/W2 single-exposure
  epochs spanning ~14 years, in the usual WISE cadence of ~1-day visit
  clusters every ~6 months — plenty of parallax cycles for the SGL track
  model.
- **W3/W4 (waste-heat-sensitive bands) exist for all five corridors**
  (2010 cryo phase only, 21–50 frames each) — the band-specific
  interpretation split in the hypothesis freeze is exercised everywhere.
- Counts are within a factor ~1.3 of each other; no corridor is
  coverage-starved, so target prioritization can rest on endpoint
  quality and background rather than raw epoch count.

## Caveats

- Proxy only: the true corridor is not a fixed point — it sweeps an
  annual-parallax ellipse (up to ~375" at 550 AU) plus secular drift and
  spans the relay-distance range. Real discovery uses sglseti envelopes
  per epoch, not this cone. At these count levels the conclusion
  (viability) is robust.
- Antipode positions computed from catalog star positions without
  propagating proper motion of the target (irrelevant at 0.55° scale for
  a viability check; Barnard's moves ~2.5' over the span).
- Counts are frame-band inventory rows, not usable-pixel coverage.
