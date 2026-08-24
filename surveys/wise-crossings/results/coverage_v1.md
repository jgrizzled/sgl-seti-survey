# WISE crossings coverage intersection v1 — results

Run 2026-08-24 against hypotheses.md freeze v1.0 (+ the coverage-stage
quality amendment). Inputs: `crossings/wise_v1` (`xng-ffc922f7bbca`,
WISE-spacecraft observer); IRSA TAP discovery snapshots under
`runs/wise-crossings/` (45 cones — the elongation-pre-gate targets;
the other 34 targets are recorded `elongation_null` without queries,
enforced by the archive's own pointing law; 71 snapshots incl.
retries; one transient malformed TAP row triggered the resume path).
A frame covers an event if its t_mid lies in t_ca ± √(r²−b²)/v⊥ and
its nominal footprint contains the per-epoch propagated star
position. "Primary" = `qual_frame > 0`, `moon_sep ≥ 15°` (the frozen
`qa_status = A` gate named a nonexistent value — archive uses
Prelim/Reviewed/Final — and was dropped by the documented pre-search
amendment).

## Channel A, 1.0 AU rung (the sole in-scope channel)

| population | events |
|---|---|
| in-era total | 1,272 |
| elongation-null (34 targets, not queried) | 602 |
| queried (45 targets) | 670 |
| covered (≥ 1 primary in-window epoch) | **91** |
| per band | W1 91 / W2 91 / W3 5 / W4 4 |

Every covered event has both W1 and W2 in-window frames (bands ride
the same scan); W3/W4 exist only for the handful of events reaching
back into the 2010 cryogenic months. Off-window baselines are rich:
211–343 primary epochs per (target, band) for the multi-window
targets — WISE's fixed-elongation revisits deliver the baseline the
optical surveys lacked, even though (thresholds.md) the same geometry
defeats the temporal-control family.

Per-(target, band) primary epoch MJD lists are in `epochs_v1.json`
(the threshold freeze evaluates pseudo-window support offline from
them). qa_status distribution across queried frames: Prelim 49,876 /
Reviewed 13,138 / Final 6,357.

Products: `coverage_v1_events.ecsv` (per event, per-band in-window
epoch/visit counts and era totals, status column), `epochs_v1.json`,
`coverage_v1_summary.json`, `coverage_partial_v1.jsonl` (resume
ledger).
