# TESS crossings coverage gate v1 — PASSED (pre-freeze, metadata only)

Run 2026-08-24 (`scripts/coverage_gate.py`); **corrected same day**: the
first pass inherited a year typo in the HEASARC sector table (s046 end
date "2022 Dec 30" for 2021), which inflated sector 46 to 393 days and
falsely claimed the 2022-spring grazes of van-maanen (b = 0.078 R☉)
and gj-1276 as covered. The script now clamps implausible sector spans
(24–32 d expected) to start + 28 d and flags them (s046, s097
clamped); every number below is from the corrected run. Inputs: the
TESS-spacecraft crossing list `crossings/tess_v1` (`xng-a943f0f3dbe4`
— mandatory: Earth-center misplaces TESS-frame grazing events by up
to 0.33 R☉ / 3 h), TESScut sector lookups per (channel, target) and
the HEASARC sector table, snapshotted under `runs/tess-crossings/`.
Day-level precision; the survey stage refines with FFI timestamps and
the ~1 d mid-sector downlink gap. No pixels touched.

## Observed coverage (sectors completed by 2026-08-24, corrected)

| channel / rung | events covered | targets | full-window |
|---|---|---|---|
| **B 1.2 R☉ (grazing)** | **2** | 2 | **2 — complete ingress→egress** |
| **B 2.5 R☉ (grazing)** | **2** | 2 | **2** |
| B 0.1 AU | 3 | 3 | 1 |
| A 0.1 AU | 5 | 4 | 2 |
| A 1.0 AU | 292 | 86 | 19 |

All narrow-rung rows come from the **Ecliptic sectors** (42–44, 71,
91) at 600 s or 200 s FFI cadence — hundreds to thousands of samples
per window, versus 1–6 epochs in every previous survey.

**The grazing table (the cell no other archive can reach):**

| target | b | t_ca | sector (cadence) | windows |
|---|---|---|---|---|
| wolf-359 | 0.55 R☉ | 2021-09-05 | s42 (600 s) | 0.6 d + 1.3 d, both full; the 11.3 d wide window full as well |
| teegarden | 1.04 R☉ | 2025-05-06 | s91 (200 s) | 0.3 d + 1.3 d, both full; wide window 7.5 d partial |

Channel A narrow rung adds two complete windows of near-grazing
on-star geometry: gj-1276 (b = 0.97 R☉, s42) and van-maanen
(b = 0.52 R☉, s43), each a full ~11.5 d light curve of the star at
600 s cadence, plus partial teegarden ×2 and a gj-908 edge. All
targets TESS-unsaturated (T ≈ 11–14).

**Genuinely uncovered, stated for the record:** the two deepest
grazes of the TESS era — van-maanen b = 0.078 R☉ (2022-04-03) and
gj-1276 b = 0.49 R☉ (2022-03-03) — fell between sectors (2022 spring,
no ecliptic pointing), as did van-maanen's b = 0.017 R☉ graze of
2019-04 (pre-ecliptic-campaign). The ecliptic campaigns are the sole
supply of narrow-rung coverage; future ones (and s92's scheduled
wolf-359 antipode coverage) extend it.

Artifact notes: TESScut returned one bogus sector id ("1751") for the
ross-154 antipode — no date entry, ignored; A 1.0 AU remains covered
everywhere (292 events) with its constraint-only status carried over.

## Verdict

**Proceed to the survey chain.** Two complete grazing-crossing light
curves (one at 200 s cadence) plus two complete near-grazing on-star
windows are data no other archive holds, and the pulse-period cell
(2 × cadence → window length) opens for the first time. Design
requirements for the freeze: TESScut-based adapter (new);
chord-profile matched filter in time as the primary statistic; ring
controls rescaled to arcminutes (21″ pixels); flare discrimination on
channel A; windows and b from the spacecraft-frame list only; depth
expectation T ≈ 17–18 per stacked window (kW-scale grazing-cone power
limits — the value is temporal structure, not depth).
