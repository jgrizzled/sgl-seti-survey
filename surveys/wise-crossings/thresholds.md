# WISE crossings threshold freeze v1.0 → v1.1

Frozen 2026-08-24 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json` (v1.1,
`sha256:620beeea…5abb`), bound to content hashes of `hypotheses.md`,
`configs/elongation_gate_v1.json`, `results/coverage_v1_events.ecsv`,
`results/epochs_v1.json`, and `results/saturation_cut_v1.ecsv`. No
pixel data was touched at any point. Single channel by scope
(hypotheses v1.0): A, 1.0 AU rung; the ZTF/PS1 final constructions
adopted (S > max(T, 0), empirical variance rescale floor 1,
WEIGHT_CAP 20×, ≤ 120-epoch baseline clip, 8 designated temporal
pseudo-window controls, seed 20260825, dev fraction 0.30).

## Amendment v1.1 (pre-pixel)

The v1.0 offset pool (±23–97 d, redraws ±113/127 — designed for the
0.6–11.6 d grazing/0.1 AU windows) cannot clear this rung's ~116 d
windows: an offset must exceed the window length yet stay inside the
~249 d inter-window gap (annual same-side recurrence), so v1.0
yielded 0 searchable units by construction. v1.1 rescales the pool to
the rung's window length — a frozen deterministic scan ±130, ±135, …
±245 in 5 d steps, first 8 offsets satisfying the *unchanged*
validity rule (no overlap with any same-rung window; ≥ 2
epoch-supported shifted windows). The amendment generalizes the
redraw mechanism; nothing about the statistic, gates, or rule
changed.

## Outcome: 0 searchable units — the survey is constraint-only

| class | n | detail |
|---|---|---|
| searchable units | **0** | — |
| single-window (constraint-only) | 18 | 10 targets, one covered window each |
| gate-blocked (constraint-only) | 12 | 6 targets × W1/W2, 2–3 covered windows, 211–343 off-window epochs — but 0–6 valid offsets of the required 8 |

Even under the rescaled pool the achievable offset count is 2–6:
epoch support fails (WISE visits are two ~5 d clumps per year, and a
shifted 116 d window catches one only part of the time), and targets
with merged double-length (~232 d) windows (van-maanen, gj-1276,
gj-518) admit **zero** geometrically valid offsets at any scale. The
8-control standard is retained rather than weakened per-survey —
lowering it after seeing which targets would become searchable is
exactly the data-driven design drift the frozen discipline exists to
prevent.

Combined with the elongation gate this is the survey's result: for a
fixed-attitude elongation-90° surveyor, every crossing channel is
either geometrically invisible (B, A-0.1) or — via the one surviving
rung's annual ~116 d windows — impossible to control with the frozen
temporal-pseudo-window family. The designated future route to a
searchable WISE crossings experiment is a v2-style null ensemble
(time scrambling / trajectory randomization with a calibrated global
error rate), pre-registered from scratch; the covered-window
inventory produced here (91 events with in-window primary frames)
feeds the joint-stage coincidence ledger regardless.
