# WISE beam-crossings survey (Pipeline B) — report

Third Pipeline-B (§3.5) survey: WISE/NEOWISE L1b single exposures
against Earth's crossings of hypothesized Sun–star relay beam axes,
on the spacecraft-observer crossing list `crossings/wise_v1` —
the first crossings survey whose event geometry and archive share the
same observer. Executed 2026-08-24. **Outcome: a structural null — 0
searchable units, 0 trials, no signal statistic was ever formed.**
The survey's product is the pair of gates that force that outcome,
plus the covered-window inventory they leave behind.

## 1. Construction and provenance chain

1. **Spacecraft crossing list** `crossings/wise_v1`
   (`xng-ffc922f7bbca`): 5,170 b(t) minima, 2010–2024, WISE
   L1b-header spacecraft observer (exact at frame epochs).
2. **Hypothesis freeze v1.0** `surveys/wise-crossings/hypotheses.md`
   (+ two documented pre-search amendments: the qa_status quality
   value, the v1.1 control-pool rescale).
3. **Elongation gate** (`configs/elongation_gate_v1.json`) —
   computed pre-freeze, scope-defining.
4. **Coverage intersection** (`results/coverage_v1.*`): 45 TAP cones,
   91 covered events, per-band epoch inventories.
5. **Saturation cut** (`results/saturation_cut_v1.*`): W1/W2
   searchable population = 11 faint targets (white dwarfs, T/Y
   dwarfs, latest-M); every ordinary M dwarf within 10 pc is
   brighter than the W1/W2 single-frame saturation limits.
6. **Threshold freeze v1.0 → v1.1**
   (`configs/threshold_freeze_v1.json`, `thresholds.md`): 0
   searchable units; 30 constraint-only unit-rows.

## 2. Gate I — the elongation theorem

WISE observes only on the great circle at solar elongation ≈ 90°.
Every beam-crossing geometry places the beam source near the
observer–Sun axis during the window (elongation ≈ 180° for
interception, ≈ 0° for the sunward combinations). Over the frozen
event list:

| channel / rung | in-era events | observable in-window |
|---|---|---|
| B 1.2 R☉ / 2.5 R☉ / 0.1 AU | 50 / 67 / 102 | **0 / 0 / 0** |
| A 0.1 AU | 101 | **0** |
| A 1.0 AU | 1,272 | 102 (45 targets) |

The workhorse channel of the ZTF and PS1 crossings surveys — downlink
pre-lens interception at the antipode — is *invisible in principle*
to a fixed-attitude terminator-orbit surveyor, at every beam radius.
Only the widest uplink rung survives, through events at large impact
parameter or high ecliptic latitude where the crossing time decouples
from opposition.

## 3. Gate II — the control-geometry theorem

The surviving rung's windows are ~116 d long, recurring annually with
~249 d gaps (and merging into ~232 d double windows in some years).
Under the frozen discipline (8 designated temporal pseudo-window
controls; a valid offset must move every window clear of every
same-rung window *and* land ≥ 2 shifted windows on primary epochs):

- the ZTF/PS1 offset pool (±23–97 d) cannot clear a 116 d window at
  all — 0 units (freeze v1.0);
- a pool rescaled to the window length (±130–245 d, amendment v1.1,
  pre-pixel) reaches only 2–6 valid offsets per unit — epoch support
  fails against WISE's two ~5 d visit clumps per year, and the
  merged double windows admit zero offsets at any scale — still 0
  units.

The 8-control standard was retained rather than weakened per-survey.
Result: 18 single-window + 12 gate-blocked unit-rows, all
constraint-only; **zero trials, zero exceedance budget, no pixel ever
searched.** This closes, from a third independent direction, the
finding both optical surveys made about the 1.0 AU rung: a d = 1
wide-beam hypothesis whose windows tile the calendar cannot be
temporally controlled by the frozen pseudo-window family. A future
searchable WISE crossings experiment needs a v2-style null ensemble
(time scrambling / trajectory randomization with a calibrated global
error rate), pre-registered from scratch.

## 4. What survives: the covered-window inventory

91 events on 43 targets have in-window primary W1+W2 frames
(W3/W4: 5/4 cryo-era events), with 211–343-epoch off-window baselines
per band — tabulated per event in `results/coverage_v1_events.ecsv`
with per-band epoch lists in `epochs_v1.json`. These rows enter the
joint-stage covered-window ledger (plan §11.2 step 6) as
*coverage-without-statistic* entries: cross-archive coincidence
bookkeeping can still use them (a candidate from another archive
whose window WISE covered can be checked by targeted forced
photometry as follow-up, outside this survey's frozen scope).

**Not constrained:** everything. No calibrated sensitivity, exclusion
or reference depth is claimed for any cell; the mid-IR (3–5 µm)
crossing regime remains untested, now with the design requirements
for testing it on record.

## 5. Lessons for the remaining crossings steps

1. **Check the attitude law before the cadence.** The plan's step-5
   gate asked whether NEOWISE's 6-month cadence supports the window
   durations; the elongation constraint is stricter and categorical.
   Fixed-elongation surveyors (WISE, and any similar IR mission)
   cannot do interception-channel crossings work; TESS-like
   anti-solar pointing can (its list exists: `crossings/tess_v1`).
2. Saturation, not depth, bounds blended IR searches on the 10 pc
   sample: only WDs and T/Y dwarfs are photometrically accessible in
   W1/W2 at L1b frame level.
3. Pre-pixel amendments are cheap and honest: both freeze defects
   (nonexistent qa_status value; mis-scaled control pool) were found
   and documented by frozen-input computation before any image was
   touched — the discipline's ordering (coverage → cut → freeze →
   pixels) did its job.

## 6. Products

`surveys/wise-crossings/{hypotheses.md, thresholds.md, configs/*,
results/*, scripts/*}`; spacecraft crossing list
`crossings/wise_v1/`; TAP snapshots under `runs/wise-crossings/`
(71 snapshots, ~70 k frame rows). No image products were fetched.
