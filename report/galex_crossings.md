# GALEX beam-crossings survey (Pipeline B) — report

**Result: 0 candidates.** Blind confirmatory run over 5 band-units /
15 frozen trials (14 effective): **1 exceedance vs 1.56 expected
control crossings**, adjudicated retained-ambiguous (wolf-359 A NUV
burst — micro-flare leading interpretation; not promotable). The
survey opens the programme's UV pulse-period cell: 5 ms photon-level
searches of the 0.1 AU-rung crossing windows in 2007–2010, with the
first in-band harmonic of the 1064 nm hypothesis family (266 nm ∈
NUV) and kW-class time-averaged limits on coherent UV pulse trains.

## 1. Construction and provenance chain

Ninth crossings survey (plan §5.12, queue item 4). Substrate: the
MAST gPhoton photon database — direct photon-event retrieval, no
images (adapter `sglsurvey/adapters/mast_gphoton.py`; all queries
snapshotted under `runs/galex-crossings/`). Chain, all 2026-08-26 on
the dev machine: reachability recon → hypothesis freeze v1.0
(D1–D8; pre-freeze contact declared, remedy D1a) → **amendment
v1.1** (aspect gate `flag % 2 == 0`, the gPhoton PhotonPipe
convention; restores the 1,637 s gj-1276 B 2010 visit from the
pervasive flag-64 state, astrometry verified at 0.37″ median) →
coverage (47 snapshots; 5 in-window units confirmed; grazing rungs
structurally uncovered) → threshold freeze v1.0 (statistics S_rate /
S_burst / S_period; 8-control max rule, S > max(T, 0)) → dev stage
(pseudo-units only: calibration ZP_eff 19.592 ± 0.088, background
0.209 cts/s per 8″ aperture; segment gate → gj-1276 A NUV/FUV and
wolf-359 A FUV constraint-only; **amendment v1.2**: S_period on
U(0, 5 ms)-jittered times — the 200 Hz tick otherwise corrupts the
H-test grid; the in-situ positive control: all three statistics
detected the real 2009-03-25 wolf-359 flare) → blind confirmatory +
injections (`results/confirmatory_v1.json`,
`completeness_v1.json`). Universal list `crossings/universal_v1`
Earth-center (LEO budget 0.010 R☉); era 2003-06-07 → 2013-05-01.

## 2. Results

**Trials (15 frozen; 14 effective — ross-128 A FUV S_period
degraded at the frozen 10-photon control gate):**

| unit | S_rate | S_burst | S_period | verdict |
|---|---|---|---|---|
| gj-1276 B NUV (2007 + 2010) | 0.275 < 0.349 | 5.37 < 5.49 | 40.9 < 47.0 | clean |
| gj-1276 B FUV (2007) | 0.018 < 0.037 | 1.41 < 2.63 | 0 (gate) | clean |
| wolf-359 A NUV (2007) | 1.61 < 3.84 | **4.68 > 4.13** | 27.9 < 32.6 | 1 exceedance |
| ross-128 A NUV (2007) | 0.455 < 0.509 | 3.92 < 4.21 | 29.1 < 35.1 | clean |
| ross-128 A FUV (2007) | 0.0455 = T (tie, not >) | 2.51 < 2.64 | constraint-only | clean |

Constraint-only lanes (computed, labeled, zero trials): gj-1276 A
NUV/FUV and wolf-359 A FUV (dev segment gate); the D1a `forced_dev`
S_rate on gj-1276 B 2010 — measured 0.235 cts/s at the locus,
consistent with the 0.209 background (its pre-freeze recon contact
is cited here as frozen).

**The exceedance, adjudicated
(`results/adjudication_wolf359_v1.json`):** a 5-photon cluster in a
0.5 s boxcar (0.80 expected), z = 4.68 vs the control threshold
4.13. The frozen ladder: FRED-morphology test — unresolved at 5
photons (the 5 s light curve is otherwise flat; contrast 1.9 at bin
scale vs the dev flare templates' 3–21.6×); two-band discriminator —
powerless at this fluence (0.067 FUV photons expected in-window;
0 observed); strict ≤ 25′ re-run — unchanged; asteroid mechanism —
excluded by timescale (an MBA cannot modulate at 0.5 s); no
detector-fixed signature. The veto cannot positively fire, and no
second covered wolf-359 A window exists in-archive →
**retained-ambiguous, not promotable** (frozen rule: promotion
requires independent recurrence, never in-window evidence alone).
Leading interpretation: an ordinary wolf-359 (CN Leo) micro-flare —
the star's own off-window 97 s control segments reached z = 4.13
and its full visits z = 14–106; a 4.68 on-star burst sits inside
its demonstrated activity. Statistically the family closes at 1
exceedance vs 1.56 expected.

### 90 %-recovery depths (`results/completeness_v1.json`; 231 nm effective in NUV, 154 nm in FUV)

Persistent (window-locked, through the identical chain; NUV via the
measured ZP_eff 19.592, FUV via nominal ZP with the NUV-measured
0.49 mag aperture term transferred — a declared ±0.3 mag budget):

| unit | r90 (cts/s) | m90 (AB) |
|---|---|---|
| gj-1276 B NUV 2007 (109 s) | 0.20 | 21.3 |
| gj-1276 B NUV 2010 (1,637 s; D1a reference, labeled) | 0.12–0.2 | ≈ 21.9 |
| gj-1276 B FUV 2007 | 0.07 | ≈ 21.2 |
| ross-128 A NUV (110 s) | 0.12 | 21.9 |
| ross-128 A FUV | 0.04 | ≈ 21.8 |
| wolf-359 A NUV (variability-limited) | 2.5 | 18.6 |

Pulses (the new cell; 90 %): **4 photons** per pulse on the 1,637 s
unit at both 0.05 s and 0.5 s widths (5 on the 100 s-class units;
2 in FUV where the background is ~0). Coherent trains (80 %, 1,637 s
unit, orbital drift sampled ± 2.56 × 10⁻⁵): **60 photons total for
P = 0.5–50 s** (0.037 cts/s time-averaged ≈ NUV 23.2-equivalent);
120 at P = 0.05 s (drift smear, honestly carried by the frozen
injection model).

## 3. Physical interpretation (0.1 AU cone; monochromatic line via band effective width)

Transmitter-power limits for an isotropic-within-cone emitter at
z ≥ 550 AU radiating through the 0.1 AU-radius cone at Earth
(the PTF/ZTF convention):

- **Persistent relays:** ≳ 20–33 kW (NUV; gj-1276 B and ross-128 A,
  both bands ≈ 16–28 kW FUV) in the covered 2007 windows; ≈ 20 kW
  as the labeled D1a reference on the 2010 window; ≈ 415 kW at the
  flare-noise-limited wolf-359 A.
- **Pulsed transmitters (first UV pulse constraints in the
  programme):** single pulses ≳ 1.7 MW peak (0.5 s) / 17 MW (0.05 s)
  — ~0.8 MJ per pulse — through the cone in the covered windows;
  **coherent pulse trains at P = 0.5–50 s excluded above ≈ 6 kW
  time-averaged** on the 2010 window — the programme's deepest
  time-averaged crossing constraint in any band, and its first at
  sub-second cadence resolution.
- **266 nm (frequency-quadrupled Nd:YAG) falls in the NUV band** —
  these are the first limits of the programme covering a harmonic
  of the 1064 nm line family; 532/1064/1550 nm remain unconstrained
  here as everywhere.

Caveats, stated as frozen: coverage is 5 windows on 3 targets in
2007–2010 (one two-window recurrence pair, gj-1276 B); the grazing
rungs (1.2 / 2.5 R☉) are structurally uncovered by GALEX's visit
cadence; teegarden's antipode is archive-empty; gj-908 A is
rim-limited; d = 1 while geometry holds, with duty-cycle scaling
and schedules avoiding Earth-crossing windows declared
unconstrained.

## 4. Lessons for the next adapters

1. **Aspect/photon flag conventions must be read from the pipeline
   source, not sampled** — the recon's flag-0 sample was a 2007
   accident; flag 64 is the normal late-mission state
   (`flag % 2 == 0` is the rule; v1.1).
2. **Tick quantization corrupts photon-level periodicity
   statistics**: any archive with discretized time stamps needs
   de-quantization jitter in the frozen statistic (v1.2), or every
   frequency commensurate with the tick rate is a false coherence.
3. **The MCAT carries near-duplicate per-visit rows** — dedupe
   (3″) before isolation cuts or calibrator selection.
4. Petal-pattern observing fragments off-window data into 60–110 s
   dwells: segment-based temporal controls work, but band coverage
   differs per visit — **control supply is a per-band property**
   (wolf-359's long visits are NUV-only → its FUV lane died at the
   gate).
5. On-star units at flare stars are variability-limited, not
   sky-limited: the control-segment construction absorbs stellar
   activity into the threshold honestly (wolf-359 T_rate 3.84 vs
   sky-limited ~0.5), and sub-fluence bursts adjudicate to
   retained-ambiguous — a UV on-star channel wants a
   simultaneous-two-band fluence gate at design time, not
   adjudication time.

## 5. Products

`surveys/galex-crossings/`: hypotheses.md (v1.0 + v1.1 + v1.2),
thresholds.md (+ dev assessment), notes/galex_recon_2026-08-26.md,
results/{era_scope_v0, coverage_v1*, nonlinearity_cut_v1, dev_v1*,
jitter_validation_v0, confirmatory_v1, completeness_v1,
adjudication_wolf359_v1}.json/.md/.ecsv, scripts/{coverage_intersect,
freeze_thresholds, dev_stage, galexlib, flare_template,
confirmatory_search, completeness_supplement,
adjudicate_wolf359_burst}.py; `runs/galex-crossings/` (snapshots +
query records); adapter `sglsurvey/adapters/mast_gphoton.py`.
Covered-window rows for the joint ledger at the next §5.7 refresh;
the wolf-359 A NUV retained-ambiguous burst is recorded for any
future UV mission (UVEX-class) that can re-cover a wolf-359 window —
no in-archive recurrence test exists.
