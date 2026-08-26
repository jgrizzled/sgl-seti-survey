# PTF crossings dev stage v1 — results

Run 2026-08-26 under threshold freeze v1.0
(`scripts/dev_search.py`; products under `runs/ptf-crossings/`,
110/110 cutout pairs fetched, zero 404s — PTF serves every listed
product, unlike ZTF's ~1-in-8). Dev scope per freeze D8: gj-1276 A
0.1 AU R (the one searched dev trial), wolf-359 B 0.1 AU R
(single-epoch class), ross-128 A R (saturation-excluded exercise).

## Realizations fixed at dev (implementation-defined at freeze)

1. **Calibrator selection**: MAST DR2 mean stars restricted to the
   transform's dwarf locus (r−i ∈ [0, 0.8] for R; g−r ∈ [0.2, 1.2]
   for g), err ≤ 0.1, nDet ≥ 5, predicted mag ∈ [15.5, 19.5] (R) /
   [16.0, 20.0] (g). The unrestricted set inflated per-frame ZP
   scatter to ~0.25 (off-locus red stars); restriction drops a test
   frame from 0.252 to 0.047.
2. **All cutouts 256 px** (258″): the first-pass 128-px off-window
   cutouts held only ~2 calibrators at these high-galactic-latitude
   fields and failed the ≥ 5 gate everywhere.
3. **Asteroid control uses a dedicated calibrator cone** at the
   asteroid's position (it roams outside the field cone).
4. Star position per epoch = linear interpolation of the per-event
   star_icrs positions (parallax amplitude ≪ pixel at 4–5 pc).

## Search results (`dev_search_v1.json`)

| unit | outcome |
|---|---|
| wolf-359 B 0.1 AU R (single-epoch class, 0 trials) | machinery clean end-to-end: 1 in-window + 27/34 off-window usable, k = 1.0 (floor, n = 26), **S = 0.101 vs T = 2.099** (8/8 ring controls) — no reportable exceedance; nearest catalogued static source 8.2″ (no annotation) |
| gj-1276 A 0.1 AU R (1 trial) | **constraint-only under the frozen offset-validity gate**: 1 of 12 pseudo-window offsets valid (campaign cadence leaves ±23–127 d unpopulated; 40 usable epochs cluster in other seasons). S = −3.84 computed and reported (one-sided positive → null); no searchable trial |

Calibration chain across both fields: 72 frames calibrated, median
ZP scatter 0.085 mag, median 6 calibrators/frame; gate attrition
~10 % (B field) to ~47 % (A field, scatter-dominated — the gate
rejecting genuinely bad photcalflag-0 frames as designed).

**Trials ledger after dev**: the single dev trial resolved
constraint-only → 0 searched dev trials, 0 exceedances. Confirmatory
population: 10 trials (6 units, all channel B), expected control
crossings 1.11.

## Saturation gate (`saturation_verify_v1.json`)

- **R band verified, no amendment**: stars at R_pred 12.9–13.0 fire
  dmask bit 8/6 in the 5-px core; 14.06 and fainter are clean — the
  empirical boundary brackets the frozen sat ≈ 14.0 within the
  0.5 mag amendment trigger. E_R = 14.5 stays (conservative).
- **g band one-sided**: no field star bright enough to fire the bits
  was available (brightest sampled 13.76, clean); the g boundary is
  unbracketed from above. No searchable consequence — no confirmatory
  g-band A unit exists.
- **ross-128 excluded-class exercise**: bit 8 + bit 6 both fire at
  the star core (R_est 9.86) — **exclusion confirmed**.

## Positive control (`asteroid_control_v1.json`)

(8971) Leucocephala, V = 18.9, three exposures of one night at the
wolf-359 antipode field, SkyBoT positions (err ≤ 0.4″), identical
chain. **Position-locked recovery in all three; internal RMS
0.025 mag** across independently calibrated frames (≤ 0.1 gate met
on the chain-precision reading). Absolute offset +0.24 vs Horizons
V (18.893, snapshotted) + assumed V−R = 0.4 — within the combined
prediction systematics (H ± 0.2, color class ± 0.1, rotation phase);
the AB flux scale is independently pinned by the PS1-DR2 star
calibration. No chain defect. For a stricter absolute control at a
later stage: an asteroid with known taxonomy/rotation, or PS1
catalog photometry of the same object.

## Consequences for the confirmatory run

1. Machinery locked: fetch → dmask → star-calibrated matched filter
   → k rescale → S/T with ring controls all validated on real data.
2. The confirmatory run is **all channel B** (van-maanen B 2.5 g,
   B 0.1 g+R; ross-128 B 2.5 R, B 0.1 g+R): fetch, measure, reduce
   once under the frozen chain with the dev-fixed realizations.
3. Expect calibration-gate attrition of order 10–50 % per field —
   coverage epoch counts thin at the gate exactly as declared;
   single-event units may drop below the 2-epoch S_event gate at
   measure time (reported as constraint-only if so, per the frozen
   rule).
4. Completeness injections (per frozen spec) run with the
   confirmatory measure, sharing the stamp-response machinery.
