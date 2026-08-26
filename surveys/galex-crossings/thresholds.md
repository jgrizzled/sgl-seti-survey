# GALEX crossings threshold freeze v1.0

Frozen 2026-08-26 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json` (`sha256:18346431…017e47`), bound
to content hashes of `hypotheses.md` (v1.0 + amendment v1.1) and the
coverage products (`coverage_v1_events.ecsv`, `_summary.json`,
`coverage_v1.md`). **No photon touched, no signal statistic
formed** — the only archive contact at this stage is the frozen
section-6 nonlinearity cut (MCAT NUV magnitudes of the three
channel-A stars; catalog annotation-class, snapshot-disciplined).

## Substrate

Direct photon-event retrieval (freeze D3): photon flag = 0, aspect
flag % 2 == 0 (v1.1), boresight ≤ 33′; r_ap = 8″ apertures (B: the
deduplicated z-grid aperture set per event at the visit epoch);
live time = gated aspect seconds. Rate → flux only via the in-visit
MCAT-star regression at the ≤ 0.2 mag gate — a dev-stage
deliverable (the recon's ×2.6 naive-rate discrepancy must close
before any depth is quoted).

## Nonlinearity cut (frozen rule §6)

`results/nonlinearity_cut_v1.json`: wolf-359 NUV 19.28 **ok**;
ross-128 NUV 21.18 **ok**; gj-1276 — no MCAT NUV source within 15″
of the event-epoch position → **ok by construction** (fainter than
any nonlinearity concern; the high-PM per-visit identification is a
dev item). No unit lost; the expectation stated at the hypothesis
freeze is confirmed by the rule, not assumed.

## Search units — 8 units × 3 statistics = 24 trials

Unit = (target, channel-rung, band); unit statistic = max over its
eligible covered events. All 0.1 AU rung.

| unit | events (usable s) | S_rate | S_burst | S_period |
|---|---|---|---|---|
| gj-1276 B NUV | 2007 (109) + 2010 (1,637) | 2007 only (D1a) | both | both |
| gj-1276 B FUV | 2007 (109) | ✓ | ✓ | ✓ |
| gj-1276 A NUV / FUV | 2007-09 (92 / 92) | ✓ | ✓ | ✓ |
| wolf-359 A NUV / FUV | 2007-03 (97 / 97) | ✓ | ✓ | ✓ |
| ross-128 A NUV / FUV | 2007-03 (110 / 110) | ✓ | ✓ | ✓ |

The D1a lane: S_rate on gj-1276 B 2010 is `forced_dev` —
computed and reported constraint-only, outside the trial family,
its pre-freeze contact cited at any adjudication. Ledger (zero
trials): grazing rungs (0/39 + 0/49 covered), gj-908 A
(rim-limited), teegarden/van-maanen/ross-154 (§ coverage_v1.md).

## Statistics (numeric, frozen)

- **S_rate** = N_aperture / t_live (counts s⁻¹); B takes the max
  over the z-family apertures, mirrored in controls.
- **S_burst** = max over widths w ∈ {0.05, 0.5, 5, 50 s} (boxcar,
  step w/2, usable seconds only) of
  (C_max(w) − λw) / √max(λw, 0.5), λ = N_tot/t_live of the same
  series.
- **S_period** = max H (de Jager H-test, m ≤ 20) over a geometric
  period grid P ∈ [0.02 s, t_span/3] at 5×-oversampled Fourier
  spacing Δ(1/P) = 1/(5 t_span). Gate: ≥ 10 in-aperture photons,
  else constraint-only. The LEO orbital phase smear is not
  corrected — injections carry it (freeze §4), so short-period
  completeness states the loss honestly.

## Controls and rule

8 controls per trial; threshold **T = max over the 8 controls**;
exceedance **S > max(T, 0)**; **24 trials → expected control
crossings 2.67** (exchangeable rate 1/9 per trial, the standing
crossings convention).

- **B (pseudo-positions):** 8 positions on the same visit at the
  locus boresight radius (same detector annulus), position angles
  locus + k·40° (k = 1..8), each carrying the translated z-family
  aperture pattern; a control within 30″ of an MCAT source brighter
  than NUV 21 or within 60″ of the locus segment rotates +5° until
  valid (frozen resolution rule).
- **A (temporal segments):** 8 contiguous usable segments of the
  unit's in-window duration, drawn from the star's off-window visits
  (outside every in-era 0.1 AU window of that target-channel),
  time-ordered / equally spaced / non-overlapping, seed 20260826.
  Fewer than 8 valid segments → constraint-only (assessed at dev:
  ross-128 21 visits, gj-1276 12, wolf-359 5 long visits — segment
  counts are ample in all three).

## Veto ladder (calibrated at search, frozen order)

(1) A-channel stellar-flare veto — FRED-morphology template test +
the two-band discriminator (a line is single-band; flares brighten
FUV+NUV together; FUV live on all three A units); (2) SkyBoT
known-object census at every burst exceedance (B antipodes in the
opposition asteroid stream; an MBA transits the 8″ aperture in
minutes); (3) detector-fixed (xi/eta) vs sky-fixed clustering +
NUV hotspot check; strict re-run at boresight ≤ 25′; (4) recurrence
— gj-1276 B holds the 2007/2010 pair; elsewhere surviving
exceedances are retained-ambiguous (promotion never on in-window
evidence alone).

## Dev stage (pseudo-units only, D8) — deliverables

Machinery on off-window visits at the unit positions: (i) formal
flag-64 astrometric verification (v1.1 item); (ii) live-time /
aperture / dead-time calibration closing the ×2.6 discrepancy, then
the MCAT-star ≤ 0.2 mag gate; (iii) control-validity census (B
rotation rule outcomes; A segment counts); (iv) flare census at the
three A stars feeding the veto template — if a flare is found,
S_burst must recover it (conditional positive control); (v)
end-to-end pseudo-unit runs of all three statistics with injection
completeness (persistent, boxcar, periodic trains with the orbital
smear model). Blind confirmatory follows only after every dev gate
passes.

## Dev-stage assessment (2026-08-26; results/dev_v1.md)

The "assessed at dev" items resolved, per the frozen rules:

- **A control-segment gate (< 8 → constraint-only):** wolf-359 A NUV
  51 ✓, ross-128 A NUV 35 ✓, ross-128 A FUV 25 ✓; **gj-1276 A NUV 4,
  gj-1276 A FUV 3, wolf-359 A FUV 0 → constraint-only** (the gj-1276
  star's only long visit, 2010-10-23, is rim-only at > 33′; wolf-359's
  long visits are NUV-only). Confirmatory family: **5 band-units ×
  3 statistics = 15 trials, expected control crossings 1.67**
  (was ≤ 24/2.67 at the freeze).
- **Amendment v1.2** (user-approved): S_period on tick-jittered times
  (the 5 ms quantization artifact; hypotheses.md v1.2).
- Gates: flag-64 astrometry **PASS** (16 source-blocks, median
  centroid offset 0.37″); calibration **PASS** (ZP_eff
  19.592 ± 0.088, 15 deduped MCAT calibrators; blank-sky background
  0.209 cts/s per 8″ aperture); B pseudo-position rule valid with
  zero rotations on 3 blocks; B pseudo-unit null-clean on all three
  statistics; **the conditional flare positive control fired** — the
  A pseudo-unit segment landed on the real 2009-03-25 wolf-359 flare
  and all three statistics detected it through the full chain; FRED
  veto templates extracted (rise ≤ 10–20 s, decay e-fold 30–40 s,
  contrast up to 21.6×; one slow-rise event noted — the two-band
  discriminator carries non-FRED flares).
- Injection machinery calibrated on the 1,381 s pseudo-unit:
  persistent 90 % recovery between 0.05–0.1 cts/s (≈ NUV 22.3);
  0.5 s pulses at ~5 photons; trains recover to P = 10 s under the
  jittered statistic. Known property, stated for the record: the
  wolf-359 A thresholds inherit stellar-variability (flare) noise
  from the control segments — conservative by construction.
