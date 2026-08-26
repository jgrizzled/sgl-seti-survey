# PTF crossings hypothesis freeze v1.0 (Pipeline B)

Drafted 2026-08-26 after the reachability recon
(`notes/ptf_recon_2026-08-26.md`) and the era scoping
(`results/era_scope_v0.json`), before any in-window pixel or
photometric data is touched at survey positions (pre-freeze contact
declared in §10). **FROZEN 2026-08-26: the user approved every §11
decision as recommended (D1–D8).** Sixth crossings survey; the construction is the PS1
crossings freeze (`surveys/ps1-crossings/hypotheses.md` v1.0 — the
warp-direct-substrate variant of the ZTF construction) with PTF
substitutions, each marked **[PTF]**.

Input: `crossings/universal_v1` (`xng-a09e2db7681d`, Earth-center,
1980→2028, both link directions, every b(t) minimum kept). PTF era:
**MJD 54891–57051** (2009-03-01 → 2015-01-28, the measured global
span of `ptf.ptf_procimg`, 2,818,277 CCD exposures; late-iPTF
2015–2017 was never publicly released). The era is concurrent with
PS1's (54985–57067) but from a different site, cadence pattern, and
filter set: covered windows here are **independent re-observations of
PS1-era windows**, including events PS1 lost to chip-gap masks — not
extensions of the covered-window record's time span.

## 1. Observable channels

Identical to the PS1 freeze: channel A (uplink interception,
`inbound`, blended with the star, near solar opposition) and channel
B (downlink pre-lens interception, `outbound`, at the star antipode,
z-parameterized parallax-reflex track 0.36–6.5″/day). Sunward
combinations out of scope here (§5.8 item 9 owns them). **[PTF]**
Channel-B antipode fields have no Pipeline-A corridor infrastructure
to reuse (PTF was never a Pipeline-A archive): field selection,
per-frame calibrator sets, and WCS handling are built fresh by the
adapter, on the `irsa_ztf.py` pattern with the recon's endpoint facts.

## 2. Beam-radius ladder

The ZTF/PS1 ladder verbatim: B 1.2 R☉ / B 2.5 R☉ / B 0.1 AU /
A 0.1 AU / A 1.0 AU, flat-chord window t_ca ± √(r²−b²)/v⊥ per event.

- **A 1.0 AU: constraint-only** (standing window≈season theorem,
  carried forward — PTF's campaign cadence concentrates epochs in the
  same season the window occupies).
- **[PTF] B 1.2 R☉: coverage-without-statistic.** The era scoping
  (§9) finds **zero in-window epochs** on the 0.3–0.6 d windows at
  any target — the closest epoch (van-maanen 2011, 0.36 d from t_ca)
  falls just outside the ±0.3 d window edge. The rung enters the
  covered-window ledger as structurally uncovered; the deep-graze
  event family is instead tested through the 2.5 R☉ rung, whose
  1.3 d windows do hold epochs (§9).

## 3. Wavelength **[PTF]**

Mould R (~577–727 nm) and SDSS-like g (~400–550 nm); occasional
Hα rows are excluded by a `fid ∈ {1, 2}` filter gate. The
frequency-doubled 532 nm line falls in g (the ZTF/PS1 test); R
carries the generic-leakage interpretation at its ~658 nm effective
center. 1064 nm and 1550 nm remain outside both bands — declared
unconstrained, as in every optical survey of the programme.
Hypothesis flux = monochromatic line converted via band effective
width. Both bands are primary where they have in-window epochs; the
scoped units are R-dominated (g is a minority of epochs at every
scoped position).

## 4. Duty cycle

The PS1 declaration verbatim: d = 1 while geometry holds (primary;
channel-B geometry enforces window-locked transience); d = 0.1 as a
flat ×d flux scaling of the 60 s exposures (secondary, no re-search);
declared unconstrained — pulse periods between the exposure time and
the window length, and schedules avoiding Earth-crossing windows.
No separate pulse statistic (PS1 precedent; PTF cadence within a
night is campaign-dependent and not frozen as a rung).

## 5. Detection constructions **[PTF: scie-direct substrate]**

The public PTF archive serves no difference images. Both channels use
**star-calibrated forced PSF photometry directly on level-1 science
cutouts** — the PS1 warp-direct estimator transferred: per-frame zero
point from field calibrators through the identical filter, with a
filter-median fallback where calibrator-sparse.

- **Calibration (measured basis).** Metadata `photcalflag` = 0 for
  ~100 % of epochs at 6 of 8 scoped positions (photometric-night
  flag): absolute header zero points (`MAGZPT`, `PHTCALEX`) exist for
  most frames but are non-photometric-night values, and `PHTCALEX=0`
  frames (3/30 in the off-window sample) have none. The per-frame
  field-star calibration is therefore **primary for all frames**;
  header MAGZPT is a cross-check diagnostic only. Frames whose
  star-calibration fails the gate (calibrator count / ZP scatter,
  thresholds set at the threshold freeze) are unusable, not
  fallback-calibrated.
- **Static sky is in the photometry.** Track nodes within 2″ of a
  catalogued static source are annotated at search time
  (`ptf_objects` snapshot, ≥ 3-detection rule); the flux-consistent
  catalogued-static test is the calibrated veto at adjudication.
  The level-2 reference image and its `depcov`/`uncert` products are
  annotation and veto inputs, not a differencing substrate.
- **A (blended):** window-locked excess vs the same star's off-window
  epochs at the per-epoch propagated star position; saturation rule
  §6; off-window-support gate frozen at the threshold freeze.
  Discriminators as PS1 (phase-lock, cross-epoch chromatic anomaly,
  chord shape/recurrence).
- **B (track):** per event, per z-grid point: star-calibrated forced
  photometry along the predicted track, shift-and-stack over the z
  family; null/thresholds from off-window epochs on the same tracks
  plus 8 offset control trajectories.
- **Veto ladder (B):** (1) SkyBoT known-object census per exceedance
  epoch (antipodes sit in the opposition asteroid stream);
  (2) rate test against the predicted 0.36–6.5″/day retrograde track
  rate; (3) **[PTF] same-night repeat test** — campaign cadence
  frequently delivers ≥ 2 exposures per night; an ordinary mover
  (≳ 15″/hr) leaves a fixed track position between same-night
  exposures while the relay track moves < 0.2″ — where same-night
  pairs exist this is the PS1 TTI-pair veto, and where they don't the
  ladder simply lacks that rung (recorded per exceedance);
  (4) recurrence on the recomputed track at any second covered
  window. Ring controls cross static sources at the same areal rate;
  the threshold absorbs the confusion floor.

## 6. Saturation rule (channel A) **[PTF]**

Frozen point-source saturation estimates for 60 s P48 exposures:
sat ≈ R 14.0 / g 14.5 (Law et al. 2009 R ≈ 14; g scaled by zero-point
difference). Rule per band, the established shape:

    excluded  m_est < E_b        E = {R 14.5, g 15.0}
    marginal  E_b ≤ m_est < E_b + 0.5
    ok        m_est ≥ E_b + 0.5

with m_est from a frozen catalog transform built at the cut stage.
The dev stage must verify E against ≥ 1 bright-star frame (dmask
bit-8 extent vs magnitude); an amendment adjusts E_b if the empirical
level disagrees by > 0.5 mag. Expected consequence, stated for the
record: van-maanen (R ≈ 12.4) and ross-128 (R ≈ 10.1) are
excluded-class in R on-star; teegarden is marginal-class in R and
ok-class in g — the A-channel unit list is decided by this rule at
dev, not silently.

## 7. Quality masks **[PTF]**

- **Exact (search-stage usability):** per-exposure `dmask`
  (ancillary `anciltype=dmask` — **select ancillaries by
  `anciltype`, never slot index**; slot order varies). Bit
  definitions pinned from Laher et al. 2014 Table 15:
  0 aircraft/satellite track, 1 object detected, 2 high dark
  current, 4 noisy, 5 ghost, 6 CCD bleed, 7 radiation hit,
  8 saturated, 9 dead/bad, 10 NaN, 11 dirt on optics, 12 halo;
  3/13–15 reserved. **Fatal template = 65533** (every bit except
  2¹ = object detected — a source at the locus is the signal; bleed
  and rad-hit are set in tandem per the paper; ghost/halo are
  bright-star circular masks of 85–450 px, exactly the artifact
  class the programme masks).
- **Listing-level (coverage stage):** none — every level-1 row at
  the position counts as covered (metadata `infobits` observed 0
  throughout the scoped sample; no archival bad-quality bit is
  documented for PTF level 1). Coverage freezes *before* dmask and
  calibration-gate attrition, as in PS1 — a "covered" window may
  still be lost at the mask.
- **Strict (exceedance re-runs):** seeing ≤ 2.5″ (metadata), fatal
  template, and the calibration gate at its strict threshold.
- **Depth accounting:** header `LIMITMAG` is present in only ~1/3 of
  frames (measured) and is **not** used. Per-epoch 5σ depths are
  computed by the photometry chain itself (PSF + background +
  star-calibrated ZP) — the same numbers the injections calibrate.
  Off-window measured expectation: R median ≈ 21.1 (range
  19.8–21.7 over seeing 1.5–3.5″); the 60 s single-epoch depth is
  comparable to PS1 warps, and the survey adds no new depth class —
  its value is *which windows* it covers, not how deep.

## 8. Injections, completeness, positive control

Stamp-response injections into the fetched cutouts spanning flux,
seeing, and detector position (the ZTF/PS1 machinery), calibrated
per unit; the C1 rule holds — no depth statement without the
injection-measured response. Positive control: one numbered
main-belt asteroid crossing a scoped antipode field, selected via
SkyBoT at dev, run through the identical
cutout→mask→calibration→statistic chain and recovered against its
predicted magnitudes (PS1 precedent: (60000)/220000 to ≤ 0.1 mag).

## 9. Era scope (from the frozen universal list + recon epoch lists)

`results/era_scope_v0.json`, flat-chord windows, per-era positions,
0.02° boxes. Units with any in-window coverage:

| unit | events covered / in era | in-window epochs |
|---|---|---|
| van-maanen B 2.5 R☉ | 1 / 12 | 2 (the 2011 b = 0.278 R☉ deep graze — the PS1-mask-lost family) |
| van-maanen B 0.1 AU | 4 / 12 | 25 |
| van-maanen A 0.1 AU | 1 / 12 | 11 |
| ross-128 B 2.5 R☉ | 2 / 12 | 9 |
| ross-128 B 0.1 AU | 2 / 12 | 67 |
| ross-128 A 0.1 AU | 1 / 12 | 4 |
| teegarden A 0.1 AU | 1 / 12 | 24 |
| gj-1276 A 0.1 AU | 1 / 12 | 3 |
| wolf-359 B 0.1 AU | 1 / 12 | 1 |

Zero-coverage ledger entries: B 1.2 R☉ everywhere (§2); teegarden,
gj-908, ross-154 antipodes and the ross-154 star (no PTF epochs at
all); every other unit-event without in-window epochs. gj-908's star
has 60 epochs, none in-window. The coverage stage re-derives this
table with fresh snapshot-disciplined queries and the exact per-event
positions; §9 numbers are scoping, not the frozen coverage record.

## 10. Pre-freeze data-contact declaration

All probe files live only in the session scratchpad and are
re-pulled fresh at the coverage stage. Contact at survey positions
before this freeze, in full:

1. **Metadata (timestamps, seeing, filters, WCS, filenames)** at all
   14 star/antipode positions — the same class every coverage stage
   touches pre-threshold-freeze; in-window timestamps were counted
   (the §9 scan) but no flux, pixel, or catalog quantity at a locus
   was formed into any window-locked statistic.
2. **Pixel cutouts, off-window only:** two cutouts at the van-maanen
   antipode (2011-01-21, 72 d from the nearest t_ca; 2013-03-25,
   9 d — outside the 5.8 d half-window of the widest searchable
   rung), fetched to verify WCS/cutout mechanics; 30 2-px header
   stamps at centers offset 0.05° from three positions (depth
   sampling, §7), all off-window.
3. **Static-catalog rows** (`ptf_objects`, era-integrated means)
   within 0.005° of the van-maanen antipode; nearest returned source
   14″ from the position — outside the 2″ static-annotation radius.

No window-locked quantity was formed from any of it; no ATLAS-style
`forced_dev` remedy is proposed (**decision D1**).

## 11. Freeze decisions — all adopted as recommended (user, 2026-08-26)

- **D1** — accept §10 as declared: no forced-dev remedy; scratchpad
  probe files never enter the repo; coverage stage re-pulls under
  snapshot discipline.
- **D2** — era MJD 54891–57051; input `crossings/universal_v1`
  (Earth-center; Palomar topocentric offset carried as the standard
  ≤ 0.8 %-of-rung budget term, as ZTF).
- **D3** — substrate: scie-direct star-calibrated forced PSF
  photometry, field-star ZP primary for all frames (measured
  photcalflag basis, §5); level-2 references as annotation/veto
  inputs only; no differencing.
- **D4** — unit = (target, channel-rung, band); statistics
  **S_event** (max single-event window-locked chord excess, primary)
  and **S_stack** (recurrence stack, only where ≥ 2 covered events:
  van-maanen B 0.1 AU, ross-128 B 2.5 R☉ and B 0.1 AU) — 1–2 trials
  per searched unit, tallied at the threshold freeze; d = 0.1 as ×d
  scaling; no pulse statistic.
- **D5** — masks and depth per §7 (dmask fatal template 65533;
  no listing-level gate at coverage; strict = seeing ≤ 2.5″ +
  fatal + strict calibration gate; LIMITMAG unused).
- **D6** — channel-A saturation rule per §6 (E_R 14.5 / E_g 15.0,
  dev verification with ± 0.5 mag amendment trigger).
- **D7** — controls: 8 offset control trajectories per event
  (PS1 construction) primary; temporal pseudo-windows reported
  secondarily only where the campaign cadence populates the offsets
  (assessed per unit at the threshold freeze, not assumed).
- **D8** — dev/confirmatory split: **dev** = wolf-359 B 0.1 AU
  (1 epoch — machinery), gj-1276 A 0.1 AU, ross-128 A 0.1 AU
  (saturation-rule exercise); **confirmatory** = van-maanen B 2.5 R☉,
  B 0.1 AU, A 0.1 AU; ross-128 B 2.5 R☉, B 0.1 AU; teegarden
  A 0.1 AU. Both headline B families (van-maanen deep-graze,
  ross-128) stay blind. Zero-coverage units enter the ledger as
  coverage-without-statistic, no split needed.
