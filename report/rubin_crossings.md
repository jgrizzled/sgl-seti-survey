---
title: "Rubin DP2 crossings survey (Pipeline B) — report"
date: 2026-08-26
status: "complete — 0 candidates (1 blind trial, 0 exceedances vs 0.11 expected)"
---

# Rubin DP2 crossings survey (plan §5.13)

Ninth Pipeline B survey; the programme's first on a **catalog-level
substrate** — the detection layer is the Rubin DP2 difference-image
pipeline (`dp2.DiaSource`), consumed through association at predicted
relay-track positions, with no pixel access (DP2 Early publishes no
visit/difference images). Chain: recon
(`surveys/rubin/notes/rubin_recon_2026-08-26.md`) → hypothesis freeze
v1.0 + amendment v1.1 → coverage → threshold freeze → dev
(pseudo-units only) → blind confirmatory, all 2026-08-26, on the dev
machine. Survey docs: `surveys/rubin-crossings/`; snapshots:
`runs/rubin-crossings/` (326 TAP + SkyBoT query-snapshot records
across coverage/dev/confirmatory, reruns included).

Era: MJD 60790.117 → 61047.354 (2025-03-26 → 2025-12-08; 28,698
LSSTCam visits, ugrizy, 30 s). Processing pinned v30 / DM-55060.
Input: `crossings/universal_v1` Earth-center (Cerro Pachón site
offset ≤ 0.009 R☉ — inside the standing 0.010 R☉ budget).

## 1. Result

**0 candidates.** The search family was **1 unit × 1 statistic =
1 blind trial**: ross-128 downlink (channel B), 0.1 AU rung, r band —
visit 2025091300612 (MJD 60932.282, Δt = −5.66 d from t_ca), all
five z-grid positions on detector 102 (magLim r = 23.366), four live
after the static-sky exclusion (z = 550 AU sits 2″ from a bright
static Object → `not_constrainable` at that rung point; searched
family z ∈ {1000, 2500, 5500, 10000} AU). **S_det = 0** — zero
DiaSources of any kind in the 60″ discovery cone (not merely zero
associations), all 8 controls valid (rotation census recorded),
T = 0, no exceedance (0 observed vs 0.11 expected control
crossings).

## 2. Constraint (threshold statement — not injection-calibrated)

DP2 (Early) carries no archive-injection population at this field
(dev census: 0 flagged injections in 11 visits) and we cannot inject
into pixels we do not have, so per the frozen completeness route the
constraint is referenced to the served per-detector 5σ point-source
depth, **explicitly not injection-calibrated** (C1 rule: this is a
threshold statement, not an exclusion):

- The visit samples the beam cylinder at **b(t) = 20.9 R☉ =
  0.0973 AU** — a rim sample of the 0.1 AU rung (beams narrower than
  0.0973 AU are unconstrained by it; the event's 1.86 R☉ core was
  never sampled, see §3).
- A relay at z ∈ {1000–10,000} AU on the ross-128 axis, broadband in
  r, illuminating a 0.1 AU-radius top-hat beam containing Earth
  during the visit, would exceed the 5σ depth above **P ≈ 1.3 kW**
  band-integrated (magLim 23.366 → F_ν = 1.64 µJy; Δν_r ≈
  1.09 × 10¹⁴ Hz → F ≈ 1.8 × 10⁻¹⁸ W m⁻²; × π (0.1 AU)² =
  7.0 × 10²⁰ m²). Equivalent EIRP at z = 1000 AU: ~5 × 10¹¹ W
  (cone gain 4 × 10⁸).
- This is the **deepest wide-rung single-epoch flux threshold in the
  covered-window record** (~4× deeper in flux than ZTF's m90 21.8
  wide-rung median, which — being injection-calibrated — remains the
  reference exclusion; the two statements are different record
  classes and are not merged).
- Duty-cycle reach: one 30 s sample at Δt = −5.66 d constrains
  persistent-during-window emission; schedules avoiding that 30 s
  are unconstrained (freeze §4).

## 3. Covered-window ledger (no statistic)

- **Grazing rungs structurally uncovered** (±0.35 d-class windows,
  0 visits): B 1.2 R☉ — wolf-359 (b = 0.70 R☉), teegarden (0.99);
  B 2.5 R☉ — + ross-128 (1.86: the same event whose 0.1 AU rim was
  searched; its photosphere-grazing core window fell between visits).
- **A 0.1 AU saturation-limited**: ross-154 (G = 9.13, frozen §6
  rule) — 16 on-detector visits including two at Δt = −0.31 d
  sampling b ≈ 3.6 R☉; recorded coverage-without-statistic. The
  catalog substrate cannot serve bright-star channel A; this cell
  waits for image-level handling or a different archive.
- **A 1.0 AU**: 85 in-era events / 84 targets, geometry-only (the
  programme-wide deferred wide-beam rung, §5.7).
- In-era 0.1 AU events with no candidate visits: B ross-154,
  teegarden, wolf-359; A gj-1276, gj-908, teegarden, van-maanen
  (b = 0.32 R☉ — the deep-graze family's 2025 window fell outside
  DP2's visited sky).

## 4. Validation summary

- **Association-rate null** (dev, 53 off-window pseudo-units through
  the identical chain, unit visit blind-guarded): 0 associations,
  0 exceedances — the sparse-DiaSource regime makes the 1/9 budget
  very conservative.
- **Positive control PASS**: SkyBoT-predicted asteroid 2006 SE393
  associated at 0.308″ / S_det = 31.3 through the frozen chain, with
  the associated `diaSourceId` matching DP2's own `ssObjectId` link
  (timing incl. TAI→UTC, astrometry, gates verified end-to-end).
- **Amendment v1.1** (found at dev, before any unit contact): control
  locus-avoidance 10″ → 2.5″ — the v1.0 value was geometrically
  impossible for inner-z offsets; statistic untouched.
- Pre-freeze recon contact at the unit field (aggregate cone counts,
  static coadd) declared in hypotheses §10; remedy D1 kept the unit
  blind with the declaration citable at adjudication — no exceedance
  arose, so no adjudication occurred.

## 5. Lessons and hand-offs

1. **Catalog-level crossings surveys work end-to-end** — association
   + pseudo-position controls + SkyBoT positive control transfer the
   v2 discipline to a substrate with no pixels; the missing piece is
   injections, which caps every constraint at threshold-statement
   class (route ii). When Rubin publishes visit/difference images
   (late 2026), the image-level re-run upgrades this unit to an
   injection-calibrated exclusion — the designated follow-up.
2. **The sparse-DiaSource null is not a free pass**: with 0–1
   DiaSources per 60″ cone, S_det = 0 on both unit and controls is
   the typical outcome and the exceedance budget is loose; any
   association at all is effectively adjudication-bound. The veto
   ladder, not the threshold, carries the false-positive control on
   this substrate.
3. **Rim-sample honesty**: a wide-rung window can be covered by a
   visit that samples only the cone's rim; the constraint statement
   must carry b(t_visit), not b_min. (Here: 20.9 R☉ vs the 1.86 R☉
   core.)
4. Standing maintenance (§5.7): fold the one searched unit + ledger
   rows into the covered-window ledger at the next refresh; re-run
   the intersect when DP2 grows to its full "January 2026" span or
   Year-1 lands (more in-window epochs, recurrence tests, and the
   van-maanen 2025 deep-graze window may become covered).
