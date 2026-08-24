# PS1 beam-crossings survey (Pipeline B) — report

Second Pipeline-B (§3.5) survey: a search of Pan-STARRS1 DR2
single-epoch warps for optical emission during Earth's crossings of
hypothesized Sun–star relay beam axes, 2009–2015 — windows disjoint
from and eight years before the ZTF crossings era, extending the
covered-window record backwards. Executed 2026-08-24. **No
candidates.**

## 1. Construction and provenance chain

1. **Universal crossing list** `crossings/universal_v1`
   (`xng-a09e2db7681d`): every b(t) minimum for 88 endpoints × 2 link
   directions, Earth-center, 1980→2028.
2. **Hypothesis freeze v1.0** `surveys/ps1-crossings/hypotheses.md`:
   the ZTF channel constructions with PS1 substitutions — era MJD
   54900–57300 (warps span 54985–57067); warp-direct star-calibrated
   substrate (PS1 has no difference images); grizy all primary
   (532 nm line in g, band-center lines in r/i/z/y; 1064/1550 nm
   outside the bands and unconstrained); the A 1.0 AU rung declared
   constraint-only at freeze (window ≈ observing season, ZTF lesson 1
   adopted rather than re-derived); frozen saturation levels with a
   dev-stage verification gate; Haleakalā-vs-Earth-center offset
   carried as a ≤ 1.5 % budget term.
3. **Coverage intersection** (`results/coverage_v1.*`; 93 discovery
   cones, skycell WCS cache seeded from Pipeline-A): B wide rung
   14/47 events covered on **all seven** narrow-rung targets (PS1's
   season centers on opposition — where channel-B windows sit),
   typically 1–3 same-night TTI visits per window; grazing rungs 1/27
   and 1/34 (epoch-starved as the cadence predicts) — the one covered
   grazing event is van-maanen at **b = 0.28 R☉** (2010-04-03, TTI
   i-band pair), the deepest graze in the programme; A 0.1 AU 18/46,
   A 1.0 AU 294/563 (coverage recorded; per-window depths deferred
   with the rung's constraint-only status).
4. **Channel-A saturation cut** (`results/saturation_cut_v1.*`,
   grizy): survivors gj-1276 (all bands), teegarden (g + marginal r),
   van-maanen (y); rule verified at dev against 32 bright-star warps
   + header `CELL.SATURATION` — true saturation is ≥ 1 mag brighter
   than every frozen level, so the levels stand.
5. **Threshold freeze v1.0** (`configs/threshold_freeze_v1.json`,
   seed 20260824): 19 search units (2 A + 17 B), 8 designated
   controls (B: 20/30/40″ rings; A: temporal pseudo-windows),
   exceedance S > max(T, 0), expected 2.1 control crossings; the ZTF
   v1.1/v1.2 amendments (variance rescale, same-rung pseudo-window
   exclusion, retired R ratio) adopted at freeze time; new TTI-pair
   mover veto; 2″ catalogued-static annotation. Dev split: B-wide
   teegarden + wolf-359; A dev empty (singleton strata).
6. **Dev search** (`results/dev_search_v1.md`): 0 exceedances, **no
   amendment** — no reference image means no PS1 analog of the ZTF
   dev A-channel failures. Findings: P1 correlated exact-mask
   attrition (a unit's 1-night epochs share pointing, so `CONV.BAD`
   chip-gap bands kill whole units), P2 partial ring controls
   recorded as anomalies, P3 the off-window variance-rescale
   realization (fixed antipode point, all-era same-band epochs).
7. **Confirmatory search** (`results/confirmatory_v1.*`): 432 epoch
   cutouts, 0 fetch failures; 12 B rows + 4 A rows.
8. **Completeness** (`results/completeness_v1.json`): exact
   stamp-response injections (Moffat β=3 at the warp `CHIP.SEEING`,
   `sglsurvey.inject.stamp_response` — linear-filter exact), 200
   bootstrap draws per unit per magnitude (grid 14–22; 22.0 =
   grid-censored lower limit), background drawn from the unit's
   off-window sample pool (ring samples where that pool is starved),
   injected z drawn from the recoverable z set with dead z declared,
   recovery against the unit's frozen threshold.

## 2. Results

**11 searched units (all channel B; channel A entirely
constraint-only by the frozen gates): 2 exceedances vs 2.1 expected
control crossings over the frozen 19-unit family — the budget rate.
Both fully adjudicated through the frozen ladder; 0 candidates.**

- **ross-128 r `evt-3f7ad50e70` (S 5.06, T 3.61): vetoed** —
  flux-consistent catalogued static. The DR2 *stack* catalog places
  static sources 1.46″ (r = 22.5–22.8) and 1.99″ (r = 23.5) from the
  in-window track nodes; the measured node fluxes (r ≈ 21.6–22.1)
  sit in their blended matched-filter response envelope. SkyBoT
  census clean; TTI-pair consistency excludes movers.
- **gj-1276 i `evt-344c12d32a` (S 3.32, T 2.30): retained-ambiguous,
  non-promotable.** A night-consistent i ≈ 22.7–23.4 signal
  (1.1–2.3 σ per epoch over 3 usable exposures spanning ~35 min) at
  the near-static z = 550 node. Census clean; main-belt rates
  excluded (no 21″ displacement); TTI-consistent; **no static
  counterpart** — the DR2 stack catalog is empty within 6″ at ~23+
  depth, so the static veto cannot fire; recurrence cannot run (the
  only covered PS1 window). Amplitude sits inside the 1/9
  control-crossing budget. Promotion requires recurrence at a second
  covered window; the designated test is the joint crossings stage —
  gj-1276 has ZTF-covered wide-rung windows 2018–2026.

**Mask attrition (the survey's structural finding).** 7 of the 18
covered-window unit-rows lost every in-window sample to `CONV.BAD`
chip-gap bands — correlated per unit because each window's epochs are
1–2 same-pointing TTI nights. **The van-maanen b = 0.28 R☉ grazing
family is among them**: the photosphere-grazing rung ends
*unconstrained*, recorded in the covered-window ledger as
nominal-covered / mask-unusable. Unit-level attrition (~40 %) far
exceeds the ~25 % per-epoch Pipeline-A rate, as the correlation
predicts.

### 90 %-recovery depths (AB mag; line hypothesis at 532 nm in g, band centers elsewhere)

| unit population | n | m90 |
|---|---|---|
| B 0.1 AU searched units, all bands | 11 | median 21.6 (4 grid-censored ≥ 22) |
| — g (532 nm hypothesis) | 4 | 19.2 / 20.6 / 21.6 / ≥22 (median 21.1) |
| — r | 3 | 20.5 / 21.1 / 21.4 |
| — i | 4 | 21.7 / ≥22 / ≥22 / ≥22 |
| A reference depths (constraint-only, threshold = max(T,0)) | 3 | gj-1276 g 17.9 / i 18.0; teegarden g 21.3 |

Two units (wolf-359 r, gj-908 g) exclude the z = 550 grid node from
their depth statement (node masked in every usable epoch; z-dead sets
declared per unit). The shallow tail (ross-154 g 19.2, T = 35) is the
warp-direct confusion floor in its near-Galactic-plane antipode field
— the ring controls cross catalogued stars and the frozen threshold
absorbs them.

## 3. Physical interpretation (532 nm line in g unless noted)

**Downlink pre-lens.** A relay at 550 AU beaming through the 2.5 R☉
solar-grazing cone (gain ≈ 9×10⁹) toward the target star would have
been detected at 90 % confidence during the searched 2010–2014
windows at transmitter power ≳ **~280 W** (g median m90 = 21.1 →
EIRP 2.5×10¹² W at 549 AU); the deepest units reach ~170 W
(all-band median 21.6) and the grid-censored i units ≲ 37 W
(752 nm leakage). The constraint applies only to the 9 searched
events on 6 targets — a per-event statement, not a continuous
monitor — and lands in windows no other survey has covered: every
one predates the ZTF era.

**Uplink interception.** No searchable channel-A unit exists (the
single-phase cadence starves the pseudo-window controls), so PS1
sets only reference depths: a 10-m-class 532 nm transmitter at
teegarden would have appeared above the star's off-window scatter at
≳ ~1 kW (reference m90 21.3, threshold 0, no calibrated
false-alarm control), gj-1276 at ~23–32 kW. These are qualified
sensitivity statements, not calibrated exclusions.

**Not constrained:** 1064/1550 nm links (outside grizy); the
photosphere/coronal grazing rungs (their one covered event is
mask-unusable); transmitters scheduled to avoid Earth crossings;
pulse periods between ~30 s and the window length; the 5 remaining
uncovered wide-rung events; δ-independent — all seven narrow-rung
targets are in-footprint, the losses here are cadence and masks, not
declination.

## 4. Lessons for the next crossings adapters

1. **Correlated mask attrition is the PS1-class failure mode**: when
   a survey visits a field 1–2 nights per window with repeated
   pointing, per-epoch mask attrition (~25 %) converts to ~40 %
   whole-unit loss, and it took the survey's single best event
   (van-maanen b = 0.28 R☉). Future single-phase archives (DECam)
   need a coverage stage that tests the exact mask at the narrow-rung
   nodes *before* the threshold freeze counts a window as covered.
2. The warp-direct (no differencing) substrate works: the confusion
   floor lands in the ring-control threshold (T up to 35 in a
   Galactic-plane field) instead of producing false candidates, and
   the deep DR2 *stack* catalog supplies the decisive
   static-vs-transient discriminator the ZTF construction got from
   difference imaging.
3. The TTI-pair veto is free and sharp: same-night pairs separate
   ordinary movers (arcsec–arcmin) from relay tracks (< 0.2″) with no
   extra data.
4. Exceedance budgeting again behaved as designed: 2 observed vs 2.1
   expected, one vetoed by a calibrated test, one retained-ambiguous
   and handed to the cross-archive recurrence test rather than
   discretionary dismissal.

## 5. Products

`surveys/ps1-crossings/{hypotheses.md, thresholds.md, configs/*,
results/*, scripts/*}`; universal list under `crossings/universal_v1/`;
snapshots, skycell WCS cache, masks and cutouts under
`runs/ps1-crossings/` (717 MB: coverage snapshots + masks 332 MB,
dev cutouts 83 MB, confirmatory cutouts 293 MB);
the retained-ambiguous exceedance and the mask-unusable
covered-window rows are the survey's inputs to the joint crossings
ledger (plan §11.2 step 6).
