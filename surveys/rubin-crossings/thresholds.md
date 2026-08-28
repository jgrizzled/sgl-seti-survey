# Rubin DP2 crossings threshold freeze v1.0 (+ amendment v1.1)

Frozen 2026-08-26 by `scripts/freeze_thresholds.py` →
`configs/threshold_freeze_v1.json`
(v1.0 `sha256:fe3cc037…21cf61`; re-frozen under amendment v1.1 as
`sha256:9384e97a…0b44da4`), bound to content hashes of
`hypotheses.md` v1.0+v1.1, the coverage products, and the frozen inputs
(`targets/universal_v2.json`, `crossings/universal_v1/events.ecsv`).
**No DiaSource row touched, no signal statistic formed** — the only
input at this stage beyond coverage metadata is the §6 saturation cut,
read from the already-frozen universal-list target file (no new
archive query).

## Saturation cut (frozen rule §6)

Ross 154 **G = 9.126** → **excluded** (< 16.5). Channel A 0.1 AU
closes **coverage-without-statistic (saturation-limited)**: its 16
on-detector visits — including two at Δt = −0.31 d sampling
b ≈ 3.6 R☉ — enter the covered-window ledger with no statistic. No
other channel-A unit had candidate visits.

## Search family — 1 unit × 1 statistic = 1 trial

| unit | event | visit | Δt | b(t_visit) | detector | magLim |
|---|---|---|---|---|---|---|
| **ross-128 B 0.1 AU r** | evt-5197f21cd67b (t_ca 60937.94, b_min 1.855 R☉) | 2025091300612 (r, 60932.284) | −5.657 d | 20.92 R☉ = 0.0973 AU | 102 | 23.366 |

**1 trial → expected control crossings 1/9 ≈ 0.111** (exchangeable
rate per trial, the standing crossings convention); FWER α = 0.05.
The visit samples the 0.1 AU cone at 97 % of its radius (rim sample);
the constraint statement is written for the wide rung, not the
1.86 R☉ core.

## Statistic (numeric, frozen)

- **S_det** = max over gate-passing associated DiaSources of
  `psfFlux`/`psfFluxErr`; 0 if no association survives.
- Association: within **r_assoc = 1.0″** of any deduplicated z-grid
  apparent position (550/1000/2500/5500/10000 AU at the visit epoch;
  positions closer than r_assoc to each other merge). At this unit
  the z-family spans offsets ~36.4″ (550 AU) → ~2.0″ (10,000 AU)
  from the axis point.
- Gates: freeze §7 exactly — on-detector with non-null magLim;
  pixelFlags exclusion set (bad / saturatedCenter / crCenter / edge /
  nodataCenter / interpolatedCenter / suspectCenter / streakCenter);
  archive-injection exclusion (`pixelFlags_injected*`);
  `psfFlux_flag`, `centroid_flag`.
- Static-sky exclusion: predicted position within **2″** of a
  `dp2.Object` with any-band psfMag < magLim(r) + 0.5 = 23.87 →
  position excluded; all positions excluded → `not_constrainable`.
- `reliability`: no absolute cut (freeze D1b); recorded as annotation
  on any association, cited at adjudication.

## Controls and rule

**8 pseudo-position controls per unit**: the full z-family offset
pattern about the axis point, rotated **k·40° (k = 1…8)** about the
axis point at the visit epoch (offset magnitudes preserved); a
control position failing the static exclusion, the footprint gate,
or the **≥ 2.5″ locus-avoidance (amendment v1.1** — the v1.0 10″
rule was geometrically impossible for the inner-z positions, found
at dev) rotates **+5°** until valid (frozen resolution rule;
resolutions recorded; a pattern unresolvable after a full turn is
recorded invalid). Threshold **T = max over the 8 control S_det**;
exceedance **S > max(T, 0)**.

## Veto ladder (frozen order, hypotheses §5)

(1) known-object census — dp2.SsSource association + SkyBoT cone at
the visit epoch; (2) static/residual check — Object +
ForcedSourceOnDiaObject history at the association position;
(3) artifact annotation — flags, `reliability`, trail fits;
(4) recurrence — no second in-era covered window exists; a surviving
exceedance is retained-ambiguous by construction (promotion never on
in-window evidence alone; recurrence falls to future releases).

## Dev stage (pseudo-units only, D8) — gates before confirmatory

(i) association-rate null over all off-window visits at the ross-128
antipode field (the identical S_det + control chain on epochs where
no relay is predicted); (ii) control-rule census (rotation
resolutions, static-exclusion rate); (iii) positive control — a
SkyBoT-predicted known object on an off-window visit must associate
a gate-passing DiaSource through the identical chain; (iv)
archive-injection census near the unit field: decides whether the
constraint is injection-calibrated (§8 route i) or a
magLim-referenced threshold statement (route ii). Blind confirmatory
runs only after every dev gate passes.

## Dev-stage assessment (2026-08-26; results/dev_v1.md)

All four gates resolved (`results/dev_v1.json`; the v1.0-rule run
preserved as `dev_v1_rule10as_superseded.json`):

- **Association-rate null — PASS.** 62 off-window pseudo-units at the
  ross-128 antipode field (blindness guard asserted the unit visit
  was never queried); 53 searched (9 off-detector/edge), **0
  exceedances, association rate S_det > 0: 0/53** — the 1/9
  exchangeable budget is very conservative for this sparse-DiaSource
  substrate (60″-cone DiaSource counts 0–1 per visit).
- **Control-rule census — PASS under v1.1.** 161 of 424 control
  patterns needed extra rotation (inner-z geometry); all resolved;
  the v1.0 10″ rule was unresolvable at small Δt epochs — the
  amendment's cause.
- **Positive control — PASS.** SkyBoT-predicted position of asteroid
  2006 SE393 (visit 2025071700513, TAI→UTC converted) landed 0.308″
  from the catalog DiaSource; the frozen chain associated it at
  S_det = 31.3 and the associated diaSourceId matches the catalog's
  own ssObjectId link. Two earlier attempts where the 3′ SkyBoT cone
  returned only a different, unrelated object (144″ away) are
  recorded as non-matches, not failures of the chain.
- **Archive-injection census — route (ii).** 0 injected-flagged
  DiaSources in 11 off-window r visits at the field: DP2 (Early)
  carries no usable injection population here. **The confirmatory
  constraint is a magLim-referenced threshold statement, explicitly
  labeled not injection-calibrated (freeze §8; C1 rule — no
  exclusion claim).**
