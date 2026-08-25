# ATLAS + ASAS-SN crossings hypothesis freeze v1.0 (Pipeline B)

Drafted 2026-08-25 after the reachability recon
(`notes/atlas_asassn_recon_2026-08-25.md`) and the era scoping
(`results/era_scope_v0.json`), before any archive data is touched at
survey positions (one recon exception, declared in §10).
**FROZEN 2026-08-25**: the user approved every §11 decision as
recommended (D1–D8; D1 includes the recommended re-pull — the recon
probe files live only in an ephemeral session scratchpad, were never
in the repo, and every coverage-stage light curve will be pulled
fresh under snapshot discipline). Execution note: the coverage stage
is **deferred until the project moves to a server** (see
`notes/server_migration_resume.md`); the ATLAS-era end date and the
ASAS-SN v2 currency re-measure are taken at coverage start on that
machine, per §2. Fifth crossings survey; ZTF/PS1/TESS constructions
inherited with ATLAS substitutions, marked **[ATLAS]**.

Input: `crossings/universal_v1` (Earth-center, 1980→2028, both link
directions). Earth-center is valid for ground observers: the
site-vs-geocenter baseline is ≤ R⊕ = 0.0092 R☉ — 0.8 % of the
tightest rung — and the four ATLAS sites differ by less; declared as
a budget term, no derivative observer list needed.

## 1. Archives and roles **[ATLAS]**

- **ATLAS forced photometry** — the search substrate for both
  channels: difference-flux (channel B, antipodes) and
  reduced/difference (channel A, on-star) tphot series at arbitrary
  positions, MJD 57227 → freeze date, c/o bands, all-sky. Access:
  token API (`scripts/atlas_api.py`), serial per-account queue.
- **ASAS-SN Sky Patrol v2** — supplementary, catalogued sources only:
  (a) the window **coverage-fraction ledger** (epoch lists of sources
  adjacent to each survey position — coverage records, never
  constraints); (b) a channel-A on-star series for unsaturated
  targets, V band reaching back to 2013 (pre-ATLAS for two years).
- **ASAS-SN Sky Patrol v1** (reCaptcha-gated) — a **manual
  instrument**: usable only for pre-registered, enumerated follow-ups
  (e.g. adjudicating an ATLAS exceedance; saturated-star ML
  photometry), never for a survey statistic.

## 2. Eras

| substrate | era (MJD) | note |
|---|---|---|
| ATLAS north of δ −50 | 57227 → coverage-start date | 02a from 57227, 01a joins 57928; end date recorded when the coverage stage begins |
| ATLAS δ < −50 | 59578 → coverage-start date | southern units only; **affects only the A 1.0 AU constraint rung** — every narrow-rung position sits at δ −24…+17 (measured, `era_scope_v0.json`) |
| ASAS-SN v2 V | ~56595 → 58386 | cameras ba–bh |
| ASAS-SN v2 g | ~58033 → DB ceiling (60841 at recon; re-measured at coverage start) | cameras bi–bt |

## 3. Channels, rungs, in-era scope

Channel definitions, beam-radius ladder, and the flat-chord window
construction are the ZTF freeze's, verbatim. Sunward combinations
remain out of scope (universal declaration; §5.8 item 9 owns the
question). In-era events (ATLAS era; per-target detail in
`era_scope_v0.json`):

| channel / rung | events / targets | note |
|---|---|---|
| B 1.2 R☉ | 88 / 4 | van-maanen (b 0.24 R☉!), wolf-359 (0.69), gj-1276 (0.72), teegarden (0.96) — 22 events each |
| B 2.5 R☉ | 110 / 5 | + ross-128 (1.83) |
| B 0.1 AU | 154 / 7 | + ross-154 (3.29), gj-908 (12.2) |
| A 0.1 AU | 154 / 7 | same 7 systems, star side |
| A 1.0 AU | 1,891 / 88 | **constraint-only, frozen at draft** — the standing window≈season theorem (4 independent confirmations) |

**The survey's new cell is recurrence.** Every prior archive saw at
most a handful of narrow-rung windows (PS1 lost the deepest to
correlated chip-gap masks; TESS's two deepest fell between sectors).
ATLAS holds **22 semiannual windows per narrow-rung target** with
nightly 4-exposure sampling — the van-maanen b = 0.24 R☉
photosphere-grazing family, unconstrained by every previous survey,
finally gets a searchable, *recurrence-stacked* test. This is the
window coverage-fraction / duty-cycle rationale the queue adopted
this item for.

## 4. Bands and wavelength **[ATLAS]**

o (560–820 nm) primary — ~3.5× the epochs; c (420–650 nm) secondary
and the chromatic check (c contains the frequency-doubled 532 nm
line; o contains 752 nm — the PS1 gj-1276 anomaly line falls in-band,
and the gj-1276 antipode *is* a B target here, but stacked depths
(~21) cannot reach that anomaly's 22.8 AB amplitude; declared, same
honesty as TESS). 1064/1550 nm remain outside every band in the
programme. ASAS-SN V/g as in §1. Hypothesis flux = monochromatic line
converted via band effective width, as always.

## 5. Duty cycle and temporal models

- **d = 1 chord — primary**, at three timescales: per-exposure (30 s),
  per-night quad (4 exposures over ~1 h — recon-confirmed), and the
  **recurrence stack across all in-era windows** of a unit (§6).
- **Pulse / flare-beacon — secondary**: max per-exposure statistic
  inside windows (vs controls) — the intra-hour quad rung sits
  between TESS's 200 s–window cell and the nightly rung; nightly
  sampling over 22 windows constrains duty cycle to roughly
  d ≳ (windows hit)/(windows sampled) for window-locked beacons.
- **Declared unconstrained**: sub-exposure pulses beyond ×d scaling;
  schedules avoiding Earth-crossing windows; single-window-only
  transmissions below per-window depth (the stack assumes recurrence).

## 6. Detection construction **[ATLAS]**

Per-event, windowed, fixed-position forced-photometry tasks
(`mjd ∈ [t_ca − 110 d, t_ca + 110 d]` — window + off-window baseline
+ pseudo-window range; windowed tasks keep queue time per task well
below the measured 25–70 min full-history cost).

**Positions.** Per event, at the tabulated apparent axis position at
t_ca (PM-propagated per event — the per-era drift of both stars
(≤ 5″/yr) and antipodes makes a single mean position invalid across
the era). Track drift within a window is 6.5″/d × window/2 at
z = 550 AU: ≤ 4″ for the grazing rungs (inside the 3.7″ FWHM —
handled as a per-epoch response factor), but ±37″ for B 0.1 AU at
z = 550. Draft proposal: grazing rungs = 1 position/event; B/A
0.1 AU = 3-position mini-track (t_ca, ±half-window at z = 550),
nested-z handled as position selection per epoch. Task budget ≈
(242 grazing-rung event-rows × 1 + 308 × 3) ≈ 1,170 windowed tasks —
days of queue drain, submitted as a background queue (§9 of the
recon). **Decision D3.**

**Statistics per unit** (unit = target × rung × band; **decision D4**):

1. **S_stack** — recurrence-stacked chord-weighted window-locked
   excess over all in-era events of the unit (primary; the new cell).
2. **S_event** — max single-event chord amplitude (secondary;
   catches non-recurring transmissions at per-window depth).
3. **S_pulse** — max per-exposure S/N in-window (tertiary; pulse
   cell).

All statistics run on **detrended** series: the TESS lesson
(undetrended chord statistics saturate their budget with
low-frequency drift) plus ATLAS's known template steps (MJD 58417,
58882) require a frozen detrending layer — draft: per-position
running median over off-window epochs, ~30 d window (**decision
D5**). Baseline and empirical variance rescale k = median(r²/v)/0.4549
(floor 1) from off-window epochs, per position and band (the
established recipe).

**Controls.** 8 temporal pseudo-windows (offsets ±23/47/71/97 d,
same-rung exclusion, the v1.1 rule): geometrically available here —
narrow-rung windows are 0.6–11.5 d against a ~183 d semiannual
recurrence, and nightly cadence populates every offset (the
theorem that killed temporal controls elsewhere applies only to
season-scale windows). Expected control crossings budget: 1/9 per
searched trial, tallied at the threshold freeze.

**Channel A (blended, 0.1 AU rung).** Difference-flux at the
PM-propagated star position; the parallax-factor systematics template
+ variance rescale constructions (ZTF v1.1/v1.2) transfer. ATLAS
difference photometry subtracts the static star by construction, but
bright stars leave saturation/PSF-wing residuals: frozen saturation
rule of the ZTF shape (excluded / marginal / ok with a 0.5 mag
margin) against an o/c-band saturation limit and per-target estimated
o/c magnitudes from a frozen catalog transform — table built at the
cut stage, exclusions reported per target (**decision D2**).

**Veto ladder (B):** (1) SkyBoT known-object census at every
exceedance epoch — antipode fields sit in the opposition asteroid
stream; (2) **intra-quad consistency** — a main-belt mover crosses
the 3.7″ PSF in ≲ 15 min, so at most ~1 of the 4 quad exposures is
contaminated; window-locked signals must be quad-consistent; (3)
template-step annotation — exceedances within ~30 d of MJD 58417 /
58882 (or a WPDATE change where images are pulled) are re-run on the
reduced-mode series; (4) recurrence — the stack demands it, and any
S_event exceedance is tested on the other 21 windows; (5) rate test
vs z-track prediction where a mini-track exists.

## 7. Quality masks **[ATLAS]**

Primary = the FAQ recipe verbatim (`atlas_api.faq_quality_mask`:
duJy < 10000, err == 0, chip-edge x/y margins, 1.6 < maj/min < 5,
−1 < apfit < −0.1, mag5sig > 17, Sky > 17; keeps 89–93 % — measured).
Strict (exceedance re-runs): additionally mag5sig > 18.5,
|MJD − template step| > 30 d, and drop H-filter rows. ASAS-SN v2:
quality == 'G' primary.

## 8. Flux scale, injections, completeness **[ATLAS: no image access]**

ATLAS fluxes arrive server-calibrated (AB, µJy). Without pixel access
the stamp-response injection machinery does not apply; draft:
**response-model injections** into the flux series — signal µJy ×
R(offset(t), maj/min(t)) with per-epoch offsets from track-vs-task
position — with the absolute scale validated by two gates before any
depth is quoted (the C1 rule in spirit: never combine a measured
scale with an unnormalized response): (a) the **MPC positive
control** — a numbered main-belt asteroid through the identical
task→mask→detrend→statistic chain, recovered against its predicted
magnitudes (the server's `mpc` mode makes this one task); (b) a
small image-request subsample if (a) leaves residuals > 0.2 mag.
**Decision D6.**

## 9. Expected sensitivity, declared before search

Measured single-exposure mag5sig ≈ 19.1; nightly quad ≈ 19.9. A
grazing window (0.6–1.3 d) holds 1–2 nights ≈ 4–8 exposures → per-
window ≈ 19.9–20.3; the 22-window recurrence stack (~100–180
exposures) ≈ 21.3–21.7 where windows are covered (weather/sun-gap
attrition measured at coverage stage — the coverage-fraction number
is itself a deliverable). Power scale: ZTF's 130 W (2.5 R☉, 21.8 AB)
rescales to **a few hundred W persistent-recurrent through the
grazing cones — including the 0.24 R☉ photosphere-grazing cone for
the first time in the programme**; per-window (non-recurrent)
limits ~kW-class. Depths must not be promised past what measured
coverage delivers; the threshold freeze sets the searchable-unit
population from the coverage stage, as always.

## 10. Pre-freeze data-contact declaration

The 2026-08-25 recon pulled two full-history light curves before any
freeze: the **wolf-359 antipode** (a prospective confirmatory B
position) and the σ Dra antipode (A 1.0 AU constraint rung only).
Only night counts, unit provenance, and quality-cut fractions were
computed — no window-locked quantity was formed. Remedy (draft,
following the DECam pilot precedent): **wolf-359 is forced into the
dev split** (`forced_dev`) at the threshold freeze; σ Dra needs no
remedy (no searchable rung). **Decision D1.**

## 11. Freeze decisions — all adopted as recommended (user, 2026-08-25)

- **D1** — accept the §10 remedy: wolf-359 `forced_dev`; van-maanen,
  gj-1276, teegarden (+ ross-128 at 2.5 R☉) stay confirmatory-
  eligible. (Alternative: also quarantine the recon probe files and
  re-pull under snapshot discipline — recommended regardless.)
- **D2** — channel A saturation rule: adopt the ZTF-shape rule with
  ATLAS o/c limits (proposed sat ≈ 12.5 + 0.5 margin) and the frozen
  catalog transform; same for ASAS-SN V/g (sat ≈ 10.5 + 0.5 —
  expected to exclude gj-908, ross-154, likely ross-128/van-maanen
  in V).
- **D3** — position construction: 1 fixed position per grazing-rung
  event; 3-position mini-track for the 0.1 AU rungs; per-epoch
  response factors carry the residual offset. (Alternative: collapse
  B 0.1 AU to high-z only with 1 position, halving the task budget.)
- **D4** — unit = (target, rung, band) with the three-statistic
  family (S_stack primary / S_event / S_pulse) — trials = 3 per
  searched unit, tallied against the 1/9 budget at threshold freeze.
- **D5** — detrending layer: per-position ~30 d off-window running
  median under all statistics (TESS lesson; template steps).
- **D6** — completeness by response-model injection with the MPC
  positive-control gate (≤ 0.2 mag) and image-request fallback.
- **D7** — eras of §2 (ATLAS end = freeze date; ASAS-SN v2 ceiling
  re-measured at freeze) and ASAS-SN roles of §1 (coverage ledger +
  unsaturated channel A only; v1 manual-only, enumerated follow-ups).
- **D8** — dev/confirmatory split: proposed dev = wolf-359 (D1) +
  ross-154 or gj-908 (wide-b, insensitive rungs — machinery
  validation); all grazing-rung targets except wolf-359 confirmatory.
  Seed declared at threshold freeze.
