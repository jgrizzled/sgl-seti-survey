# Solar Orbiter / SoloHI sunward-channels crossings survey

**0 candidates.** Fourth heliospheric substrate for the sunward
(direction × side) beam combinations first searched with SOHO/LASCO,
STEREO-A HI-1 and PSP/WISPR, using SoloHI level-2 inner-tile frames
of Solar Orbiter orbits P04–P12 (2022-01 → 2026-04) on a
Solar-Orbiter-observer crossing list. Blind confirmatory: **54
searched trials over 18 units, 8 exceedances vs 6.0 expected control
crossings, all adjudicated — 4 control-interpolation artifacts
(gj-1111 S2 ×3, gj-251 S2), 2 uncatalogued moving objects through the
patch (teegarden S2), 2 extended level offsets of alternating sign
(gj-876 S1); every recurrence stack of the LASCO/HI-1/WISPR family
null.** Recurrence-stack depth V_eq 9.1–13.9 (median 10.5) →
**downlink 14 MW – 1.2 GW through the 0.1 AU cone (median 340 MW),
10-m-class uplink 33–350 MW (median 89 MW)** — LASCO-class, ~1 mag
deeper than WISPR's stacks, as declared at the freeze; the survey's
additions are the **pre-t_ca half of every window** (WISPR saw the
post-t_ca half), the outer cone b_e 7.5–21.5 R☉ from 0.3–0.9 AU, the
first star-fixed baseline on an inner-heliosphere sunward substrate,
and first sunward constraints on 13 targets beyond the seven-system
family (20 targets, 22 units searched over dev + confirmatory).
Survey docs `surveys/solohi-crossings/`; plan §5.23.

## 1. Hypotheses and geometry

The two sunward combinations excluded by declaration from every
night-sky survey (plan §5.11): **S1 — downlink post-lens** (outbound,
target side; apparent source = the star's antipode, a fixed ICRS
point) and **S2 — uplink past the Sun** (inbound, anti-target side;
apparent source = the target star). Observer: `crossings/solo_v1`
(Horizons −144 at 10-min sampling, 4,624 events; Solar Orbiter crosses
every axis twice per 150–230-d Venus-resonant orbit at a per-target
repeating heliocentric radius, stepping only at the Venus flybys).
The decisive geometry (§5.23 pass,
`notes/solo_solohi_geometry_pass_2026-09-06.md`): the source's
elongation is ε = atan(b_e / r_along), so the 0.1 AU rung subtends
6°–20° from 0.28–0.9 AU and its in-beam arc runs through the two inner
SoloHI tiles (λ −5.2° → −45°, |β| ≤ 20°, **anti-ram side** — the field
model was measured from 48 L2 headers: the tiles are fixed in the
orbit-plane frame at every epoch, with 0.49° seams and the inner seam
on the orbital plane) for a median 42 h per event, always in the
**pre-t_ca half** of the window. The grazing cones never approach the
5.2° edge (≤ 1.1° / 2.4°); the S2 1 AU rung is out of scope (D7).
Freeze: `hypotheses.md` D0–D9 (approved 2026-09-06), `thresholds.md`
v1.0 with amendments v1.1 (saturation mask) and v1.2 (structure-noise
epoch gate, colour calibrators V ≥ 4.5), both dev-driven and frozen
before any confirmatory pixel, all sha-pinned.

## 2. Chain

SoloHI L2 inner-tile frames (`1ft` 1024 × 960, `2ft` 960 × 1024;
F-corona in, no background-subtracted science product exists),
fetched alternately from the NRL day tree and the SOAR data service
(the same files; 4 MB each, size- and signature-checked), header
celestial WCS refined by a Hipparcos match (translation + affine;
0.56 px rms median), **saturation mask** (pixels ≥ 0.85 × DSATVAL:
the F-corona plateau along the sunward edge in the ≥ 45-s regime),
25-px median high-pass, top-hat aperture (r 3 px = 3.7′, PSF FWHM
1.8 px) at ICRS-fixed patches, per-frame 2-D colour-corrected star
ZP (+0.371 mag/(B−V) from 6.9e5 calibrator records; uncertainty gate
0.10 mag), **star-fixed differential** for the static content
(E′ = F − median F over the same patch on the same-tile off-beam
frames before the arc, ≥ 10 baseline epochs; the Gaia DR3 template as
a secondary annotation), 8 parallel-track control patches at
orbit-latitude offsets ±1.5°–±6° (one-sided toward the tile interior
for 125 of 178 unit-events), same-frame differential with the frozen
quadratic-in-latitude control interpolation, **structure-noise gate**
(epochs whose annulus noise exceeds 1.5 units — the inner ~120 px of
every arc, ε ≲ 8°, where the unresolved corona scatter is 6–11 units
per frame — are dropped: 49,388 of 193,902 confirmatory arc records, 25 %),
per-event z, **S_stack** (Σz/√n), **S_event** (max z), **S_pulse**
(max over consecutive frame pairs of the smaller standardised
excess), T = max over the 8 controls, S > max(T, 0). Vetoes:
Mercury/Venus within 10°, Earth/Jupiter within 5° of any patch, 12-px
proximity, the 96 recon frame-event contacts. Per-frame footprint
from the frame's own WCS inside DSTART/DSTOP (12-px margin).

Volume: 46,682 frames fetched (9,313 dev v1.0 + 9,313 dev v1.1 +
28,056 confirmatory; ~185 GB transient, purged per day; no fetch
failed), 27,973 confirmatory frames usable (99.7 %), stamp response
0.95 ± 0.07 (12,198 stamps).

## 3. Results (blind confirmatory, 2026-09-06)

`results/confirmatory_v1.md`, `confirmatory_search_v1.json`. 27
confirmatory units: **18 searched** (75 included events over 3–8
orbits each), 9 `constraint_only` with < 3 events surviving the edge
gates and the baseline rule (ez-aqr S2, gj-1002 S2, gj-1111 S1,
gj-1276 S2, gj-876 S2, ross-128 S1, ross-154 S2, teegarden S1,
wolf-359 S1 — their arcs lie mostly inside the sunward-edge band).
**54 trials, 8 exceedances vs 6.0 expected**, adjudicated under the
frozen ladder with re-fetched pixel stamps, the Horizons `@-144`
known-object census, aligned arc stacks, and one post-blind
diagnostic reduce with a robust median control interpolation (the
frozen statistics stand; `confirmatory_v1.md` §Adjudication):

- **Control-interpolation artifacts (4)** — gj-1111 S2 (S_stack 19.2
  vs 13.3, S_event 16.9 vs 12.8, S_pulse 38 vs 30) and gj-251 S2
  (S_pulse 5.48 vs 5.42): a bright star in the +1.5° control patch
  loses its core to the saturation mask on single frames (E′ −47 to
  −54 units), and the frozen quadratic, extrapolated to the source
  offset on a one-sided ladder, turns that into +50–60 units on the
  source; the source pixels are flat at those epochs (−3.6 … +0.4
  units). Diagnostic median interpolation: gj-1111 S2 S_stack −3.5 vs
  5.6, S_pulse 2.4 vs 5.4.
- **Uncatalogued moving objects (2)** — teegarden S2 S_event (P12,
  16-epoch arc) and S_pulse: compact peaks displaced 4 px between
  consecutive frames and absent before and after (2026-02-23 20:02 →
  20:26; 2025-03-24 05:13 → 05:37, diffuse); no H < 7 asteroid within
  1°; the star's own baseline 0.14–0.19 units.
- **Extended level offsets (2)** — gj-876 S1 S_stack 4.6 vs 3.7 and
  S_event 7.9 vs 4.85: per-orbit z +5.9 (P09), −5.4 (P10), +7.9 (P11)
  on a starless antipode patch; aligned mean stacks of 169 / 79 / 86
  frames are uniform ±0.01–0.04 units/px over 15 × 15 px with no
  point-source peak; marginal against 8 controls spanning −4.2 …
  +3.7; no consistent-sign recurrence.

The LASCO/HI-1/WISPR cross-substrate family: van-maanen S1 (S_stack
−1.9 vs 7.8) and S2 (0.8 vs 8.9), wolf-359 S2 (−4.3 vs 1.1), ross-128
S2 (0.8 vs 8.4), gj-1276 S1 (−7.6 vs 14.2), teegarden S2 (1.8 vs 4.4)
— null; wolf-359 S1, ross-128 S1, teegarden S1, gj-1276 S2 and
ross-154 S2 fell below three events.

**Dev** (`results/dev_v1.md`, `notes/dev_machinery_log_2026-09-06.md`):
5 units, 12 trials, 0 exceedances vs 1.33 expected; gj-908 S2 (V 9.0)
baseline null (z median −1.2 over 4 orbits), gj-908 S1 z median 0.00,
61-vir S2 `bright_star_saturated` (the V 4.7 star's aperture fails
the mask/validity gates in 5 of 6 events).

## 4. Completeness and controls

`results/completeness_v1.json`, `power_limits_v1.json` — injections
into the real null series against the fixed confirmatory thresholds
(stamp response 0.95 ± 0.07), C1-normalised; the pulse cell injects
two consecutive frames. Recovery threshold = max(T, observed null S,
0); units whose null exceeds T are degenerate for that statistic
(gj-1111 S2 all three, gj-876 S1 stack/event, teegarden S2 event/
pulse).

| statistic | m90 V_eq (16 units) | S1 downlink through the 0.1 AU cone | S2 uplink, 10-m class |
|---|---|---|---|
| S_stack (recurrence, 3–8 orbits) | 9.1–13.9, median 10.5 | 14 MW – 1.2 GW, median 340 MW | 33–350 MW, median 89 MW |
| S_event | 8.4–11.4, median 9.7 | ~2–3× shallower | |
| S_pulse (≥ 2 frames, 12–48 min) | 4.8–7.9, median 6.8 | ~10–100 GW | |

Band conversion at λ_c 610 nm, W_eff 250 nm (532 nm in band;
1064/1550 nm as broadband leakage), ±0.3 mag band systematic, and a
declared **±0.3 mag flux-scale systematic** (V 4–5 calibrators read
0.06–0.08 mag faint, V 3–4 0.30 mag: the SoloHI bright-star
nonlinearity; the ensemble ZP is fitted on V ≥ 4.5). Family limits
(S_stack): van-maanen S1 690 MW / S2 730 MW (10-m uplink 320 MW),
wolf-359 S2 330 MW (uplink 45 MW), teegarden S2 95 MW (uplink 33 MW),
ross-128 S2 330 MW (uplink 89 MW), gj-1276 S1 1.2 GW. Deepest units:
gj-783 S1 (V_eq 13.9, 14 MW), gj-674 S1 (13.4, 23 MW), ez-aqr S1
(12.6, 48 MW).

**Positive control (D5)**: (4) Vesta through tile 1 on 2022-08-20 →
09-14 (118 frames, r 0.84 → 0.61 AU, ε 6°–26°, V 7.5–7.9 from Solar
Orbiter) recovered at **−0.05 ± 0.14 mag** per frame by the identical
chain along the Horizons track; (1) Ceres 2024-06-30 (ε 5°–15°) yields
0 usable frames — its passage lies entirely in the edge band the
gates exclude. In situ: gj-908 S2 baseline null (dev); the Hipparcos
ensemble ZP −19.87 … −20.02 across 16–145 s exposures (recon).

## 5. What this survey adds

- **A fourth instrument on the 0.1 AU sunward cells**, at LASCO-class
  depth (recurrence stacks ~0.3 GW, ~1 mag deeper than WISPR's) from
  0.3–0.9 AU, covering the **pre-t_ca half of every sunward window**
  that no other substrate sees (LASCO/HI-1 sample the Earth-side
  windows, WISPR the post-t_ca ram-side half) and the outer cone
  b_e 7.5–21.5 R☉. Its own limit is the sunward edge: inside ~120 px
  (ε ≲ 8°) the unresolved corona structure and, in the long-exposure
  regime, the saturation plateau make the frames unusable for point
  sources, which removed the inner quarter of every arc and 9 of 27
  confirmatory units.
- **First sunward constraints on 13 targets** beyond the seven-system
  family (20 targets, 22 units searched over dev + confirmatory), and
  the first heliospheric-imager substrate with an **off-beam
  star-fixed baseline** for every arc — the differential that removed
  the template systematics WISPR had to adjudicate.
- **The 12–48-min sunward pulse cell is systematics-limited** (V_eq
  4.8–7.9): uncatalogued movers through the fixed patches and
  single-frame control excursions set the thresholds.
- Honest framing: this survey does not deepen the cells HI-1 set at
  12–61 MW; proximity to the beam buys visibility, not flux.

## 6. Open items and hand-offs

1. **v2 construction lessons** (recorded, not applied post-blind): a
   robust (median or clipped) control interpolation — the frozen
   quadratic extrapolates single-control excursions into the source on
   one-sided ladders (4 of the 8 exceedances; the diagnostic reduce
   gives 4 vs 6.0 expected with lower thresholds); control patches
   masked for Gaia G ≤ 8 wherever the saturation mask can bite; a
   local (not global-median) fill under the mask before the
   high-pass; an inverse-variance per-epoch weighting in place of the
   hard 1.5-unit gate.
2. The sunward-edge band (ε ≲ 8°, ~120 px) is unusable at this
   resolution; a K/F-corona structure model or the 16-s perihelion
   frames only would be needed to recover it.
3. P13+ (perihelion 2026-08-19) unreleased at the snapshot
   (2026-04-10); `solo_v1` and the survey extend at the yearly
   refresh (SPK to 2030-11-20).
4. The seam rows (ross-128 P05–P09, in-plane) and the 9 constraint-
   only units are ledger entries; wolf-359 S1 / ross-128 S1 /
   teegarden S1 / gj-1276 S2 / ross-154 S2 need the 16-s perihelion
   geometry (P04–P05-like exposures) to become searchable.
5. The outer tiles (ε 25°–45°) were used only as baseline; they hold
   the S2 1 AU transits (36 stars, 214 covered star-orbit pairs),
   out of scope.
