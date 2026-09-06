# PSP/WISPR sunward-channels crossings survey

**0 candidates.** Third heliospheric substrate for the sunward
(direction × side) beam combinations first searched with SOHO/LASCO
and STEREO-A HI-1, using Parker Solar Probe WISPR-I level-3 frames of
encounters E1–E27 (2018-11 → 2026-03) on a PSP-observer crossing list.
Blind confirmatory: **90 searched trials over 30 units, 15
exceedances vs 10.0 expected control crossings, all adjudicated —
9 static-content / latitude-curvature systematics (Gaia-matched),
3 bright-star-class, 2 arc-exit ramps, 1 extended coronal front;
every recurrence stack of the LASCO/HI-1 family null.**
Recurrence-stack depth V_eq 7.8–12.0 (median 9.3) → **downlink
0.08–3.9 GW through the 0.1 AU cone (median 1.2 GW), 10-m-class
uplink 0.07–4.4 GW (median 0.95 GW)** — LASCO-class, as declared at
the freeze; the survey's additions are the inner cone (b_e 2.6–14 R☉)
from a new vantage, first sunward constraints on 16 targets beyond
the seven-system family (23 targets searched), the first in-band 532 nm sunward substrate,
and a 5-min pulse cell that turned out to be systematics-limited.
Survey docs `surveys/wispr-crossings/`; plan §5.18–§5.19.

## 1. Hypotheses and geometry

The two sunward combinations excluded by declaration from every
night-sky survey (plan §5.11): **S1 — downlink post-lens** (outbound,
target side; apparent source = the star's antipode, a fixed ICRS
point) and **S2 — uplink past the Sun** (inbound, anti-target side;
apparent source = the target star). Observer: `crossings/psp_v1`
(Horizons −96 at 10-min sampling, 7,448 events; PSP crosses every
axis twice per 88–150-d orbit at a per-target repeating heliocentric
radius, so each target's crossing geometry recurs every encounter).
The decisive geometry (O5 pass, §5.18): the source's elongation is
ε = atan(b_e / r_along), so the 0.1 AU rung that subtends ≤ 5.7° from
1 AU subtends 14°–65° from 0.05–0.25 AU and its in-beam arc runs
through WISPR-I (13.5°–53.5°, ram side) for a median 23 h per event;
the post-t_ca half of every window is the ram-side half. The grazing
cones never enter the field (≤ 6.9° / 13.7°); the S2 1 AU rung is out
of scope (D7). Freeze: `hypotheses.md` D0–D10 (approved 2026-09-05),
`thresholds.md` v1.0 with amendments v1.1 (pre-dev) and v1.2–v1.5
(dev-driven, pre-confirmatory), all sha-pinned.

## 2. Chain

NRL level-3 synoptic WISPR-I frames (F-corona "grand minimum" model
subtracted, rescaled to L2 units by the header (r/0.2 AU)^2.3 factor;
E2's unbinned frames 2 × 2-averaged on load), header celestial WCS
refined by a Hipparcos match (translation + affine; 0.65 px rms
median), 15-px median high-pass, top-hat aperture (r 3 px = 7.6′,
PSF FWHM 1.6 px) at ICRS-fixed patches, per-frame 2-D colour-corrected
star ZP (colour coefficient +0.286 mag/(B−V) from 588k calibrator
records; uncertainty gate 0.10 mag), **stellar template** for the
static content (Tycho-2 ∪ Hipparcos + the S2 target star, with a
measured template response E = F − (0.07 + 0.73 T) from 5,371 control
patches, self-calibrated for S2 stars with 3 ≤ V < 8), 8 parallel-
track control patches at orbit-latitude offsets ±1.5°–±6° (same
frames), same-frame differential, per-event z, **S_stack** (Σz/√n),
**S_event** (max z), **S_pulse** (max over consecutive frame pairs of
the smaller standardised excess), T = max over the 8 controls,
S > max(T, 0). Vetoes: Mercury/Venus within 10°, Earth/Jupiter within
5° of any patch, 12-px proximity to any body, the four recon frame
epochs. Per-frame footprint from the frame's own WCS (12-px margin).

Volume: 62,915 synoptic WISPR-I frames fetched (24,645 dev + 38,270
confirmatory; every failed fetch recovered on rerun; ~250 GB
transient), 56,159 usable (89 %; the rest fail the ZP gate — mostly
inside 0.07 AU where the per-star scatter reaches ~1 mag), 0.65 px
astrometry, stamp response 0.98 ± 0.04.

## 3. Results (blind confirmatory, 2026-09-06)

`results/confirmatory_v1.md`, `confirmatory_search_v1.json`. 36
confirmatory units: **30 searched** (588 included events), 5 not
covered (procyon-a/b, luyten-star, gj-54, gj-1087 S1 — sources at
orbit latitude +17.5° to +19.6°, above the real field's +15° north
edge; the planning model's ±20° was optimistic there), fomalhaut S2
`bright_star_unsearchable` (V 1.2). **90 trials, 15 exceedances vs
10.0 expected**, adjudicated under the frozen ladder with the Gaia
DR3 static-content check and the latitude-curvature test
(`results/confirmatory_v1.md` §Adjudication):

- **Static content / latitude curvature (9)** — gj-1087 S2 (×3),
  gj-667-c S1 (×2), ross-128 S1, gj-581 S1, van-maanen S1 (×2): the
  source-patch excess is present at the same level in every
  encounter and equals, to ~20 %, the Gaia-predicted stellar flux
  missing from or mis-weighted in the template (Tycho-incomplete
  G 10–12.5 stars; encircled-fraction errors for bright neighbours
  2–3 px off-centre), and/or vanishes when the controls are
  interpolated quadratically in orbit latitude (sources on the
  brightness ridge near the orbital plane: van-maanen S1 S_stack
  17 → 5 < T; gj-581 S1 25 → 11 < T).
- **Bright-star class (3)** — gj-783 S2 (V 5.3): the declared v1.3
  limit; its template had counted the 1.6″/yr star more than once
  (Tycho positions not proper-motion propagated), the per-encounter
  calibration scatter of a S/N-100 star gives S_stack/S_event, and
  the "pulse" is a 2-frame spike of the star itself.
- **Arc-exit ramps (2)** — gj-514 S2, gj-1087 S2 pulses: the last
  2–4 frames of an arc at the 53.5° field exit.
- **Extended coronal front (1)** — wolf-1061 S1 pulse (E22, 20 min,
  6–7σ): the pixels show a uniform 0.5–1 unit/px elevation over the
  whole 11 × 11 stamp with no point-source peak, within 3° of the
  inner edge.
- gj-667-c S2 S_event: a single-encounter level offset on a patch
  dominated by the HD 156384 AB pair (V 6.7, 35″), in the
  V3-processed E01; the unit's other 20 events are null.

The LASCO/HI-1 cross-substrate family: van-maanen S2, wolf-359 S1,
teegarden S1/S2, gj-1276 S2, ross-128 S2 — S_stack −7.2 … 12.0 vs T
15–24, S_event 6–13 vs T 16–48, all null; van-maanen S1 and ross-128
S1 are the curvature/static cases above.

**Dev** (`results/dev_v1.md`): 5 units, 15 trials, 1 exceedance vs
1.67 expected (61-vir S2, a monotonic 0.95 → 0.88 calibration drift
of the V 4.7 star; retained, non-promotable). gj-908 S2 (V 9.0)
template null: z median +0.3 over 24 events.

## 4. Completeness and controls

`results/completeness_v1.json`, `power_limits_v1.json` — injections
into the real null series against the fixed confirmatory thresholds
(stamp response 0.98 ± 0.04 from 14,570 stamps), C1-normalised; the
pulse cell injects two consecutive frames (v1.5). Recovery threshold
= max(T, observed null S, 0); units whose null already exceeds T are
degenerate for that statistic (ross-128 S1, gj-581 S1, van-maanen S1
S_stack; gj-667-c S2 S_event).

| statistic | m90 V_eq (24–26 units) | S1 downlink through the 0.1 AU cone | S2 uplink, 10-m class |
|---|---|---|---|
| S_stack (recurrence) | 7.8–12.0, median 9.3 | 0.08–3.9 GW, median 1.2 GW | 0.07–4.4 GW, median 0.95 GW |
| S_event | 6.6–12.5, median 7.9 | ~3–10× shallower | |
| S_pulse (≥ 2 frames, 5–15 min) | 3.2–6.3 (9 units) | ~10–100 GW | |

Band conversion at λ_c 615 nm, W_eff 250 nm (532 nm in band;
1064/1550 nm as broadband leakage), ±0.3 mag band systematic, and a
declared **±0.3 mag flux-scale systematic** (the ZP surface reads
V 3–6 calibrators 0.23–0.33 mag faint because S/N-selected faint
calibrators are noise-boosted). Family limits (S_stack): wolf-359 S1
1.3 GW, teegarden S1 1.2 GW / S2 0.42 GW, gj-1276 S2 4.4 GW,
van-maanen S2 1.1 GW, ross-128 S2 69 MW (10-m uplink).

**Positive control (D5)**: (4) Vesta through WISPR-I at E22 (52
frames, r 0.10–0.19 AU, V 7.05–7.15 from PSP) recovered at **−0.19 ±
0.32 mag** per frame by the identical chain along the Horizons
track (source/control flux ratio 6.2); the E10 perihelion passage
(r 0.06 AU, 15-s exposures) is not recovered per frame, consistent
with the V 5.4–6.3 single-frame depth there. In situ: gj-908 S2
template null (dev); the Hipparcos ensemble ZP stable to 0.04 mag
across a 3× exposure change (recon).

## 5. What this survey adds

- **A third instrument on the 0.1 AU sunward cells**, at LASCO-class
  depth (recurrence stacks ~1 GW) from a vantage 5–20× closer to the
  Sun, covering the **inner cone** (b_e 2.6–14 R☉ on the perihelion-
  side crossings) that LASCO C3 saw only from 1 AU and HI-1 could not
  see at all (its arcs sit at b_e 14–22 R☉). Depth, not geometry, is
  the limit: single frames reach V 5.4–7.9, and the template
  systematics cap the stacks at V_eq ~9.
- **First sunward constraints on 16 targets** beyond the seven-system
  family (23 targets, 35 units searched over dev + confirmatory, all
  with ≥ 10 encounters), and the first in-band 532 nm sunward
  substrate.
- **The 5-min sunward pulse cell is opened but systematics-limited**:
  particle hits (5-exposure sums) and multi-frame coronal transients
  put the control thresholds at 4–70σ even with two-frame persistence.
- Honest framing: this survey does not deepen the cells HI-1 set at
  12–61 MW; proximity to the beam buys visibility, not flux.

## 6. Open items and hand-offs

1. **v2 construction lessons** (recorded, not applied post-blind): a
   Gaia DR3-based template (G < 15, proper-motion propagated,
   measured encircled-fraction curve) in place of Tycho/Hipparcos —
   9 of 15 exceedances were static-content residuals Gaia accounts
   for; quadratic-in-orbit-latitude control interpolation in place of
   the median (the ridge bias for near-plane sources); a field model
   with the measured +15°/−25° vertical extent (5 units were indexed
   outside the real field).
2. Frames inside ~0.07 AU are photometrically structure-dominated
   (per-star scatter ~1 mag); the deepest perihelion arcs contribute
   little — a K-corona/streamer model (L3 already removes the
   F-corona) would be needed to use them.
3. E28/E29 (2026-06, 2026-09) unreleased at the snapshot; `psp_v1`
   and the survey extend at the yearly refresh (Horizons SPK end
   2026-12-01).
4. WISPR-O (ε 50°–108°) carried as a ledger substrate; its arcs cover
   the outer cone only.
5. Solar Orbiter SoloHI (5°–45°, 0.28–0.9 AU) remains the natural
   fourth substrate on the same machinery (`--observer solo`).
