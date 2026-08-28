# LASCO sunward-channels crossings survey

**0 candidates.** First search of the sunward (direction × side) beam
combinations — excluded by declaration from every previous survey in
the programme — using SOHO/LASCO C2/C3 coronagraph imaging,
1996–2026. Blind confirmatory: 30 searched trials over 13 units,
2 exceedances vs 3.3 expected control crossings, both adjudicated
(one cosmic-ray single frame, one CME-period window); every
recurrence statistic null. Survey docs
`surveys/heliospheric-crossings/`; execution record
`notes/project_history.md` §9; plan §5.11.

## 1. Hypotheses and geometry

The sunward geometry study (2026-08-25) showed the two out-of-scope
combinations are confined to solar elongation ≲ 6° for every frozen
beam-radius rung — coronagraphs are the unique archival substrate:

- **S1 — downlink post-lens** (outbound, target side): Earth inside
  the relay's beam after it grazed the Sun; the apparent source is
  the solar-limb graze point at radius b_e(t), at the position angle
  of the observer's transverse axis offset. The S1 sky position is a
  **fixed ICRS point** (the antipode direction); the Sun sweeps past
  it.
- **S2 — uplink past the Sun** (inbound, anti-target side): the
  target star near conjunction, ε ≈ arcsin(b_e/1 AU).

Observer: `crossings/soho_v1` (Horizons −21; SOHO's L1 halo
cross-track ≈ 0.9 R☉ was measured directly in the recon star check —
Earth-center geometry is invalid at grazing precision). Era
1996-01 → 2026-06 (NRL currency). No z-grid exists: both channels'
positions are z-independent. Rungs: S1 2.5 R☉ (C2, wing regime
[2.2, 2.5] R☉ — the 1.2 R☉ rung is `not_constrainable` behind the
C2 occulter), S1/S2 0.1 AU (C3, usable annulus [4.4, 29] R☉).
S2 1.0 AU declared out of scope (non-exclusive twilight skirt).

## 2. Chain

Level-0.5 synoptic frames (SDAC + NRL trees, anonymous), per frame:
azimuthal-median radial-profile detrend → synthetic celestial WCS
(Sun-from-SOHO + P-angle + roll auto-resolution) → Hipparcos star fit
(2-star minimum on C2, 5 on C3; radius-resolved colour-cut star ZP) →
forced photometry at the per-epoch geometry positions.
Statistic (amendments v1.1/v1.2): ring-differential (source minus the
median of 8 same-radius PA-ring controls on the same frame), night
medians, baseline-median centring, empirical night-level σ;
S_stack (recurrence, primary) / S_event / S_pulse (C2 full-cadence
only); rule S > max(T, 0), T = max over rings. Masks: planet
proximity, bright-planet in-FOV veto (Venus/Jupiter — found via
ross-154's exact 8-yr Venus-synodic contamination cycle in dev),
Tycho-2 VT ≤ 11 star proximity (with the S2 source-star exemption),
occulted-plateau epoch gate, aperture validity. S2 units with
V ≤ 12.5 targets carry a stellar-flux template (measured colour
system c_C3 = 0.429 mag/(B−V), 8,414 calibrators).

## 3. Results (blind confirmatory, 2026-08-26)

**0 candidates** — 30 trials, 2 exceedances vs 3.3 expected; all
S_stack null (`surveys/heliospheric-crossings/results/confirmatory_v1.md`
has the full unit table). Adjudications under the frozen ladder:

- gj-1276 S1 2.5 R☉ **S_pulse 18.7/17.4** → one 12-min frame
  (2014-09-04 07:12, +2078 with negative neighbours): single-frame /
  cosmic-ray class, **vetoed** (persistence rule).
- van-maanen S1 2.5 R☉ **S_event 4.51/1.19** → one window
  (2020-10-05/06): a ~3 h all-position-angle annulus disturbance
  whose onset follows a CDAW-catalogued C2 CME by 48 min; both ±25°
  rings swing by hundreds in both signs; no recurrence in 27 sibling
  windows. **Adjudicated CME-period systematic**, retained,
  non-promotable. (The catalogued CPAs differ from the source PA —
  the azimuthal-median detrend couples all PAs during corona-wide
  transients; a design note for any v2.)
- ross-128 S1 0.1 AU: **constraint-only — permanently blended**: the
  fixed S1 sky position lies 121″ (2.2 px) from a VT 7.2 star.
  Ledger: nominal-covered / resolution-blended.

Dev stage (5 units, 11 trials, 1 exceedance vs 1.2): ross-154's
Venus-cycle stray light vetoed; ross-128's marginal S_event on the
2015 St. Patrick's Day CME storm window adjudicated; **gj-908
(V 8.98) detected window-locked at the predicted S2 position and
nulled by the stellar template — the survey's in-situ positive
control** (detect → model → null).

## 4. Completeness and controls

Stamp-measured response (Gaussian stamps through the identical
detrend+photometry chain; R ≈ measured per camera) + injections into
the real confirmatory null series against the fixed confirmatory
thresholds, 100 draws/flux (seed 20260825;
`results/completeness_v1.json`). **The coronal night-to-night
systematics, not photon noise, set the floor** — measured m90 depths
are far shallower than naive stack scaling:

| cell | m90 (V) | power through cone |
|---|---|---|
| S1 2.5 R☉ recurrence stacks (C2) | 4.5–7.3 | **40–530 MW downlink** |
| S1 2.5 R☉ pulse (gj-1276; ≥ 1 frame/12 min) | 9.3 | **6.7 MW·(12 min)** |
| S1 0.1 AU recurrence stacks (C3) | 6.2–8.4 | 3–26 GW downlink |
| S2 0.1 AU recurrence stacks (C3) | 6.2–7.6 | 10-m-class uplink 0.7–1.1 GW |

Full per-unit table incl. S_event: `results/power_limits_v1.json`
(line-equivalent conversion, W_eff 300 nm C3 Clear / 100 nm C2
Orange, declared; the 532 nm doubled line is in-band for C3 only).

**Flux-scale validation:** 217 single-frame S/N > 8 star measurements
unbiased (−0.02 ± 0.39 mag); ensemble colour system scatter
0.53 → 0.38 after correction; era stability — year-to-year scatter of
the chain's scale 0.064 mag, pre-2017 vs post-2018 offset −0.066 mag
(the level-1-era boundary; `results/zp_validation_v1.json`).
**Uranus control** (`results/uranus_control_v1.json`): recovered at
the predicted moving position at ~20σ in a 250-frame stack, but
1.65 mag fainter than its V — consistent with Uranus's
methane-absorbed red spectrum in the red-weighted Clear band; the
≤ 0.2 mag flux gate is carried by the stellar validation and the
gj-908 red-dwarf control, and a ± ~0.3 mag band-conversion systematic
is declared on the depths.

## 5. What this survey adds

The sunward cell — archivally virgin before this programme — is now
searched at MW-class power: the deepest statement is the **pulse
cell** (≥ 12-min pulses ≳ 6.7 MW through the 2.5 R☉ cone) and the
grazing-wing recurrence stacks (≳ 40–70 MW persistent-recurrent for
wolf-359/gj-1276). These are shallow next to the night-sky channels
(ZTF's 130 W) — the corona is a hostile background — but they are the
only constraints of any kind on the Sun-side beam combinations, and
the survey validated the full coronagraph substrate chain
(astrometry, roll, photometry, observer geometry) for the STEREO/WISPR
follow-ons.

## 6. Open items and hand-offs

- STEREO HI-1 (own drifting-observer list; V ~ 13 depth would improve
  the sunward cell by ~5 mag where its FOV applies) and WISPR (PSP
  observer-geometry pass) — plan §5.8 queue.
- v2 design notes: remove in-window days from baseline picks (~5 %
  conservative contamination, recorded); a CME-robust detrend
  (per-PA-sector profile instead of full-azimuth median); inner-field
  (r ≲ 13 R☉) scale carries larger uncertainty (sparse calibrators).
- The 1.2 R☉ rung stays open at every archive (C2-occulted here;
  mask/sector-gap-lost elsewhere) — a future-observation note.
- Yearly refresh: extend `soho_v1` past the 2026-10 Horizons SPK end;
  the 2026-07+ NRL-currency window tail (7 events) is recoverable.
