# PS1 crossings hypothesis freeze v1.0 (Pipeline B)

Frozen 2026-08-24, before any coverage intersection or signal search.
Input: `crossings/universal_v1` (`xng-a09e2db7681d`, Earth-center,
1980→2028, both link directions, every b(t) minimum kept). PS1 era for
this survey: MJD 54900 → 57300 (generous brackets; the DR2 warp
archive actually spans MJD ≈ 54985–57067, 2009-06-03 → 2015-02-13,
measured from the 24,852 warps in `runs/panstarrs/coarse_v1`), MAST
DR2 warp products only. Second crossings survey; construction is the
ZTF crossings freeze (`surveys/ztf-crossings/hypotheses.md` v1.0)
with PS1-specific substitutions, each marked **[PS1]** below. The PS1
era is disjoint from ZTF's (2018–2026): every covered window here is
a new entry in the covered-window record, extending it back to
2009–2015.

## Observable channels

Of the four (link direction × side-of-axis) combinations, exactly two
are PS1-observable; the other two arrive from the solar direction
(the lensed/unlensed beam past the Sun) and are declared out of scope.

**Channel A — uplink interception** (`inbound`, `axis_distance_au > 0`).
Earth between Sun and star; beam arrives from the star's direction near
solar opposition. Source: blended with the star (a transmitter using the
star's own SGL, or any station within the seeing disc), plus a station
annulus — a direct transmitter offset a AU transverse appears at
a/d ≈ up to ~10″ for a ≤ 30 AU at the nearest targets. Crossing times
are insensitive to station offset (axis shift at Earth ~ a × 1 AU/d).

**Channel B — downlink pre-lens interception** (`outbound`,
`axis_distance_au < 0`). Earth between relay and Sun; beam arrives from
the relay = star antipode, near the anti-solar point. The relay at
z ∈ [550, 10 000] AU has apparent parallax-reflex motion v⊕/z
(6.5″/day at 550 AU, 0.36″/day at 10 000 AU) along a predicted track;
position and rate are functions of the single parameter z. **[PS1]**
The channel-B antipode locus is by construction the Pipeline-A corridor
of the same target: the corridor skycells, calibrator-density grades
(`surveys/panstarrs/targets/overlay_v1.md`), skycell WCS cache, and
per-frame star-calibration machinery are reused directly.

## Beam-radius ladder (assumed effective radius at Earth)

Identical to the ZTF freeze:

| channel | radius | AU | rationale |
|---|---|---|---|
| B | 1.2 R☉ | 0.0056 | photosphere-grazing downlink; beam radius at 1 AU ≈ solar impact parameter |
| B | 2.5 R☉ | 0.0116 | coronal-avoidance grazing annulus |
| B | 0.1 AU | 0.1 | wide-beam / pointing-margin downlink |
| A | 0.1 AU | 0.1 | diffraction-limited ~10 m-class transmitter at 3–5 pc |
| A | 1.0 AU | 1.0 | ~1 m-class transmitter; effectively any listed minimum |

Window per event and radius r > b: t_ca ± √(r² − b_min²)/v⊥ (flat-chord
approximation; b_min, v⊥ from the events table). Representative full
durations at v⊥ = 30 km/s, b ≪ r: 0.6 d (1.2 R☉), 1.3 d (2.5 R☉),
11.5 d (0.1 AU), 115 d (1 AU).

**[PS1] The A 1.0 AU rung is declared constraint-only at freeze time.**
ZTF dev lesson 1 (carried forward, not re-derived): its window (~115 d)
is comparable to the observing season, the crossing geometry
phase-locks to the sidereal year, and PS1's single-phase cadence
concentrates all epochs in that same season — temporal controls fail by
construction. The rung is reported with coverage and depth only; it
carries no discovery statistic and enters no trials count.

## Wavelength **[PS1]**

PS1 grizy constrains ~400–1000 nm: g 394–559, r 550–689, i 690–819,
z 819–922, y 918–1001 nm (half-max passbands). The canonical SGL link
wavelengths remain outside: 1064 nm falls redward of the y half-max
cutoff and 1550 nm is far outside — both declared unconstrained. PS1
tests the frequency-doubled 532 nm line (mid-g, exactly as ZTF) and
generic in-band leakage in r/i/z/y. Hypothesis flux is a monochromatic
line converted to per-band AB magnitude via the band effective width;
the frozen line wavelengths are 532 nm (g) and the band effective
centers 617/752/866/962 nm (r/i/z/y, generic-leakage interpretation).
All five bands are primary (they share MAST availability and cadence;
there is no ZTF-style sparse-band caveat).

## Duty cycle

- **d = 1 while geometry holds** — primary. In channel B the *geometry*
  enforces window-locked transience regardless of transmitter schedule
  (relay radiates toward the Sun; Earth is in the beam only in-window).
- **d = 0.1** — secondary; pulse period ≪ the exposure (30–45 s
  nominal: g 43, r 40, i 45, z 30, y 30 s) integrates to a flat ×d flux
  scaling, reported without re-search.
- **Declared unconstrained:** pulse periods between ~30 s and the window
  length (missed-epoch probability per epoch ≈ 1 − d), and transmitters
  scheduled to avoid Earth-crossing windows.

## Detection constructions **[PS1: warp-direct substrate]**

PS1 publishes no difference images. Both channels use star-calibrated
matched-filter forced photometry directly on warp cutouts (the
Pipeline-A PS1 estimator: per-frame zero point from DR2 `mean` stars
through the identical filter, filter-median fallback where the field is
calibrator-sparse), not on a differenced substrate. Consequences,
frozen here:

- The static sky is *in* the photometry. Channel-B track nodes within
  2″ of a DR2 catalogued static source (≥ 3 detections, the
  `sglsurvey/vetting.py` snapshot rule) are annotated at search time;
  the flux-consistent catalogued-static test is the corresponding
  calibrated veto at adjudication. Ring controls cross static sources
  at the same areal rate, so the threshold absorbs the confusion floor
  (v1 lesson: confusion, not per-pixel noise, sets the floor).
- Channel A measures the star itself, so the statistic is a
  window-locked excess against the same star's off-window epochs — but
  with no reference image there is no PM-dipole/reference systematic,
  and the ZTF v1.1 parallax-factor template is *not* imported. The
  off-window baseline is thin (~10–15 epochs per band over the era);
  units without enough off-window support are constraint-only (gate
  frozen in the threshold freeze).

**A (blended):** star-calibrated forced photometry at the per-epoch
propagated star position over all epochs; statistic = window-locked
excess vs the same star's off-window distribution. Discriminators in
order of power: (1) phase-lock to the crossing ephemeris; (2) chromatic
anomaly — single-band excess with a flat other band (laser signature;
stellar variability is broadband-correlated; PS1 bands are not
simultaneous, so this discriminator is cross-epoch and weaker than
ZTF's same-night version); (3) chord light-curve shape and recurrence
at successive crossings with b-modulated duration. Station annulus
≤ 10″: point-source search, same construction as B without the
z-track. Feasibility: saturation cut per band (rule below); sensitivity
is contrast-limited, not sky-limited; this channel is supplementary
to B.

**B (track):** per event, per z-grid point: predicted track over the
window; star-calibrated forced photometry on warp cutouts along the
track, shift-and-stack over the z family. Null/thresholds from
off-window epochs on the same tracks plus 8 offset control trajectories.
Veto ladder: (1) MPC/known-object census per epoch — low-b events
belong to near-ecliptic stars whose antipodes sit in the opposition
asteroid swarm; (2) rate test — candidate must move at the predicted
~0.3–6.5″/day retrograde rate for its z; (3) **[PS1] TTI-pair motion
test** — 3π warps come in same-night pairs ~15–40 min apart, over which
an ordinary mover (≳ 15″/hr) displaces ≥ several arcsec while the relay
track moves < 0.2″: a same-night pair showing consistent flux at a
fixed track position is an automatic ordinary-mover veto ZTF lacked;
(4) recurrence at a second covered window on the recomputed track. No
reference-image contamination check is needed (no differencing);
in-window epochs of the *calibration* star sample are unaffected by
construction (calibrators are field stars, not the hypothesis source).

## Saturation rule (channel A) **[PS1]**

Frozen per-band point-source saturation estimates for 3π warp
exposures, derived from GPC1 full well ≈ 65 ke⁻, peak-pixel fraction
≈ 4 % at 1.2″ FWHM / 0.258″ pixels, per-band e⁻/s zero points and
nominal exposure times: sat ≈ g 13.3 / r 13.3 / i 13.4 / z 12.5 /
y 11.5. Rule per band, mirroring the ZTF cut (sat + ~0.5 margin,
+0.5 marginal, transform accuracy ±0.3 absorbed):

    excluded  m_est < E_b        E = {g 14.0, r 14.0, i 14.0, z 13.0, y 12.0}
    marginal  E_b ≤ m_est < E_b + 0.5
    ok        m_est ≥ E_b + 0.5

The estimate is deliberately conservative (~0.5 mag beyond the derived
sat). Warp headers carry `CELL.SATURATION` and masks carry SAT/STARCORE
bits, so the dev stage must verify the rule against ≥ 1 bright-star
warp before any confirmatory A search; an amendment tightens or relaxes
E_b if the empirical level disagrees by > 0.5 mag. Optically-dark
targets (late-T/Y) cannot saturate and carry the empty-field
sensitivity statement.

## Quality masks **[PS1]**

Listing-level: `badflag == 0` (uniformly true across the 24,852
corridor warps; retained as a formal gate). Coverage counts every
badflag-0 warp ("primary" at coverage stage); seeing is unavailable
until header/cutout fetch, so — unlike ZTF — the primary mask carries
no seeing term at coverage time. Search-stage usability is the exact
warp mask (fatal template 16255, `surveys/panstarrs/hypotheses.md`
v1.0 §9): Pipeline-A experience is ~75 % of nominal-footprint hits
usable, and coverage counts freeze *before* that attrition — a window
"covered" by a single epoch may still be lost at the mask. Strict mask
(exceedance robustness re-runs): PSF FWHM ≤ 2.5″ from the warp header,
plus the fatal-template test.

## Observer

Universal list is Earth-center; PS1 observes from Haleakalā. The
topocentric offset ≤ R⊕ = 4.3 × 10⁻⁵ AU is 0.8 % of the smallest rung
radius (1.2 R☉ = 0.0056 AU) and shifts window edges by ≤ 2 min at
v⊥ = 30 km/s: carried as a budget term, exactly as Palomar was for
ZTF. Apparent relay positions (channel B) are evaluated per epoch with
the astropy builtin ephemeris, Earth center.

## Era event counts (from the frozen universal list, pre-coverage)

Channel A: 563 events / 86 targets (b < 1 AU); 46 / 7 at b < 0.1 AU
(gj-1276, gj-908, ross-128, ross-154, teegarden, van-maanen,
wolf-359). Channel B: 554 events in era; 27 / 4 targets at b < 1.2 R☉,
34 / 5 at 2.5 R☉, 47 / 7 at 0.1 AU (same seven targets; all antipodes
at δ > −30°, inside the PS1 footprint — no ZTF-style southern loss on
the narrow rungs).

Coverage intersection (which windows contain usable exposures at the
predicted positions) is the next stage and precedes any threshold
freeze: thresholds must not promise a depth the cadence cannot deliver.
Expected shape, stated for the record: ~110–130 warps per corridor
skycell clustered in the observing season centered on opposition —
which is where channel-B windows sit by geometry — so wide-rung events
should typically catch 1–3 epochs (single-epoch class dominant) and the
grazing rungs (~1 d windows) a handful of events at best.
