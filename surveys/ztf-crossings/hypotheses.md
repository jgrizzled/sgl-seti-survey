# ZTF crossings hypothesis freeze v1.0 (Pipeline B)

Frozen 2026-08-23, before any coverage intersection or signal search.
Input: `crossings/universal_v1` (`xng-a09e2db7681d`, Earth-center,
1980→2028, both link directions, every b(t) minimum kept). ZTF era for
this survey: MJD 58178 (2018-03-01) → 61275 (2026-08-23), public
products only (`ipac_gid=1`).

## Observable channels

Of the four (link direction × side-of-axis) combinations, exactly two
are ZTF-observable; the other two arrive from the solar direction
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
position and rate are functions of the single parameter z.

## Beam-radius ladder (assumed effective radius at Earth)

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

## Wavelength

ZTF constrains ~400–900 nm: zg 408–552, zr 562–695 primary; zi
opportunistic (sparse coverage), matching the v2 band policy. The
canonical SGL link wavelengths (1064, 1550 nm) are outside ZTF's
response — ZTF tests the frequency-doubled 532 nm line (mid-g) and
generic in-band leakage. Hypothesis flux is a monochromatic line
converted to per-band AB magnitude via the band effective width.

## Duty cycle

- **d = 1 while geometry holds** — primary. In channel B the *geometry*
  enforces window-locked transience regardless of transmitter schedule
  (relay radiates toward the Sun; Earth is in the beam only in-window).
- **d = 0.1** — secondary; pulse period ≪ 30 s exposure integrates to a
  flat ×d flux scaling, reported without re-search.
- **Declared unconstrained:** pulse periods between ~30 s and the window
  length (missed-epoch probability per epoch ≈ 1 − d), and transmitters
  scheduled to avoid Earth-crossing windows.

## Detection constructions

**A (blended):** forced photometry on the star position over all epochs;
statistic = window-locked excess vs the same star's off-window
distribution (~2 windows/yr, ≥ 17 windows in-era per target). Discriminators
in order of power: (1) phase-lock to the crossing ephemeris; (2)
chromatic anomaly — single-band excess with a flat simultaneous other
band (laser signature; stellar variability is broadband-correlated);
(3) chord light-curve shape and recurrence at successive crossings with
b-modulated duration. Station annulus ≤ 10″: point-source search on
difference images, same construction as B without the z-track.
Feasibility: ZTF saturates at ~12.5–13 mag — endpoints brighter than the
per-band saturation limit are excluded from the blended search and the
exclusion is reported per target, not silently dropped. Sensitivity is
contrast-limited (~1–2 % of stellar flux), not sky-limited; this channel
is supplementary to B.

**B (track):** per event, per z-grid point: predicted track over the
window; forced photometry on ZTF difference images along the track,
shift-and-stack over the z family. Null/thresholds from off-window
epochs on the same tracks plus 8 offset control trajectories (v1
discipline: confusion, not per-pixel noise, sets the floor). Veto
ladder: (1) MPC/known-object census per epoch — critical: low-b events
belong to near-ecliptic stars whose antipodes sit in the opposition
asteroid swarm; (2) rate test — candidate must move at the predicted
~0.3–6.5″/day retrograde rate for its z (asteroids ~30″/hr fail fast,
static sources fail slow); note these rates are *below* intra-night
tracklet thresholds, so the census cannot veto by matching the
candidate itself, only ordinary movers; (3) semiannual recurrence on
the recomputed track (analog of the parallax-phase veto). Reference
images must be checked for in-window contamination and masked if any
reference epoch falls inside a window.

## Quality masks

Reuse ztf-v2 masks: primary = archival bad-quality false, seeing ≤ 4″;
strict adds seeing ≤ 2.5″, maglimit ≥ 19.5, moon illumination ≤ 0.8.

## Era event counts (from the frozen universal list, pre-coverage)

Channel A: 720 events / 86 targets (b < 1 AU), 59 / 7 at b < 0.1 AU.
Channel B: 722 events in era; 35 / 4 targets at b < 1.2 R☉, 43 / 5 at
2.5 R☉, 60 / 7 at 0.1 AU.

Coverage intersection (which windows contain usable exposures at the
predicted positions) is the next stage and precedes any threshold
freeze: thresholds must not promise a depth the cadence cannot deliver.
