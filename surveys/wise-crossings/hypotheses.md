# WISE crossings hypothesis freeze v1.0 (Pipeline B)

Frozen 2026-08-24, before any coverage intersection or signal search.
Input: `crossings/wise_v1` (`xng-ffc922f7bbca`, **WISE-spacecraft
observer** — the first crossings survey whose event list and archive
share the same observer; L1b-header table sha-pinned, Earth-center
fallback in > 2 d gaps), era MJD 55203–60523 (2010-01-07 →
2024-08-01) minus the hibernation gap MJD 55593–56639 (2011-02-01 →
2013-12-13, no frames). Third crossings survey; construction follows
the ZTF/PS1 chain with WISE substitutions, each marked **[WISE]**.

## The elongation gate (scope-defining, computed pre-freeze)

WISE observes only on the great circle at solar elongation ≈ 90°
(sun-synchronous terminator orbit). Every beam-crossing geometry puts
the beam source near the observer–Sun axis during the window
(elongation ≈ 180° for the interception channels; ≈ 0° for the
sunward ones ZTF declared out of scope). Computed over the frozen
`wise_v1` list (window sampling, |elongation − 90°| < 4°, hibernation
excluded):

| channel / rung | in-era events | with any in-window elongation-90 epoch |
|---|---|---|
| B 1.2 R☉ | 50 | **0** |
| B 2.5 R☉ | 67 | **0** |
| B 0.1 AU | 102 | **0** |
| A 0.1 AU | 101 | **0** |
| A 1.0 AU | 1,272 | **102** (45 targets) |

**Channel B (the ZTF/PS1 workhorse) and the A 0.1 AU rung are
geometrically null for an elongation-90° observer and are declared out
of scope** — a scope reduction fixed by spacecraft attitude law, not by
cadence or depth; the report must state it as such, never as searched
coverage. The survey's sole channel is **A (uplink interception) at
the 1.0 AU rung**: wide windows (≈ 115 d at b ≪ 1 AU) and
large-b/high-ecliptic-latitude geometry (where t_ca decouples from
opposition, and the ecliptic-pole family lives at elongation ~90°
year-round) put 102 events onto WISE's scan circle (frozen gate: `configs/elongation_gate_v1.json`).

**[WISE] The A 1.0 AU rung is *searchable here*, unlike ZTF/PS1.**
Their constraint-only declaration rested on the season-lock degeneracy
(window ≈ the observing season, all epochs in-window). WISE visits
every position at fixed elongation for ~1–10 d every ~6 months across
14 years: in-window epochs are the 1–2 visits intersecting the
window, and the off-window baseline is the ~20+ other visits — no
degeneracy. Whether pseudo-window controls have epoch support is a
per-unit gate (threshold freeze), not a rung-wide failure.

## Channel definition

**Channel A — uplink interception** (`inbound`,
`axis_distance_au > 0`), 1.0 AU rung (≈ 1 m-class transmitter at
3–5 pc; effectively any listed minimum). Source blended with the star
plus the ≤ 10″ station annulus. Window per event:
t_ca ± √(r² − b²)/v⊥ from the `wise_v1` events table (b, v⊥ computed
with the spacecraft observer).

## Bands and wavelength **[WISE]**

W1 (3.35 µm) and W2 (4.60 µm) primary — full era. W3 (11.56 µm) and
W4 (22.1 µm) opportunistic: cryogenic era only (MJD ≤ 55414 full-cryo,
W3 to ≈ 55469), which intersects at most a handful of events.
Hypothesis flux is a monochromatic line at the band effective
wavelength converted to per-band Vega magnitude via the band width —
the generic in-band-leakage interpretation. The canonical SGL link
wavelengths (1064, 1550 nm) are *not* covered by any WISE band and
remain unconstrained; what WISE adds is the first crossings test of
mid-IR leakage (3–5 µm), a regime no optical survey touches.

## Duty cycle

- **d = 1 while geometry holds** — primary (7.7 s W1/W2 exposures;
  pulse periods ≪ 8 s integrate to a flat ×d scaling, reported without
  re-search).
- **Declared unconstrained:** pulse periods between ~8 s and the
  window length; transmitters scheduled to avoid Earth-crossing
  windows; and — because only one rung survives the elongation gate —
  every beam radius except ~1 AU-class at Earth.

## Detection construction **[WISE: blended A, L1b substrate]**

Matched-filter forced photometry on L1b `int` cutouts at the per-epoch
propagated star position (linear PM fit through the era events'
tabulated star positions), the v2 total-flux kernel and fatal-mask
template (`WiseExactFootprint.FATAL_MASK`), per-frame `magzp` flux
scale. Statistic = window-locked excess vs the same star's off-window
baseline (robust mean, 3×3σ clip), weighted stack over covered
windows, one-sided positive; empirical variance rescale
k = median(r²/v)/0.4549 (floor 1) from the clipped off-window epochs
(the v1 confusion-floor rule — WISE is where it was learned);
WEIGHT_CAP 20× (v1 lesson 4: one frame can dominate a stack).
Off-window baseline capped at a deterministic ≤ 120-epoch
uniform-in-time subsample per (target, band) (the ZTF v1.1 clip
precedent — WISE baselines would otherwise run to thousands of
frames). Discriminators: (1) phase-lock to the crossing ephemeris;
(2) W1:W2 color anomaly (single-band excess with flat other band —
same-visit W1/W2 frames are near-simultaneous, so unlike PS1 this
discriminator is same-epoch); (3) recurrence at successive covered
windows. Station annulus ≤ 10″ ≈ 1.5 W1 PSF: inside the blended core
at 6″ resolution — subsumed into the blended statistic, not a separate
search.

## Saturation rule (frozen; single-frame L1b saturation W1 ≈ 8.1 /
W2 ≈ 7.0 / W3 ≈ 3.8 / W4 ≈ −0.4 Vega)

    excluded  m_est < E_b        E = {W1 8.6, W2 7.5, W3 4.3, W4 0.1}
    marginal  E_b ≤ m_est < E_b + 0.5
    ok        m_est ≥ E_b + 0.5

m_est: census `W1mag` where present (249/386 CNS5 rows); fallback
W1 ≈ Ks − 0.1 (M-dwarf sequence); W2 ≈ W1 − 0.2; W3 ≈ W2 − 0.1;
declared ±0.3, absorbed in the margin. Nearby M dwarfs are bright in
the mid-IR: the cut is expected to remove *most* of the 45 viable
targets, leaving white dwarfs, latest-M, and T/Y dwarfs — the
searchable population is decided by the cut, not assumed.
Optically-dark objects (late-T/Y) are IR-*bright* for WISE: wise-0855
(W2 ≈ 14) is an ordinary searchable target here, not an empty field.

## Quality masks

Primary: `qual_frame > 0`, `moon_sep ≥ 15°`, and the exact
fatal-mask template at the star position (search stage). Strict
(exceedance re-runs): adds `saa_sep ≥ 5°` and `qual_frame ≥ 10`.
*Amendment 2026-08-24 (coverage stage, pre-search): the originally
frozen `qa_status = A` gate named a value that does not exist in the
merge table (the archive uses Prelim/Reviewed/Final); qa_status is
dropped from the primary mask and its distribution is recorded in the
coverage summary. No signal statistic had been formed.*

## Era event counts (from the frozen `wise_v1` list, pre-coverage)

A 1.0 AU: 1,272 in-era events / 88 targets; 102 / 45 targets pass the
elongation pre-gate (`configs/elongation_gate_v1.json`). True coverage (frames in-window at the star) is
the next stage and is the authoritative gate — the archive's own
pointing law enforces elongation, so coverage intersection against
the L1b frame inventory subsumes the geometric estimate. Thresholds
must not promise a depth the cadence cannot deliver.
