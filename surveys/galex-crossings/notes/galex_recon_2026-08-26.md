---
title: "GALEX/gPhoton reachability recon for a Pipeline B crossings survey"
date: 2026-08-26
status: "recon complete — all routes probed live, anonymous end-to-end; freeze not started"
---

# GALEX/gPhoton recon (plan §5.8 item 4)

Goal: verify the access route for the UV pulse-period crossings
survey — time-tagged GALEX photons (5 ms stamps) over the universal
crossing windows (`crossings/universal_v1/`), channels A (on-star)
and B (antipode). The plan row said "gPhoton photon database via
MAST" and was marked _unprobed_; this recon replaces that with
endpoint facts.

Headline: **fully reachable, anonymous, end-to-end** — the gPhoton
photon database at MAST answers arbitrary-position, arbitrary-time
photon-event queries with millisecond (5 ms-tick) timestamps, and a
real in-window photon pull succeeded on the survey's best unit
(gj-1276 B 2010: 1,695 NUV photons from a 30″ box over the 1.6 ks
visit). Two structural findings reshape the survey relative to the
plan row: (1) the **grazing rungs are structurally uncovered** — the
±0.35 d flat-chord windows at 1.2/2.5 R☉ never coincide with GALEX's
sparse visit cadence, at any of the 10 grazing-family positions —
so this becomes a coverage-without-statistic ledger entry like PTF's
B 1.2 R☉ rung; (2) the searchable cell is the **0.1 AU rung**
(±5.8 d windows): 5 usable in-window units at 4 target-channels,
4 of them with simultaneous FUV — the programme's first two-band
photon-level units, in a 2007–2010 era that predates ZTF and TESS.

## Service probes (2026-08-26, dev machine, all anonymous)

| Probe | Status | Notes |
| --- | --- | --- |
| MAST CAOM cone (`mast.stsci.edu/api/v0/invoke`, `Mast.Caom.Cone`) | **up** | 232 GALEX rows at the van-maanen antipode (AIS tiles, a 1,681 s GII visit, GIS grism spectra); obs metadata route for depth/exposure annotation |
| gPhoton Mashup query service `mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest` | **up** | free-form SQL over the photon DB, `format=extjs` JSON. **Host matters: `mast.stsci.edu` 404s for this path — `mastcomp` only** |
| `fGetNearbyAspectEq(ra, dec, radius_arcmin, t0_ms, t1_ms)` | **works** | returns `htmID, time, band, distance` per aspect second; `band` ∈ {`NUV`, `FUV`, `FUV/NUV`} (**substring-match it** — `'FUV/NUV'` means both detectors on); `distance` = boresight offset in arcmin |
| `NUVPhotonsV` / `FUVPhotonsV` photon views | **works** | columns `zoneID, time, cx, cy, cz, x, y, xa, ya, q, xi, eta, ra, dec, flag`; time-range + ra/dec box WHERE clauses answer in seconds; verified count + row pulls |
| `fGetTimeRanges` | **does not exist** | name from older gPhoton docs — "Invalid object name"; use `fGetNearbyAspectEq` |
| MCAT `GR6Plus7.dbo.photoobjall` / `visitphotoobjall` | **works** | visit-level NUV/FUV mags (probed to NUV 23.8 near ross-128); `imgrun` exists but its column names differ from the docs — pin at freeze |

## Photon-database facts

- **Era (measured, `select min(time), max(time) from aspect`)**:
  2003-06-07 05:02 → 2013-05-01 18:02 UTC — essentially the full
  mission. FUV detector rows end in the 2009-05 failure era (the
  2010 in-window visit is NUV-only, as expected).
- **Time system**: "GALEX time" = UNIX − 315,964,800 s (GPS epoch);
  DB `time` columns are **milliseconds** of GALEX time. Photon
  stamps land on **5 ms ticks** (verified on pulled rows).
- **Aspect cadence** 1 s; a "visit" reconstructs as a contiguous
  aspect-time run (gap > 2 s splits). AIS visits ≈ 100 s; MIS/GII
  ≈ 1.5–1.7 ks (one eclipse).
- **Photon pull scale (measured)**: 1,695 NUV photons in a 30″×30″
  box over the 1,637 s gj-1276 B 2010 visit ≈ **1.0 ct/s of
  sky+source background** in an aperture-scale box — Poisson
  counting statistics for a pulse search are cheap to compute
  directly from DB pulls; no image products needed.
- `flag = 0` on all pulled photons; flag semantics not yet pinned
  (open item).
- Detector geometry: FOV ~1.25° ∅; boresight `distance` ≤ 33′
  adopted at recon as the usable-detector cut (rim artifacts beyond;
  freeze revisits).

## Coverage probed (fGetNearbyAspectEq, 37.5′ discovery radius)

Full GALEX-era scan of the universal list: grazing family
(≤ 2.5 R☉: van-maanen, wolf-359, teegarden, gj-1276, ross-128 —
19–20 events per target-channel in-era, including the van-maanen
b = 0.27 R☉ April/October families) plus the two extra 0.1 AU-rung
targets (gj-908 b_min 0.057 AU, ross-154 0.015 AU). Aspect visits
per position:

| Position | Visits | Total s | Years |
| --- | --- | --- | --- |
| van-maanen star / antipode | 2 / 6 | 216 / 1,795 | 2006 / 2004+2008 |
| wolf-359 star / antipode | 5 / 4 | 5,193 / 373 | 2006–10 / 2004+2007 |
| teegarden star / antipode | 2 / **0** | 204 / 0 | 2007 / — |
| gj-1276 star / antipode | 12 / 115 | 2,131 / 10,824 | 2003–11 / 2006–10 |
| ross-128 star / antipode | 21 / 7 | 5,138 / 3,238 | 2004–09 / 2003–08 |
| gj-908 star / antipode | 10 / 5 | 1,942 / 2,196 | 2006–08 / 2004+2011 |
| ross-154 star / antipode | 3 / 4 | 285 / 308 | 2006 / 2006–07 |

### In-window scan (flat-chord windows, strict ≤ 33′ cut)

**1.2 R☉ and 2.5 R☉ rungs: zero in-window seconds, all targets,
both channels.** The ±0.3–0.4 d grazing windows never meet the visit
cadence — the van-maanen deep-graze family in particular has visits
(2004 AIS, 2008 GII) but none in any window. Grazing rungs enter the
ledger as coverage-without-statistic.

**0.1 AU rung (windows ±5.8 d): 5 usable units at 4
target-channels** (`results/era_scope_v0.json`):

| Unit | Visit (UTC) | Offset from t_ca | NUV s | FUV s | Boresight dist |
| --- | --- | --- | --- | --- | --- |
| wolf-359 A 2007-03-03 | 2007-02-26 11:21 | −5.3 d | 97 | 97 | 22.1′ |
| gj-1276 A 2007-09-05 | 2007-09-01 01:52 | −4.1 d | 92 | 92 | 24.1′ |
| gj-1276 B 2007-03-03 | 2007-02-26 11:19 | −4.8 d | 109 | 109 | 30.2′ |
| gj-1276 B 2010-03-03 | 2010-03-06 21:57 | +3.9 d | **1,637** | 0 | 23.8′ |
| ross-128 A 2007-03-17 | 2007-03-22 15:40 | +4.7 d | 110 | 110 | 17.6′ |

Edge case: gj-908 A 2007-09-21 has a 109 s visit in-window but at
36.3′ minimum boresight distance — outside the 33′ cut (detector
rim). The freeze decides whether a rim-response argument admits it;
default is exclusion. teegarden and van-maanen contribute nothing;
ross-154 nothing.

## What this means for the survey design

1. **The cell is the UV pulse-period cell at the 0.1 AU rung** —
   5 ms photon stamps over 92–1,637 s in-window exposures, in
   2007–2010 (only PS1 and PTF otherwise cover this era, neither at
   sub-second cadence). Four units carry simultaneous FUV: a
   two-band coincidence test on any pulse candidate, new to the
   programme.
2. **Substrate = direct photon-event counting from the DB** — no
   images, no cutouts. A pulse/burst statistic on binned photon
   arrival times with local (annulus or off-position) background is
   the natural construction; calibrated flux limits need the
   effective-area + flat chain (gPhoton's gAperture or gPhoton2's
   local raw6 pipeline are the fallback routes if absolute
   calibration is required — decide at freeze, the pulse statistic
   may need only relative counts).
3. **Small survey**: ~5 units / 4 target-channels — PTF-scale
   (6 units / 10 trials), so the PTF confirmatory machinery
   (per-unit trials, frozen thresholds, pseudo-window controls
   drawn from the same visit off-window... n.b. most visits are
   only ~100 s ≈ window ≪ visit is impossible — **pseudo-position
   controls on the same visit** are the natural control family,
   not pseudo-windows) transfers with one structural change.
4. **On-star channel A units (ross-128, wolf-359, gj-1276, all
   M dwarfs) have a known astrophysical false-positive: NUV/FUV
   flares.** The stellar PSF sits in the same aperture-scale region
   as the predicted axis position. The PTF gj-1276 lesson (on-star
   temporal controls defeated by campaign cadence) recurs here in
   sharper form — the freeze must declare a stellar-flare veto
   (e.g. flare morphology + position centroid vs the axis offset)
   before any data contact.
5. **Observer**: GALEX is ~700 km LEO — `universal_v1` Earth-center
   events are valid under the standing 0.010 R☉ / ~4 min budget
   (the WISE/SPHEREx convention). No new observer list needed.
6. Runs on the dev machine fine — DB queries answer in seconds,
   total in-window photon volume is ~10⁴ events. No
   server-migration dependency.

## Open items before the freeze

- Pin photon `flag` semantics and the usable-photon rule; pin the
  boresight-radius cut (33′ adopted ad hoc) and the rim question for
  gj-908 A.
- Per-event probe positions: recon used first-event axis positions
  per target-channel (drift over the era is ~0.01–0.05°, inside the
  discovery radius); recompute exactly per event at freeze.
- Decide the statistic's calibration route: relative counting
  (Poisson, local background) vs gAperture/gPhoton2 absolute chain;
  and the pulse-period grid (5 ms floor → window length) with the
  trials accounting.
- Stellar-flare veto design for the three on-star units; check MCAT
  variability records at those stars.
- Dead-time and hotspot masks (NUV hotspot positions), aspect
  `flag` column (not yet examined), and whether photon `ra/dec` are
  aspect-corrected (they appear to be — pulled photons cluster at
  sky positions, not detector coordinates).
- Pin `imgrun`/exposure-table column names for visit metadata;
  snapshot the Mashup responses and the gPhoton docs at freeze.
- Barycentric timing: spacecraft-UTC vs TDB drift is irrelevant for
  window membership (±5.8 d) and ≲ ms-scale over a 1.6 ks visit —
  record as a declared budget, not a correction.
