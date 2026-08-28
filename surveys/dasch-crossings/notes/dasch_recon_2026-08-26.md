---
title: "DASCH DR7 reachability recon for a Pipeline B crossings survey"
date: 2026-08-26
status: "recon complete — all routes probed live, anonymous end-to-end; crossings-list backward extension built the same day; freeze not started"
---

# DASCH recon (plan §5.8 item 5)

Goal: verify the access route for the pre-1980 century of
annually-recurring crossing windows — Harvard scanned plates,
1885–1992, DR7 — and execute the plan row's prerequisite: extend the
universal crossing list backward from its 1980 start
(`crossings/universal_1885_v1/`, built 2026-08-26; §"Crossings
backward extension" below), exercising the §7 model-accuracy budget
at old epochs. The plan row said "DASCH lightcurve + cutout services
at Harvard" and was marked _unprobed_; this recon replaces that with
endpoint facts.

Headline: **fully reachable, anonymous, all-sky, and far richer than
the plan row assumed.** The DR7 access layer is the Starglass REST API
(five JSON POST endpoints — exposure lists, catalog sources,
century-long lightcurves *with per-epoch non-detection limits*,
FITS cutouts, and per-plate subregion photometry including
uncatalogued detections). Every probed grazing-family position (5
stars + 5 antipodes) has ~9,000–15,000 overlapping exposures,
~2,300–3,100 of them photometrically calibrated, 150–340 reaching
B ≥ 15. And the science case is stronger than expected: **every
annual van-maanen crossing of the DASCH era is a photosphere-class
graze (b ≈ 0.59–0.69 R☉ through the century)**, with ~18–28 windows
per channel holding calibrated exposures within ±1 d — the largest
grazing-window recurrence sample in the programme by an order of
magnitude (PTF: 6 events at ±3 d; TESS: 2 resolved grazes).

## Service probes (2026-08-26, dev machine, all anonymous)

Base URL `https://api.starglass.cfa.harvard.edu/public/` (no key;
"lower" rate limits — never hit during ~20 recon queries at ~10 s
each). A registered tier `/full/` exists (free Starglass account →
40-char `x-api-key`) with higher limits — the ATLAS-style `.env`
token pattern if scale demands it. **WAF note:** requests with
python-urllib's default User-Agent get HTTP 403; any curl-like UA
passes.

| Probe | Status | Notes |
| --- | --- | --- |
| `GET /health` | **up** | `{"status": "ready"}` |
| `POST /dasch/dr7/queryexps` `{ra_deg, dec_deg}` | **works** | all exposures whose footprint covers the point; 27-column CSV-in-JSON (list of strings, first = header); 9,000–15,000 rows/position, ~10 s, ~2 MB |
| `POST /dasch/dr7/querycat` `{ra_deg, dec_deg, radius_arcsec, refcat}` | **works** | refcat sources in box; carries `gsc_bin_index` + `ref_number` (the lightcurve key), PM columns, `pos_epoch` |
| `POST /dasch/dr7/lightcurve` `{gsc_bin_index, ref_number, refcat}` | **works** | van Maanen: 3,544 rows, 1890–1989, 64 columns incl. `magcal_magdep`, `magcal_local_rms` (median 0.168 — matches the documented ~0.15), `limiting_mag_local`, `time_accuracy_days`, AFLAGS/BFLAGS quality bits. 1,552 detections / 1,992 non-detections; **every non-detection row carries `limiting_mag_local`** (median 11.6 — mostly shallow plates; 266 deeper than the B = 12.6 star, the "missing points" caveat set) but **none carries `time_accuracy_days`** — timing for non-detections joins from the exposure list. Every detection has nonzero `aflags` (bit semantics needed — not a reject boolean) |
| `POST /dasch/dr7/cutout` `{plate_id, solution_number, center_ra_deg, center_dec_deg}` | **works** | base64(gzip(FITS)); 835×835 px ≈ 20′, 1.44″/px, `RA---TAN` WCS, target lands 0.9″ from centre. Header is minimal: no DATE-OBS/exptime (metadata rides with the exposure row); CD-diagonal-only suggests a resampled north-up product — pin at freeze |
| `POST /dasch/dr7/platephot` `{plate_id, solution_number, center, refcat}` | **works** | calibrated SExtractor detections in a plate subregion, lightcurve-row schema; **includes uncatalogued detections** (5/50 blank `ref_number` on the probe) — a catalogue-level antipode transient search needs no pixel work. Probe returned exactly 50 rows — cap/subregion size unverified (open item) |
| `POST /dasch/dr7/mosaic_package` `{plate_id, binning ∈ {1,16}}` | **works** | pre-signed S3 URLs (15-min expiry) for the fpack'd full-plate mosaic; binned-16 available for cheap whole-plate QA. Full mosaics average 750 MB — cutouts are the working substrate |

OpenAPI specs (public + registered) at
`https://starglass-api-documentation.s3.amazonaws.com/prod/{public,registered}-api.json`
— snapshot at freeze. Docs: `dasch.cfa.harvard.edu/dr7/` (web-apis,
lightcurve-columns, exposurelist-columns, known-issues pages fetched
this session). `daschlab` 1.0.0 on PyPI is the official client; the
adapter will use the REST endpoints directly per house pattern
(snapshots, content hashes).

## DR7 facts

- **Corpus**: ~430,000 plates scanned (of >550,000 in the
  collection), ~1880–1990 (measured at our positions: 1885–1989);
  mosaics at 11 µm; ~400 TB total. ~97% of mosaics have astrometric
  solutions (Astrometry.Net + distortion refinement); ~89% of those
  have photometric calibrations.
- **Photometry**: calibrated against two refcats. **APASS DR8 → B
  band, "excellent long-term stability" — the science choice**;
  ATLAS-refcat2 (g) has *documented false long-term trends*
  (chromatic mismatch) and is preferred only for astrometry. ~24
  billion magnitudes / ~250 million sources; typical lightcurve RMS
  ~0.15 mag; typical depths B 14–16, deep plates to B ≈ 18.3
  (probed: plate a26884, 1949, lim 18.3).
- **Exposure-list columns of use**: `series platenum solnum expnum`
  (plates can carry multiple exposures — the row unit is the
  *exposure*), `expdate` (UTC, real HH:MM resolution — 1,120
  distinct times at one position), `exptime` **in minutes** (median
  60, p10 6, p90 113; multi-hour and multi-night exposures exist,
  max 3,167 min), `limMagApass/limMagAtlas` (per-exposure local
  limiting mag; present only where photometrically calibrated),
  `medianColorterm*`, `wcssource` (`imwcs` = solved / `logbook` =
  position from logbooks only — not photometry-usable), `centerdist
  edgedist` (deg). The `epoch` column is the coordinate equinox
  (2000.0), **not** the date.
- **Timing**: lightcurve `time_accuracy_days` mode 0.0007 d ≈ 60 s
  (1,384/3,544 van Maanen epochs); a tail of 1.0 d (logbook date
  only). Times are plate UTC — barycentric/light-time handling is
  ours, and the docs warn the correction varies across wide plates.
  **Window-overlap must use the [start, start+exptime] interval, not
  the start epoch** — a 60-min exposure is a meaningful fraction of a
  grazing chord window (±0.3–0.8 d).
- **Known issues** (dr7/ki/, all directly relevant to a
  single-detection transient search): missing lightcurve points
  (non-detection claims need image checks), source splitting
  (merge close catalog entries), undetected plate defects ("science
  results hinging on a single detection should include detailed
  vetting … potentially up to physical examination of the plate"),
  undetected blends (low-resolution series), incorrect WCS on
  meteor/patrol wide-field series. These set the vetting/annotation
  design for any exceedance.
- High-PM handling: querycat rows carry PM and `pos_epoch`; van
  Maanen (µ ≈ 3″/yr → ~5′ over the century) resolves to a single
  catalog source with a 3,544-point lightcurve, so DASCH's
  PM-aware matching works at our µ range — verify per target at
  coverage (source splitting is the known failure mode).

## Coverage probed (grazing family, stars + antipodes)

APASS-calibrated exposures per position (counts; `lim ≥ 15` in
parentheses):

| Target | A = star | B = antipode |
| --- | --- | --- |
| van-maanen | 2,868 (290) | 2,991 (334) |
| wolf-359 | 2,717 (213) | 2,440 (255) |
| teegarden | 3,018 (276) | 2,600 (337) |
| gj-1276 | 2,436 (260) | 2,732 (199) |
| ross-128 | 3,119 (231) | 2,259 (146) |

Per-decade profile (all positions alike): ramp from the 1890s,
peak 1930s–40s (~60–85 calibrated exposures/yr/position), the
**Menzel gap 1953–1968** (1960s: 3–52 per decade), partial recovery
1970s–80s, end by 1990. DASCH is genuinely all-sky (Arequipa /
Bloemfontein southern stations) — both hemispheres covered at every
probed position.

## In-era window scan (van-maanen prototype; full family below)

Trial single-target crossings run 1885→1993 (432 events; identical
construction to `universal_1885_v1`): **every annual event of both
channels sits at b ≈ 0.59–0.69 R☉** — van-maanen's axis grazes the
photosphere continuously through the DASCH century (the modern-era
deep-graze family, b 0.078–0.37 R☉, is the tail end of a secular
deepening). Against the calibrated exposure lists:

| Channel | era events | covered ±1 d | ±3 d | exposures in ±3 d |
| --- | --- | --- | --- | --- |
| van-maanen A (Oct family) | 216 | 18 | 48 | 118 |
| van-maanen B (Apr family) | 216 | 28 | 46 | 115 |

The ±1 d numbers are the realistic proxy for the 1.2 R☉ flat-chord
window (±0.28 d at b ≈ 0.6); the 0.1 AU rung (±5.8 d) will hold
nearly all ±3 d windows. Even at one-third attrition this is a
**~10–30-window recurrence stack per channel at photosphere-grazing
impact parameters** — the cell the PS1 masks and TESS sector gaps
destroyed, now with a century of independent recurrences.

### Addendum: full-family scan over `universal_1885_v1` (same day)

Build landed (`xng-e2d1063af9d0`: 37,096 events, 0 invalid, 91 MB;
13,053 degraded — dominated by the pre-1941 `long_propagation_span`
annotation). **Overlap consistency check passed**: all 4,104 events
in the shared era (1980-06 → 1992-06) match `universal_v1` 1:1 with
max |Δt_ca| = 23 s and max |Δb| = 3.4×10⁻⁷ R☉ — the backward
extension reproduces the canonical list exactly where they overlap.

Era census (1885–1993): the ≤ 1.2 R☉ rung exists for **four
targets** — van-maanen (A and B: *all* 216 events per channel),
gj-1276 (likewise all 216 per channel), wolf-359 (A 152 / B 184),
teegarden (B only, 38) — and ≤ 2.5 R☉ adds ross-128 (all events).
Probed-position intersection (calibrated exposures, ±1 d / ±3 d
proxies):

| Target-channel | 1.2 R☉ events | ±1 d | ±3 d |
| --- | --- | --- | --- |
| van-maanen A | 216 | 18 | 48 |
| van-maanen B | 216 | 28 | 46 |
| gj-1276 A | 216 | 21 | 34 |
| gj-1276 B | 216 | 21 | 45 |
| wolf-359 A | 152 | 15 | 32 |
| wolf-359 B | 184 | 18 | 35 |
| teegarden B | 38 | 0 | 2 |
| (2.5 R☉ rung adds) teegarden A / ross-128 A / ross-128 B | 216 each | 20 / 21 / 17 | 41 / 40 / 35 |

**Three targets carry persistent photosphere-grazing annual families
through the entire DASCH century**, with ~15–28 covered windows per
target-channel at ±1 d (~120 across the family) — roughly an
order of magnitude more grazing-rung windows than every previous
survey combined. The 0.1 AU rung tracks the ±3 d numbers
(~35–48 per channel). Definitive flat-chord windows, exposure-
interval overlap, and the timing-accuracy gate move to the coverage
stage.

## Crossings backward extension (`crossings/universal_1885_v1`)

Built 2026-08-26 per the plan row's prerequisite:
`python -m sglsurvey.crossings --observer earth --start 1885-01-01
--stop 1993-01-01 --out crossings/universal_1885_v1` (Earth-center —
valid for ground plates; same registry, model
`tusay2022_eq5_7_v1`, every-minimum construction, no beam-radius
cut). The window overlaps `universal_v1` over 1980–1992 by design:
the shared era is a consistency check on the extension (see
addendum), and the DASCH survey intersects this one self-contained
product rather than stitching two lists at a seam inside the archive
era.

Old-epoch validity: events before ~1941 carry sglseti's
`long_propagation_span` annotation (`validity = degraded`) — the
declared 75-yr linear-motion bound
(`LINEAR_PROPAGATION_SPAN_YEARS`, providers.py) measured from the
2016.0 catalog epoch, direction-dependent via the light-time
retardation. This is conservative bookkeeping, not a physical
limit; the measured budget is below.

## §7 model-accuracy budget at 1885 (measured this session)

- **Earth ephemeris**: astropy `builtin` (ERFA epv00) vs JPL DE440S
  (which starts 1849): heliocentric Earth position differs by
  **≤ 6 km ≈ 1×10⁻⁵ R☉** across 1885–1992 (the ~130 km difference in
  the SSB frame cancels in the Sun-relative geometry the axis uses).
  Negligible.
- **Target astrometry**: propagation is rigorous space motion
  (`apply_space_motion` → ERFA pmsafe), so perspective acceleration
  is modelled, not an error term. Numerical sensitivity at 1885 for
  the grazing family: a generous 5 km/s RV error moves the star
  direction 26–171 mas (worst wolf-359); a 0.1 mas/yr µ error adds
  ~13 mas. Combined axis displacement at Earth: **≤ 0.0002 R☉**,
  i.e. 50× below the standing 0.010 R☉ LEO observer budget and
  ~3,000× below the 1.2 R☉ rung; t_ca shifts ≲ 5 s. Gaia DR3 errors
  are far smaller than these test perturbations.
- **Unmodelled orbital motion** stays the one open annotation for
  century baselines (registry rationale "negligible at
  WISE-baseline" was not vetted for 131 yr); the grazing family is
  isolated (van Maanen isolated WD, Gaia RUWE 1.21), so this is a
  per-target check at freeze for any candidate unit outside that
  family.
- **DASCH-side timing** (`time_accuracy_days` up to 1.0 d for
  logbook-dated plates) dominates every geometry term — a coverage-
  stage gate (e.g. require time accuracy ≪ window half-width),
  not a geometry-model problem.

## What this means for the survey design

1. **The cell is the century-scale grazing recurrence stack** —
   van-maanen A and B are permanent photosphere-grazing families
   (~10–30 covered windows per channel), with the other grazing
   targets' rungs to be enumerated by the full scan. This is the
   deepest extension of the §5.10/§5.9 van-maanen story: PTF added
   2009/2011/2013; DASCH adds ~a century of independent annual
   recurrences ending where the 1980 list begins.
2. **Substrate = lightcurve/platephot-first, cutouts for vetting.**
   Channel A (on-star) rides the lightcurve endpoint directly
   (calibrated mags + per-epoch limits at the star). Channel B
   (antipode) rides platephot (calibrated detections incl.
   uncatalogued ones + `limiting_mag_local`) — image-level forced
   photometry only if the dev stage shows the catalogue level is
   leaving depth on the table. This is a *new substrate shape* for
   the programme (no prior survey searched at catalogue level as
   primary), and the frozen statistic must be designed for it.
3. **B ~ 12–16 typical / 18 best depth** frames the power cell:
   photosphere-grazing (≤ 2.5 R☉ cone) at B ~ 15 is ~kW-class
   relay power sensitivity era-independent of anything else in the
   programme — and the *pre-1941* windows are the first coverage of
   any kind before the PS1 era, closing backward to 1890.
4. **APASS refcat for science, ATLAS for astrometry** (documented
   false trends in ATLAS-refcat2 photometry). Colorterm columns
   exist per exposure — emulsion heterogeneity is annotatable.
5. **Exposure-interval overlap, not epoch matching**: 60-min median
   exposures vs ±0.3–0.8 d chord windows means in-window integrated
   fractions are computable per exposure — a duty-cycle-friendly
   statistic (and multi-hour exposures can contain a *complete*
   ingress→egress at the 1.2 R☉ rung).
6. **Runs on the dev machine.** ~10 s/position queries, MB-scale
   pulls, anonymous. ~176 positions × a handful of endpoint calls is
   an evening, not a queue drain (contrast ATLAS §5.9). The
   registered API key is an optional accelerant.
7. **Single-detection vetting is the statistic's hard part**: DASCH's
   own docs put the burden (defects, blends, missed points) on the
   user. The exceedance-adjudication design must budget for
   cutout-level (and worst-case plate-level) vetting of every
   over-threshold window, and the "missing lightcurve points" issue
   means *non-detection-based* coverage claims need spot image
   checks — fold into the coverage-stage QA sample.

## Open items before the freeze

- Pin AFLAGS/BFLAGS/`local_bin_reject_flag`/`plate_quality_flag` bit
  definitions and a fatal-bit template (dr7/catalog-columns and
  flags docs); decide the quality gate for "calibrated exposure".
- Resolve the platephot 50-row question (cap? subregion size?) and
  whether it serves arbitrary positions on solved-but-uncalibrated
  exposures.
- Pin the cutout product's resampling status (CD-diagonal TAN — is
  it interpolated? affects any pixel-level photometry and defect
  vetting).
- Decide the series/class policy: exclude `logbook`-positioned
  exposures; decide on meteor/patrol wide-field series (known bad
  WCS risk) and multi-exposure plates (`expnum > 0` semantics).
  Probed-position series census (calibrated fraction varies 1–77%):
  `ai` dominates raw counts (47k, 12% calibrated), then `fa` (15k,
  10%), `ac` (11k, 49%), `am` (6.5k, 60%), `dnb` (1.9k, 71%) — the
  per-series depth/quality split belongs in the coverage stage.
- Full-family flat-chord window scan over `universal_1885_v1` →
  the searched-unit list and rung structure (the PTF §9 pattern).
- Timing gate: choose the `time_accuracy_days` threshold per rung
  (1.0 d logbook dates are useless for 1.2 R☉ chords but fine for
  the 0.1 AU rung); verify barycentric convention end-to-end at
  one epoch.
- Positive control: pick a known variable (or the RY Cnc tutorial
  star) through the full chain; for astrometry/motion, a bright
  asteroid with a DASCH-era Horizons ephemeris if platephot serves
  it.
- Snapshot the OpenAPI specs, docs pages, and a full per-position
  queryexps response set at freeze (content hashes per house rule).
- Registered-API key: create the Starglass account only if public
  rate limits bite at coverage scale (`.env` pattern per ATLAS).
