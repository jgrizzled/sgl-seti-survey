---
title: "Rubin Science Platform (DP2 Early) recon for both pipelines"
date: 2026-08-26
status: "recon complete — service probes + geometry intersects; no survey started"
---

# Rubin DP2 recon (plan §6 → unblocked; access granted, RSP_API_TOKEN in .env)

Goal: DECam-style reachability recon of the Rubin Science Platform now
that data-platform access is granted, and a geometry-only intersect of
DP2 (Early) against the universal corridor list (Pipeline A) and the
universal crossing list (Pipeline B), to decide what Rubin survey(s)
to define and when. No pixel or per-row catalog data were used for any
hypothesis-relevant statistic — boresight cones, aggregate counts and
service probes only.

Scripts: `surveys/rubin/scripts/recon_{fetch_visits,intersect,probe_catalogs,probe_images}.py`;
results: `surveys/rubin/results/recon_*_v0.json`; bulk snapshots in
`runs/rubin/recon/` (28,698-row `dp2.Visit` snapshot + probe cutouts).

## Service probes (2026-08-26, dev machine; all authenticated)

| Service | Status | Notes |
| --- | --- | --- |
| TAP sync (`data.lsst.cloud/api/tap/sync`) | **works** | Bearer header on raw curl, or pyvo `AuthSession` with `x-oauth-basic` basic auth. ~1 min server timeout; aggregate GROUP BY over 5.1M-row VisitDetector in 7.5 s |
| TAP async (UWS `/async`) | **works** | pyvo `submit_job` → COMPLETED; 28,698-row full-table pull, snapshot-able |
| Spatial ADQL (Qserv) | **works** | `CONTAINS(POINT('ICRS',…), CIRCLE('ICRS',…)) = 1`; `dp2.Object` uses `coord_ra/coord_dec`, Visit/VisitDetector/DiaSource use `ra/dec` |
| SIA2 (`/api/sia/dp2/query`) | **works** | POS/DPSUBTYPE/MAXREC verified; returns ObsCore rows with DataLink `access_url` |
| DataLink (`/api/datalink/links`) | **works** | per-dataset links: `#this` = signed Google Cloud Storage URL for the full file (time-limited signature — fetch promptly or re-resolve; full deep_coadd is cell-based `lsst_cells_v2`) + 3 cutout service descriptors |
| SODA cutout (`/api/cutout/sync`) | **works** | `ID` + `CIRCLE`/`POLYGON` + `CUTOUTDETAIL` ∈ {Image, MaskedImage, Exposure}; sync = 303 redirect to result (follow redirects). 2′ Exposure cutout: 4.3 MB in 2.0 s |
| Schemas visible | — | `dp2`, `dp1`, `dp02_dc2_catalogs` (simulated DC2), `ivoa` (ObsCore) |

Token scopes exercised: `read:tap`, `read:image`. Quotas/rate limits
undocumented; nothing throttled at recon volumes.

**The Exposure-detail cutout is the headline product fact**: one call
returns IMAGE + MASK + VARIANCE + the spatially-varying PSF model grid
+ aperture corrections + background model + PROVENANCE (the input
visit list). That is the complete v2-engine substrate (pixels, masks,
variance, PSF) in a single request — no separate `-msk`/dqmask/wtmap
fetches as in WISE/DECam, and no empirical-PSF fallback needed.

## DP2 (Early) contents

- LSSTCam, Rubin pipelines v30 (processing run `DM-55060`), obs era
  **MJD 60790.1 → 61047.4** (2025-03-26 → 2025-12-08); 28,698 visits
  (u 1,964 / g 4,166 / r 4,856 / i 7,959 / z 5,809 / y 3,944), exposure
  30 s (93 %) or 38 s. Coadds + catalogs over ~3,000 deg²; single-epoch
  catalogs over ~15,000 deg².
- **Images in ObsCore: deep coadds only** (925,460). No visit or
  difference images until late 2026 (plan §6 expectation unchanged).
  DP1 (ComCam, 7 fields, late 2024) has the full image suite
  (raw/visit/difference/template/coadd) — but its fields are ≥ 7° from
  every corridor antipode and no in-era crossing event position lands
  within 2° of any field: **DP1 = adapter-development sandbox only**.
- Catalogs (key columns for us):
  - `Visit` (boresight, `expMidptMJD`) — snapshot taken;
  - `VisitDetector` (5.15M rows) — **per-detector corner coordinates**
    (`llcra…urcdec`: exact detector footprints without touching
    pixels), per-detector `magLim`, `seeing`, `zeroPoint`;
  - `DiaSource` (difference-image detections): real/bogus
    `reliability`, `scienceFlux`/`templateFlux`, trail fits,
    `pixelFlags_injected*` (archive-side injections are flagged),
    `forced_PsfFlux`;
  - `ForcedSource`: per-visit forced PSF photometry at every Object
    position on **both** the visit image (`psfFlux`) and the
    difference image (`psfDiffFlux`) — verified 48/48 epochs with
    diff flux on a sample ross-128-antipode object;
  - `DiaObject`, `ForcedSourceOnDiaObject` (same, keyed to DiaObjects);
  - `SSObject`/`SSSource` + replicated MPC tables — 299,343 known
    solar-system objects with per-detection links = ready-made
    positive-control machinery;
  - `IsolatedStarStellarMotions` (PM/parallax fits vs Gaia DR3),
    `Object` (1,248 cols), `ShearObject`, `CoaddPatches`.
- Depth (per-detector `magLim` averages, whole DP2): u 23.4 / g 24.2 /
  r 23.6 / i 23.3 / z 22.7 / y 21.7; mean seeing 1.17–1.36″. Single
  epochs ~1–1.5 mag deeper than ZTF; comparable to DECam long
  exposures but survey-uniform cadence and calibration.
- Sanity note: catalog probes at the ross-128 antipode 6′ cone —
  6,786 Objects, 436 DiaObjects, 1,190 DiaSources (mean reliability
  0.06–0.21: bogus-dominated, cuts required as with ZTF).

## Pipeline A intersect (visit boresight within 1.5° of antipode)

19/76 universal-list systems covered (`recon_intersect_v0.json` has
the full table). Headliners:

| corridor | N visits | 30-d bins | bands |
| --- | --- | --- | --- |
| Ross 128 | 45 | 3 | ugrizy |
| Wolf 359 | 33 | 3 | ugrizy |
| GJ 1111 | 29 | 2 | ugrizy |
| GJ 625 | 21 | 2 | g i r y z |
| GJ 687 | 19 | 3 | g i r y z |
| Struve 2398 AB | 16 | 3 | g i r y z |
| σ Dra | 14 | 3 | g i r y z |
| GJ 1221 | 12 | 3 | g i r y z |
| GJ 251 | 11 | 2 | ugrizy |
| Wolf 1069 / GJ 13157 | 9 / 8 | 2 | — |

plus 8 more at 1–2 visits. The southern DECam family (GJ 625, GJ 687,
Struve 2398, σ Dra, GJ 1221, Wolf 1069, GJ 13157, 61 Cyg) is already
accumulating Rubin epochs at DECam-class depth with survey-uniform
processing. Calendar spread is only 1–3 30-day bins per corridor —
Early DP2 alone gives limited parallax-phase leverage; that arrives
with Year-1/DR1.

## Pipeline B intersect (in-window visits within 1.5° of event position)

Windows = per-event `beam_radius / v_perp` half-widths (reproduces the
frozen ±0.35 d grazing / ±5.8 d 0.1 AU values); DP2 era as above.

| rung | in-era events | events with in-window in-cone visits |
| --- | --- | --- |
| B 1.2 R☉ | 5 | 0 |
| B 2.5 R☉ | 6 | 0 |
| B 0.1 AU | 9 | **1** — ross-128, t_ca MJD 60937.94, **b = 1.86 R☉**, 1 r visit in ±5.8 d |
| A 0.1 AU | 9 | **1** — ross-154, t_ca 60859.41, b = 3.4 R☉, 9 visits (i r y z) |
| A 1.0 AU | 178 | 36 (one event each, 36 targets) |

The ross-128 hit is notable: a **grazing-family event** (b < 2.5 R☉)
with an in-window Rubin visit — the ±0.35 d grazing windows themselves
caught nothing (expected at ~30 s/visit scheduler cadence), but the
wider rung window contains a difference-image epoch at r ~ 24 depth on
a b = 1.86 R☉ crossing. Nothing this deep exists in any covered-window
ledger row.

## Assessment — what surveys this recon supports

Two projects, staged:

1. **Now — catalog-level DP2 survey (both pipelines), TAP-only
   adapter.** Substrate = `DiaSource`/`DiaObject` screening +
   `ForcedSource`/`ForcedSourceOnDiaObject` per-visit forced
   difference photometry; depth r ~ 24 single-epoch. Pipeline B first
   (the ross-128 b = 1.86 R☉ + ross-154 windows are frozen-rung
   events; small unit count fits the v2 discipline), Pipeline A
   screening second (19 corridors, phase-limited). Structural caveat
   to freeze around: **no arbitrary-position forced photometry** —
   ForcedSource exists only at Object/DiaObject positions, so the
   §3.4 catalog-absence rule applies; the decision rule must be
   built on association within a matching radius, with the
   completeness cost measured (or the unit declared
   coverage-without-statistic where no Object anchor exists).
   Archive-side `pixelFlags_injected` must be cut explicitly.
2. **Late 2026 — image-level v2 survey** when visit + difference
   images publish: the Exposure-cutout product shape means the full
   v2 chain (own injections, exact masks, PSF in hand) with the
   simplest fetch layer in the programme. DP1's difference images are
   the adapter-development target in the meantime.

Also unlocked for existing maintenance items: `SSObject`-based
positive controls, and `IsolatedStarStellarMotions` as an external
astrometric cross-check.

## Open items before any freeze

- Decide the association radius + completeness treatment for the
  no-arbitrary-position-forced-photometry caveat (drives the whole
  catalog-level decision rule).
- `reliability` cut calibration: bogus-dominated cones (mean 0.06–0.2)
  need a frozen threshold with a measured efficiency — check whether
  DP2 publishes the real/bogus training characterization.
- Visit-level saturation/bright-limit for the on-star channel A units
  (ross-154 is V ~ 10.4: likely saturated in 30 s — the DiaSource
  route may be dead on-star; check `pixelFlags_saturatedCenter`
  behavior near bright stars before promising channel A).
- Quotas/fair-use for a ~10⁴-query screening pass (nothing throttled
  at recon volume; ask/observe before scale-up).
- Whether Early DP2 will be superseded in-place by later DP2
  processings (freeze must pin `DM-55060`/v30 provenance; re-check at
  freeze time).
- Era end: visit table stops 2025-12-08 but the release notes say
  "to January 2026" — confirm whether more visits appear in this
  release before freezing the era bounds.
