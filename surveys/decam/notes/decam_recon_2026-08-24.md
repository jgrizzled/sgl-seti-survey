---
title: "DECam/NOIRLab recon for a Pipeline A southern survey"
date: 2026-08-24
status: "recon complete — service probes + coverage sweep; pilot not started"
---

# DECam recon (for the survey now at plan §4.6; complete — `report/decam_survey.md`)

Goal: assess what a DECam adapter + 3-corridor pilot (Lalande 21185,
σ Dra, HD 219134) on the ZTF/PS1 pattern needs, for the 15 southern
corridors (δ < −30°) unreachable by PS1/ZTF.

## Service probes (2026-08-24, dev machine)

| Service | Status | Notes |
| --- | --- | --- |
| Astro Archive API (`astroarchive.noirlab.edu/api`, v7.1) | **works, anonymous** | `adv_search/find` POST JSON (snapshot-able), `limit=20000` returned 12,652 rows in one page |
| `api/retrieve/<md5>/?hdus=0,N` | **works** | per-CCD retrieval: dqmask HDU ~0.3 MB, image CCD ~5.4 MB in 4.5 s, wtmap CCD ~3 MB. Plain HTTP `Range` is ignored (200 + full file). This is the cutout mechanism. |
| Data Lab SIA (`datalab.noirlab.edu/sia/nsa`) | works, anonymous | 754 rows at σ Dra antipode incl. pre-DECam Mosaic II epochs (back to ≥2005 — possible baseline extension) |
| Data Lab `/svc/cutout` | works, anonymous | 30×30 pix FITS in 0.45 s, **but WCS comes back TAN — TPV distortion terms stripped**. Use `?hdus=` CCD fetches for anything astrometric. |
| Data Lab TAP (`/tap/sync`, `/tap/async`) | **down/unusable today** | sync: 0 bytes after 300 s; async: job stuck QUEUED. Blocks NSC DR2 (`nsc_dr2.object`/`meas`) screening. Retry later. |
| Data Lab query manager (`/query/query`) | rejects anonymous | needs a (free) Data Lab account — user action if NSC is wanted |
| VizieR fallbacks | work | ATLAS RefCat2 (`J/ApJ/867/105/refcat2`) has rows at δ −69.7 (all-sky, star calibration); Gaia DR3 via `I/355` for the flux-consistent static veto; DES DR2 = `II/371` (first cone at Lalande antipode returned empty — recheck footprint/query before relying on it) |

## Product facts (instcal, community pipeline)

- Every instcal exposure has matching `image` (med 321 MB), `dqmask`
  (med 6 MB fpacked), `wtmap` (med 154 MB) — 185/185/185 at σ Dra.
- dqmask: 62 HDUs (61 CCDs), per-CCD 2046×4094, full **TPV WCS (22 PV
  terms) in every CCD header**, `MAGZERO`/`MJD-OBS`/`EXPTIME`/`FILTER`/
  `PHOTFLAG` in the primary → the WISE `-msk` / PS1 skycell-mask trick
  transfers: dqmask alone = exact astrometry + usable pixels.
- Image CCD HDU: CompImageHDU float32, verified readable after
  `?hdus=0,25` fetch.
- File identity: archive `md5sum` is per *file*; an `?hdus=` subset has
  different bytes → record archive md5 + local sha256 of the subset
  (same pattern as a cutout).

## Coverage sweep (instcal images, center within 1.1° of antipode)

| corridor | antipode | N | span | cal. months | med exp |
| --- | --- | --- | --- | --- | --- |
| lalande-21185 | 345.83, −35.95 | 57 | 2015→2024 | 5 | 90 s |
| ross-248 | 175.48, −44.17 | 79 | 2013→2023 | 8 | 60 |
| 61-cyg | 136.75, −38.76 | 62 | 2015→2026 | 4 | 30 |
| struve-2398 | 100.68, −59.64 | 279 | 2016→2025 | 6 | 90 |
| groombridge-34 | 184.61, −44.02 | 89 | 2013→2026 | 8 | 90 |
| gj-1221 | 87.02, −70.88 | **12,652** | 2012→2026 | 8 | 20 |
| gj-338 | 318.58, −52.68 | 237 | 2013→2023 | 4 | 90 |
| gj-625 | 66.36, −54.30 | 318 | 2012→2026 | 6 | 90 |
| gj-687 | 84.10, −68.33 | **1,158** | 2012→2025 | 12 | 100 |
| gj-251 | 283.70, −33.27 | 52 | 2014→2024 | 6 | 90 |
| sigma-dra | 113.10, −69.65 | 185 | 2012→2024 | 7 | 90 |
| hd-219134 | 168.34, −57.17 | 81 | 2013→2026 | 8 | 60 |
| wolf-1069 | 126.52, −58.58 | 35 | 2016→2023 | 6 | 90 |
| gj-3512 | 310.33, −59.49 | 205 | 2013→2023 | 5 | 90 |
| gj-13157 | 147.92, −59.29 | 151 | 2014→2026 | 7 | 30 |

All 15 southern corridors covered; ~15.6k exposures total. Geometric
hit fraction on an active CCD for center-within-1.1° is roughly
3.18 deg²/3.8 deg² ≈ 0.8 before the precise pass. Filters are grizY
(+u, + occasional narrowband — cut in the freeze). gj-1221 and gj-687
antipodes sit near the LMC: data-rich but heavily crowded (SMASH/DELVE
depth); treat as a crowding class like SPHEREx's Galactic-plane
corridors, and expect to cap/select the 12.6k gj-1221 exposures.

Unlike PS1 (97:3), calendar-month spread (4–12 months) means most
corridors get **both parallax phases** → the phase veto works
single-archive.

## What the survey needs (estimate)

1. **Adapter** `sglsurvey/adapters/noirlab_decam.py` (~PS1-size,
   ~350–500 lines): discover = adv_search POST (snapshot raw JSON);
   nominal footprint = static 61-CCD focal-plane layout (derive corner
   offsets once from one reference dqmask) around `ra/dec_center`;
   exact footprint = dqmask CCD TPV WCS + bitmask; fetch = `?hdus=`
   image/wtmap CCD subsets. No interface change to `base.py` expected.
2. **Survey dir** `surveys/decam/`: hypotheses freeze (CTIO observer —
   new sglseti site entry; grizY bands; exposure/quality cuts incl.
   PHOTFLAG, EXPTIME ≥ ~30 s, obs_type=object; heterogeneous-PI data
   is the norm, not survey-uniform cadence), overlay builder for the
   15 southern corridors, profile.py (~200 lines) for the v2 engine.
   First post-v2-programme survey → v2 discipline from day one
   (dev/confirmatory split, null ensemble, image-level injections,
   frozen rule) — no v1-style exploratory pass.
3. **Flux scale**: per-frame star calibration (PS1 lesson 1) against
   ATLAS RefCat2 / Gaia DR3 via VizieR; header MAGZERO as cross-check.
4. **PSF**: no published PSF product → per-CCD empirical PSF from
   field stars (PS1 per-skycell pattern).
5. **Screening/statics**: Gaia DR3 via VizieR for the flux-consistent
   catalogued-static veto (all-sky, works at −70°); NSC DR2 meas/object
   when Data Lab TAP recovers (or free account for query manager).
6. **Positive control**: known asteroid through the full chain, as ZTF/PS1.
7. **Data budget**: pilot (3 corridors, 323 exposures) ≈ 5–10 GB.
   Full southern overlay ~15.6k exposures × ~10 MB fetched ≈ 150 GB
   transient with PS1-style per-batch purge; gj-1221 selection policy
   decides the real number.

## Open items before the pilot freeze

- Retry Data Lab TAP (NSC DR2) on another day; decide NSC vs
  Gaia/RefCat2-only screening.
- Verify DES DR2 (II/371) VizieR footprint at the δ > −65 corridors.
- Decide whether Mosaic II epochs (SIA `nsa`, pre-2012) are in scope
  (separate instrument = separate adapter or deferred).
- Confirm expected single-epoch depth from wtmap-based noise on a real
  CCD (90 s g ~ 23.5 AB expected → stacks deeper than ZTF's 22.5).
- Pilot corridor choice: plan names lalande-21185 (57 exp), sigma-dra
  (185), hd-219134 (81); if 57 proves too thin after the precise pass,
  gj-625 (318) or struve-2398 (279) are the natural substitutes.

## Addendum 2026-08-24 (later) — `astro-datalab` client + account question

Probed the `astro-datalab` PyPI package (v2.22.1, GitHub
astro-datalab/datalab) as the NSC DR2 route, while TAP remained dead a
second day (connection-level timeout, http 000).

**Findings (all anonymous, no account):**

- `dl.queryClient` sync queries **work anonymously** against
  `datalab.noirlab.edu/query` — including while TAP is down. Yesterday's
  "not logged in" error was a malformed raw-HTTP call, not an auth wall.
- `nsc_dr2.object` cone queries are fast: 277 / 802 / 3,997 objects in
  3′ cones at the lalande / sigma-dra / hd-219134 antipodes, 0.4–0.9 s.
- `nsc_dr2.meas` (per-exposure detections — the ZTF-objects-style
  screening layer): **direct q3c cone times out (>180 s)**; query by
  `objectid IN (…)` from an object-cone first — 0.4 s, returns
  mjd/filter/mag_auto/exposure-name rows. The exposure names map to
  Astro Archive files.
- `nsc_dr2.exposure` within 1.1° of sigma-dra: 112 rows, MJD max 57903
  (June 2017) → **NSC DR2 is time-partial** vs the 2012–2026 instcal
  archive (DR2 ingested through ~2019 at best here). Catalog screening
  covers only part of the epoch range; forced photometry covers all.
- Async queries and MyDB **require an account** ("provided security
  token is invalid" for anonymous).
- Dependency footprint: pulls astropy/numpy/pandas/pyvo/matplotlib/
  specutils etc. (41 pkgs in a clean venv) — all but specutils/httplib2/
  pycurl-requests already in this repo's env; cost is small.

**Assessment:**

- *Package: yes, worth using* — as the screening/calibration query
  route in the DECam adapter (pin it; snapshot the SQL + CSV response
  verbatim, same discipline as TAP/IBE snapshots). It is the only
  currently-working programmatic route to NSC DR2, and sidesteps the
  query-manager REST details. Discovery/fetch stays on the Astro
  Archive API (no dep). Alternative if the dep is unwanted: extract the
  ~30-line REST call (GET /query/query, X-DL-AuthToken header) — but
  the pinned package is less brittle against service changes.
- *Account: not needed for the pilot.* Sync anonymous covers object
  cones + meas-by-objectid at corridor scale. An account adds, in
  usefulness order: (1) **async queries** — needed only for full-cone
  meas sweeps or survey-wide control-corridor pulls that exceed the
  sync timeout; (2) **MyDB** — upload predicted track tables and
  crossmatch server-side, attractive for the LMC corridors
  (gj-1221: 12.6k exposures) and scale-up; (3) VOSpace + Data Lab
  Jupyter — compute near data, but against this project's local
  snapshot/reproducibility model. Decision: start anonymous; register
  a free account at scale-up if sync limits bite.

## Addendum 2026-08-24 (build day) — archive failure modes + validations

Found while fetching all 214 pilot usable exposures (each mode ~0.5-1%
of files; both recoveries now in the adapter):

- **`?hdus=` can 500** for specific files. Fallback: full-file fetch +
  local HDU extraction.
- **HDU order can differ between siblings** of one exposure (seen on a
  wtmap vs its dqmask). The fetch-time EXTNAME assertion catches it;
  recovery reads the per-file EXTNAME order from the `api/header/<md5>/`
  page (HTML, parseable; returns `BADFFILE` JSON for corrupt files).
- **Served bytes can mismatch the archive's own md5sum** (repaired
  files served under the original md5 record). The fallback accepts
  them only after parsing as FITS with matching EXPNUM
  (content-verified), with a printed note.
- **One pilot image file is unrecoverably corrupt server-side**
  (`de55324b…`, EXPNUM-file at hd219134): served bytes are not FITS at
  all (no fpack/gzip magic). The exposure is recorded missing —
  archive-side data loss, 1 of 214.
- Threading note: astropy's warning logger crashes
  (`dictionary changed size during iteration`) with concurrent imports
  in ThreadPoolExecutor workers — first fetch pass should warm imports
  or run `--workers 1` on retries.

**Validations (pre-freeze, `surveys/decam/results/`):**

- `observer_validation.json`: sglseti CTIO site vs independent astropy
  computation at z = 550 AU: site effect 7–16 mas, residual ≤ 1.7 mas
  (typical 0.3 mas) over 3 endpoints × 5 epochs — PASS (target 20 mas).
- `flux_scale_check.json`: 24 frames, 6 bands, NSC stars 16–19.5.
  Star-ZP MAD 0.014–0.10 mag per frame (crowded fields at the high
  end) → per-frame star calibration is well-fed. **MAGZERO convention
  = counts (not counts/s)**; header MAGZERO is per-frame unreliable
  (most within ±0.2 mag of star ZP; outliers +3.3 (i), −1.2 (r); u
  systematically +2.7) → star calibration mandatory, header ZP only a
  cross-check.
- `surveys/decam/targets/overlay_v1.{json,md}`: all 15 southern
  corridors graded — 10 ok / 5 crowded / 0 sparse after the draft
  quality cut (EXPTIME ≥ 30 s, grizY, obs_type object): quality-exposure
  counts 48 (wolf1069) – 2,033 (gj1221); non-pilot queue ordered
  ok-by-depth first.
