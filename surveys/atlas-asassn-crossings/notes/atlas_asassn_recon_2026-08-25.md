---
title: "ATLAS + ASAS-SN reachability recon for a Pipeline B crossings survey"
date: 2026-08-25
status: "recon complete incl. post-account ATLAS probes (token live, two full-history tasks pulled and measured); no blocking actions left; freeze not started"
---

# ATLAS + ASAS-SN recon (plan §5.8 item 2)

Goal: verify the access routes for the window coverage-fraction /
duty-cycle-beacon rung — nightly-cadence forced photometry over the
universal crossing windows (`crossings/universal_v1/`: 16,586 events,
88 endpoints), channels A (on-star) and B (antipode). Both routes were
marked _unprobed_ at adoption; this recon replaces general knowledge
with endpoint facts.

Headline: **both services are reachable and alive; the plan's central
assumption holds for ATLAS but only half-holds for ASAS-SN.** ATLAS
serves forced photometry at arbitrary positions (both channels) but
needs a free account (user action). ASAS-SN's scriptable service (Sky
Patrol v2) serves **pre-extracted light curves of catalogued sources
only** — no arbitrary-position photometry — so it cannot do channel B
point photometry; the OSU v1 service that can is behind reCaptcha
(human-in-the-loop only).

## Service probes (2026-08-25, dev machine)

| Service | Status | Notes |
| --- | --- | --- |
| ATLAS forced phot (`fallingstar-data.com/forcedphot`) | **up, API live** | REST API; anonymous `/queue/` GET → 403 (auth wall works); `/api-token-auth/` responds (needs username+password). Site banner: full outage recovered 2026-08-21 |
| ATLAS registration (`/forcedphot/register/`) | up | Free account: username, email, password only — **user action, then token via API** |
| ATLAS docs (`/apiguide/`, `/faq/`, `/resultdesc/`, `/api/schema/`) | up, complete | Full REST walkthrough incl. throttling protocol; OpenAPI 3 schema served — snapshot it at freeze |
| ASAS-SN Sky Patrol v2 (`asas-sn.ifa.hawaii.edu/skypatrol/`) | **up, HTTP only** | Port 443 refused — docs and API are plain HTTP. Client API host is `http://asassn-lb01.ifa.hawaii.edu:9006` |
| Sky Patrol v2 queries (client `skypatrol` 0.6.21) | **work, anonymous** | cone_search 0.7–3.3 s; light-curve download (threads=1) ~2 s/source; `adql_query` works with simple predicates (`BETWEEN`; `DISTANCE(POINT…)` → 400) |
| ASAS-SN Sky Patrol v1 (`asas-sn.osu.edu`) | up | Arbitrary-coordinate on-demand photometry form — **behind Google reCaptcha**, not scriptable; manual use only |

## ATLAS facts (home/apiguide/faq/resultdesc pages)

- **Coverage: now all-sky.** Four 0.5 m units — Haleakala, Mauna Loa,
  El Sauce (Chile), Sutherland (South Africa) — "whole sky with a
  cadence of 1 day between −50 and +50 and 2 days in the polar
  regions". The plan's "δ > −50" rationale is stale in our favour:
  southern corridors are covered too (southern units later-era;
  per-position epoch start must be measured post-account). o ~ 19.5,
  c/o filters, AB; 4 × 30 s exposures per night spaced over ~1 h
  (quad structure — an intra-hour sampling rung for the duty-cycle
  cell). Occasional H-alpha ("H") rows appear in output.
- **API flow** (all under token auth): POST `/queue/` with
  `ra, dec, mjd_min[, mjd_max, use_reduced, radec list, callback_url]`
  → 201 + task URL (429 throttle returns machine-parseable wait time)
  → poll task (or HTTPS `callback_url`, fired once, no retries) →
  `result_url` (short-lived) → space-separated text → DELETE task.
  `queuepositions.json` gives queue positions without listing tasks.
- **Photometry**: tphot PSF forced fit (pixels 1.86″, FWHM ~2 px), on
  **difference images** (default — signed µJy flux; the exact
  channel-B substrate, ZTF-crossings-style) or **reduced/target
  images** (`use_reduced` — for on-star channel A; beware tphot sky
  over-subtraction in crowded fields, it fits one PSF only).
- **Output columns** (fixed schema): MJD (exposure start, **not
  barycentric**), m/dm, uJy/duJy, F, err, chi/N, RA/Dec, x/y,
  maj/min/phi, apfit, mag5sig, PA_DEG, Sky, Obs (data-file id →
  provenance). Standard cleaning recipe (FAQ, after Rest et al.):
  `duJy<10000, err==0, 100<x,y<10460, 1.6<maj,min<5, −1<apfit<−0.1,
  mag5sig>17, sky>17`.
- **Template ("wallpaper") changes** cause step discontinuities in
  long difference light curves: v1→2 near MJD 58417 (2018-10-26),
  v2→3 near MJD 58882 (2020-02-03); exact template in WPDATE/WPDIR
  FITS headers since early 2021 (image requests exist as a separate
  request type). These are frozen-systematics-template material for
  the on-star channel (cf. `notes/learnings.md` §8).
- **Built-in MPC forced photometry** (`mpc (12267)` in the RA/Dec
  field, packed designations for unnumbered) — the asteroid
  positive-control chain comes for free, best for MPC uncertainty
  parameter 0 objects.

## ASAS-SN facts

**Sky Patrol v2** (Hart et al.; server deployment 0.6.20, 2025-09-11):

- Official client is PyPI **`skypatrol`** (0.6.21, author ASAS-SN,
  import name `pyasassn`); the PyPI package literally named `pyasassn`
  (0.6.4) is a stale twin with a broken `pyarrow==4.0.1` pin — do not
  use. Client methods: `cone_search`, `query_list`, `adql_query`,
  `random_sample`, `simbad_lookup`, `solar_system_object`; 20 input
  catalogs (`master_list` ~100M sources, `stellar_main`, external-ID
  catalogs incl. morx, aavsovsx, asteroids/comets).
- **Anonymous end-to-end works.** Light-curve columns: asas_sn_id, jd,
  flux, flux_err, mag, mag_err, limit, fwhm, image_id, camera,
  quality (G/B), phot_filter (V/g). V = era 2013–2018 (cameras
  ba–bh), g = 2017– (bi–bt). Probed cadence: g ~1.4 d/epoch mean,
  V ~2.5 d.
- Southern corridors covered: σ Dra antipode (δ −69.77) has 3,525
  epochs from 2014-04.
- **Catalogued sources only** — the structural limit. Nearest
  `master_list` source to the wolf-359 antipode is ~1′ off (8″
  pixels, ~15″ FWHM: not a substitute measurement). Channel B point
  photometry is **not servable from v2**. On-star channel A is:
  both probed targets (teegarden, σ Dra) have sources within 1′
  (saturation caveat below).
- **Currency gotcha**: max JD ≈ 2460841 (2025-06-15) in both probed
  fields — ~14 months stale today, consistent with the server's
  "reprocessing a number of lightcurves" banner (2025-05-22).
  Re-measure currency at freeze; windows after mid-2025 are not yet
  servable from v2.
- Client gotcha: `download=True` uses multiprocessing (spawn) — code
  must run from a script file under `if __name__ == '__main__':`
  (stdin/exec use spawn-loops forever); `threads=` is settable.
  Transport is plain HTTP on port 9006 — snapshot discipline
  unaffected (public data), but note it in provenance.
- ADQL dialect is restricted: simple column predicates work
  (`BETWEEN`), `DISTANCE(POINT(...), POINT(...))` → 400. Use
  `cone_search`/box predicates.

**Sky Patrol v1** (OSU, real-time arbitrary-position photometry):

- Form: RA/Dec (J2000), "days to go back" (default 20), proper motion
  + epoch (Gaia DR3 default), method = aperture / image-subtraction
  (**with or without reference flux** — the no-ref-flux mode is
  difference-flux, signed) / saturated-star ML (Winecki & Kochanek
  2024 — relevant because ASAS-SN saturates near V ≈ 10–11 and many
  of our 88 targets are brighter). V+g, APASS-calibrated zero points.
- POST target is a plain Rails endpoint (`/sky-patrol/computation/check`)
  but the form carries **Google reCaptcha** — scripting it is off the
  table (ToS + practicality). Treat v1 as a manual instrument: fine
  for a handful of designated follow-ups or spot checks, not for the
  survey statistic.

## What this means for the survey design

1. **ATLAS is the workhorse for both channels.** Arbitrary-position
   difference-flux forced photometry over the full archive at every
   star and antipode, all-sky, 1–2 d cadence, with a built-in
   asteroid positive control. 88 endpoints × 2 roles = 176 positions
   (plus control corridors); each is one full-history task
   (`mjd_min` ≈ 57000). Throughput/throttle limits are per-account
   and unmeasurable until registration.
2. **ASAS-SN role narrows to what v2 can serve**: (a) on-star channel
   A light curves for unsaturated targets, extending cadence coverage
   back to 2013–2017 (pre-ATLAS-south era in the south, V band);
   (b) window **coverage-fraction measurement** — the epoch list
   (image_id/camera/jd) of catalogued sources adjacent to any
   position is a direct record of when the cameras observed that
   field, usable for coverage accounting at antipodes even where no
   photometry at the exact point exists. Frame (b) as coverage, not
   constraint — a nearby-source epoch list is not a null measurement
   at the antipode.
3. The reCaptcha wall means any v1 arbitrary-position measurements
   (e.g. adjudicating an ATLAS exceedance at an antipode, saturated
   on-star curves via the ML method) are **manual, enumerated,
   follow-up actions** — pre-registered as such, like the APF/VLASS
   hand-offs in §5.5.
4. Cross-check at freeze: whether the ATLAS 4-exposure/night quad
   (~1 h span) supports a pulse-cell statement between the TESS
   cadence rung and the nightly rung.

## Addendum 2026-08-25 (later) — post-account ATLAS probes

Account registered (token in `.env` as `ATLAS_API_TOKEN`, gitignored).
Two full-history difference-flux tasks (`mjd_min=55000`) run to
completion; raw light curves in the session scratchpad (regenerable —
identical queries at freeze go under `runs/` snapshot discipline).

| probe position | rows | MJD span | distinct good nights | wall time |
| --- | --- | --- | --- | --- |
| wolf-359 antipode (δ −7.0) | 4,190 | 57227 → 61275 (2015-07-24 → **2026-08-22**) | 940 | 26.5 min |
| σ Dra antipode (δ −69.8) | 3,024 | 59578 → 61275 (**2022-01-29** → 2026-08) | 602 | 68.5 min (queued behind the first) |

Measured facts:

- **Archive is current to within ~3 days of real time** (last epochs
  MJD 61275–61276) — no ingestion lag; recent crossing windows are
  servable, unlike Sky Patrol v2's 2025-06 ceiling.
- **Southern history starts 2022.** Obs-file prefixes date the units:
  01a Haleakala from MJD 57928 (2017-06), 02a Mauna Loa from 57227
  (2015-07, archive start), 03a El Sauce from 59725 (2022-05), 04a
  Sutherland from 59683 (2022-04). At δ −69.8 only 03a/04a contribute
  → **δ < −50 corridors have no ATLAS epochs before 2022** — frame
  southern window-coverage claims accordingly (2022+ windows only).
- **Quad structure confirmed**: median 4 rows/night at both
  positions — the intra-hour pulse rung is real. Night density: ~940
  nights/11.1 yr (north, incl. seasonal sun gaps), ~600/4.6 yr
  (south).
- **FAQ quality recipe keeps 93% / 89%** of rows; median per-exposure
  mag5sig 19.19 / 19.07 — o≈19.5 nominal holds; nightly quad stacks
  reach ~19.9–20.
- **Throughput**: two rapid POSTs accepted without throttling, but
  execution is **serial per account** (second task sat queued 26 min
  until the first finished). At ~25–70 min/full-history position,
  176 positions ≈ 3–6 days of wall clock — plan the production pull
  as a days-long background queue drain (callback_url or slow poll),
  or probe `radec`-list batching at freeze. Result URLs are
  short-lived; fetch-then-delete (done here) is the pattern.
- Result text carries `###` comment prefixes (strip before parse);
  `Obs` id (e.g. `01a58384o0596o`) is the per-exposure provenance key.

## Open items before the freeze

- Probe `radec` list batch submission (size limit, whether it
  parallelizes within one account) before choosing the production
  pull pattern; confirm `use_reduced=true` returns the same Obs ids
  as difference mode for the on-star channel.
- Sky Patrol v2: re-measure DB currency (is the 2025-06 ceiling
  reprocessing lag or ingestion lag?); decide the saturation cut for
  channel A from the actual target magnitudes; verify camera-era
  coverage per corridor (V vs g).
- Decide the coverage-fraction estimator for antipodes from v2 epoch
  lists (nearest-source vs field-level union) before any window is
  scored.
- Naming/carpentry: `surveys/atlas-asassn-crossings/` scaffold,
  observer = Earth-center vs site (four ATLAS sites; ASAS-SN five
  stations) — the universal list is Earth-center, and site parallax
  at 550+ AU is far below the beam scales, but record the decision.
- Snapshot the ATLAS OpenAPI schema and both docs pages at freeze;
  pin `skypatrol` 0.6.21.
