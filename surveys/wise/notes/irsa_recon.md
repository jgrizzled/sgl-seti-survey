---
title: "IRSA reconnaissance — WISE/NEOWISE data products and services"
date: 2026-08-18
status: "verified against live IRSA services on 2026-08-18 unless marked TO VERIFY"
---

# IRSA reconnaissance for the WISE shakedown

Everything below was confirmed by live queries against
`https://irsa.ipac.caltech.edu` on 2026-08-18 from this machine (no
authentication required), except items marked **TO VERIFY**.

## Services

- **TAP** — `https://irsa.ipac.caltech.edu/TAP/sync` (and `/async`).
  ADQL over all tables below; CSV/VOTable output. This is the primary
  discovery interface: frame inventory tables carry footprint corners and
  WCS, so coarse discovery can run server-side.
- **IBE** — `https://irsa.ipac.caltech.edu/ibe/data/wise/<dataset>/p1bm_frm/`
  for image retrieval. Verified dataset directories: `allsky/`, `allwise/`,
  `cryo_3band/`, `merge/`, `neowiser/`, `postcryo/`, `prelim/`,
  `prelim_postcryo/`. **Per-frame path template (verified 2026-08-18):**
  `<scangrp>/<scan_id>/<frame_num %03d>/` containing, per band,
  `<scan_id><frame_num %03d>-w<band>-int-1b.fits`, `-w<band>-msk-1b.fits.gz`,
  `-w<band>-unc-1b.fits.gz`, artifact tables `-art-w<band>-{D,H,P}.tbl`,
  frame flags `-fflag-1b.tbl`, and a `.md5` sidecar for every product
  (checksums for free — use them in Observation records). **Cutouts
  verified:** append `?center=<ra>,<dec>&size=<N>pix` to the FITS URL
  (a 120-pix W1 cutout returned HTTP 200, 61 KB, vs ~4 MB full frame).
- SIA also exists at IRSA but TAP + IBE appears sufficient; decide in the
  adapter whether SIA adds anything (plan §3.3 allows either).

## Frame inventory tables (L1b images) — discovery targets

Verified row counts, 2026-08-18:

| Mission phase | Image inventory table | Rows | Bands |
| --- | --- | --- | --- |
| 4-band cryo (2010) | `allsky_4band_p1bm_frm` | 5,964,417 | W1–W4 |
| 3-band cryo | `allsky_3band_p1bm_frm` | 1,165,452 | W1–W3 |
| Post-cryo (2010–11) | `allsky_2band_p1bm_frm` | 1,802,540 | W1–W2 |
| NEOWISE-R (2013–2024) | `neowiser_p1bm_frm` | 53,776,224 | W1–W2 |
| NEOWISE-R merged | `wise.neowiser_merge_p1bm_frm` | 62,708,633 | W1–W2 |

`prelim*` tables are marked superseded at IRSA — exclude them.
**Verified 2026-08-18:** `wise.neowiser_merge_p1bm_frm` is the union of
all four mission phases — it contains band-4 rows at MJD 55203 (Jan 2010,
4-band cryo), and the row counts add exactly
(5,964,417 + 1,165,452 + 1,802,540 + 53,776,224 = 62,708,633). This is
the "WISE merged L1b" of the plan and the natural single discovery
target; the corresponding IBE dataset is `merge/`.

Each `p1bm_frm` table has one row per band per frame. Key columns
(verified on `neowiser_p1bm_frm`):

- Identity: `scan_id`, `scangrp`, `frame_num`, `band`, `wrelease`
- Timing: `mjd_obs` (mid-point, UTC), `date_obs`, `exptime` (7.7 s W1/W2
  nominal — so plan §3.3's midpoint-is-sufficient note applies, but keep
  the start/stop interval model anyway)
- Footprint: corner positions `ra1/dec1 … ra4/dec4`, plus full WCS
  (`crval*`, `crpix*`, `ctype*` with distortion, `cd*_*`, `pxscal*`, `pa`)
  and astrometric uncertainties (`crder1/2`, `csdradec`)
- Quality/context: `moon_sep`, `magzp`, `magzpunc`, and more (frame-level
  quality columns like `qual_frame` live in the `p1bs_frm` metadata
  tables — join on `scan_id`,`frame_num`)

Discovery-rate sanity check: ~63M NEOWISE frame-band rows over the whole
sky means a small corridor intersects a very manageable number; a
positional TAP query with corner-polygon or `CONTAINS(POINT(...))`
constraints per epoch chunk should return thousands, not millions, of rows
per pilot target. **Verified 2026-08-18:** IRSA TAP supports
`CONTAINS(POINT('ICRS',crval1,crval2), CIRCLE('ICRS',ra,dec,r))` on the
62.7M-row merge table; a 0.55° cone count takes roughly 1–2 minutes via
`/TAP/sync`. Fine for pilot-scale discovery; use `/TAP/async` for
production runs and consider corner-column bounding-box prefilters if it
becomes a bottleneck.

## Source-level and catalog products (screening layer, plan §4.5)

Verified tables and row counts:

- `allsky_4band_p1bs_psd` (9.5e9), `allsky_3band_p1bs_psd` (3.7e9),
  `allsky_2band_p1bs_psd` (7.3e9), `neowiser_p1bs_psd` (2.0e11) — single
  exposure source tables
- `catwise_2020` (1.9e9) and `catwise_2020_reject` (3.4e8) — motion-aware
  catalog + reject table (plan: reject tables are first-class)
- `unwise_2019` (2.2e9) catalog; `unwise.unwise_neo3_images` coadd
  metadata exists at IRSA but time-resolved unWISE epoch coadds may need
  outside hosting — **TO VERIFY** where time-resolved unWISE lives
- `neowiser_p1ba_mch` (2.3e8) — known-solar-system-object association
  list (veto/positive-control input; IRSA flags caveats in its README)
- **TO VERIFY:** single-exposure *reject* table names for each phase, and
  VarWISE availability via TAP

## Practical notes for the adapter

- Snapshot policy: store the exact ADQL, HTTP request, response bytes,
  and `wrelease` values with every discovery run (plan §4.3).
- Row-level product identity: (`scan_id`, `frame_num`, `band`) + mission
  phase table + `wrelease` is the natural immutable Observation key.
- Full-frame L1b FITS are ~4–8 MB each; a pilot corridor with a few
  thousand frames is a few tens of GB — feasible locally, but IBE cutouts
  along the locus would cut this by orders of magnitude. Decide retention
  policy per plan §6 before bulk pulls.
