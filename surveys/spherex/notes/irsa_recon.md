# IRSA SPHEREx reconnaissance (2026-08-20)

What this sub-project queries, and the facts verified from the dev
machine before the adapter was written. Reference: SPHEREx Explanatory
Supplement (IRSA), `https://irsa.ipac.caltech.edu/data/SPHEREx/docs/`.

## Discovery: IRSA TAP

Tables: `spherex.obscore` (1,324,071 rows), CAOM views
`spherex.observation`, `spherex.plane`, `spherex.artifact`,
`spherex.provenance_input`. Collections in ObsCore:

| obs_collection | dataproduct_type | calib_level | rows |
|---|---|---|---|
| `spherex_qr2` | image | 2 | 1,079,504 |
| `spherex_qr2_deep` | image | 2 | 265,805 |
| `spherex_qr2_cal` | measurements | 1 / 3 | 74 |

QR1 is not served separately (superseded by the QR2 reprocessing). All
rows `data_rights='public'`. Epoch range 2025-05-24 → 2026-08-11
(MJD 60789–61233) at the pilot corridors.

- `obs_id` = `YYYYWww_nX_NNNN_k` (week, pointing, exposure); one row per
  (obs_id, detector). `energy_bandpassname` = `SPHEREx-D1..D6`.
- `s_region` = 4-vertex ICRS polygon (3.5° field); `s_ra/s_dec` centre;
  `t_min/t_max` bracket the ~113 s exposure; `s_pixel_scale` 6.15";
  `s_resolution` 3.08 (listed; actual PSF FWHM ≈ 5.3").
- `access_url` is a DataLink document (science file, IBE cutout
  endpoint, calibration files, S3 location).
- ADQL: `CONTAINS(POINT, CIRCLE)` and `CONTAINS(POINT, s_region)`
  work; `INTERSECTS(CIRCLE, s_region)` is rejected. Discovery therefore
  selects field centres within (cone radius + 2.6° half-diagonal) and
  applies the polygon footprint client-side. Verified: 349 polygon hits
  at the Lalande antipode = TAP's own `CONTAINS(POINT, s_region)` count.
- Product URI and **published MD5** come from
  `spherex.artifact` (`uri`, `contentchecksum`, `contentlength`,
  `producttype='science'`) joined through `spherex.plane` on
  `obs_publisher_did`. The join runs in ~20 s per corridor without
  `TOP`; with `TOP n` the planner degrades to ~90 s.
- Level-2 pipeline version string (`l2b-v20-YYYY-DDD`) varies within a
  week — URIs cannot be constructed from `obs_id`; always use the join.

## Products

Full Level-2 file: 71,634,240 bytes, MEF:

| HDU | content | notes |
|---|---|---|
| 0 PRIMARY | `VERSION` only | |
| 1 IMAGE | 2040×2040 float32, MJy/sr | `RA---TAN-SIP` (A/B/AP/BP order 3), `MJD-AVG`, `PSF_FWHM`, `OMEGA_MEDIAN` (arcsec²), `X_SC…` geocentric state; alternate WCS `A` (raw pixels) and `W` (`WAVE-TAB` → WCS-WAVE) |
| 2 FLAGS | int32 bitmask | bits in `MP_*` header cards (see hypotheses §9) |
| 3 VARIANCE | MJy²/sr² | |
| 4 ZODI | MJy/sr | zodiacal model |
| 5 PSF | 121 × 101 × 101 | 11×11 detector zones (`XCTR_k/YCTR_k`), 10× oversampled (0.615"/px) |
| 6 WCS-WAVE | 1 row: X[9], Y[9], VALUES[9,9,2] | wavelength and bandwidth (µm) on a 9×9 grid of raw pixel coordinates |

- **S3 mirror:** `s3://nasa-irsa-spherex/` (us-east-1), anonymous HTTPS
  at `https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/<key>` with
  byte-range support. Measured 20 MB/s for 8 MB ranges, ~0.15 s
  latency; IBE cutouts 4 MB/s. The adapter assembles slim cutouts by
  range reads (9 requests, ~8.7 MB per 80-px cutout, ~3 s).
- **IBE cutouts:** `<ibe url>?center=ra,dec&size=Npix` return the full
  MEF (all extensions cut identically, PSF cube intact: 5.15 MB for
  any size; extension selection not supported; never gzipped).
- No per-exposure source catalogs in QR2; no difference images.
- Cutout WCS: IBE shifts `CRPIX1/2`, `CRPIX1A/2A` and `CRPIX1W/2W`;
  the adapter does the same.
- Empirical check on one D1 cutout: usable fraction 97.8 %; robust
  residual scatter 1.18× the VARIANCE model; matched-filter 5σ depth
  ≈ 18.6 AB at 0.75 µm; PSF central-pixel fraction 0.66.
