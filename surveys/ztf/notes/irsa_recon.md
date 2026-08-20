# IRSA ZTF reconnaissance (2026-08-20)

What this sub-project queries, and the facts verified from the dev
machine before the adapter was written. Authoritative reference: ZTF
Science Data System (ZSDS) Explanatory Supplement,
<https://irsa.ipac.caltech.edu/data/ZTF/docs/ztf_pipelines_deliverables.pdf>
(§7 file products, §10.3 mask bits, §10.4 INFOBITS, §13 cautionary
notes). Section numbers below refer to that document.

## Discovery: IBE metadata search

`https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?POS=ra,dec&SIZE=deg&ct=csv&mcen&WHERE=...`

One row per CCD-quadrant exposure ("readout channel", `rcid` 0–63).
Columns used: `field ccdid qid rcid fid filtercode pid nid expid
imgtype obsjd exptime filefracday seeing airmass moonillf maglimit
infobits crpix1 crpix2 crval1 crval2 cd11 cd12 cd21 cd22 ra1..dec4
ipac_pub_date ipac_gid`. The IBE `POS/SIZE` search returns quadrants
whose footprint overlaps the box (verified: a 0.01° box at the Ross 128
corridor returns 335 rows from two field/quadrant combinations).

- **Access groups:** `ipac_gid=1` public (MSIP), 2 partnership,
  3 Caltech. Only gid 1 is fetchable anonymously; the adapter filters
  `ipac_gid=1` in the WHERE clause. Ross 128 point: 164 public frames
  2018-06-01 → 2026-06-01 (g 71 / r 92 / i 1); all 30 s.
- **Quality:** metadata `infobits` bit 25 (value 33554432) is the
  archival "bad-quality" flag (§13.4 note 2); header INFOBITS does not
  carry it. Header `STATUS=1` always in archive (fatal bits already
  removed, §10.4). Adapter keeps bit-25 rows as Observations but marks
  them `quality_flags.bad_quality=True`; the pipeline excludes them.
- **Timing:** `obsjd` = exposure **start** (§13.1 note 2);
  `t_mid = obsjd + exptime/2`. Filename `filefracday` is ~1 s earlier
  than OBSJD (§13.1 note 5) — never derive time from the filename.

## Products: IBE data tree

`https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/YYYY/MMDD/dddddd/ztf_YYYYMMDDdddddd_FFFFFF_<filtercode>_cCC_o_qQ_<suffix>`

| suffix | content | cutout-able |
|---|---|---|
| `sciimg.fits` | calibrated quadrant image, TPV WCS, MAGZP/SEEING/MAGLIM in header | yes |
| `mskimg.fits` | 16-bit mask, same WCS | yes |
| `scimrefdiffimg.fits.fz` | PSF-matched science − reference difference (fpacked) | yes (returns 2-HDU file) |
| `diffimgpsf.fits` | PSF of the difference image | n/a |
| `psfcat.fits` | PSF-fit source catalog (cols defined in §10.6) | n/a |
| `sexcat.fits` | SExtractor aperture catalog | n/a |

Cutouts: append `?center=ra,dec&size=Npix`; the returned header keeps
the full TPV WCS with CRPIX shifted (verified: CRPIX1=-407.5 for a
50-pix cutout). WCS is `RA---TPV` — astropy handles it. Pixel scale
1.012"/pix; quadrant 3072×3080 pix ≈ 52'×52'.

- **Missing products:** ~1 in 8 public metadata rows have no served
  product directory (404 on the `dddddd/` directory itself). Sample of
  8: 7 present, 1 missing (2023-10-28 g-band). The adapter records
  `products.available=false` rather than dropping the Observation.
- **No checksums** are published alongside products (unlike WISE
  `.md5`); the adapter records sha256 of the bytes it received.
- **Difference images** only exist where a reference image predated
  the epoch (§13.3 note 3); early-2018 epochs may lack them. The
  forced-photometry stage must fall back to `sciimg` with its own
  background model when `scimrefdiffimg` is absent.

## Mask bits (§10.3)

| bit | meaning | fatal for us |
|---|---|---|
| 0 | aircraft/satellite track | yes |
| 1 | contains SExtractor detection | **no** (a source at the locus is the signal) |
| 2 / 3 | low / high responsivity | yes |
| 4 | noisy | yes |
| 5 | ghost from bright source | yes |
| 6 | ghost from charge spillage (OBSMJD ≤ 58779) | yes |
| 7 | pixel spike / rad hit | yes |
| 8 | saturated | yes |
| 9 | dead | yes |
| 10 | NaN | yes |
| 11 | contains PSF-extracted source | **no** |
| 12 | halo from bright source | yes |
| 13–15 | reserved | no |

Fatal template = 6141 (the ZSDS "uncontaminated" template).

## Catalogs for screening

- IRSA TAP `ztf_objects_dr24` (DR24 lightcurve objects; per-object
  `ra dec oid nobs ngoodobs medianmag ...`); lightcurve epochs via the
  IRSA lightcurve service. Objects are per filter/field/quadrant.
- Per-frame `psfcat.fits` (§10.6: `ra dec flux sigflux mag snr chi
  sharp flags`) — the deepest single-epoch list; `flags` is the OR of
  mask bits in a 5×5 box.
- ZTF alerts are not archived at IRSA; difference-image forced
  photometry here replaces them.

## Observer

Palomar P48: lon −116.8650°, lat +33.3563°, h 1712 m. Topocentric vs
Earth-center locus shift at 550 AU ≤ 0.016" (sglseti `Observer.from_geodetic`,
verified against 6371 km / 550 AU). Negligible against padding; the
site observer is still declared so the run config is honest about the
observer state.
