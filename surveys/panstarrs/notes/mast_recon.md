---
title: "Pan-STARRS1 DR2 at MAST — reconnaissance for the PS1 adapter"
date: 2026-08-20
---

# PS1 DR2 access facts (verified 2026-08-20 from the dev machine)

All services below are at STScI/MAST — **no IRSA dependency**.

## Single-epoch images ("warps")

- Listing: `https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=&dec=&filters=grizy&type=warp`
  returns one row per warp covering the *skycell that contains the query
  point*: `projcell subcell ra dec filter mjd type filename shortname badflag`.
  `type=warp.mask` and `type=warp.wt` list the mask and weight planes
  (same filename stem + `.mask.fits` / `.wt.fits`). ~1 s per call.
  Counts at the three pilot antipodes: ross128 117, epsind 131, proxima 90
  warps (≈20 per filter; MJD 55000–57000, i.e. 2009–2014).
- Files: `https://ps1images.stsci.edu` + `filename`. Full warp image
  ≈30 MB, full **mask ≈3.3 MB (fpack uint16, HDU 1 `COMPRESSED_IMAGE`,
  6240 x 6243)** — the WISE `-msk` trick works again: masks alone give
  exact astrometry and valid pixels.
- Cutouts: `https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra=&dec=&size=<pix>&format=fits&red=<filename>`
  — 1200² in ~1.5 s, 2400² (23 MB) in ~3 s (≈7.5 MB/s per stream).
  Output is a single primary HDU, BITPIX -32, with the **full skycell
  WCS and CRPIX shifted for the cutout** ("Reference pixel shifted for
  cutout"). fitscut adds `RADESYS=FK5`, `TIMESYS=TAI` (cosmetic).
- **Gotcha: fitscut renders mask value 0 as NaN** (the full mask is 78 %
  zeros; the cutout mask has no zeros and 87 % NaN). Use full masks, or
  treat NaN-in-cutout-mask as 0.
- All warps of one skycell share one pixel grid / WCS (TAN, CDELT
  0.25"/pix, CRVAL = projection-cell centre, PC = [[-1,0],[0,1]]). A
  warp covers only the part of the skycell the exposure touched; GPC1
  chip/cell gaps appear as NaN (≈29 % of pixels in the probe cutout).
- Header facts: `MJD-OBS` (exposure start, TAI per fitscut), `EXPTIME`
  (60 s in the 2009 probe; nominal 3π: g 43, r 40, i 45, z 30, y 30 s),
  `HIERARCH FPA.ZP` (24.64), `HIERARCH CHIP.SEEING` (FWHM pixels),
  `AIRMASS`, `HIERARCH CELL.GAIN` (0.998 e/DN), `HIERARCH CELL.SATURATION`,
  observatory `FPA.LONGITUDE 10.417 h W / FPA.LATITUDE 20.7071 / FPA.ELEVATION 3048 m`
  (Haleakalā), per-chip MD5 of the full skycell image
  (`MD5_SkyChip_SkyCell_0`).
- Mask bits (header MSKNAMnn/MSKVALnn): DETECTOR 1, FLAT 2, DARK 4,
  BLANK 8, CTE 16, SAT 32, LOW 64, SUSPECT 128, CR 256, SPIKE 512,
  GHOST 1024, STREAK 2048, STARCORE 4096, CONV.BAD 8192, CONV.POOR
  16384, BURNTOOL 128 (alias), MASK.VALUE 8575 (IPP's own "bad" template
  = DETECTOR|FLAT|DARK|BLANK|CTE|SAT|LOW|CR|CONV.BAD), MARK.VALUE 32768.
  Most common values in the probe: 0, 8192 (CONV.BAD, 9.5 % of the
  skycell, image finite), 8232, 8230, 1024 (GHOST), 8193.
- Weight plane: positive floats (median 160 in the probe); verify
  variance convention empirically against the image background scatter
  before trusting it (done in `sample_tensor.py`).

## Catalogs (MAST catalogs API, CSV)

- Per-epoch detections: `https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/detection.csv?ra=&dec=&radius=&pagesize=&page=&columns=[...]`
  — 116k detections within 0.2° of the ross128 antipode in 7 s.
  Columns used: objID, detectID, filterID (1–5 = grizy), obsTime (MJD),
  ra, dec, psfFlux/psfFluxErr (Jy), psfQfPerfect, infoFlag/2/3, zp,
  expTime, imageID.
- Mean objects: `.../dr2/mean.csv` (raMean, decMean, nDetections,
  {g..y}MeanPSFMag/Err, qualityFlag) — 44.7k objects within 0.2° in 4.5 s.
  Used for the per-warp empirical zero point (star recovery) and the
  recurrence context.

## Not used

- MAST CAOM TAP lists PS1 stacks, not warps. Rubin RSP needs a token
  (401). No checksum service for warps beyond the in-header MD5 of the
  full skycell image (not applicable to cutouts) — recorded checksums are
  sha256 of received bytes, as for ZTF.

## Measured during the pilot

- `detection.obsTime` is **~50–60 s later than the warp `MJD-OBS`**
  (consistent with exposure end / start + EXPTIME + overhead), so
  warp-to-detection matching uses a ±120 s window. About 15 % of warp
  epochs have no detections in the DR2 table at all (ross128 skycells:
  126 detection epochs vs 155 warp epochs).
