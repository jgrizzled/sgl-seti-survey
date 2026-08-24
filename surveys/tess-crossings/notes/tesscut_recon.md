# TESScut recon (2026-08-24, pre-freeze)

One real cutout fetched and inspected before freezing: the wolf-359
grazing-crossing antipode (`evt-e316050714db`, b = 0.55 R☉,
t_ca MJD 59462.23), sector 42, 15×15 px
(`tess-s0042-1-1_344.080150_-6.986940_15x15_astrocut.fits`, 16 MB
zipped).

Facts the freeze relies on:

- **Endpoint**: `mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=&dec=&x=&y=&units=px&sector=`
  → zip of one FITS per sector/camera/ccd. Sector lookup:
  `/api/v0.1/sector?ra=&dec=` (used by the coverage gate; one bogus id
  "1751" observed once and ignorable). Cutouts are cheap: a 15×15
  full-sector cube is ~16 MB.
- **Structure**: PRIMARY + `PIXELS` bintable + `APERTURE` image ext.
  PIXELS columns: TIME, TIMECORR, CADENCENO, RAW_CNTS, **FLUX (e-/s)**,
  FLUX_ERR, **FLUX_BKG, FLUX_BKG_ERR** (SPOC background is provided —
  no self-built background model needed), QUALITY, POS_CORR1/2,
  FFI_FILE. `APERTURE` carries the TAN WCS; the requested position
  landed at pix (7.53, 7.87) of 15×15 — centered, sub-pixel accurate.
- **Time**: TIME is BTJD (TIMESYS=TDB, BJDREF 2457000), barycentre-
  corrected (TIMECORR column carries the correction). MJD(TDB) ≈
  BTJD + 56999.5. s42 span 59447.20–59472.66 → the grazing t_ca sits
  mid-sector with ~13 d of in-sector off-window baseline either side.
  UTC-vs-TDB and barycentre-vs-spacecraft light time (≤ ~8 min
  combined) are ≪ the shortest window (0.3 d) — carried as a budget
  term.
- **Cadence**: exactly 600.0 s in s42 (10 min, sectors 27–55); 3,534
  cadences per sector; s56+ are 200 s, s1–26 are 30 min.
- **Quality**: QUALITY == 0 for 57.6 % of s42 cadences; the dominant
  nonzero flags are 2048 (straylight) and 18432 (straylight
  combination) in contiguous blocks, plus a handful of bit-5 desat
  cadences. QUALITY == 0 is the natural primary mask; strict adds a
  straylight-adjacency margin.
- **Fluxes**: e-/s; background level ~946 e-/s/px with ~20 e-/s frame
  scatter at this ecliptic field. TESS magnitude scale
  T ≈ 20.44 − 2.5 log₁₀(flux e-/s) as the header-level convention;
  per-cutout star calibration against TIC magnitudes through the
  identical filter (the PS1 star-ZP lesson) is the planned primary
  flux scale with 20.44 as fallback.
- **Scale**: 21″ pixels; track motion at z = 550 over an 11.6 d window
  ≈ 3.6 px; a 25×25 px cutout holds track + arcminute ring controls +
  margin. PSF FWHM ≈ 1.3–2 px (undersampled — the SPHEREx sub-pixel
  lesson applies; SPOC PRF models exist per camera/CCD for
  injections).
