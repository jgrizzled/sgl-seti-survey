# PTF/iPTF corridor survey — hypotheses v1.0 (FROZEN 2026-09-10)

Status: **frozen** (user approved the §8 decisions 2026-09-10; the
hash of this document is pinned in `configs/v2_freeze.json`). Version
history: v0.1 draft (2026-09-10) accompanied the geometric stages;
v1.0 freezes the §8 decisions and the positive-control selection. The v2 discipline applies
from the first analysis (DECam precedent, `surveys/decam/hypotheses.md`):
stratified dev/confirmatory split with a declared seed and the
pilot corridors forced into the development set, null-ensemble
thresholds, FWER α = 0.05 with BH q-values informational, image-level
injections. There is no v1-style exploratory pass.

Why this survey (plan §4.13): Pipeline A has PS1 (2009–2014, one
parallax phase per corridor) and ZTF (2018–) on the northern sky; PTF
is the only other 1″-class multi-epoch optical archive of the PS1 era,
from a different site, cadence and filter set, and its campaign
cadence revisits fields across the year — it can supply the parallax
phases PS1 lacks before 2018 and an independent re-observation of
the PS1-era corridors. Depth is PS1-class or shallower (R ~ 21 per
epoch), so the cell it adds is phase coverage and era-independent
persistence, not depth.

## 1. Physics cell (inherited)

Same baseline local-artifact hypothesis as WISE v2.1 / ZTF v2 /
PS1 v2 / DECam v1.0: compact station-kept relay on the anti-star
focal ray, log-uniform relay-distance prior 550–10,000 AU, Rx and Tx
as separate hypotheses, duty cycle ≥ 0.5 instantiated by the four
temporal injection families, |µ_resid| ≤ 1″/yr (L∞ box), unresolved
source. Bands g (fid 1, SDSS-like) and R (fid 2, Mould R) carry the
reflected-sunlight and self-luminous-optical interpretations with
flat-Fν injections; Hα frames (fid ≥ 3) are excluded by the adapter's
discovery gate and never enter the cell.

## 2. Geometry

- Model/ephemeris: Tusay2022Eq57V1 + astropy built-in (as all v2
  surveys); 99 % confidence, +10″ padding, 5″ tolerance at discovery;
  1″ precise tolerance, 10″ seed step.
- Observer: Palomar P48 (`palomar-p48`, lon −116.8650°, lat 33.3563°,
  h 1712 m) via `GeometryContext.ztf_default()` — PTF and ZTF are the
  same telescope and site, so the ZTF pilot's validated site model
  (16 mas bar) carries over unchanged. Mid-exposure evaluation: 60 s
  frames move the locus < 5 mas.
- Era: MJD 54850–57200 discovery window over the public level-1
  archive (2009-03-01 → 2015-01-28; late iPTF never released).
- Covariance MC (N = 2000, seed 20260910, epochs 2009.5 / 2012.0 /
  2014.5, z = 550 and 10,000 AU); cross-track tensor dimension where
  σ_xt(99 %) > 0.5 × nominal FWHM (g 2.3″, R 2.1″) — the v2 common
  rule.

## 3. Archive and products

IRSA IBE `ptf/images/level1`: one row per CCD exposure (12 CCDs,
CCD 3 dead; 2048 × 4096 px at 1.01″/px), explicit-column discovery
(the default column set omits fid, WCS, checksums), `scie` image +
`dmask` ancillary selected by `anciltype` (slot order varies), MD5
published for every product (prefix-verified — some `achecksumN` are
truncated). Cutouts (`center`/`size=Npix`) preserve the full
TAN-SIP + PV WCS with shifted CRPIX (verified at the crossings recon
and again here on 600-px stamps). **No difference images exist**:
the search image is the sky-subtracted science image (PS1
warp-direct pattern), so static sources are in the photometry and
the catalogued-static flux-consistency test applies at full weight
(no per-epoch difference-fraction scaling).

Products retained per usable exposure: one 600-px `scie` cutout and
its `dmask` cutout on the same pixel grid (centred on the union of the
Rx and Tx loci at mid-exposure; 600 px = 606″ encloses the ≤ 375″
z = 550 AU locus with margin, and exceeds the 384-px calibrator-density
floor measured at the crossings survey).

Data-quality convention (draft): fatal dmask template 65533 (Laher
et al. 2014 Table 15 — every bit except 2¹ object-detected) as
primary; strict = primary + seeing ≤ 2.5″ + |moonillf| ≤ 0.8; loose =
no listing-level cut (the fatal template is part of the substrate,
not a mask). `photcalflag` and `infobits` are recorded, not gated
(`photcalflag` ≈ 0 for nearly every frame; `infobits` observed 0
throughout the crossings sample).

## 4. Flux scale

Per-frame star calibration through the identical matched filter (PS1
lesson 1; DECam §4): calibrators are PS1 DR2 `mean` objects from the
corridor's PS1 corridor-screen cone (`runs/panstarrs/screen_v1`, 0.2°,
the same 62 corridors), predicted into the PTF bands by the frozen
transforms of the crossings survey — g = gMeanPSFMag; R = r −
0.153 (r − i) − 0.117 (Jordi et al. 2006 Cousins R) **+ 0.21 mag**
(Cousins R Vega → AB, Blanton & Roweis 2007) so that both bands are
reported in AB — with the dwarf-locus colour restriction (r − i ∈
[0, 0.8] for R, g − r ∈ [0.2, 1.2] for g), nDetections ≥ 5, errors
≤ 0.1, predicted mag 15.5–19.5 (R) / 16.0–20.0 (g). zp_star =
median(m_pred + 2.5 log10 F_mf) over calibrators with good_frac ≥ 0.7
and S/N ≥ 5; gate ≥ 5 survivors and robust scatter ≤ 0.2 mag, else
the frame is **unusable** (never header-MAGZPT fallback: MAGZPT is a
non-photometric-night value and blank on PHTCALEX = 0 frames). Header
MAGZPT is recorded as a cross-check annotation. The Mould-R vs
Cousins-R band mismatch is declared as a ±0.3 mag systematic on the
R flux scale (crossings declaration), not propagated into the
statistic. Verified against catalogued stars ≤ 0.2 mag before any
depth is quoted (SPHEREx rule).

## 5. Screening and vetoes

Static-test catalogue: the same PS1 DR2 mean objects, per band with
the transformed magnitude (objects lacking i fall back to r − 0.117),
≥ 3 detections — the flux-consistent catalogued-static veto of
`sglsurvey/vetting.py` at search radius 30″, factor 2. No archive
psfcat screen (PTF level-1 ships a SExtractor catalogue that the v2
engine does not consume; the layer-1 catalogue screen is not part of
the decision rule in any v2 survey). Parallax-phase test: PTF's
campaign cadence gives 1–12 calendar months per corridor — the
overlay records the 30-day day-of-year bins covered; the held-out
epoch test is not available (closed archive) and is not a veto in
any v2 survey. All other rules are annotations unless
injection-calibrated (v2 common rules).

## 6. Controls

Positive controls through the full chain (384-px sci + dmask cutouts
at the Horizons position, site 675; PTF matched filter; in-frame PS1
star calibration; `sglsurvey.control.rescore` under the frozen rule
and without the clip), predicted magnitudes = Horizons V + solar
colours (g = V + 0.25, R_AB = V − 0.15). The ZTF/PS1/DECam control
(60000) Miminko has only 10 PTF frames over 2009–2015 (4 surviving the
ZP gate) — PTF's coverage of any one asteroid's track is too sparse —
so the controls were found inside the pilot-corridor exposures by
SkyBoT (`scripts/asteroid_control_recon.py`, ross128 + vanmaanen
fields on the ecliptic): **(798452) 2012 QR36**, V 20.3–20.5, 46
frames / 11 nights in both bands (catalogue regime: median
single-frame S/N 8.5), and **(388125) 2005 UP482**, V 21.6, 35 frames
/ 11 nights (stack regime: below the single-frame limit). Results
(pre-freeze, `runs/ptf/v2/control/`): 798452 recovered at R̃ 6.7 (g,
27 epochs) / 3.1 (R, 13), rank p = 0.020 (0/48 ring controls above),
no-clip recovered 20.50 / 20.22 vs predicted 20.39 / 20.13 (+0.11 /
+0.09 mag — the flux scale confirmed against Horizons + solar
colours); 388125 recovered at R̃ 4.8 (R, 16 epochs), rank p = 0.020,
21.52 vs 21.32 predicted (+0.20). Negative controls and thresholds:
v2 null ensemble (48-offset ring primary; time scrambling and
trajectory randomisation as annotations), survey-wide max-R̃ null,
FWER ≤ 0.05. Negative controls
and thresholds: v2 null ensemble (48-offset ring primary; time
scrambling and trajectory randomisation as annotations), survey-wide
max-R̃ null, FWER ≤ 0.05.

## 7. Known survey-specific risks (to resolve at freeze)

- **Coverage is campaign-driven and very heterogeneous** (0 to ~670
  CCD exposures per corridor at discovery; ~half the corridors have
  fewer than 20). The epoch floor (min 5 epochs per cell) decides
  which corridors are searchable; the overlay grades them
  (`targets/overlay_v1.md`) and only searchable corridors enter the
  freeze — the rest are ledger-only coverage records.
- Single-phase corridors (all epochs in one 30-day day-of-year bin)
  cannot carry the parallax-phase annotation; they stay in the cell
  (the statistic does not require both phases) and are flagged.
- Depth: per-epoch R ~ 20–21.5 depending on seeing (median 2.4″);
  default m90 21.0 in both bands for the injection window when no
  prior depth exists.
- WCS sampling: the level-1 WCS is TAN-SIP + PV; the engine's
  linear-Jacobian sampling error across a 600-px stamp, measured on
  237 mask cutouts (`results/wcs_linearity_check.json`), is 0.08″
  median / 0.39″ max within ±200 px (the z = 550 AU half-span) and up
  to 1.5″ at the ±300-px corners — above half the 0.5″ locus
  tolerance, so the profile samples through the full SIP inverse
  (`linear_wcs = False`; 0.1 s per 475k points, no cost).
- Flux scale measured on the 620 pilot frames
  (`results/flux_scale_check.json`): 616 calibrated (4 scatter-gate
  failures), ZP scatter median 0.052 (g) / 0.065 (R) mag with 25 / 45
  calibrators per frame; star ZP − header MAGZPT = +3.85 (a different
  convention — the header value is never used).
- No archive difference images and no epoch hold-out (closed
  archive): the persistence test rests on the parallax-phase
  annotation and the cross-archive joint stage (a future PS1 + PTF +
  ZTF joint on a common µ reference epoch is a rebuild, not a
  relabel — lesson §10).
- T0 = MJD 56000 (PTF mid-baseline) for the standalone survey;
  the joint stage's T0 = 59800 would under-sample the µ family at
  PTF's epochs (learnings §4).

## 8. Frozen decisions (approved 2026-09-10)

1. **Dev/confirmatory split**: stratified by the WISE v2 confusion
   class (shared across surveys), seed 20260910, dev fraction 0.30,
   with **ross128, kapteyn, ltt1445 and vanmaanen forced into the
   development set** — the first three shaped the calibration and the
   tensor-build shakedown, vanmaanen's exposures fed the asteroid
   control recon. The engine's `forced_dev` mechanism records them.
2. **Searchable set**: the 40 corridors graded searchable by the
   overlay (best (endpoint, role, band) cell ≥ 5 usable precise-pass
   exposures, `targets/overlay_v1.json`); the 22 below-floor or empty
   corridors enter the ledger as coverage-without-statistic and are
   not part of the family.
3. **Flux scale**: per-frame PS1 DR2 mean-star zero points through the
   identical matched filter, g = PS1 g, R = Jordi (2006) Cousins R
   + 0.21 mag (AB), dwarf-locus calibrator rule, gate ≥ 5 stars and
   robust scatter ≤ 0.2 mag else the frame is unusable; header
   MAGZPT never used. ±0.3 mag Mould-vs-Cousins systematic declared on
   the R scale, not propagated.
4. **WCS sampling**: full TAN-SIP inverse (`linear_wcs = False`); the
   linear-Jacobian error (0.39″ max inside the locus span, 1.5″ at
   stamp corners) exceeds half the 0.5″ locus tolerance.
5. **Grid and epochs**: NZ = 192 uniform in 1/z, 5 × 5 µ on ±1″/yr,
   T0 = MJD 56000, ZP 25 AB, single-epoch clip |S_e| ≤ 5, per-frame
   weight cap 20×, min 5 epochs per cell, no epoch hold-out split
   (closed archive; the post-hoc test is unavailable, not an
   annotation).
6. **Quality masks**: primary = fatal dmask template 65533 + metadata
   seeing ≤ 4″; strict = primary + seeing ≤ 2.5″ + |moonillf| ≤ 0.8;
   loose = fatal template only. `photcalflag` and `infobits` recorded,
   not gated. All frames are 60 s object exposures in g or R (adapter
   gate fid ≤ 2).
7. **Static veto**: PS1 DR2 mean objects of the corridor's PS1 screen
   cone, per band with the transformed magnitude, ≥ 3 detections,
   search 30″, factor 2, at full weight (scie-direct substrate, no
   difference-fraction scaling); the only calibrated veto.
8. **Positive controls**: (798452) 2012 QR36 (catalogue regime, both
   bands, `configs/asteroid_control_v1.json`) and (388125) 2005 UP482
   (stack regime, R, `configs/asteroid_control_388125.json`), both run
   before the freeze and recovered (§6).
9. **Injection window**: default m90 21.5 ± 2 mag in both bands (no
   prior depth for this archive); 400 injections per cell over the
   four temporal families; completeness reported at threshold and
   final-candidate level.
