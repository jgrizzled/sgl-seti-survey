# PS1 pilot v1.0 — stage summary (2026-08-20)

Corridors `ross128`, `epsind` (A and B endpoints), `proxima`; hypotheses
`ps1-hypotheses-v1.0`; registry `pilot_wise_2026.yaml` (shared, hash as
recorded in every run config). All runs under `runs/panstarrs/` (22 GB,
regenerable; gitignored).

## Coarse discovery (`coarse_v1`)

| corridor | warps listed | skycells | coarse hits (rx / tx) |
|---|---|---|---|
| proxima | 302 | 3 | 193 / 193 of 302 |
| ross128 | 455 | 4 | 249 / 249 of 455 |
| epsind | 594 | 4 | A 280 / 138; B 316 / 200 of 594 |

1,351 distinct Observation records (warps, MJD 55086–57060 = 2009-09 →
2015-01), 3,890 coarse evaluations, 1,818 hits. Listing service
`ps1filenames.py` sampled on a 0.12° hex grid over each discovery cone
(13–19 s per endpoint); one full mask per skycell cached for its WCS.

## Precise pass (`precise_v1`)

1,818 evaluations on full skycell masks (3.3 MB each; ~1,000 fetched in
5.6 min with 8 workers): **892 usable, 463 partial, 457 unusable, 6
unknown** (mask not served). Mean usable fraction along the locus
0.5–0.8; 70 % of usable evaluations have disjoint covered z-intervals
(GPC1 cell/OTA gaps cross the arc). 617 distinct usable warps.

Mask lesson: `CONV.BAD` (8192) alone marks the diagonal OTA-gap bands
where the resampler interpolated across missing data; the image is
finite there, so the mask (not image NaNs) must define usability.

## Catalog screening (`screen_v1`)

DR2 `detection` cone queries (0.19–0.26° radius): 472,777 (proxima),
109,330 (ross128), 191,388 (epsind) detections; `mean` tables 40,938 /
44,783 / 64,201 objects. `obsTime` is ~50–60 s after warp `MJD-OBS`
(match window ±120 s). **6,040 ScreenMatch records** within 10″.
Recurrence: bin occupancy saturated (background 0.93–1.00 per bin per
visit); fixed-z point filter (1.5″) best cases 2–6 of 29–33 visits with
support scattered along the arc (different field stars), one static
star at Ross 128 rx (z ≈ 577 AU) and one at Proxima tx. Nothing
track-following. `runs/panstarrs/screen_v1/recurrence_report.json`.

## Cutouts and flux maps

617 image + weight fitscut cutouts, 400–3,600 px (19 GB). `.wt` is
variance (median 160 vs robust background variance 161 in the probe).
Matched-filter S distribution is ~1.6× wider than unit Gaussian
(resampling-correlated noise) — handled by the empirical control
thresholds, never by the formal variance.

## Flux scale (star calibration)

Per-warp zero point from DR2 `mean` stars through the same matched
filter: 345 of 416 maps with ≥ 5 calibrators (Ross 128's small
high-latitude cutouts often have < 5 → global per-filter fallback;
7 fallback epochs dropped for header ZPs > 0.5 mag off the filter
median, i.e. non-photometric):

| filter | ZP_star − FPA.ZP | MAD | N |
|---|---|---|---|
| g | +3.573 | 0.069 | 40 |
| r | +3.484 | 0.089 | 53 |
| i | +3.580 | 0.081 | 114 |
| z | +3.152 | 0.069 | 73 |
| y | +3.129 | 0.123 | 65 |

= 2.5 log10(EXPTIME) (4.08/4.00/4.13/3.69/3.69) minus ≈ 0.5 mag of
Gaussian-kernel throughput — the same loss the ZTF asteroid control
measured, now calibrated per frame. Frame-to-frame scatter (~0.15 mag)
tracks seeing, which is why the per-frame calibration matters.

Bug found and fixed during the pilot: a per-corridor fallback silently
used the raw header ZP where a corridor had no calibrated frames in a
filter (Ross 128 g/r came out 3.5 mag too shallow in the first
calibration run; that run's records were discarded before any report).

## Sample tensors (`calib_v1/tensors`)

360-node 1/z grid × 5×5 µ grid × 9 trajectories; skycell duplicates
collapsed (one exposure on two overlapping skycells). Epochs per
endpoint-role: ε Ind A 114/109, ε Ind B 110/104, Proxima 68/68,
Ross 128 93/93. **Parallax phase split is ~97:3** everywhere (ε Ind
110:4, Proxima 67:1, Ross 128 90:3): the 3π cadence revisits a field at
the same season each year.

## Injection calibration (`calib_v1`, AnalysisRun `run-d6ff91459819`)

320 Constraint records (302 recovery curves, 18 not-constrainable);
6 Candidate records (threshold exceedances), **all vetoed as
single-phase data**; no survivor. Thresholds T = 4.8–9.8 (8 offset
controls). Median 90 %-recovery depths (duty ≥ 0.5, |µ| ≤ 1″/yr, AB):

| filter | median m90 | range over endpoint-roles |
|---|---|---|
| g | 21.13 | 20.69–21.38 |
| r | 21.05 | 20.34–21.29 |
| i | 20.97 | 19.78–21.20 |
| z | 19.88 | 19.63–20.02 |
| y | 19.01 | 18.44–19.16 |

Proxima (Galactic plane) is 0.5–1 mag shallower in i; ε Ind A rx g has
only 7 % of cells defined (14 epochs, few survive the N ≥ 5 floor).

## Positive control (`control_v1`)

Asteroid (60000), Horizons site F51, 431 candidate nights (19.8 ≤ V ≤
21.3), 30 warps found by per-night skycell listing, 12 usable frames.
Recovered at 0″ offset: g S = 30.7 (T 3.0), z S = 15.5, y S = 11.4
(T 7.0). Recovered magnitudes on the star-calibrated scale: g 20.45
(Horizons V 20.14 + solar g−V 0.25 = 20.39), z 19.87 (V 20.32 − 0.45 =
19.87), y 19.80 (V 20.32 − 0.55 = 19.77): **flux scale verified to
≤ 0.1 mag** end to end. (z control threshold 36.7: one offset control
trajectory crossed a star; single-frame S/N 9.8 in z.)

## Throughput

Listing 1 s/call; full mask 3.3 MB ≈ 1 s; fitscut 1200² 1.5 s, 2400²
3 s (≈ 7.5 MB/s per stream, 6 streams sustained); catalog cone of
0.2° ≈ 10–25 s. Whole pilot ≈ 1.5 h wall-clock including one
re-fetch. No MAST rate limiting encountered.
