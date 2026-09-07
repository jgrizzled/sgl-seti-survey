---
title: "High-energy archives (Chandra / XMM / eROSITA / Swift / Fermi-LAT) reachability recon"
date: 2026-09-07
status: "recon complete — every route probed live, anonymous end-to-end; freeze not started"
---

# High-energy archives recon (plan §5.26)

Goal: verify access routes and coverage for the two uses named in the
plan row — (i) opportunistic high-energy coincidence with the universal
crossing windows (`crossings/universal_v1/`, channels A on-star and B
antipode, rungs 1.2 R☉ / 2.5 R☉ / 0.1 AU, flat-chord windows as in
every crossings survey) and (ii) persistent-source tests on the 88
anti-star corridors (Pipeline A step 1, catalogue level). The row was
"not adopted; no recon". This note replaces it with endpoint facts.

Headline: **the crossing windows are a wide-field-monitor cell, not a
pointed-telescope cell.** Two of the five archives are blind to the
windows *in principle*: XMM-Newton (pointed + slew) must keep the
solar aspect angle in 70°–110° and eROSITA scans the great circle at
~90° elongation, while both channels sit at elongation 180° ± 6° inside
every window — the same Gate I elongation theorem that closed the WISE
and SPHEREx crossings (§5.3). The pointed archives that *can* point at
opposition contributed **one** in-window observation in 27 years
(Swift XRT+UVOT on wolf-359, 2018-03-06, +2.9 d in the 0.1 AU window,
4.2 ks, with a 4.2 ks UVOT UV-grism spectrum). The substrate that
actually covers the windows is the two all-sky monitors: **Fermi-LAT**
(≥ 100 MeV, every grazing window 2008–2018-03 has 1–19 ks of in-FoV
livetime; after the 2018-03 solar-array anomaly one window in three
has none) and **Swift/BAT** (15–150 keV coded mask: half of the
grazing windows and 95 % of the 0.1 AU windows carry ≥ 1 ks within 30°
of the boresight). Both are photon-event archives at arbitrary
positions with ms timing, so the cell that opens is a **high-energy
pulse/burst coincidence cell** on the windows — a photon-starved,
nearly background-free LAT count (7 photons > 100 MeV within 1° in 3
days at the van-maanen antipode) and a BAT rate/imaging test. For the
corridors, eROSITA-DE now serves **arbitrary-position upper limits
(DR1 eRASS1 and the 2026-07-31 DR2 eRASS:3)** by API for the western
Galactic hemisphere, and 5XMM-DR15, CSC 2.1, 2SXPS/LSXPS, BAT 157-month
and 4FGL-DR4 answer 7′ cones — a catalogue-level persistent-source
screen is one afternoon of work with no new adapter.

## Service probes (2026-09-07, server, all anonymous)

| Probe | Status | Notes |
| --- | --- | --- |
| HEASARC Xamin TAP `heasarc.gsfc.nasa.gov/xamin/vo/tap/sync` | **up** | `FORMAT=csv` unsupported on sync — use `FORMAT=text` (pipe-separated + trailing count lines). `chanmaster`, `xmmmaster`, `swiftmastr`, `fermilweek`, `erassmastr` (sky tiles, no times), `xmmssc`/`xmmstack` (5XMM-DR15), `swift2sxps`/`swiftlsxps`, `erass1main`/`erass1hard`, `swbat157m`, `fermilpsc` (4FGL-DR4) all answer ADQL cones in 1–3 s |
| XSA TAP `nxsa.esac.esa.int/tap-server/tap/sync` | **up** | `xsa.v_public_observations`; **`xsa.v_slew_exposure` has a per-exposure footprint polygon — `1=CONTAINS(POINT(...), slew_exposure_fov_spoly)` answers "which slews crossed this position"** (the only slew-coverage route found; HEASARC has slew *sources* only) |
| CXC `cda.harvard.edu/cxctap` and `cda.cfa.harvard.edu/csc21tap` (CSC 2.1) | **up** (303 → results) | cxctap has no `ivoa.obscore` (tables `cxc.observation`, `ivoa.obsplan`; HEASARC `chanmaster` used instead); CSC TAP rejects ADQL geometry (`CIRCLE`/`CONTAINS`/`cone_distance` "not found") and ignores `FORMAT=csv` — RA/Dec box + VOTable parse + local cut works. Products: `cxc.cfa.harvard.edu/cdaftp/byobsid/<d>/<obsid>/primary/` (evt2 listing verified; the `cda.harvard.edu/cdaftp` host 404s) |
| eRODat `erosita.mpe.mpg.de/erodat/` | **up** | `/catalogue/SCS?CAT={DR1_Main,DR1_Hard,DR1_Supp,DR1_Clusters,DR2_Main,DR2_Hard}` cone; `/upperlimit/service` and `/upperlimit/service_multi` (≤ 10⁴ positions per POST; `DR1_eRASS1` bands 021/022/023/024/02e, `DR2_eRASSc3` bands 023/024); `/skyview/skytile_search_api/` (tile + `de_sky`); DR1 event lists at `/data/download/{tile[3:]}/{tile[:3]}/EXP_010/em01_{tile}_020_EventList_c010.fits.gz` (14.7 MB, verified; the old `/dr1/erodat/...` paths 301 → 404) |
| eROSITA-DE DR2 (2026-07-31) | read | **catalogue-only** release: eRASS:3 cumulative (2019-12-12 → 2021-06-16, 556 survey days), western Galactic hemisphere; no event data beyond DR1's eRASS1 |
| Fermi LAT weekly files `heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/{photon,spacecraft}/` | **up** | 944 weeks w009–w953 (2008-08-04 → 2026-09-07, current to the day); spacecraft files 2.6–2.9 MB (30-s attitude + livetime), photon files ~180 MB/week (`p305`) |
| Fermi LAT data server `fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi` | **works** | POST arbitrary position / radius / time / energy; query `L260907093219AD225ECD61` (van-maanen antipode, 5°, 2020-04-01 → 04, ≥ 100 MeV) served `_PH00.fits` (150 photons, 46 GTIs = 218 ks) + `_SC00.fits` in ~2 min; photon `TIME` at µs resolution |
| Swift obs tree `heasarc.gsfc.nasa.gov/FTP/swift/data/obs/YYYY_MM/<obsid>/{xrt/event,bat/survey,uvot}` | **up** | the in-window obsid 00010119025 has cleaned PC and WT event lists (`sw…xpcw3po_cl.evt.gz`) and BAT survey DPH |
| UKSSDC `swift.ac.uk/API/` (swifttools: LSXPS, XRT products, upper-limit server) | up (docs) | the `swifttools` pip package is the route to arbitrary-position XRT upper limits and on-demand products; not exercised |
| MAST gPhoton | n/a | the UV half of this row ran as §5.12 |

All 2,106 responses are snapshotted under
`runs/highenergy-crossings/recon/` (`records/query_snapshot.jsonl` +
`snapshots/`), plus the 448 LAT weekly spacecraft files
(`fermi_spacecraft/`), the LAT query products and the eRASS1 tile
194096 event list (1.2 GB in all). Scripts:
`scripts/recon_scan.py` (stages `masters slew bat lat lat01 catalogs
routes` → `results/recon_*_v0.json`), `scripts/recon_erodat_rows.py`,
`scripts/recon_summary.py` → `results/recon_summary_v0.md`.

## Geometry: which archives can see a crossing at all

At `t_ca` the Earth is on the Sun–star axis, so channel A (star) and
channel B (anti-star relay direction) are both at solar elongation
180°; over the ±5.8 d 0.1 AU window the elongation stays ≥ 174°, over
the ±0.3–0.7 d grazing windows ≥ 179°. Pointing constraints then decide
everything:

| Archive | Pointing law | Windows |
| --- | --- | --- |
| XMM-Newton (pointed and slew) | solar aspect angle 70°–110° | **invisible in principle** — measured: nearest slew over any of the 14 channel positions is 72–105 d from a `t_ca`; nearest pointed observation of wolf-359/teegarden/ross-154 sits at ±77–110 d from the windows (May/December campaigns on a March star) |
| eROSITA (SRG survey) | scan great circle ⊥ Sun (elongation ~90°) | **invisible in principle** — measured on the DR1 tile 194096 event list: the van-maanen antipode was scanned 2019-12-23 14:29 → 12-24 22:30 (9 passes of ~40 s every 4 h), 100 d before the 2020-04-02 crossing |
| Chandra | Sun > 46.4°, anti-Sun allowed (thermal limits) | allowed but **never happened**: 11 observations of the 7 stars/antipodes in 27 yr, nearest +7.0 d from a 0.1 AU `t_ca` (ross-128 field, 2023-03-25, ACIS-S), all others ≥ 14 d |
| Swift XRT/UVOT | Sun > 46°, anti-Sun allowed | **1 unit** (below) out of 177 pointed observations of the 14 positions |
| Swift BAT | same, 1.4 sr half-coded FoV | window coverage by geometry, quantified below |
| Fermi LAT | survey mode, 2.4 sr FoV, all sky every ~3 h (pre-2018) | window coverage quantified below |

So for the crossing cell XMM and eROSITA join WISE/SPHEREx under Gate I
(coverage-without-statistic is not even reachable — there are no
in-window data to ledger), Chandra is an empty pointed record, and the
survey substrate is Swift + Fermi.

## Pointed X-ray coverage (HEASARC master tables, era-mean positions; `results/recon_masters_v0.json`)

Cones of 25′ (Chandra), 20′ (XMM) and 15′ (Swift) at the 14 channel
positions over each mission's whole life, then the flat-chord window
test per event and rung:

| Channel | Target | Chandra | XMM pointed | Swift |
|---|---|---|---|---|
| A | gj-1276 | 0 | 0 | 0 |
| A | gj-908 | 0 | 3 (CGCG 381-051) | 1 |
| A | ross-128 | 1 | 0 | 4 (2026 "DragPointing", 0 s) |
| A | ross-154 | 2 (2002 ACIS-S 61 ks, 2007 HRC-I 49 ks) | 2 | 16 |
| A | teegarden | 1 (2019 HRC-I 48 ks) | 1 (2021, 33 ks) | 17 (2017-11 → 2018-01 campaign) |
| A | van-maanen | 0 | 0 | 0 |
| A | wolf-359 | 7 (2022-06 ACIS-S, 2023-11 HRC-S) | 10 (2004–06, 2021-12) | 113 (2017, 2018, 2021-12, 2022-06) |
| B | gj-1276 | 0 | 0 | 0 |
| B | gj-908 | 0 | 0 | 1 (AT2025aiad) |
| B | ross-128 | 0 | 0 | 1 (0 s) |
| B | ross-154 | 0 | 0 | 0 |
| B | teegarden | 0 | 0 | 2 (GRB230116c) |
| B | van-maanen | 0 | 4 (Mars, 2014-06) | 2 (SGWGS-194) |
| B | wolf-359 | 0 | 0 | 0 |

**In-window: one unit.** wolf-359 A 0.1 AU, `t_ca` 2018-03-03 13:24
(b = 0.79 R☉; the 1.2/2.5 R☉ windows of the same event are ±2.9 h /
±7.6 h and not covered): Swift obsid 00010119025 "CNLeo", 2018-03-06
10:02 → 23:43 UTC (+2.9 d), XRT PC 4,167 s, **UVOT UV grism 4,161 s**
(the only UVOT mode used), BAT survey 4,198 s, pointing 3.0′ from the
star. Products verified in the obs tree. The UV grism (1700–2900 Å,
R ~ 75–100) is an on-star in-window *spectrum* — it belongs to the
spectral-archive family (§5.17, channel A 0.1 AU rung) as much as to
this row; the XRT event list is the X-ray flare context for the same
window. Nothing at any grazing rung, nothing on channel B.

The teegarden Swift campaign (17 visits 2017-11-19 → 2018-01-03, 36 ks)
starts 11 days after the 2017-11-08 `t_ca` — a near miss of the kind the
GALEX and PTF recons also recorded; the XMM 2021-12 and Chandra 2022-06
/ 2023-11 wolf-359 campaigns are 90–100 d from the March windows,
exactly the XMM constraint and, for Chandra, presumably scheduling.

XMM slew survey (`results/recon_slew_v0.json`): 33 slew exposures cross
the 14 positions (0–6 each; none over ross-128 B, ross-154 B, wolf-359
B); nearest to a `t_ca` +72 d. **0 in-window.**

## Swift/BAT coverage of the windows (`results/recon_bat_v0.json`, 308 per-event queries)

For every event in the Swift era (2004-11-20 →) the `swiftmastr` rows
overlapping the 0.1 AU window with non-zero BAT survey/event exposure
within 40° of the channel position, exposure pro-rated to the window
overlap, binned by pointing offset (≤ 20° ≈ ≥ 50 % coding, ≤ 30°
partial, ≤ 40° edge):

| Channel | Rung | Windows (t_ca ≤ today) | ≥ 1 ks within 30° | ≥ 1 ks within 20° | median s (30°) | p90 s (30°) |
|---|---|---|---|---|---|---|
| A | 1.2 R☉ | 86 | 45 | 19 | 1,068 | 3,964 |
| A | 2.5 R☉ | 108 | 79 | 51 | 3,673 | 9,768 |
| A | 0.1 AU | 151 | 147 | 142 | 40,166 | 77,728 |
| B | 1.2 R☉ | 88 | 40 | 21 | 924 | 4,665 |
| B | 2.5 R☉ | 109 | 76 | 43 | 2,706 | 8,046 |
| B | 0.1 AU | 153 | 146 | 135 | 32,815 | 73,314 |

Per target-channel, grazing windows (≤ 2.5 R☉) with ≥ 1 ks inside 20°:
7–11 of 21–22 for every one of the ten grazing-family target-channels
(van-maanen A and B 11/21 and 11/22). Example: van-maanen B 2020-04-02
(b 0.25 R☉, 1.2 R☉ window ±7.6 h): 6 pointings, 6.8 ks within 30°,
5.5 ks within 20°, best 16.9° off-axis. BAT's 15–150 keV sensitivity at
a few ks partial coding is ~ 10–30 mCrab (≈ 2–6 × 10⁻¹⁰ erg cm⁻² s⁻¹) —
a burst/flare rate cell, not a persistent one. The survey-mode DPH
(5-min bins) and the event-mode data (when a GRB trigger happened to be
running) are in the obs tree per obsid; per-position BAT light curves
at arbitrary sky points need the `batsurvey` mosaicking chain (HEASoft),
not a web service — the one real adapter cost on this row.

## Fermi-LAT coverage of the windows (`results/recon_lat_v0.json`)

In-FoV livetime from the weekly spacecraft files (30-s bins: boresight
angle ≤ 60°, zenith angle ≤ 100°, `DATA_QUAL > 0`, `LAT_CONFIG = 1`,
pro-rated to the window) for every grazing-family window of both
channels in the Fermi era (2008-08-04 →) plus the van-maanen 0.1 AU
windows:

| Channel | Rung | Windows | livetime min / median / max (ks) | duty | windows < 1 ks |
|---|---|---|---|---|---|
| A | 1.2 R☉ | 73 | 0.0 / 5.3 / 13.6 | 0.14 | 12 |
| A | 2.5 R☉ | 91 | 0.0 / 14.3 / 27.6 | 0.14 | 12 |
| A | 0.1 AU (all 7 targets, `recon_lat01_v0.json`) | 127 | 0.0 / 132.6 / 280.7 | 0.13 | 1 |
| B | 1.2 R☉ | 73 | 0.0 / 5.5 / 18.9 | 0.13 | 12 |
| B | 2.5 R☉ | 91 | 0.0 / 15.3 / 45.8 | 0.14 | 14 |
| B | 0.1 AU (all 7 targets) | 127 | 0.0 / 127.5 / 269.8 | 0.13 | 7 |

The zeros are all after **2018-03-16** (the solar-array-drive anomaly,
after which the rocking profile changed): pre-anomaly every one of the
97 grazing-family windows (78 at 1.2 R☉) has ≥ 1.2 ks (median 5.7–6.6 ks at 1.2 R☉,
15–16 ks at 2.5 R☉); post-anomaly 12/33 (A) and 12/35 (B) of the
1.2 R☉ windows have none — the LAT boresight stayed > 80° from the
opposition point for the whole window (e.g. gj-1276 A 2019-09-05:
θ_min 82°) — and eight 0.1 AU windows are empty for all 11.6 days
(teegarden B 2023–2026, ross-154 B 2020–2022, teegarden A 2018) —
all post-2018, all with the anti-Sun point held > 60° off the
boresight for the whole window; pre-2018 the 0.1 AU windows carry
37–281 ks (median 148 ks). So "Fermi sees everything every three hours" holds
for 2008–2018 only; each post-2018 window needs its own livetime
number (this table supplies it). The van-maanen deep-graze family
(b 0.24–0.28 R☉, 1.2 R☉ window ±7.6 h) has 4.6–18.9 ks in 17 of 18
windows.

Photon route verified: 150 photons ≥ 100 MeV within 5° over 3 days at
the van-maanen antipode (2020-04-01 → 04), 7 within 1°, 2 within 0.5°,
none > 1 GeV within 1° — the on-position background is ~2 photons
d⁻¹ deg⁻² above 100 MeV at |b| = 58°, so a burst statistic is a
near-zero-background Poisson count with the LAT PSF (68 % ~ 1° at
1 GeV, 5° at 100 MeV) as the aperture; no LAT source within 30′ of any
of the 14 positions in 4FGL-DR4 (see catalogue table). The relevant
prior art is LAT's own GRB/flare pipelines; nothing has ever been
searched at these positions and times.

## Catalogue screens on the corridors (`results/recon_catalogs_v0.json`)

Cones at the 88 anti-star corridor positions (7′ — the 550 AU relay
parallax is 375″, so 7′ contains every corridor track) and the 7
deep-family stars (2′), 4FGL-DR4 at 30′ and BAT-157m at 12′
(`results/recon_summary_v0.md` has the full match list):

| Catalogue | Corridors with ≥ 1 source | Sources (corridors) | Deep-family stars detected |
|---|---|---|---|
| eRODat DR2_Main (eRASS:3, 2026-07-31) | 34/88 | 187 | ross-128, wolf-359 |
| eRODat DR1_Main = HEASARC `erass1main` | 26/88 | 91 | ross-128, wolf-359 |
| `erass1hard` | 0/88 | 0 | — |
| Swift LSXPS / 2SXPS | 11/88 / 6/88 | 18 / 8 | ross-154, teegarden, wolf-359 |
| 5XMM-DR15 `xmmssc` / `xmmstack` | 3/88 | 52 / 169 | gj-908, ross-154, teegarden, wolf-359 |
| CSC 2.1 | 2/88 | 18 | ross-154, teegarden |
| BAT 157-month | 0/88 | 0 | — |
| 4FGL-DR4 (30′) | 9/88 | 10 | — (none within 30′ of any of the 14 channel positions) |

- The eROSITA-DE half: **36 of the 88 corridors** (and ross-128,
  wolf-359 of the deep family) lie in the western Galactic hemisphere.
  The upper-limit API returned limits for all 36 in one POST: DR1
  eRASS1 band 024 exposure median 132 s (73–1,593), UL median
  7.1 × 10⁻¹⁴ erg cm⁻² s⁻¹ (1.1–15 × 10⁻¹⁴), counts median 0 (max 30,
  the LMC corridors); DR2 eRASS:3 exposure median 347 s (196–4,378),
  UL median 3.4 × 10⁻¹⁴ (0.55–6.3 × 10⁻¹⁴). This is the persistent
  soft-X-ray constraint on every DE corridor, ready-made.
- Crowding: gj-687 (221 matches) and gj-1221 (49) have their anti-star
  points in the LMC; sigma-dra, struve-2398, van-maanen (32),
  gj-625, wolf-1069 have 13–30. A 7′ presence test is meaningless
  there; the SGL-track test (a 20–375″ annual parallax that no
  catalogued source shows) is the actual screen.
- Closest matches to a corridor centre: **teegarden's antipode has a
  persistent X-ray source at 1.0′** (LSXPS J145306.5−165245 and
  3eRASS J145306.7−165249, positions agreeing to ~5″ across 2005–2025
  Swift stacks and the 2019–21 eRASS:3) — the first concrete
  track-test case: a fixed position over 20 years excludes the
  z ≤ 10,000 AU relay track (≥ 20″ parallax) once the per-epoch LSXPS
  detections are pulled; luyten-star (2SXPS/LSXPS 1.4′), gj-687
  (5XMM 1.0′ ×3, eRASS 1.2′), gj-783 (eRASS 1.2′, det_like 235),
  gj-1221 (eRASS 1.3′, det_like 129), luhman16 (CSC 1.05′), the
  groombridge-34 / struve-2398 / gj-338 pairs (1.2–1.9′). Nothing
  within 1′ of any corridor centre.
- On-star (channel A context): wolf-359 (5XMM 11 detections / 121
  stacked rows, LSXPS 3, eRASS1 23 ct in 70 s), teegarden (5XMM 5,
  CSC 2, LSXPS 1), ross-154 (5XMM 6, CSC 4, 2SXPS 2), gj-908 (5XMM 6),
  ross-128 (eRASS1 14 ct in 96 s); gj-1276 and van-maanen undetected
  everywhere. The quiescent + flaring X-ray star is the on-star
  false-positive population for any channel-A X-ray unit, exactly as
  the GALEX recon found for NUV/FUV flares.

## What this means for the survey design

1. **Reframe the row.** "Chandra / XMM / eROSITA / Swift / Fermi
   event products for coincidence with the windows" is, after the
   recon, "Fermi-LAT + Swift/BAT all-sky photon-event coincidence with
   the windows, plus one Swift XRT/UVOT unit". XMM and eROSITA are
   structurally blind to every window (Gate I), Chandra is empty, and
   the pointed X-ray archives contribute exactly one 0.1 AU unit in 27
   years.
2. **The cell is a high-energy pulse/burst cell** (LAT ≥ 100 MeV;
   BAT 15–150 keV), nearly background-free at the LAT and
   systematics-limited at the BAT, on 97 pre-2018 + 85 post-2018
   grazing-family windows (2.5 R☉; 78 + 68 at 1.2 R☉) and every 0.1 AU
   window. It is the high-energy
   twin of the TESS/GALEX pulse cells (§5.6, §5.12): a photon-count
   excess inside the window versus the same position's out-of-window
   count — the pseudo-window family is *available* here (unlike WISE
   Gate II) because the monitors observe the position continuously.
   Order of magnitude (not injection-calibrated): a LAT non-detection
   at ~3 photons in 5 ks of in-FoV livetime with ~5,000 cm² effective
   area is F ≲ 10⁻⁷ ph cm⁻² s⁻¹ (> 100 MeV) ≈ 6 × 10⁻¹¹ erg cm⁻² s⁻¹
   — orders of magnitude above the optical cells in energy flux, but a
   band nobody has constrained at these positions and times; the
   power convention (flat-top footprint of radius b_rung, P = F·πb²,
   as in the §5.17 freeze) is fixed at the freeze, not here.
3. **Units are cheap to build and free of adapters on the LAT side**:
   the data server query (or the weekly photon files + `gtselect`)
   returns photons at any position/time; the livetime per window is
   already computed here from the spacecraft files. The BAT side needs
   HEASoft `batsurvey` for arbitrary-position light curves (the
   transient monitor covers catalogued sources only) — that is the
   adapter cost of the row.
4. **Corridors**: the persistent-source screen is catalogue-level and
   a one-day job: eRODat upper limits at all 88 anti-star positions
   (DR1 eRASS1 and DR2 eRASS:3 for the western half), 5XMM-DR15 /
   CSC 2.1 / 2SXPS-LSXPS / BAT-157m / 4FGL cones, and a per-hit
   SGL-track test (a 550 AU relay moves 375″ yr⁻¹ peak-to-peak; the
   eRASS1 vs eRASS:3 epoch pair and multi-epoch XMM/Chandra
   detections resolve that trivially). Sources inside 7′ are common
   at eRASS:3 depth (2–4 per cone), so the test is a track test, not
   a presence test.
5. **Two hand-offs to other rows**: the Swift UVOT UV grism of
   wolf-359 in the 2018 0.1 AU window → spectral-archive family
   (§5.17 channel A, a new instrument and band for that ledger); the
   XRT event list of the same visit is the X-ray flare context for
   the §5.12 GALEX stellar-flare veto design.
6. Observer: Swift and Fermi are ~550 km LEO — `universal_v1`
   Earth-center events are valid under the standing 0.010 R☉ / ~4 min
   budget. No new observer list.
7. Runs on the server as-is: ~1 GB of weekly spacecraft files are
   cached; per-window photon pulls are ~MB; a full BAT survey chain
   for ~200 windows is HEASoft-bound, hours not days.

## Recommended freeze items (for the user; not started)

1. **Adopt the row as a Fermi-LAT (+ Swift/BAT) crossings survey**,
   channels A and B, rungs 1.2 R☉ / 2.5 R☉ / 0.1 AU, era 2008-08-04 →
   2026-09, on the seven deep-family targets; retire XMM and eROSITA
   from the crossings half of the row as Gate-I structural nulls (no
   ledger rows possible) and Chandra as an empty pointed record.
2. **Primary statistic** = LAT in-window photon count within an
   energy-dependent PSF aperture (≥ 100 MeV, 1° at 1 GeV) against a
   pseudo-window ensemble drawn from the same position ±60 d with
   matched livetime (the exposure is known per pseudo-window from the
   spacecraft files); secondary = a burst statistic on sub-window
   timescales (the TESS pulse-cell construction). Declare the
   livetime floor (1 ks) below which a window is
   coverage-without-statistic (the post-2018 zeros).
3. **BAT arm as a second family**: `batsurvey` light curves at the
   14 positions over the windows with ≥ 1 ks inside 20°; rate excess
   vs the same position's out-of-window survey points; systematics
   controls from the 157-month mosaics. Decide at freeze whether the
   BAT arm is confirmatory or `constraint_only` given the mCrab floor.
4. **The single pointed unit** (wolf-359 A 0.1 AU 2018, XRT PC
   4.2 ks) as a pre-registered look: on-star X-ray light curve in the
   window with the stellar-flare veto from §5.12; UVOT grism handed
   to §5.17.
5. **Corridor screen** as a separate one-day Pipeline A catalogue
   task (item 4 above) with its own tiny hypothesis file: eRODat UL
   at 88 positions, cone catalogues, SGL-track test on every hit,
   ledger rows `catalogue_screen_only`.
6. **Positive controls**: a LAT-detected GRB or a known flaring
   blazar epoch processed through the same window/pseudo-window
   machinery; for BAT, a catalogued transient at the window date.
7. **Before any data are touched**: pin the LAT event class/type
   (`P8R3_SOURCE`, `evtype 3`), the zenith cut, the energy range and
   the PSF aperture per energy; snapshot the data-server caveats
   page; decide whether the per-window photon pulls go through the
   data server (async, ~2 min each) or the local weekly photon files
   (180 MB/week, ~40 GB for the windows' weeks).

## Plan corrections

- §5.26 description: the row's crossing half is Fermi-LAT + Swift/BAT
  (+ one Swift XRT unit); XMM-Newton and eROSITA cannot observe any
  crossing window (solar-aspect 70–110° / scan at 90° elongation) and
  belong to the corridor half only; Chandra has no in-window
  observation.
- The corridor half is catalogue-level and includes a service that
  did not exist when the row was written: eROSITA-DE upper limits at
  arbitrary positions (DR1 eRASS1; DR2 eRASS:3 released 2026-07-31,
  catalogue-only).
- Gate I of §5.3 now has four members: WISE, SPHEREx, XMM-Newton,
  eROSITA.
- Fermi's "continuous all-sky" coverage is true only before
  2018-03-16; window-by-window livetime is required afterwards.

## Open items (recon-level, cheap)

- BAT: confirm the partial-coding fraction vs offset angle
  (the 20°/30°/40° bins are geometric proxies) from the BAT
  instrument log or a `batsurvey` run on one in-window obsid.
- The teegarden-antipode source (1.0′): pull the per-observation
  LSXPS detections and the eRASS1 vs eRASS:3 positions — the first
  SGL-track test, ten minutes.
- The 2026-06 Swift "DragPointing" rows (0 s exposure) at ross-128
  are an unfamiliar mode — check the Swift status pages before
  counting 2026 BAT exposure.
- MAXI/GSC (2–20 keV, every 92 min, 2009–) and INTEGRAL are the
  obvious extensions of the monitor cell; not in the row, not probed.
