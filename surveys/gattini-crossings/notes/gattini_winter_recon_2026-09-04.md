---
title: "Palomar Gattini-IR + WINTER reachability recon for a Pipeline B crossings survey"
date: 2026-09-04
status: "recon complete — Gattini-IR DR1 catalog probed live (anonymous); WINTER closed; no public image substrate in either — item 6 BLOCKED, no freeze"
---

# Gattini-IR + WINTER recon (plan §5.8 item 6)

Goal: verify the access route for the >900 nm time-domain cell —
epochal near-IR imaging over the universal crossing windows
(`crossings/universal_v1/`), channels A (on-star) and B (antipode).
The plan row said "Palomar Gattini-IR (J, 2018–, J ~16) + WINTER
(Y/J/Hs, 2023–) … archive access uncertain — recon first
(_unprobed_)". This note replaces that with endpoint facts.

Headline: **neither archive can carry a crossings search today.**
Palomar Gattini-IR (PGIR) has a public DR1, but it is a light-curve
*catalog of 2MASS point sources* (NOIRLab Astro Data Lab), not an
image archive, and two properties of the catalog defeat both
channels: photometry is forced at 2MASS-epoch (1998–2001) positions,
so every in-era target star — all nearby high-proper-motion M dwarfs
and van Maanen — reads as empty sky in its own entry, and no antipode
has a catalogued source within the PSF; times are stored as float32
Julian dates quantised to 0.25 d. The PGIR epochal stacks (the real
substrate) are not served anywhere public. WINTER has no public data
release, no IRSA holding, and a login-only collaboration portal.
Three plan corrections follow (§"What this means"). Item 6 is
**blocked on a non-public substrate**, not on adapter cost.

## Service probes (2026-09-04, dev machine, all anonymous)

| Probe | Status | Notes |
| --- | --- | --- |
| Data Lab TAP `datalab.noirlab.edu/tap` (sync, ADQL) | **up** | schema `pgir_dr1`: `exposures` (4,133,007 rows), `photometry` (47.6 G rows), `sources` (148,860,045 rows) + 5 `x1p5__` cross-match tables (Gaia DR3, AllWISE, unWISE, NSC DR2, SDSS specobj). ADQL geometry (`CONTAINS/POINT`) fails on the REAL-typed `tmcra/tmcdec`; `COUNT(DISTINCT)` unsupported |
| Data Lab `queryClient` (SQL, `astro-datalab` 2.22.1, already a repo dependency) | **works** | `q3c_radial_query` cone on `sources` 0.2 s; `photometry` by `pts_key` ~0.6 s (indexed); JOIN to `exposures` on `stackquadid` fine. All queries snapshotted under `runs/gattini-crossings/recon/` |
| Data Lab SIA `sia/pgir`, `sia/pgir_dr1` | **404** | the same SIA pattern answers for `nsc_dr2` → no image service exists for PGIR |
| Data Lab storage `pgir://` | 500 / 401 | no public file tree |
| IRSA TAP + IBE | **no holding** | no `winter`, `gattini` or `pgir` table or image set; Gator `mission=WINTER` → ORA-00942 |
| `winter.caltech.edu` | login-only | instrument/spec pages + username/password form; no data, archive, or release statement |
| WINTER papers (Lourie+ 2020; Frostig+ 2022; instrument paper arXiv:2512.16753, Dec 2025) | read | pipeline `mirar` open-source; alerts to Kowalski/Fritz; **no public-data or proprietary-period statement anywhere** |

## PGIR DR1 facts (measured)

- **Era (global)**: `obsjd` 2458407.0 → 2459868.0 =
  **2018-10-15 → 2022-10-15** (MJD 58406.5–59867.5), 1,130 nights
  (`nightid` 13–1474). PGIR keeps observing; DR1 is a 4-year cut.
- **Footprint**: δ > −28.5°, 30,470 deg² (Data Lab page; paper:
  1,329 fields of 4.96° × 4.96°, 16 sub-quadrants of 1.24° each,
  21,264 sub-quadrant "matchfiles"). All 14 in-era target-channels
  are covered (80–224 epochs each; table below).
- **Depth**: `limmag` (5σ, per quadrant stack) median 14.1 J Vega
  globally; **14.3–14.7 at the target positions** (13.5 at ross-154,
  Galactic-plane). Plan said "J ~16": that is the AB-equivalent of
  the *nominal* extragalactic depth (14.9 Vega ≈ 15.8 AB) — the
  measured per-epoch Vega depth is ~14.5. Saturation `saturmag`
  J ≈ 8.5 (2018–2020) → ≈ 6.0 (2020+), mean 6.6.
- **Exposures**: stacks of dithered frames; `exptime` 71.3 s (60 %),
  72 s (28 %), 45 s (10 %), with 28.8/90/172.8 s minorities. Pixels
  8.7″ native, 4.35″ drizzled; PSF `numnoisepixels` mean 30
  (effective PSF footprint radius ≳ 15–25″ depending on which pixel
  scale the column uses — not pinned).
- **Photometry model**: PSF photometry on each quadrant stack at the
  **2MASS point-source position** (`pts_key` = 2MASS key; J < 15.5
  outside / < 13 inside the plane). Non-detections are kept as rows
  (`magpsf` NaN, `magpsflim` = 3σ limit) → the epoch list of any
  catalogued source is the complete visit list of its sub-quadrant.
  Six flag bits (paper): F1 west of meridian, F2 above saturation,
  F3 outside reliable range, F4 airmass > 2, F5 ZP off by > 0.75 mag,
  F6 contaminant brighter than m+2 in the PSF aperture. Magnitudes
  Vega, 2MASS-calibrated, with a per-exposure J−H colour term.
- **Timing** — the first defect. `obsjd` is float32 *at the source*
  (paper table types; the HDF5 matchfiles carry the same), so at
  JD 2.458×10⁶ the resolution is **0.25 d**: one night (`nightid`
  18, 5,040 quadrant stacks) has exactly two distinct values; the
  whole DR1 has 2,596 distinct `obsjd` over 1,130 nights. No
  mid-time, no start/mid convention beyond "Julian date of exposure
  start", no per-frame UT anywhere in DR1. The `filename` carries a
  per-night sequence number (`stack_fs<N>_fid<field>_q<quad>`) that
  could in principle be interpolated to ~minutes using the bin
  boundaries as clocks, but that is an unvalidated reconstruction.
- **Proper motion** — the second defect. The forced position is the
  2MASS-epoch one; every target with an in-era event at ≤ 0.1 AU is
  a nearby high-PM star, and the PGIR-era star position (events
  table) sits 15–120″ from the 2MASS entry. Measured on-star
  entries (`meanmag` vs 2MASS J): van-maanen 16.6 vs 11.7 (star's
  own entry 61″ away reads sky; a J 13.9 neighbour 14.6″ from the
  moving star reads 12.7 — the star's flux, blended and drifting),
  wolf-359 17.1 vs 7.1 (neighbours at 21″/34″ read 9.3/12.1),
  teegarden 16.6 vs 8.4, gj-1276 17.0 vs 14.0, ross-128 10.8 vs 6.5
  (28″, partial + saturated pre-2020), gj-908 9.8 vs 5.8 and
  ross-154 8.3 vs 6.2 (saturated, `psfcontam` 4). **Channel A has no
  valid on-star light curve for any in-era target.**
- **Antipodes**: nearest catalogued 2MASS source 33″ (ross-154 B) to
  117″ (wolf-359 B) from the relay position — outside any plausible
  PSF footprint. **Channel B has no photometry at all in DR1**, only
  the visit list.

## Coverage probed (visit lists from the 3 nearest 2MASS sources, 6′ cone)

Era-scope (geometry only, sunward axis crossings excluded, δ > −28.5°):
B 1.2 R☉ 16 events / 4 targets; B 2.5 R☉ 20 / 5; B 0.1 AU 28 / 7;
A 0.1 AU 28 / 7 — 4 events (one per year 2019–2022) per
target-channel. In-window epochs counted under the 0.25 d
quantisation as strict / nominal / loose (bin inside / bin centre
inside / bin overlaps the flat-chord window):

| Unit | Epochs | Nights | limmag | nearest 2MASS | 1.2 R☉ s/n/l | 2.5 R☉ s/n/l | 0.1 AU s/n/l | closest epoch (d) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| van-maanen A | 93 | 93 | 14.49 | 14.6″ | — | — | 6/6/6 | 0.94 |
| van-maanen B | 118 | 98 | 14.43 | 62.0″ | 0/0/0 | 0/0/0 | 8/8/8 | 1.06 |
| teegarden A | 142 | 138 | 14.65 | 89.5″ | — | — | 18/19/20 | 0.03 |
| teegarden B | 110 | 110 | 14.49 | 42.9″ | 0/0/0 | 0/0/0 | 11/11/11 | 1.15 |
| wolf-359 A | 111 | 109 | 14.43 | 21.5″ | — | — | 4/4/4 | 1.49 |
| wolf-359 B | 224 | 115 | 14.55 | 117.4″ | **2/2/4** | **4/4/4** | 16/18/18 | 0.006 |
| gj-1276 A | 110 | 110 | 14.59 | 58.3″ | — | — | 9/9/9 | 0.49 |
| gj-1276 B | 109 | 107 | 14.54 | 87.8″ | 0/0/0 | 0/0/0 | 4/4/4 | 1.94 |
| ross-128 A | 91 | 91 | 14.26 | 28.3″ | — | — | 3/3/3 | 0.27 |
| ross-128 B | 80 | 80 | 14.68 | 84.4″ | — | 0/0/0 | 6/6/6 | 1.11 |
| gj-908 A | 98 | 98 | 14.61 | 28.6″ | — | — | 7/8/8 | 0.25 |
| gj-908 B | 168 | 86 | 14.55 | 47.7″ | — | — | 6/6/6 | 0.13 |
| ross-154 A | 103 | 103 | 13.52 | 15.0″ | — | — | 11/11/11 | 0.43 |
| ross-154 B | 163 | 154 | 14.28 | 33.3″ | — | — | 25/26/26 | 0.37 |

Grazing family: the **van-maanen April deep-graze family
(b 0.247–0.256 R☉, 2019–2022) has zero in-window epochs** (closest
1.06 d against a ±0.32 d / ±0.67 d window); only **wolf-359 B**
(b 0.72 R☉, 2020-09-05 and 2021-09-05) is covered at 2.5 R☉ (2
epochs each, closest 0.006 d) and marginally at 1.2 R☉. The 0.1 AU
rungs (±5–6 d windows) are covered at 3–26 epochs per unit — the
cadence (4–5 d extragalactic, nightly in the plane) is fine for the
wide rung, and that is where a PGIR-image survey would live.
Full per-event numbers: `results/recon_scan_v0.json`
(`scripts/recon_scan.py`; snapshots `runs/gattini-crossings/recon/`).

## What this means for the survey design

1. **DR1 is not a search substrate.** A catalog of forced photometry
   at catalog-epoch positions cannot deliver channel A for
   high-proper-motion targets (the very targets that have small b),
   nor channel B at empty-sky antipodes; and 0.25 d time bins blur
   the grazing windows by ~40 % of their width. What DR1 *does*
   deliver — and this note records — is the **visit list**: PGIR
   imaged every in-era target-channel 80–224 times at J ~14.5, with
   wolf-359 B inside its 2.5 R☉ window twice.
2. **The substrate is the epochal quadrant stacks**, which are not
   public (no Data Lab SIA, no IRSA holding, no storage tree). A PGIR
   survey therefore needs one of: (a) a data request to the PGIR
   team (Caltech/IPAC; De, Kasliwal) for stack cutouts — or their
   forced-photometry run at PM-propagated star positions and at the
   14 channel positions over the 14 × 4 windows, a small ask; (b) a
   future image release. Either is a collaboration/queue decision,
   not an adapter task — surfaced for the user.
3. **Plan corrections.** (i) Depth: per-epoch J ≈ 14.5 Vega
   (~15.4 AB), not "J ~16"; a PGIR survey is ~5 mag shallower than
   ZTF/ATLAS, so the 0.1 AU cone-power floor scales from the ATLAS
   26–55 kW to **MW class** (×50–100 in flux) — "opening the cell",
   as the plan already framed it, not deep exclusion. (ii) Wavelength: PGIR is J only
   (~1.17–1.33 µm); **1064 nm is not in J**. The 1064 nm cell needs
   WINTER Y (0.97–1.07 µm) or a Y/z-band archive; Gattini opens the
   ~1.25 µm broadband cell only. (iii) Era: DR1 ends 2022-10-15;
   2022–2026 exists on disk at Caltech but not publicly.
4. **WINTER**: Y/J/Hs, 1.1″ pixels, 2.3–3.6″ FWHM, 120 s × 8 dithers,
   single-visit J 17.6–19.3 per board, operating since June 2023 —
   exactly the >900 nm substrate (Y contains 1064 nm) at a useful
   depth, but **closed**: no public release, no proprietary-period
   statement, no archive. A re-check belongs in the §5.7 standing
   maintenance list (IRSA holdings + `winter.caltech.edu` + the
   instrument-paper series) rather than in the queue.
5. **No adapter is written** — nothing to build against. The recon
   scan script and its snapshots are the only artefacts; the
   surveys directory stays a stub until a substrate exists.

## Open items (blocking; not freeze items)

- Decide whether to ask the PGIR team for stack cutouts / a
  PM-propagated forced-photometry run over the 56 windows (14
  units × 4 years); if yes, the ask is specified above and the
  visit list in `results/recon_scan_v0.json` is the attachment.
- WINTER: add to §5.7 maintenance (re-probe IRSA, Data Lab and the
  WINTER portal quarterly; watch arXiv for a DR announcement).
- If a PGIR image route opens: pin the PSF-footprint pixel scale
  (`numnoisepixels` units), the stack timing keyword (per-frame UT
  is presumably in the stack headers), and the flag-bit → fatal
  template, then run the standard era scope and freeze.
