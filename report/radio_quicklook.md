# Radio quick-look v1 — in-window epochs at the relay and star positions

Plan §5.15 item O1: hand-offs (a) and (b) of
`report/radio_crossings_ext.md` §4 plus the VLASS six of
`report/radio_crossings.md` §4. Executed 2026-09-04 (anonymous stage),
completed the same day with an OPAL login for the CASDA image stage.
**A look, not a search** — the radio decision (§5.5, geometry-only) is
unchanged; no threshold was frozen and no statistic is claimed. One
question per epoch: is there a catalogued or visible broadband
continuum source at the relay (channel B, antipode) or star (channel
A) position in the observation that fell inside a 0.1 AU crossing
window?

**Answer: no.** 17 epochs (10 ASKAP, 1 LoTSS, 6 VLASS): 0 coincident
catalogue components in the 10 epochs with a catalogue (7 of them
epoch-resolved), 0 pixels above 3σ in the 12 ASKAP and 26 VLASS
cutouts. The largest value at any position is 2.7σ (gj-1276's star,
VAST 2023-09), which a 16-epoch check of the same field shows to be
the maximum of a noise series (§4). One ASKAP epoch (a REJECTED FLASH
image) is inaccessible; its catalogue is too, and it failed the
survey's own validation.

## 1. What was looked at

Positions: the per-event ICRS positions of `crossings/universal_v1`
(channel A: the star, proper-motion-propagated to t_ca; channel B: the
beam-axis antipode at t_ca — its apparent motion over a ±6 d window is
under 1′). Coincident = a catalogue component within 20″ (about one
ASKAP restoring beam); the nearest component within the cone is
reported regardless. Cutouts: 3′ radius; `value` = the pixel at the
position, `peak` = maximum within one beam (BMAJ), rms = MAD of the
30–150″ (ASKAP) or 20–120″ (VLASS) annulus. Catalogue rms = median
`rms_image`/`noise` of the components within 0.5°.

Access (2026-09-04): VAST full-survey epochs have per-SBID component
catalogues in CASDA TAP
(`AS207.vast_extragal_dr1_<field>_sb<sbid>_components_v01`),
anonymous and epoch-resolved. RACS epochs were matched to the release
catalogues (`AS110.racs_low2_v01`, `racs_high_components_v01`,
`racs_mid_components_v01`); low2/high carry no SBID column, so the
match rests on the field and date being the release's single epoch
(consistent in both cases); the mid table's SBID column shows the
2024-11 teegarden observation (SB67840) is not in the release, so its
own Selavy catalogue was fetched instead. Per-SBID Selavy catalogues
of the VAST pilot, FLASH and re-observed RACS-mid SBIDs are level-5
files, not TAP tables — downloaded through CASDA DataLink (basic auth
with the OPAL account) and cone-searched locally. ASKAP image cutouts:
DataLink (auth) → per-product cutout token → SODA async job at
`casda_data_access/data/async` (no further auth; the astroquery.casda
flow). Access is per product: the VAST-pilot v1 products and the
REJECTED FLASH SB83234 return "no permission"; the pilot's v2
reprocessing is open. LoTSS: `lotss_dr3.main_sources` (ASTRON VO) is a
mosaic of every run on the pointing — not epoch-resolved. VLASS: CADC
CAOM2 lists every quick-look tile containing the position with exact
time bounds (VLASS4.1 planes lack them; tile header DATE-OBS used) and
CADC SODA serves cutouts anonymously.

Scripts `surveys/radio-crossings/scripts/quicklook_v1.py` and
`quicklook_gj1276_epochs.py`; results
`surveys/radio-crossings/results/quicklook_v1.json` and
`quicklook_gj1276_epochs.json`; raw TAP/DataLink/CAOM responses, the
Selavy catalogues and the 54 cutouts sha-pinned under
`runs/radio-crossings/quicklook/`.

## 2. ASKAP — catalogue and cutout

| ch | target | b | epoch (SBID, field) | Δt | catalogue | nearest | cutout value / peak / rms (mJy) | beam | verdict |
|---|---|---|---|---|---|---|---|---|---|
| B | **wolf-359** | 0.71 R☉ | VAST SB52549 VAST_2257-06, 2023-09-03 | −2.1 d | per-SBID (epoch-resolved), 77 in 0.5° | 288″, 1.4 mJy | +0.03 (0.1σ) / 0.29 / 0.22 | 12.8″ | empty |
| B | **wolf-359** | 0.71 R☉ | VAST SB65727 VAST_2257-06, 2024-09-10 | +5.6 d | per-SBID (epoch-resolved), 107 | 292″, 1.8 mJy | −0.12 (−0.7σ) / 0.13 / 0.17 | 13.1″ | empty |
| B | ross-154 | 3.30 R☉ | RACS-high SB34957 RACS_0651+23, 2021-12-29 | −3.4 d | RACS-high release (single epoch), 77 | 183″, 8.1 mJy | −0.06 (−0.4σ) / 0.20 / 0.16 | 12.3″ | empty |
| B | ross-128 | 1.85 R☉ | VAST pilot SB32330 VAST_2338+00, 2021-09-21 | +1.7 d | Selavy file (v2 reprocessing; v1 no permission), 49 | 219″, 9.1 mJy | −0.17 (−0.7σ) / 0.22 / 0.26 | 12.7″ | empty |
| B | van-maanen | 0.25 R☉ | RACS-low SB38682 RACS_1237-06, 2022-03-29 | −4.8 d | RACS-low2 release (single epoch), 83; Stokes V table empty | 29″, 3.1 mJy | +0.18 (0.8σ) / 0.34 / 0.23 | 13.5″ | empty at the position (§4) |
| B | teegarden | 0.98 R☉ | FLASH SB84179 FLASH_389, 2026-05-03 | −3.2 d | Selavy file, 234 (2-h FLASH integration, rms 0.08) | 163″, 0.7 mJy | +0.03 (0.4σ) / 0.10 / 0.07 | 24.3″ | empty — the deepest antipode epoch (5σ ≈ 0.36 mJy) |
| B | van-maanen | 0.24 R☉ | FLASH SB83234 FLASH_497, 2026-03-29, **REJECTED** | −5.0 d | no permission (failed validation) | — | no permission | — | not checked; the survey rejected the image |
| A | gj-1276 | 0.83 R☉ | VAST SB52549 VAST_2257-06, 2023-09-03 | −1.6 d | per-SBID (epoch-resolved), 70 | 219″, 10.7 mJy | **+0.64 (2.7σ)** / 0.75 / 0.24 | 12.8″ | noise maximum of 16 epochs (§4) |
| A | van-maanen | 0.33 R☉ | VAST SB66449 VAST_0037+06, 2024-10-07 | +1.2 d | per-SBID (epoch-resolved), 59 | 266″, 27.8 mJy | +0.33 (0.9σ) / 0.63 / 0.38 | 14.6″ | empty |
| A | teegarden | 1.08 R☉ | RACS-mid SB67840 RACS_0248+18, 2024-11-11 | +3.6 d | Selavy file (release = 2021 epoch), 87 | 156″, 2.5 mJy | +0.12 (0.7σ) / 0.37 / 0.18 | 11.9″ | empty |

Nothing within 20″ of any position in any catalogue, and no cutout
pixel above 3σ. Seven of the nine checked epochs are genuinely
epoch-resolved (per-SBID catalogue plus the epoch's own image). The
antipode channel's first content: in the two epochs where wolf-359's
relay position sat 0.8° from a VAST field centre the pixel at the
position is 0.03 and −0.12 mJy against 0.2 mJy noise, and the
teegarden 2026 antipode is empty to 0.07 mJy rms in a 2-h FLASH
image. A 3σ point source in these epochs is 0.2–0.7 mJy at
856–1656 MHz.

## 3. VLASS — cutout stage (3 GHz quick-look, 1″ pixels)

For each of the six v1 rows, every quick-look tile covering the
position (all epochs) was cut out; the in-window and tolerance-edge
tiles are listed.

| ch | target | b | epoch, tile | Δt (exact) | window | value | peak | rms | verdict |
|---|---|---|---|---|---|---|---|---|---|
| A | gj-908 | 12.3 R☉ | VLASS1.1 T11t36.J235001+023000, 2017-09-24 | +2.7 d | in | +0.29 (2.3σ) | 0.30 | 0.13 | sub-threshold; 1.2σ at the same star in 4.1 |
| A | gj-908 | 12.3 R☉ | VLASS4.1 T11t36.J235001+023000, 2025-09-24 | +2.8 d | in | +0.15 (1.2σ) | 0.30 | 0.13 | noise |
| A | ross-154 | 3.4 R☉ | VLASS1.2 T05t29.J185051-233000, 2019-07-02 | −1.5 d | in | +0.10 (0.6σ) | 0.22 | 0.16 | noise |
| A | teegarden | 1.1 R☉ | VLASS2.2 T15t04.J025431+163000, 2021-11-06 | −2.0 d | in | −0.06 | −0.02 | 0.11 | noise |
| A | gj-1276 | 0.83 R☉ | VLASS2.1 T09t35.J225404-063000, 2020-09-13 | +8.8 d | **outside** (v1 tolerance edge; half-window 5.9 d) | +0.07 | 0.30 | 0.14 | noise, out of window |
| B | wolf-359 | 0.72 R☉ | VLASS2.1 T09t35 (two overlapping tiles), 2020-09-13 | +8.3 d | **outside** (tolerance edge) | −0.15 / −0.16 | 0.19 / 0.22 | 0.14 | noise, out of window |

The exact CAOM tile times settle the two tolerance-edge rows: both
fall 2.5–3 d outside the strict window, so the strict VLASS in-window
set is four channel-A epochs, none at the antipode. Out-of-window
tiles (19 more, 26 in all) are in the JSON for reference; the largest
pixel value at any position in any epoch is the gj-908 2.3σ above,
and the largest one-beam peak is 2.9σ (gj-1276, VLASS3.1, 2023) — the
expected maximum of ~26 noise samples. A 3σ point source in these
tiles is 0.35–0.5 mJy.

## 4. Notes

- **gj-1276 star, VAST 2023-09-03 (2.7σ pixel, 0.64 mJy).** The
  in-window epoch is the largest value at that position among the 16
  VAST_2257-06 epochs CASDA holds (2023-07 → 2026-06,
  `quicklook_gj1276_epochs.json`): mean +0.09 mJy, scatter 0.25 mJy,
  16-epoch mean 0.09 ± 0.06 mJy (1.5σ) — no persistent source; no
  per-SBID catalogue component within 60″ in any of the 10 epochs
  that have one; the other two September epochs (2024, 2025) read
  +0.48 and −0.04 mJy. The maximum of 16 Gaussian samples reaches
  2.7σ about 5 % of the time. Not a detection; not a flare candidate
  either (a 12-min VAST epoch would need ≳ 1 mJy to clear 5σ).
- **van-maanen B 2022 nearest component (29″, 3.1 mJy).** Outside the
  20″ coincidence radius and, in the cutout, the pixel at the
  position is 0.18 mJy (0.8σ) — the 13.5″ beam separates them. The
  antipode is a fixed sky point (b = 0.25 R☉, the deep-graze family),
  so this is a background source that will sit 29″ away in every
  epoch; no action.
- **No power constraint is derived.** The cell is broadband
  continuum, the surveys are not injection-calibrated for it here, and
  the decision rule of §5.5 stands. The numbers above are the
  epoch-by-epoch noise floors a future analysis would start from.
- **LoTSS** (teegarden A, P043+19, 2023-11-04, −3.9 d): DR3 mosaic
  catalogue, 503 sources in 0.5°, nearest 71″ at 0.28 mJy, rms
  0.08 mJy — empty, but a mosaic of all runs, never epoch-resolved.

## 5. Access notes for the yearly re-run

CASDA DataLink is HTTP 401 anonymously; with `OPAL_USERNAME` /
`OPAL_PASSWORD` in the environment or the repo `.env` the script runs
end-to-end (basic auth on `data.csiro.au/casda_vo_proxy/vo/datalink`,
token-only SODA jobs at `casda.csiro.au/casda_data_access/data/async`,
~25 s per cutout, results served from `casda.csiro.au/download/web/`).
Per-product permissions vary (pilot v1 products and REJECTED images
are closed). The VAST cadence over VAST_2257-06 (every ~2 months)
makes the wolf-359 / gj-1276 September windows a standing re-run:
one command, ~10 min.
