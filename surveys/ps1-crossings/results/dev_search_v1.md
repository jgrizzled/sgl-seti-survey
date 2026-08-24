# PS1 crossings dev search v1 — results

Run 2026-08-24 under threshold freeze v1.0 (`dev_search.py`; plan
`configs/dev_plan_v1.json`, statistics `results/dev_search_v1.json`,
saturation gate `results/saturation_verify_v1.json`; cutouts
`runs/ps1-crossings/products/dev/`, 83 MB, 108 epoch-cutouts, 0 fetch
failures). Dev split: channel B wide rung, targets teegarden +
wolf-359 (6 units); channel A dev is empty by freeze — only the
saturation gate ran.

## Channel B — 6 dev units, 0 exceedances

| unit | status | S | T | margin | k |
|---|---|---|---|---|---|
| wolf-359 g evt-0639d9d72c4b | ok | 0.29 | 1.65 | −1.36 | 3.86 |
| wolf-359 i evt-0639d9d72c4b | ok | 1.13 | 1.91 | −0.78 | 1.99 |
| wolf-359 r evt-0639d9d72c4b | ok | 1.11 | 5.21 | −4.10 | 1.00 |
| wolf-359 i evt-a2c4ffdd0d1d | controls 7/8 | 0.43 | 2.34 | −1.92 | 1.00 |
| teegarden r evt-c349f1a1e62c | **track_masked** | — | (1.07) | — | 2.41 |
| wolf-359 i evt-f80401623d4d | ok | 1.15 | 1.84 | −0.69 | 1.00 |

Machinery validated end-to-end on the warp-direct substrate:
star-calibrated per-frame ZP (zp_star hit for ~75 % of epochs, per-band
median fallback otherwise), skycell-duplicate collapse at 0.0002 d by
best good_frac (no duplicates present in the dev epochs — all are
genuine TTI exposures 12–17 min apart), empirical variance rescale
k = 1.0–3.9 from 6–21 off-window samples per unit, nested-z max, ring
controls, catalogued-static annotation (nearest static source 3.1–19.4″
from any track node; none within the 2″ annotation radius).

## Dev findings (recorded; no freeze amendment required)

**P1 — correlated exact-mask attrition can kill a whole unit.** Every
track node of the teegarden unit, in all three of its same-night TTI
exposures, lands on mask bit 8192 (`CONV.BAD`, the resampled OTA
chip-gap band; Pipeline-A lesson 3). Same-night visits repeat the
pointing, so the ~25 % per-epoch nominal→exact attrition arrives
*correlated per unit*, not independently. The unit is reported
`track_masked` (no usable in-window data); this is a data-availability
outcome, not a statistic failure — the frozen constructions are
untouched. Confirmatory consequence: units are classified
`track_masked` when the exact mask removes all in-window samples, and
the covered-window ledger must record them as nominal-covered /
mask-unusable.

**P2 — partial ring controls.** One wolf-359 unit has 7/8 ring
trajectories with valid samples (one ring node masked); T is the max
over available controls and the deficit is recorded, following the ZTF
confirmatory precedent (its retained 5/8-control row).

**P3 — off-window realization.** The freeze left the channel-B
off-window sample implementation-defined; the dev realization —
same-band primary era epochs at the event's fixed tabulated antipode
position — yields 6–21 clipped samples per unit and k values (1.0–3.9)
consistent with the confusion-floor expectation. Adopted for the
confirmatory run. The teegarden unit sits at the n = 6 floor; units
below 6 off-samples would run at k = 1 and be flagged.

## Saturation gate (channel A) — PASS, frozen levels stand

32 bright DR2 stars (13.4–14.9 mag, both dev corridors, all five
bands) measured on in-band warp cutouts, plus the header
`CELL.SATURATION` point-source conversion:

| band | m_sat (header, median) | frozen E | empirical mask behaviour |
|---|---|---|---|
| g | 12.4 | 14.0 | STARCORE cores at 13.4–13.5; clean by 14.6–14.9 (one stray SAT pixel at 14.86) |
| r | 12.1 | 14.0 | mixed SAT cores 13.4–13.6; mostly clean at 14.7 |
| i | 12.9 | 14.0 | STARCORE at 13.4–13.6; sporadic SAT pixels to 14.8 |
| z | 12.0 | 13.0 | clean at 13.4+ |
| y | 11.4 | 12.0 | clean at 13.4+ |

True photometric saturation is ≥ 1 mag *brighter* than every frozen
exclusion level, so the levels stand (the > 0.5 mag amendment trigger
is in the unsafe direction only). The binding effect near 13.5–15 is
fatal-template mask attrition (STARCORE 4096 + sporadic single-pixel
SAT 32), which acts automatically through the good_frac ≥ 0.7 gate and
in the conservative direction: over-bright targets lose epochs, never
gain false ones. No amendment.

## Verdict

Channel-B machinery is validated clean under freeze v1.0 with no
amendment (the ZTF v1.1/v1.2 constructions imported at freeze time did
their job — no PS1 analog of the ZTF dev A-channel failures exists,
because there is no reference image). The confirmatory run may proceed:
B under v1.0 (with the P1 `track_masked` classification and the P3
off-window realization), A under the frozen gates with the saturation
levels verified.
