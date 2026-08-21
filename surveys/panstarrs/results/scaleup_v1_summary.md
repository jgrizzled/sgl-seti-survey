# PS1 scale-up v1 — 62 corridors (2026-08-20 → 21)

Driver `scripts/run_scaleup.sh` over the overlay queue (59 corridors +
3 pilot), hypotheses `ps1-hypotheses-v1.0`, registry
`pilot_wise_2026.yaml`. Wall-clock ≈ 5 h (coarse 19 min, precise 80 min,
screen 30 min, 15 batches × ≈12 min, calibration 3 min). Retained
products: records, snapshots, cutout index, manifests, 138 sample
tensors (18 GB); image products purged per batch (≈150 GB transient).

## Overlay (`targets/overlay_v1.md`)

62 of 77 corridors have Dec > −30°. Every corridor is single-phase
dominated (minor-phase fraction of warps ≤ 12 %); 30 corridors are
calibrator-sparse (< 40 DR2 stars with 15 < r < 20.5 within 0.05°).

## Pipeline A

| stage | count |
|---|---|
| Observations (warps) | 24,852 |
| coarse evaluations / hits | 55,736 / 23,184 |
| precise evaluations | 23,184: 10,419 usable, 6,108 partial, 6,616 unusable, 41 mask not served |
| distinct usable warps | 7,817 |
| ScreenMatch (DR2 detections ≤ 10″) | 134,874 |

Recurrence: bin occupancy saturated everywhere; fixed-z point filter
best cases 5–6 of ~30 visits (uncorrected P ~10⁻⁶ over ~5×10⁴ trials)
with support spread along the arc; 28 static stars flagged. Nothing
track-following.

## Flux scale

5,516 flux maps; 97.6 % with ≥ 5 calibrators (5,158 from the locus
cutout, 358 from the dedicated 1,600-px calibration cutout introduced
after the pilot); the rest use the global per-filter fallback. Per
endpoint-role the star-calibrated fraction is ≥ 86 % (median 100 %).

## Calibration (AnalysisRun `run-eaa6d89a1ec9`)

138 endpoint-roles × 5 filters = 690 cells; median 80 epochs per
endpoint-role; minor parallax phase median 6.5 % of epochs (max 14 %).

**5,520 Constraint records** (5,422 recovery curves, 98
not-constrainable). Median 90 %-recovery depth (duty ≥ 0.5,
|µ| ≤ 1″/yr, AB), 10–90 % range over endpoint-roles:

| filter | median m90 | 10–90 % |
|---|---|---|
| g | 21.06 | 20.75–21.40 |
| r | 20.86 | 20.41–21.30 |
| i | 20.72 | 20.12–21.13 |
| z | 19.88 | 19.45–20.35 |
| y | 18.86 | 18.48–19.29 |

**78 Candidate records** = threshold exceedances in 78 of 690 cells
(11.3 %; the 8-control threshold implies 12.5 % by construction, so the
census is at the chance rate). 72 vetoed automatically (single-phase
data or phase-split disagreement); 6 retained → stage-7 adjudication
(`scripts/adjudicate_candidates.py`, `runs/panstarrs/calib_v1/adjudication.json`):

| cell | S / T | split-half | verdict |
|---|---|---|---|
| gj-783 tx z, z≈9543 | 6.89 / 6.89 | 5.0 / 4.7 | **vetoed** — DR2 star 1.2″ from the major-phase position; i-band S = 10.6 at the same cell |
| lhs-1723 rx z, z≈829 | 7.46 / 6.08 | 5.7 / 4.9 | **vetoed** — DR2 star 1.7″; i-band S = 5.2 |
| gj65-a rx y, z≈829 | 6.04 / 5.72 | 1.0 / 6.6 | **vetoed** — not persistent (late half only) |
| 82-eri rx y, z≈7231, µ=(−0.5,+1) | 9.37 / 8.19 | 8.1 / 5.3 | **marginal** — y-only (g 2.0, i 1.8, z 2.7 σ); no CatWISE/unWISE/2MASS/Gaia counterpart within 4″ (a y≈19.8 brown dwarf would be W1≈15–16); µ at grid edge |
| fomalhaut rx y, z≈574, µ=(+1,0) | 5.42 / 5.15 | 2.4 / 5.1 | **marginal** — at threshold; z and µ at grid edge; no counterpart |
| gj-526 tx z, z≈559, µ=(−1,−0.5) | 5.15 / 5.13 | 4.6 / 3.2 | **marginal** — at threshold; z and µ at grid edge; other bands negative; no counterpart |

**No candidate survives as a detection.** The three marginal cells were
subsequently vetoed by the ZTF cross-archive test
(`surveys/joint/results/joint_v1_summary.md`: ZTF S = 0.19 / 1.42 / −0.11
along the PS1-fitted trajectories over 897–1,206 frames); they were the
3-in-690 tail expected at this threshold.

## Lessons

- Skycell-WCS cache must tolerate unserved masks (~1 %): try several
  warps per skycell (fixed after the first launch).
- Calibration cutouts must be centred on the locus cutout, not the
  corridor antipode (fitscut 400 when the antipode lies in a neighbouring
  skycell); 358 frames rescued by the calibration cutout.
- Per-batch purge keeps the footprint at ≈ 20 GB live; full-product
  regeneration is `precise_pass.py` + `fetch_cutouts.py`.
