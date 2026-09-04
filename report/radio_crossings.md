# Radio crossings coverage & scope decision (Pipeline B, step 7)

Closes the crossings programme (plan §11.2 step 7). Executed
2026-08-24 as a metadata-only coverage intersection — no signal
statistic was formed. Question decided: do the crossing/corridor
geometry products feed radio signal searches directly, or stay
geometry-only?

## 1. What was intersected

- **Breakthrough Listen Open Data archive** (public API; per-file
  target/position/MJD/telescope/frequency records; raw responses
  snapshotted under `runs/radio-crossings/`): one 0.5° position cone
  per registry target at the star (channel A) and at the antipode
  (channel B) — 93 cones with in-era events, 0 query failures. Files
  deduplicated to observations in 5-minute bins; a pointing counts as
  on-source when its separation is under half the primary beam at its
  center frequency (GBT/Parkes/MeerKAT dish sizes; 0.15° default when
  the frequency is absent).
- **VLASS** (four epochs 2018–2026, per-tile observing dates from the
  NRAO dynamic-schedule summary, snapshotted): star and antipode
  positions per event vs tile epoch dates, ±3 d tile-scheduling
  tolerance (declared).
- Windows from `crossings/universal_v1`, the frozen programme ladder
  (A 0.1 / 1.0 AU; B 1.2 / 2.5 R☉ / 0.1 AU), era 2016–2028 for BL.
- Scope caveat: the BL open-data archive is not all BL data (MeerKAT
  holdings are largely absent); the statements below are about the
  public archive as snapshotted.

## 2. Results

| cell | result |
|---|---|
| **BL radio, on-star, narrow rungs (A 0.1 AU)** | **zero** GBT/Parkes/MeerKAT pointings inside any window |
| BL, on-star, narrow rungs — non-radio | 3 **APF optical spectra** of teegarden (as SO0253), 2016-11-13, 5.3 d inside an 11.5 d window |
| BL, antipodes (channel B, all rungs) | **zero** pointings within even 0.5° of any antipode in-window — nobody has ever looked |
| BL, on-star, wide rung (A 1.0 AU) | abundant: 143 in-beam observations, 22 events, 18 targets (incl. proxima-cen, barnard-era targets, teegarden) |
| VLASS, on-star, A 0.1 AU | 4 strict in-window tile epochs (gj-908 ×2, ross-154, teegarden) + 1 tolerance-edge (gj-1276) |
| VLASS, antipode, B 0.1 AU | 1 tolerance-edge tile epoch (wolf-359, VLASS2.1) |
| VLASS, grazing B rungs | zero |

## 3. Decision (recorded per plan §11.2 step 7)

**The crossings/corridor geometry products stay geometry-only; this
repository runs no radio signal search.** The empirical grounds, cell
by cell:

1. **Narrow-rung on-star radio re-analysis — moot.** The one cell
   where a targeted re-analysis would have beaten Breakthrough
   Listen's published blind searches (pre-registered short windows,
   predicted drift, lower threshold) contains no radio data at all.
   Nothing to analyze.
2. **Wide-rung on-star — already answered.** The 22 in-window BL
   events sit in windows that tile ~⅓–⅔ of the calendar; a targeted
   search adds nothing over BL's published blind nulls, which are
   hereby cited as the constraint for those cells (with their own
   thresholds/drift-range caveats).
3. **Antipode radio — confirmed virgin, and confirmed empty of
   archival coverage.** The channel our optical surveys found to be
   the workhorse has zero targeted radio coverage; the "someone would
   have noticed" objection never applied, but neither does any
   archive. This is a *future-observation* recommendation (the
   geometry products specify exactly where and when to point), not an
   archival-analysis opportunity.
4. **Identified micro-follow-ups, handed off, not executed here:**
   (a) the 4 strict VLASS narrow-rung on-star epochs + 2
   tolerance-edge rows — a six-cutout quick-look transient check
   (broadband emission only; continuum imaging dilutes CW lines by
   ~10⁹); (b) the 3 teegarden APF spectra — an optical laser-line
   inspection during a crossing window, belonging to the
   spectral-archive family (plan §3.6 row 10), not radio; (c)
   re-running this intersection if/when BL's MeerKAT holdings reach
   the public archive.

The unified covered-window ledger (step 6) plus this coverage product
is the hand-off for any future radio effort: per target × rung, which
windows exist, which have any radio-adjacent coverage, and where a
dedicated observation (antipode pointings during predicted windows)
would open cells no archive can.

## 4. The crossings programme is closed

Steps 1–7 complete: universal + spacecraft-observer crossing lists;
ZTF, PS1, WISE surveys (0 candidates; one structural null with its
two theorems); the joint stage (593-row ledger; the gj-1276 anomaly
narrowed to red/line SEDs or non-persistence and refuted for flat-SED
persistence at 90 %); and this scope decision. Standing maintenance:
yearly crossing-list refresh (2028 window end); revisit SPHEREx
crossings at 3+ annual windows; the v2-style null-ensemble redesign
if the d = 1 wide-beam rung is ever to be searched with calibrated
error rates.

## 5. Extension (2026-09-04)

`report/radio_crossings_ext.md` repeats this construction on CASDA
(RACS / VAST / EMU and every other ASKAP continuum collection) and
LoTSS DR3 pointing dates. Decision unchanged. The antipode statement
in §2–3 ("nobody has ever looked") now carries a qualifier: it holds
for targeted radio (BL) and for the grazing rungs, but six validated
wide-field ASKAP epochs sit inside 0.1 AU antipode windows (wolf-359
2023 and 2024 inside the VAST primary beam).

## 6. Products

`surveys/radio-crossings/{scripts/coverage_intersect.py, results/*}`
(`bl_inwindow_v1.ecsv`, `bl_targets_v1.json`,
`vlass_inwindow_v1.ecsv`, `coverage_v1_summary.json`); snapshots and
the VLASS tile table under `runs/radio-crossings/`.
