# PS1 crossings confirmatory run v1 — results

Run 2026-08-24 under threshold freeze v1.0, constructions exactly as
dev-validated (`results/dev_search_v1.md`). Plan
`configs/conf_plan_v1.json` (12 channel-B rows: 11 confirmatory units
+ 1 single-epoch; 4 channel-A rows), 432 epoch-cutouts fetched with 0
failures (`runs/ps1-crossings/products/conf/`), statistics
`results/confirmatory_v1.json`. **No candidates.**

## Channel B

| unit | rung | status | S | T | margin |
|---|---|---|---|---|---|
| ross-154 r evt-14662f1752 | 0.1 AU | ok | 2.17 | 15.81 | −13.6 |
| ross-154 g evt-321a36d0a3 | 0.1 AU | ok | 0.66 | 35.14 | −34.5 |
| gj-1276 i evt-344c12d32a | 0.1 AU | ok | 3.32 | 2.30 | **+1.01** |
| ross-128 r evt-3f7ad50e70 | 0.1 AU | ok | 5.06 | 3.61 | **+1.44** |
| van-maanen i evt-4c4ea2c397 | 1.2 R☉ / 2.5 R☉ / 0.1 AU | **track_masked** | — | — | — |
| ross-154 g evt-6f8116b6cb | 0.1 AU | track_masked | — | — | — |
| van-maanen g evt-9184328045 | 0.1 AU | ok | 1.63 | 2.47 | −0.85 |
| ross-128 r evt-d36be8d1ff | 0.1 AU | track_masked | — | — | — |
| gj-908 g evt-ee99e18ab1 | 0.1 AU | ok | 0.99 | 9.51 | −8.52 |
| gj-908 i evt-eaafbff8ea (single-epoch) | 0.1 AU | track_masked | — | — | — |

With the 5 clean dev units, **11 units produced statistics; 2
exceedances vs 2.1 expected control crossings over the frozen 19-unit
family** (~1.2 over the 11 actually searched) — the budget behaving as
designed. The huge ross-154 control values (T = 16–35) are the
warp-direct confusion floor at work: its antipode field sits near the
Galactic anticenter plane (b ≈ +7°) and the ring controls cross
catalogued stars; the threshold absorbs them exactly as frozen.

### Exceedance adjudication (frozen veto ladder, both fully walked)

**ross-128 r evt-3f7ad50e70 (S 5.06, T 3.61) — VETOED,
flux-consistent catalogued static (rung 4).** The track is annotated
(static_min 1.55″): DR2 *stack* catalog places static sources 1.46″
(r = 22.5–22.8, primary) and 1.99″ (r = 23.5) from the in-window
nodes; the measured node fluxes (r_AB ≈ 21.6–22.1 at good_frac
0.3–1.0) sit within the blended matched-filter response envelope of
those neighbours under the epoch seeing. SkyBoT census (F51, 36″
cone): no known solar-system object. Ordinary movers additionally
excluded by TTI-pair consistency.

**gj-1276 i evt-344c12d32a (S 3.32, T 2.30) — RETAINED-AMBIGUOUS,
non-promotable.** 9 distinct exposures over ~2 h on one night
(3 usable through the exact mask; 5 on the neighbouring skycell are
CONV.BAD-masked, plus 3 skycell duplicates); the three usable epochs
show a night-consistent i_AB ≈ 22.7–23.4 signal (1.1–2.3 σ per epoch)
at the near-static z = 550 track node. Ladder: (1) SkyBoT census
clean; (2) rate test cannot discriminate within a single night for
z ≥ 550 AU but excludes main-belt rates (21″ motion over the span —
absent); (3) TTI-pair consistency excludes ordinary movers; (4) the
static test *cannot veto*: the DR2 stack catalog has **no source
within 6″** to its ~23+ depth, so no flux-consistent static
counterpart exists; (5) the recurrence test cannot run — this is
gj-1276's only covered PS1 window. Under the v2 discipline the
exceedance is retained-ambiguous, not vetoed by discretion; it is
**not promotable** (promotion requires recurrence at a second covered
window) and its amplitude sits inside the frozen 1/9 control-crossing
budget. Hand-off recorded: gj-1276 has ZTF-covered wide-rung windows
(2018–2026), so the joint crossings stage (plan §11.2 step 6) is the
designated recurrence test.

### track_masked units (finding P1 at confirmatory scale)

6 of 12 confirmatory rows lost all in-window samples to the exact
warp mask — every one to CONV.BAD chip-gap bands (verified for the
van-maanen grazing family: bit 8192 at every z node in both epochs),
correlated within each unit because a window's epochs are 1–2
same-pointing TTI nights. **The van-maanen b = 0.28 R☉ grazing event
— the survey's distinctive unit — is among them**: it is recorded in
the covered-window ledger as nominal-covered / mask-unusable, and the
photosphere-grazing downlink rung ends the survey unconstrained.
Attrition at unit level (6/12) is far above the ~25 % per-epoch
Pipeline-A rate, exactly as the correlated-attrition finding
predicted.

## Channel A — all four rows constraint-only by the frozen gates

| unit | gates (windows / off-epochs / valid offsets) | S | disposition |
|---|---|---|---|
| gj-1276 i | 2 / 41 / **1** | 2.27 | constraint-only (offset support) |
| teegarden g | 3 / **5** / 0 | −5.47 | constraint-only (baseline + offsets) |
| gj-1276 g (single window) | 1 / 11 / 0 | 0.27 | constraint-only |
| gj-1276 r (single window) | **0** / 23 / 1 | — | constraint-only (in-window epochs masked) |

The sparse single-phase cadence starves the pseudo-window controls
(≤ 1 valid offset of the required 8) — the PS1 analog of the ZTF A-1.0
finding, arriving through the epoch-support half of the frozen
validity rule. No searchable A unit exists; channel A contributes
reference depths only.

## Bottom line

0 candidates. 1 vetoed exceedance, 1 retained-ambiguous single-window
exceedance inside the control budget, 6 mask-lost units including the
grazing family, channel A entirely constraint-only. Depths:
`results/completeness_v1.json` / the survey report.
