---
title: "screen_v1 — catalog screening over all usable corridors"
date: 2026-08-18
run_dir: "runs/wise/screen_v1 (gitignored; regenerate with scripts/catalog_screen.py + screen_recurrence.py)"
---

# screen_v1 results

Layer-1 catalog screening (plan §4.5 / §3.4 item 1) over every usable
precise-pass epoch: per corridor × mission-phase × visit, one tight TAP
query against the phase-appropriate single-exposure source table over
the union locus arc (raw responses snapshotted), each returned detection
matched against every member endpoint × role locus polyline at its
epoch. Context pulls per corridor: CatWISE2020, CatWISE2020 reject
table, NEOWISE known-SSO association list.

## Pins

- Sources: coarse_v1 + precise_v1 (usable/partial epochs only)
- Registry v1.0 `sha256:82743c09…`, hypotheses `wise-hypotheses-v1.0`
- Screen radius 10" (frozen padding), locus tolerance 2", 0.5-day locus
  bins, 5-day visit gap clustering; phase-selected tables
  (`allsky_{4band,3band,2band}_p1bs_psd`, `neowiser_p1bs_psd`)
- 140 snapshotted queries, 85 MB run dir

## Match volume (43,915 ScreenMatch records)

| Endpoint | rx matches / visits | tx matches / visits | SSO-flagged |
| --- | --- | --- | --- |
| barnard-star | 2,157 / 26 | 2,214 / 26 | 0 |
| ross-154 | 2,237 / 26 | 2,613 / 26 | 6 |
| lalande-21185 | 1,853 / 25 | 1,688 / 25 | 1 |
| alpha-cen-a | 4,080 / 25 | 3,323 / 25 | 0 |
| alpha-cen-b | 4,168 / 25 | 3,728 / 25 | 0 |
| sirius-a | 3,994 / 23 | 4,061 / 23 | 0 |
| sirius-b | 3,923 / 23 | 3,876 / 23 | 0 |

~90–200 matches per visit per endpoint-role — consistent with chance
coincidence of field sources under a 10" screen along a ~6' arc
(α Cen's near-Galactic-plane corridor roughly 2× the others). Single-
visit proximity carries no signal, as expected.

## Recurrence analysis (`recurrence_report.json`)

1. **z-bin × visit occupancy saturates.** At this background density
   every log-z bin has ≥1 chance match in nearly every visit (mean bin
   coverage ≈ 1.0); the statistic has no discriminating power and the
   apparent preference for the highest-z bin is pure locus-convergence
   geometry.
2. **Fixed-z point filter** (matches within 3" of the interpolated
   locus point per trial z, 120-point log grid): mean chance support
   2.5–4.1 visits; best peaks 7–13 of 23–26 visits with uncorrected
   binomial tails down to 6×10⁻⁷ (≈10⁻³ after ~1,700 correlated
   trials).
3. **All strong peaks resolved as static background stars.** The
   dominant support clusters have (a) tightly clustered W1 magnitudes
   (one physical source), (b) sky rms 0.17–0.28", and (c) detections
   confined to a single ~5-day day-of-year window each year — i.e. one
   parallax phase only. WISE visits alternate between opposite parallax
   phases; a static star near the fixed-z point recurs only at the
   returning phase (≈ half the visits — exactly the observed 11–13 of
   23–26), while a genuine relay at fixed z would be carried to BOTH
   phases by the locus. No track-consistent candidate survives.

**Result: a defensible layer-1 null** — and a validated cheap veto (the
parallax-phase test) to formalize in the next stage's model comparison.

## Context tables (snapshotted per corridor)

CatWISE2020: 3,788–5,964 rows per corridor; reject table: 648–764 rows
where sampled; SSO association list: 171–4,748 rows. Held for candidate
vetting; not yet analyzed.

## Caveats

- Catalog absence is not a null result (plan §3.4): sources below
  single-exposure thresholds, split entries, or reject-only sources are
  invisible to this layer. Stage 2 (forced photometry + trajectory
  coaddition) addresses exactly that.
- Binomial tails assume per-visit independence and a flat per-z chance
  rate; both are approximations, quoted uncorrected, and used only for
  triage.
- Screen radius 10" and point radius 3" are declared thresholds
  (recorded in the run config), not optimized.

## Next

Stage 2 per plan §4.6: forced photometry at predicted per-epoch
positions on the L1b images (cutouts along usable tracks), joint
continuous-z + bounded-residual-motion track fitting, and explicit
model comparison (SGL track vs static star vs ordinary PM source vs
Keplerian) with the parallax-phase structure built into the
likelihood — then injection calibration (§4.7) to turn all of this
into completeness-qualified constraints.
