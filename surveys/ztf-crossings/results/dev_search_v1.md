# Dev-target search v1 — results and findings

Run 2026-08-23 under threshold freeze v1.0 (dev split only: 6 channel-A
targets, ross-154 for channel B wide rung). Products: 2,389/2,660
epochs fetched complete (sci+msk±diff cutouts, 202 MB,
`runs/ztf-crossings/products/dev/`), 0 hard failures; statistics in
`dev_search_v1.json`.

## Channel B — runs clean per the freeze

All 12 ross-154 event-band rows (10 search units, 2 single-epoch)
computed with the full 8-control ring and the max-over-z family:
**no exceedances**; R spans −0.05…0.73. S is O(1) as expected for the
normalized statistic on empty fields; T = 1.1–8.8. The track machinery,
z-max, ring controls, and single-epoch classification all behaved.
Channel B needs no amendment.

## Channel A — two frozen constructions fail (the dev split did its job)

**Finding A1 — temporal controls undefined for every unit (0/8).**
No dev unit reaches 8 valid pseudo-window offsets (best 4; ez-aqr and
gj-1111 get 2). Structural: all A dev units are on the 1.0-AU rung,
whose windows span up to ±58 d of each semiannual crossing; every
frozen offset (±23…97, redraws ±113/127) lands a shifted window inside
the unit's own or the neighboring crossing's window under the
either-rung exclusion. Note the annual phase-lock is intrinsic: the
same-side crossing recurs at ~1 sidereal year, so *any*
parallax-phase-matched control epoch is itself in-window — controls
must either accept phase mismatch or the statistic must remove
phase-dependent systematics first.

**Finding A2 — the blended statistic is systematics-dominated.**
Diagnostic S values (partial controls, no thresholds implied) are
10²–10⁴ in |S| and persist in controls: ez-aqr g S ≈ −5.5×10⁴ with
controls −1.3/−2.1×10⁴; gj-11068 S ≈ +96 (g), +596 (r) with controls
+41…+543 and spurious partial-R ≈ 1.1. Interpretation: ZTF references
are multi-year stacks, so a high-PM star sits displaced from its
reference position and the difference image carries a growing PM dipole
(positive lobe at the propagated position — gj-11068), while
brighter/marginal stars leave large negative PSF-mismatch residuals
(ez-aqr, gj-1087). The frozen raw forced-flux excess never approaches
the noise floor; a per-target systematics model (PM-dipole + static
residual fitted on off-window epochs) must be subtracted before the
window statistic is formed.

## Consequences

1. Channel B wide-rung is ready for the confirmatory run as frozen.
2. Channel A requires a freeze amendment (v1.1) before any confirmatory
   search; candidate elements, to be decided and frozen *before*
   touching confirmatory targets:
   (a) subtract a per-target PM-dipole + static template fitted on
   off-window epochs only; (b) with systematics removed, relax the
   pseudo-window exclusion to same-rung-only (a wide-beam signal leaking
   into controls only makes them conservative); (c) define R only for
   T > 0, with T ≤ 0 units reported as null-model failures.
3. The A 0.1-AU rung has no dev units (all dev targets are off-ecliptic)
   — the amendment cannot be tuned on that rung; it inherits.
4. gj-518 correctly produced no search unit (single-window,
   constraint-only), and 271 epochs lacked complete products (archive
   404s), counted and excluded, not silently dropped.

## Amendment v1.1 re-run (same day; `dev_search_v11.json`)

Template fit + same-rung offsets + S > max(T, 0) rule applied to the 8
channel-A dev units (249 additional epochs fetched; template fitted on
≤150 off-window epochs per unit).

**Offsets:** same-rung-only still leaves every dev unit at 2–4 valid
offsets — the 1.0-AU rung's own windows block the short offsets for
small-b events. Under v1.1's disposition all 8 units are
**constraint-only**; no exceedances (rule applied: none).

**Template performance:** off-window RMS drops 2–4.5× everywhere
(ez-aqr 43k→10k; gj-1087 g 16k→8.8k; gj-11068 g 227→94), i.e. the
PM-drift + parallax-factor + seeing model captures the bulk of the
deterministic residual, but bright/marginal targets retain
10²–10³-scale structure: residual scatter is far above the
matched-filter (background) variance, and pseudo-windows at ±113/127 d
retain phase-dependent offsets not shared by the real windows
(gj-1087 g: S = +1202 vs controls −172…−399). Two residual causes:
(1) per-epoch variance is background noise, not the empirical
bright-star residual scatter — the WISE v1 lesson ("confusion, not
noise, sets depth") in a new guise; (2) nonlinear seeing-dependence of
bright-star subtraction residuals survives the linear seeing term.

**Faint limit behaves:** gj-11068 g (g ≈ 16.9), the faintest unit, lands
at S = −1.41 with controls {−0.26…+4.06}, median 1.49 — order unity,
narrowly failing the |median control| ≤ 1 gate. gj-9193 r fell below
the 30-epoch fit floor (29 kept) and was correctly gated out.

**Verdict:** the v1.1 machinery works as designed — the gate correctly
passes nothing in the bright regime and the construction converges
toward validity in the faint regime. Channel A's blended search is
viable only for faint (g ≳ 15) targets even after amendment. Candidate
v1.2 refinement (user decision, not applied): rescale per-epoch
variances by the off-window residual variance per (target, band) — the
established v1 precedent — which would properly standardize S and
likely qualify the faintest units; the bright majority of channel A
remains constraint-only regardless.

## Amendment v1.2 re-run (same day; `dev_search_v12.json`)

Empirical variance rescale applied on top of v1.1 (no new data). The
k factors quantify how far the matched-filter background variance
understated the real residual scatter: k ≈ 5.8 for the faintest unit
(gj-11068 g) up to k ≈ 8.9×10⁴ (gj-1087 r). Standardized S is now
O(1–36) everywhere the template fits, and the |median control| ≤ 1 gate
is meaningful: **gj-11068 g passes it outright** (S = −0.59, controls
−0.11…+1.69, median 0.62); gj-11068 r misses narrowly (median 1.10).

All 8 units remain **constraint-only** — the binding constraint is now
purely the 2–4 valid pseudo-window offsets, i.e. the 1.0-AU rung's
window geometry, as v1.1 anticipated. No exceedances.

One new insight from the standardized numbers: bright targets show
in-window S systematically above their controls (ez-aqr g +35.9 vs ~3;
gj-9193 g +9.3 vs ≤2; gj-1087 +5.6/+7.1 vs negative controls). This is
not signal: for the 1.0-AU rung the crossing window ≈ the observing
season (the star is observable at night precisely when it is near
opposition, which is when Earth crosses the axis), so off-window fit
epochs cluster at the season edges and the template *extrapolates* into
the season core, leaving an airmass/seeing-correlated bias. This
selection coupling is intrinsic to the 1.0-AU rung and independently
justifies its constraint-only disposition. The 0.1-AU rung (windows
±≤5.75 d ≪ season, off-window epochs bracketing each window, all 8
offsets valid) does not suffer this coupling — the confirmatory
channel-A search rests there, as the coverage stage already indicated.

**v1.2 verdict: adopted. Machinery final for confirmatory** (channel B
under v1.0; channel A under v1.0+v1.1+v1.2).
