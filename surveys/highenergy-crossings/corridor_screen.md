# Corridor catalogue screen v1.0 (Pipeline A step 1; high-energy archives)

Frozen 2026-09-07 with hypotheses.md v1.0 (D11). A catalogue-level
persistent-source screen on the 88 anti-star corridors of the universal
target list — the relay corridor is the anti-star direction at
550–10,000 AU, whose annual parallax is 375″ (550 AU) to 20.6″
(10,000 AU) — using the catalogues and services the recon reached. No
flux search, no image, no statistic: the outputs are ledger rows
`catalogue_screen_only` and a per-hit **SGL-track test**.

## Inputs (annotation-class, snapshotted at the recon)

- eRODat cones DR1_Main (eRASS1, 2019-12-12 → 2020-06-11) and DR2_Main
  (eRASS:3 cumulative, → 2021-06-16) at 7′; eRODat upper limits
  (DR1 eRASS1 bands 021/024; DR2 eRASS:3 band 024) at every corridor.
- HEASARC 5XMM-DR15 per-detection rows (`xmmssc`: srcid, obsid,
  position, error, detection epoch), `xmmstack`, 2SXPS, LSXPS
  (first/last observation epochs), eRASS1 main/hard, BAT-157m, CSC 2.1
  master sources, 4FGL-DR4 (30′).

## Track test (frozen rule)

A catalogued source within 7′ of a corridor centre is on the relay
grid only if its position follows the corridor parallax track. The
smallest parallax on the grid is 20.6″ (10,000 AU). Rule, per hit:

1. **Two-epoch test.** If the hit has positions at two epochs t₁, t₂
   with predicted relay displacement Δ(z = 10,000 AU) ≥ 15″ between
   them (Δ = 206265″/z × |P⊥(E(t₁) − E(t₂))|, E the Earth's
   heliocentric position in AU, P⊥ the projection perpendicular to
   the corridor direction), and the two positions agree to < 10″ →
   **`fixed_sky`** (excluded as a relay at every z ≤ 10,000 AU… and
   beyond, up to 206265″ × |ΔP⊥|/10″).
2. **Stack test.** If the hit is a compact detection in a catalogue
   built from epochs spanning ≥ 180 d that straddle opposite parallax
   phases (eRASS:3 = three scans 6 months apart; LSXPS/2SXPS with
   first→last ≥ 180 d), then a relay at z ≤ 10,000 AU would have moved
   ≥ 2 × 20.6″ × |sin| between the contributing epochs — larger than
   the catalogue PSF (eROSITA 30″ HEW; XRT 18″ HEW) — and could not
   appear as a single compact source: **`fixed_sky`** (inference
   labelled `stack`, weaker than rule 1).
3. Otherwise **`single_epoch_unresolved`** — recorded for the yearly
   refresh.

Rule 1 uses 5XMM per-detection epochs and the eRASS1-vs-eRASS:3 pair
(eRASS1 = scan 1 alone; the eRASS:3 position is the scan 1–3 mean, so
a moving source would separate from its eRASS1 position by ≥ ⅔ of the
scan-1→scan-2 displacement). Rule 2 uses LSXPS `time`/`end_time` and
the eRASS:3 `ext` = 0 compactness.

## Upper-limit ledger

For every DE-sky corridor the eRODat DR1 eRASS1 (0.2–2.3 keV, band 024;
0.2–0.5 keV band 021) and DR2 eRASS:3 (band 024) 3σ upper limits are
ledger rows `upper_limit_only` with the service's exposure and
background; for the 52 eastern (Russian) corridors, `no_public_data`.
Persistent-source power through the corridor is not computed here
(the Pipeline-A depth statement needs a source model; the flux limits
are the record).

## Outputs

`results/corridor_screen_v1.json` / `.md`: per corridor — DE sky
yes/no, upper limits, hits per catalogue with separation and the
track-test disposition; census by disposition.
