# SPHEREx pilot v1.0 — stage summaries (2026-08-20)

Hypotheses `spherex-hypotheses-v1.0`; registry `pilot_wise_2026.yaml`
v1.5; corridors `lalande` (Lalande 21185), `gj687` (GJ 687), `sigmadra`
(σ Dra). All products regenerable from the tracked scripts
(`scripts/run_pilot.sh`); bulky data under `runs/spherex/` (3.2 GB,
8,505 slim cutouts).

## Coarse discovery (`coarse_v1`)

IRSA TAP `spherex.obscore ⋈ plane ⋈ artifact`, public Level-2 images,
MJD 60700–61400. 13,488 detector-exposures discovered (3,718 all-sky +
9,770 deep-field), epochs 2025-05-24 → 2026-08-11. 26,976 coarse
evaluations (polygon footprint, 25.5" pad), 16,978 hits.

| endpoint | exposures | rx hits D1…D6 | tx hits |
|---|---|---|---|
| lalande-21185 | 653 | 54 / 50 / 74 / 55 / 47 / 81 | 56 / 50 / 76 / 55 / 47 / 80 |
| gj-687 | 11,670 (9,770 deep) | 989 / 252 / 2,539 / 962 / 240 / 2,501 | 988 / 254 / 2,539 / 962 / 243 / 2,502 |
| sigma-dra | 1,165 | 95 / 118 / 113 / 94 / 116 / 109 | 95 / 118 / 109 / 94 / 115 / 106 |

Footprint validation: 349 polygon hits at the Lalande antipode = TAP's
own `CONTAINS(POINT, s_region)` count.

## Precise pass (`precise_v1`)

One slim cutout per hit exposure (100 px = 615", S3 byte-range reads,
~9 requests / 8.7 MB fetched / 260 kB kept, 1.3 exposures/s with 10
threads; **0 of 8,505 products missing**). SIP WCS + FLAGS template
708343, tolerance 1", seed step 10". 16,978 evaluations: 15,948
usable, 813 partial, 217 unusable (no usable-pixel crossing).

| endpoint / role | usable | covered z [AU] | epochs (MJD) | visit clusters |
|---|---|---|---|---|
| lalande-21185 rx / tx | 345 / 340 | 550–10,000 | 60809–61215 | 6 |
| gj-687 rx / tx | 7,397 / 7,412 | 550–10,000 | 60790–61227 | 4 |
| sigma-dra rx / tx | 634 / 633 | 550–10,000 | 60837–61232 | 7 |

Every corridor covers the full relay-distance range at every visit
(no chip-gap geometry: a 3.5° field with 2,040 px and a 6' locus).

## Flux-scale positive control (`control_v1`)

129 isolated 2MASS AAA stars (11 ≤ J ≤ 14.5) with CatWISE W1/W2 and
Gaia RP, 13,657 measurements through the track estimator on 1,170
cutouts. Δm = −2.5 log(measured/predicted SED):

| detector | N | median Δm | robust scatter |
|---|---|---|---|
| D1 | 2,318 | −0.20 | 0.08 |
| D2 | 1,979 | 0.00 | 0.08 |
| D3 | 5,782 | +0.03 | 0.10 |
| D4 | 1,637 | +0.10 | 0.08 |
| D5 | 928 | +0.07 | 0.07 |
| D6 | 1,013 | +0.17 | 0.12 |

All within the ±0.2 mag criterion (D1 and D6 sit where the catalogue
SED is extrapolated/interpolated across the RP–J and Ks–W1 gaps). The
0.05 mag all-band median is applied to every depth. **Before the
sub-pixel-phase estimator the same test gave −0.35 to −0.50 mag in
every detector** — the under-sampled PSF (5.3" on 6.15" pixels) loses
up to 0.5 mag when a pixel-centred matched-filter map is interpolated
at a half-pixel offset; evaluating the filter at 2 × 2 kernel phases
removed it.

## Static-sky template (`calib_v1/templates`)

Per corridor × detector, 3" tangent-plane grid, weighted linear fit in
λ across epochs, two passes with 3σ clipping. Defined on 33–56 % of
grid nodes (the rest are outside the cutout union). Median reduced χ²
0.7–1.2 in the all-sky corridors; 0.15–0.3 in the deep field, where
the per-cutout MAD noise (used for the weights) is dominated by static
confusion that the template removes. Median per-cutout variance
rescale: 1.06–1.3 (all-sky) vs 11–18 (deep field) relative to the
pipeline VARIANCE plane.

Effect on the stack: without the template, the 8 offset-control
thresholds were T = 15–70 (static stars below the single-epoch clip
add coherently; stack depth 16–17 AB, shallower than one exposure);
with it, T = 1.6–5.9 and m90 = 18.4–20.8 AB.

## Layer-1 track screen (`screen_v1`)

Single-epoch S_e > 5 peaks within 10" of the track (any z), clustered
at 6", CatWISE2020 context snapshotted:

| corridor | usable cutouts | peaks | within 10" | clusters | static (≥3 epochs or CatWISE) | unresolved |
|---|---|---|---|---|---|---|
| lalande | 347 | 15,710 | 539 | 49 | 35 | 14 |
| gj687 | 7,416 | 1,167,466 | 47,353 | 985 | 796 | 189 |
| sigmadra | 634 | 79,303 | 2,954 | 256 | 201 | 55 |

Every unresolved cluster is a 1–3-detection singleton at one or two
epochs (cosmic-ray / outlier residue at the 5σ single-epoch level);
no cluster follows the parallax track across visits.

## Injection calibration (`calib_v1`, AnalysisRun `run-55623bb1127a`)

96-node 1/z grid × 3 × 3 µ grid, 8 offset controls (±40/60/80" RA,
±50" Dec), epoch floor 5, clip |S_e| ≤ 5, weight cap 20×, 32
duty-randomised analytic injections per node, 90 % recovery, exposure
PSF, per-epoch FWHM mismatch factors. 36 (endpoint, role, detector)
searches; T = 1.6–5.9.

**8 exceedances → 8 Candidate records, all vetoed, 0 retained.**
Seven fail the parallax-phase split; σ Dra rx D4 (S = 4.3, T = 3.7,
z ≈ 2,170 AU) passes the phase split marginally (3.8 / 2.4) but the Tx
stack at the same cell — the Rx and Tx loci coincide to < 1" for
linear endpoints — gives S = 0.9 and no other detector shows anything
(|S| < 1.3): the role-coincidence veto added at this stage removes
it. Eight exceedances against an expectation of ≈ 4 (36 searches,
FAR < 1/8 each) is within the Poisson range but suggests the controls
under-sample the static residual; the scale-up should raise the
control count.

288 Constraint records (every z interval constrainable). Median m90
(duty ≥ 0.5, 90 % recovery, AB mag, 0.05 mag control offset applied):

| endpoint | D1 0.74–1.11 µm | D2 1.09–1.64 | D3 1.63–2.41 | D4 2.43–3.81 | D5 3.80–4.41 | D6 4.41–4.99 |
|---|---|---|---|---|---|---|
| lalande-21185 (rx / tx) | 19.5 / 19.6 | 19.3 / 19.8 | 20.3 / 20.0 | 20.5 / 20.1 | 19.5 / 19.6 | 19.3 / 19.0 |
| gj-687 | 19.4 / 20.2 | 19.7 / 19.1 | 19.7 / 19.9 | 20.3 / 20.8 | 19.5 / 19.5 | 20.2 / 20.7 |
| sigma-dra | 19.8 / 19.9 | 19.6 / 19.5 | 20.3 / 20.1 | 19.9 / 19.9 | 20.0 / 19.9 | 19.3 / 18.9 |

Physical conversions at the inner (550–720 AU) interval: reflected
light (albedo 0.1, solar colour) D ≳ 1.3–2.4 × 10⁵ km in D1–D4;
blackbody emitters at 700 K D ≳ 11–17 km (D4–D6), at 1000 K D ≳ 4–6 km
(D4–D6), at 400 K D ≳ 40–90 km (D5–D6). At 10,000 AU every limit
scales by ×(10,000/600)² ≈ 280.

## Findings that change the scale-up

1. **Under-sampled PSF:** pixel-centred matched-filter maps cost
   0.35–0.5 mag; sub-pixel-phase evaluation is mandatory (done).
2. **No difference images:** the static-sky template is the SPHEREx
   substitute and is worth 3 mag of stack depth. v2 should weight
   epochs by the template's own residual variance (the per-cutout MAD
   is confusion-dominated in the deep field, rescale 11–18×).
3. **Deep-field seasons are detector-dependent:** at the GJ 687
   corridor D1/D4 epochs are almost all one parallax phase and D3/D6
   the other (phase-1 counts 4–15 in D3/D6), so the per-detector phase
   split is weak there; a cross-detector phase test (D1/D4 vs D3/D6 at
   the same cell) should be added.
4. **Role coincidence** is a free, strong veto for linear endpoints
   (added at calibration, 2026-08-20).
5. **Control count:** 8 exceedances vs ≈ 4 expected; raise to 16
   offsets or use the empirical null distribution across cells.
