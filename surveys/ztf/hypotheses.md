---
title: "ZTF survey — hypothesis freeze v2.0 (decision rule and geometry specification)"
status: "v2.0 — frozen 2026-08-22 before any v2 script ran on the ZTF confirmatory set; hash in configs/v2_freeze.json"
date: 2026-08-22
---

# ZTF hypothesis freeze v2.0

The carried-forward physical-source hypothesis is reproduced directly
in §0. The decision rule, null model, completeness definition and
hold-out are those of the WISE freeze **v2.1**
(`surveys/wise/hypotheses.md` §3–§5 with its §8 amendment), applied
through the survey-agnostic engine `sglsurvey/` with the ZTF
bindings in `profile.py`.

> **Documentation restoration (2026-08-24).** This section restores the
> still-active v1 physics text after the v1 file was retired. It changes
> no search parameter, data product, candidate decision or result. The
> `hypotheses_hash` in `configs/v2_freeze.json` remains the historical
> pre-run hash and is not a hash of this post-run documentation copy.

## 0. Physical-source hypothesis carried forward from v1

1. **Targets and endpoint models.** This archive searches the
   ZTF-visible subset of registry v1.5: 69 endpoint hypotheses in 62
   corridors. Each endpoint is the registry-defined stellar component
   or explicitly modelled photocentre with fixed target-state provenance
   and covariance; barycentric and planetary endpoints are separate,
   out-of-scope hypotheses.
2. **Relay and role.** The source is a compact, actively station-kept
   artifact on or near the Sun's focal line for the endpoint. Receive
   (`Rx`) and transmit (`Tx`) geometries are searched separately.
3. **Distance and motion.** Heliocentric distance is 550–10,000 AU with
   a log-uniform physical prior. Residual motion is bounded
   component-wise by |µ*α| ≤ 1″/yr and |µ*δ| ≤ 1″/yr and is fitted with
   distance. The source is unresolved at ZTF's roughly 2″ seeing;
   extended or visibly trailed sources are different cells.
4. **Observer and timing.** The observer is the Palomar P48 site
   (longitude −116.8650°, latitude +33.3563°, height 1712 m). Loci are
   evaluated at exposure midpoint; topocentric displacement and motion
   during a 30 s exposure are negligible relative to the searched PSF
   but remain part of the geometry provenance.
5. **Bands and physical interpretation.** The searched bands are ZTF
   g, r and i (approximately 0.47, 0.64 and 0.79 µm), with constraints
   reported per filter. The primary optical interpretation is reflected
   sunlight. For a Lambertian sphere of geometric albedo p and diameter
   D at heliocentric distance z, observed near opposition with observer
   distance Δ ≈ z, the full-phase flux ratio is approximately

       F_object / F_sun ≈ p Φ(α) (D / 2 AU)^2 / (z^2 Δ^2),

   where z and Δ are in AU and Φ(α) is the phase function. A
   reflected-light limit therefore constrains p^(1/2) D only after an
   albedo and phase law are stated. The same flux search is also
   sensitive to self-luminous optical emission such as a beacon,
   sufficiently broad or recurrent leakage, or a component hotter than
   roughly 2000 K;
   those cases are reported as flux, not reflector size. The v2
   injections use a flat-Fν spectrum. A passive solar-equilibrium body
   at 550 AU is about 12 K and invisible in g/r/i, so this is not thermal
   coverage.
6. **Visibility in time.** Duty cycle is at least 0.5, instantiated by
   the persistent, exposure-flicker, visit and long-block injection
   families in item 5 below. Rare glints, low-duty pulses and emission
   deliberately directed away from Earth are separate cells.
7. **Outside this cell.** Relay swarms, off-axis infrastructure,
   inactive or dark relics, distance or residual motion outside the
   stated bounds, extended/trailed morphologies and duty < 0.5 are not
   constrained. A detection would establish an SGL-consistent moving
   point source, not its emission mechanism; a null makes no population
   or network-architecture inference.

The executed decision rule was:

1. **Statistic and null.** S(z, µ) over 192 × 5 × 5 nodes (uniform in
   1/z, 1.85″ spacing, T0 = MJD 59800), single-epoch clip |S_e| ≤ 5
   (the v1 layered-search rule), per-frame weight cap at 20 × the
   band's median frame weight, nodes with < 5 epochs undefined. Local
   threshold T from the 8 designated controls of a 48-offset ring
   (20/30/40″ × 16 angles), R = S_max/T, R̃ = R/q95(ring); the ring is
   the exchangeable null; the phase-coherence scramble and the
   trajectory randomisation are per-cell annotations; cells with a
   heavy-tailed (max > 2.5 q95) or radius-dependent ring are void.
2. **Candidate rule.** R̃ ≥ R̃_FWER, the 95th percentile of 10,000
   pseudo-experiment maxima over the family (α = 0.05), computed from
   the family's null before any real R̃ is examined. BH q-values and
   rank statements (48 controls, Wilson intervals) are reported.
3. **Calibrated veto.** The flux-consistent catalogued static-source
   / halo test against the ZTF DR objects snapshotted by v1 screening
   (≥ 3 good observations, Moffat-vs-Gaussian radial response, factor 2
   overall and in the dominant phase). The held-out-epoch prediction
   test (refit on epochs ≤ MJD 60554 = 2024-09-01, forced photometry on
   the last observing year) is an annotation with its measured rates,
   not a veto (WISE v2.1 §8). All other rules are annotations.
4. **Quality masks.** Primary: not `bad_quality` (infobits bit 25),
   seeing ≤ 4″; strict: seeing ≤ 2.5″, maglimit ≥ 19.5, Moon
   illumination ≤ 0.8; loose: v1 (all usable).
5. **Completeness.** 400 image-level injections per cell (100 per
   temporal model: persistent, flicker, visit-scale with a 5-day gap,
   long block), Moffat(β = 3) at the frame's seeing added to the v1
   hybrid search image before the matched filter (the science PSF in
   the difference image is an approximation to injecting before
   differencing; the asteroid control measures the chain on a real
   mover), continuous z (log-uniform), µ on the L∞ box, cross-track
   offsets from the propagated 99 % envelope, AB magnitude uniform on
   [m90_v1 − 2, m90_v1 + 2], flat-Fν spectrum. Threshold and
   final-candidate completeness with bootstrap 68 % intervals on 8
   gap-free reciprocal-distance intervals; worst-of-four coverage.
6. **Geometry.** sglseti Monte Carlo envelopes (N = 2,000, seed
   20260822, 2018.5 / 2022.0 / 2025.5, z = 550 / 10,000 AU); the
   cross-track dimension (3 or 5 offsets) where σ_xt,99 > 0.5 × 1.9″.
7. **Hold-out.** Random corridor split stratified by the WISE
   confusion class, seed 20260822, ≈ 30 % development / 70 %
   confirmatory; the confirmatory set is analysed once.
8. **Scope.** Targeted coverage of the 69-endpoint ZTF-visible subset
   of the frozen portfolio; no population inference.
