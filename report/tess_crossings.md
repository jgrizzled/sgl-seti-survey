# TESS beam-crossings survey (Pipeline B) — report

Sixth crossings survey and first item of the archive-expansion queue
(plan §5.8): a search of TESS Ecliptic-sector FFI cutouts for optical
emission during Earth's crossings of hypothesized Sun–star relay beam
axes — the first archive able to resolve complete ingress→egress
crossing light curves, and the first opening of the pulse-period cell
(200 s → window length). Executed 2026-08-24. **No candidates.**

## 1. Construction and provenance chain

1. **Spacecraft-frame crossing list** `crossings/tess_v1`
   (`xng-a943f0f3dbe4`): 3,390 events 2018-07 → 2026-08 over the
   88-endpoint registry, computed with the TESS spacecraft observer
   (JPL Horizons −95 SSB vectors, 6 h sampling, sha-pinned). The
   spacecraft frame is mandatory here: at grazing impact parameters the
   Earth-center list is off by up to 0.33 R☉ / 3.2 h (HEO apogee
   0.54 R☉).
2. **Coverage gate** (`surveys/tess-crossings/results/coverage_gate_v1.*`):
   TESScut sector lookups × HEASARC sector dates. The first pass was
   corrupted by a HEASARC sector-table year typo (s046 end dated 2022
   for 2021, creating a bogus 393-d "sector"); the script now clamps
   implausible sector spans to start + 28 d. Corrected coverage: two
   grazing-rung events with complete ingress→egress light curves —
   **wolf-359 b = 0.55 R☉ (s42, 600 s)** and **teegarden b = 1.04 R☉
   (s91, 200 s)** — plus two complete near-grazing channel-A windows
   (gj-1276 b = 0.97 R☉ s42, van-maanen b = 0.52 R☉ s43). The era's
   two deepest grazes (van-maanen 0.078 R☉, gj-1276 0.49 R☉, 2022
   spring) fell between sectors and remain open.
3. **Hypothesis freeze v1.0** (`hypotheses.md` + `notes/tesscut_recon.md`):
   channel A constraint-only (≥9 d windows in 27 d sectors defeat the
   8-offset temporal-control family — 4th appearance of the theorem);
   **channel B is the discovery channel** (spatial ring controls are
   window-length independent); two temporal statistics — chord LSQ
   amplitude S_c and per-cadence pulse maximum S_p; TESS band contains
   752 nm (first archive over the PS1-anomaly line wavelength, at far
   shallower depth); SkyBoT census mandatory for pulse exceedances
   (21″ pixels in the asteroid stream).
4. **Coverage refinement + TIC cut** (`results/coverage_refined_v1.*`,
   `tic_cut_v1.json`): 8 TESScut cubes (sha manifest), spacecraft-UTC
   window intersection, real cadence counts (grazing windows hold
   81–405 clean cadences); saturation excludes only gj-908 (T = 7.10;
   survey targets T ≈ 10.6–14.7 vs saturation T ≈ 6.8).
5. **Threshold freeze v1.0 → v1.1**
   (`configs/threshold_freeze_v1.json`, sha256:5bc0197e…): the
   science-array usability screen (finding T1) caught TESScut serving
   collateral non-science CCD pixels — 3/8 cutouts off-science,
   including the v1.0 dev unit (ross-128 B) → re-freeze v1.1: **6 B
   units (wolf-359 × 3 rungs + teegarden × 3 rungs, all confirmatory)
   × 2 statistics = 12 trials, 1.3 expected control crossings**; 8
   ring controls at 126/189/252″; exceedance S > max(T, 0); machinery
   validation reassigned to the zero-trial channel-A reference rows.
   First survey in the programme where no covered narrow-rung unit is
   lost to saturation or unit gates.
6. **Dev-stage validation** (`results/dev_validation_v1.*`): teegarden
   s71 star recovered to +0.07 mag of its TIC magnitude through the
   identical kernel chain; wolf-359 B field ZP scatter 0.066 (pass);
   blended channel-A series systematics-dominated (k up to 292 —
   finding T2, validating the zero-trial A design); van-maanen A
   self-calibration −0.31 mag at the FWHM-fit ceiling → declared
   ±0.3 mag scale caveat on its reference depths.
7. **Confirmatory search** (`results/confirmatory_v1.*`,
   `adjudication_v1.json`): §2.
8. **Completeness** (`results/completeness_v1.json`): exact
   stamp-response injections with the **SPOC per-camera/CCD PRF**
   (archive.stsci.edu `prf_fitsfiles`, nearest grid point to each
   cutout's CCD position, sha-recorded) through the identical kernel
   chain; 200 draws per unit per temporal model (freeze ≥ 100);
   magnitude grid T 10–18 (18.0 = grid-censored); temporal models
   d = 1 chord and d = 0.1 boxcar on the frozen 6-period log grids;
   recovery against each unit's actual frozen thresholds. Finding C1:
   the dev-stage ZP helper used a hardcoded 1.5 px kernel where the
   freeze requires the identical per-cube kernel — the flux scale was
   re-measured through the unit kernel with injection responses
   normalized at the calibration reference (flux-scale-only
   correction; no search statistic touched). Wolf-359 unit-kernel ZP
   20.80 (scatter 0.112, pass); teegarden s91 ZP 20.03 (scatter
   0.493, FAILS the 0.2 gate → teegarden depths carry a declared
   ±0.5 mag scale caveat).

## 2. Results

**6 confirmatory B units × 2 statistics: the pulse statistic — the
survey's genuinely new cell — is null in all six units (0
exceedances).** No pulse of ≥ 1 cadence (200–600 s) stands above the
ring-control maxima anywhere in the two fully resolved grazing
crossings or the wide windows — the first constraint of its kind.

The chord statistic is systematics-dominated (sector-scale
scattered-light drift; S and T in the tens; one star-contaminated
control at 927 renders wolf-359 0.1 AU deeply insensitive but null):
3 chord exceedances vs 1.3 expected (P(≥3 | 1.3) ≈ 0.14), all
adjudicated under the frozen ladder:

| unit | b (R☉) | S_c | T_c | S_p | T_p | disposition |
|---|---|---|---|---|---|---|
| wolf-359 1.2 R☉ | 0.55 | 12.5 | 15.3 | 2.5 | 4.5 | searched_null |
| wolf-359 2.5 R☉ | 0.55 | 21.4 | 20.9 | 2.9 | 5.5 | adjudicated control-crossing |
| wolf-359 0.1 AU | 0.55 | −19.4 | 927.6 | 3.9 | 137.2 | searched_null (insensitive) |
| teegarden 1.2 R☉ | 1.04 | 17.2 | 17.3 | 3.0 | 4.8 | searched_null |
| teegarden 2.5 R☉ | 1.04 | 48.3 | 24.9 | 5.2 | 7.8 | **vetoed** (chord-shape) |
| teegarden 0.1 AU | 1.04 | 79.7 | 52.1 | 6.9 | 7.9 | **retained-ambiguous** |

- **teegarden 2.5 R☉: vetoed** by the frozen chord-shape/timing test —
  split-half amplitude ratio 0.31 (a monotonic ramp, not a chord) and
  background anti-correlation −0.785 (the background-subtraction
  residual signature); the window sits in the sector's final day.
- **wolf-359 2.5 R☉: adjudicated control-crossing** — margin +0.46 is
  2 % of the statistic's scale, and the nested-window cross-check
  contradicts a real source (the 0.1 AU superset has S = −19.4; the
  1.2 R☉ subset does not exceed; a relay at b = 0.55 R☉ is in-beam
  for all three rungs simultaneously).
- **teegarden 0.1 AU: retained-ambiguous, non-promotable** — the shape
  test is weakened by window truncation at the sector end (no egress
  to test), and it shares the final-day cadences with its vetoed
  sibling; per the v2 discipline it is not dismissed by discretion.
  Promotion requires independent recurrence; teegarden's next covered
  window needs a future ecliptic sector — **the designated follow-up**
  (the programme's second retained-ambiguous row, after PS1's
  gj-1276).

Channel A produced its 4 zero-trial reference light curves (teegarden
s71: 3,424 + 460 in-window cadences; gj-1276 1,084; van-maanen 1,475)
— the only resolved on-star crossing light curves at b ≈ 0.5–1.2 R☉
in any archive; 3 rows are off-science-array (nominal-covered /
unusable) per finding T1.

## 3. Completeness (90 %-recovery depths, TESS T mag)

SPOC-PRF stamp-response injections; chord model → S_c threshold,
pulse model → S_p threshold (full 2 × 2 recovery matrices in
`completeness_v1.json`). Pulse recovery plateaus at 0.93–0.98 even
for bright injections: a d = 0.1 pulse train can fall entirely in
straylight-masked cadences — real duty-cycle incompleteness, included
in the curves.

| unit | m90 chord | m90 pulse | notes |
|---|---|---|---|
| wolf-359 1.2 R☉ | ≥ 18.0 (censored) | 16.5 | |
| wolf-359 2.5 R☉ | ≥ 18.0 (censored) | 16.3 | |
| wolf-359 0.1 AU | 15.2 | 12.1 | star-contaminated control |
| teegarden 1.2 R☉ | 17.5 | 14.8 | ±0.5 mag scale caveat |
| teegarden 2.5 R☉ | 17.6 | 14.7 | ±0.5 mag scale caveat |
| teegarden 0.1 AU | ≥ 18.0 (censored) | 14.9 | ±0.5 mag scale caveat |
| A teegarden s71 (reference) | ≥ 18.0 | 14.1 | S = 5 reference, threshold-free |
| A van-maanen s43 (reference) | ≥ 18.0 | 13.2 | S = 5 reference; ±0.3 mag |

## 4. Physical interpretation (TESS band 600–1000 nm, in-band)

**Downlink pre-lens, resolved grazing crossings (the survey's
distinctive constraint).** A relay beaming through the solar-grazing
cone toward the target star during the covered windows is excluded at
90 % above (persistent d = 1 / chord):

- wolf-359 (b = 0.55 R☉, s42): ≲ **560 W** through the 1.2 R☉ cone
  and ≲ 2.4 kW through 2.5 R☉ (both grid-censored — true limits are
  deeper); 0.1 AU rung ~2.3 MW (insensitive unit).
- teegarden (b = 1.04 R☉, s91): ~**840 W** (1.2 R☉), ~3.7 kW
  (2.5 R☉), ≲ 180 kW (0.1 AU, censored); ×~1.6 scale uncertainty from
  the failed ZP gate.

**Pulse-period cell (opened for the first time).** Any single pulse of
≥ 1 cadence during the resolved crossings is excluded at 90 % above
~2.2 kW / 1.3 MJ per pulse (wolf-359, 1.2 R☉ cone, 600 s) and ~11 kW /
2.2 MJ (teegarden, 200 s); duty-cycle miss probability 2–7 % is in the
curves. Sub-cadence pulses scale as ×(cadence/duration); pulse trains
timed to avoid clean cadences are declared unconstrained.

**Uplink reference (channel A, threshold-free).** On-star pulses at
the two resolved near-grazing windows would have stood at the S = 5
reference above T ≈ 14.1 (teegarden) / 13.2 (van-maanen) — a 10-m
diffraction-limited transmitter in the target system at MW-class
pulse power (~1–3 MW at 3.8–4.3 pc). Reference depths, not
calibrated exclusions: zero-trial rows have no controls.

**Not constrained:** the 2022-spring deepest grazes (between sectors);
transmitters scheduled to avoid Earth-crossing windows; sub-cadence
pulse structure beyond ×d scaling; the A 1.0 AU rung
(programme-wide constraint-only); depths beyond T ≈ 17–18 (four
orders of magnitude shallower than ZTF's per-epoch limits — the value
here is temporal structure, not depth).

## 5. Design lessons

1. **The undetrended chord statistic saturates its error budget with
   sector-scale scattered-light drift** — k rescales variance, not
   low-frequency structure; the 1/9 control-crossing budget is
   approximate at best for it. Any TESS-crossings v2 needs a
   differential/detrending layer under the chord filter. The pulse
   statistic, differential by construction, behaved exactly to
   budget.
2. TESScut serves collateral (non-science) CCD pixels with
   aperture = 1 and ~zero flux — a usability screen on science-array
   columns (45–2092) must run before any freeze counts a cutout
   usable (finding T1; cost this survey its only dev unit).
3. Flux calibration through a mismatched kernel is not the freeze's
   "identical kernel": measure the star ZP through the per-cube
   fitted kernel and normalize injection responses at the calibration
   reference, or depths inherit a hidden throughput factor
   (finding C1).
4. Sector-boundary truncation is the archive's structural weakness:
   both the vetoed and the retained-ambiguous exceedances live in a
   sector's final day, and the era's two deepest grazes fell between
   sectors entirely.

## 6. Products and ledger states

`surveys/tess-crossings/{hypotheses.md, thresholds.md, configs/*,
results/*, notes/tesscut_recon.md}`; crossing list
`crossings/tess_v1/`; cubes + PRFs + snapshots under
`runs/tess-crossings/` (8 cubes ~420 MB, 4 SPOC PRF files,
sha manifests); scripts under `surveys/tess-crossings/scripts/`.

For the next covered-window ledger refresh (plan §5.7): 6 searched B
unit-rows (4 searched_null, 1 vetoed, 1 retained_ambiguous — the
follow-up row), 4 channel-A constraint-only reference rows, 3
off_science_array rows (nominal-covered/unusable), and the A 1.0 AU
rung's 292 gate-level coverage events.
