# Scientific Review of the WISE SGL Survey

**Review date:** 2026-08-21  
**Recommendation:** Major revision before the survey is presented as a statistically calibrated exclusion experiment

## Executive assessment

The WISE SGL survey is a scientifically motivated and unusually transparent exploratory search. Its strongest features are its explicit hypothesis cell, frozen target registry, extensive provenance, empirical use of local controls, and candid discussion of instrumental limitations. WISE is also a sensible archival data set for searching for persistent, warm, unresolved emitters near predicted solar gravitational lens trajectories.

The current implementation does not, however, support all of the quantitative conclusions in the [survey report](../../report/wise_survey_v3.md). In particular, the reported 90% recovery limits are not end-to-end completeness estimates, the eight-control threshold does not provide a useful survey-wide false-alarm rate, and most W3/W4 cells automatically veto any threshold-crossing source because of their expected cryogenic-mission cadence. The claimed 99% trajectory confidence region is recorded but not actually propagated through the production search.

The following conclusion is scientifically defensible in the current version:

> No compelling candidate remained after application of the survey's current heuristic review rules.

The following stronger conclusions are not yet established:

- A 90%-complete exclusion throughout the entire stated distance, motion, phase, and duty-cycle search cell.
- A controlled global false-alarm probability for the null result.
- End-to-end W3/W4 constraints on approximately 100 km, 300 K structures.
- A population-level conclusion about the prevalence of SGL relays or the Picky Network hypothesis.

## Scientifically strong elements

1. **Narrow, explicit scope.** The survey specifies distance, residual motion, unresolved morphology, duty cycle, and target-role hypotheses. It also states that cold, passive, faint, resolved, or strongly intermittent technologies can remain undetected.

2. **Documented target selection.** Endpoint selection, role assignment, and target provenance are inspectable rather than reconstructed after candidate examination.

3. **Recognition of non-Gaussian backgrounds.** Local offset controls are more appropriate than treating formal matched-filter signal-to-noise as Gaussian in crowded and artifact-rich WISE images.

4. **Protection against single-frame dominance.** The weight cap is a reasonable safeguard against individual exposures determining an entire track score.

5. **Transparent correction of the flux-scale error.** The revised normalization is mathematically consistent for the assumed unit-sum Gaussian kernel, and affected products were recalculated rather than corrected only in prose.

6. **Appropriately limited population interpretation.** The report does not convert the non-detection into an occurrence-rate limit for all SGL relays.

## Major methodological findings

### 1. The detection threshold does not control the survey-wide false-alarm rate

**Severity: Critical**

The calibration sets each threshold to the maximum score found among eight offset controls ([implementation](scripts/injection_calibrate.py#L139-L170)). If the target and eight controls were independent, identically distributed, and continuously scored, a noise-only target would exceed all controls with probability

\[
P(\mathrm{target} > \max(8\ \mathrm{controls})) = \frac{1}{9}.
\]

This is a per-search rank probability, not a small false-alarm probability. Across 704 searches, the null expectation is approximately 78.2 threshold crossings. The 70 observed crossings are therefore consistent with the behavior expected from the threshold itself; they do not constitute an excess over a noise-only expectation.

The situation is less secure than the ideal rank calculation because nearby controls can be spatially correlated and can occupy systematically different backgrounds from the nominal trajectory. No family-wise error rate or false-discovery rate is calibrated across endpoints, roles, bands, distances, motions, and signal phases.

Several veto rules were also introduced or refined after examining survey candidates. Those rules may be useful for exploratory classification, but their performance on the same data cannot be interpreted as a pre-registered confirmatory false-alarm test.

**Recommendations**

- Freeze the complete detection and adjudication rule before evaluating a confirmation data set.
- Generate many more matched null trajectories using spatial offsets, time scrambling, or trajectory randomization while preserving local coverage and contamination.
- Calibrate the distribution of the maximum statistic across the entire survey, or state and control an explicit family-wise error rate or false-discovery rate.
- Report the target's rank among exchangeable controls with an uncertainty statement instead of describing the present eight-control maximum as an empirical false-alarm rate.
- Reserve held-out epochs, targets, or sky positions for confirmatory validation after exploratory rule development.

### 2. The reported completeness does not include candidate vetoes

**Severity: Critical**

The injection calibration adds an analytic signal to a precomputed sampled time series and counts recovery when the resulting score exceeds the threshold ([injection loop](scripts/injection_calibrate.py#L172-L200)). It does not pass recovered injections through the rules in [candidate adjudication](scripts/adjudicate_exceedances.py#L199-L232). The reported curves therefore measure threshold-crossing sensitivity, not the probability that a source would survive the full survey and be classified as a candidate.

This distinction matters for phase, static-source proximity, bright-star halo, marginal-excess, and cryogenic-visit vetoes. A source satisfying the stated duty-cycle hypothesis can cross the threshold but still be rejected because it appears in only one survey phase or because its projected position is close to a cataloged static source. Such outcomes are counted as recovered in the calibration but rejected in the real search.

The W3/W4 cadence veto is especially consequential. The adjudicator rejects a cryogenic cell with fewer than 60 usable samples as a "single-visit cryo cell." An audit of the generated threshold report finds:

- W3: 162 of 176 cells have fewer than 60 samples.
- W4: 164 of 176 cells have fewer than 60 samples.
- Combined: 326 of 352 W3/W4 cells are subject to this automatic veto.

This is expected from the mission design. W3/W4 observations were confined to the cryogenic mission, and typical sky locations received roughly 12--13 exposures during a visit; no new W3/W4 imaging was acquired during NEOWISE. See the official [NEOWISE Explanatory Supplement](https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec2_2.html).

**Impact**

The W3/W4 curves can be described as raw stack threshold sensitivity under the analytic source model. They cannot currently support an end-to-end non-detection constraint or the report's approximately 100 km warm-structure interpretation.

**Recommendations**

- Inject sources before photometry and run every injection through the complete candidate-vetting path.
- Measure both threshold-crossing completeness and final-candidate completeness, labeling them separately.
- Do not automatically interpret expected W3/W4 cadence as evidence against a real source.
- For W3/W4, define an independent confirmation strategy using other WISE bands, external catalogs, archival imaging, or targeted follow-up. If confirmation is unavailable, report sensitivity without an exclusion claim.

### 3. The stated 99% trajectory confidence region is not propagated

**Severity: High**

The [frozen hypothesis](hypotheses.md#L100-L106) calls for seeded Monte Carlo propagation of target covariance into a 99% locus confidence region. The production search instead evaluates the nominal adaptive locus and stores `confidence_level = 0.99` as metadata ([precise-pass implementation](scripts/precise_pass.py#L157-L205)). No production call to the covariance-propagation routine was found.

The 10 arcsec padding used during coarse frame discovery helps ensure that relevant images are retrieved, but it does not cause the final photometry to sample a 99% cross-track confidence region. The residual-motion grid likewise does not automatically represent target-catalog covariance.

The observer-position discussion contains a smaller numerical error: approximately 500 km is WISE's altitude, not its geocentric orbital radius. A roughly 6,900 km observer displacement corresponds to about 17 milliarcseconds at 550 AU rather than 2.5 milliarcseconds. WISE's approximately 525 km altitude is documented in the [JPL WISE press kit](https://www.jpl.nasa.gov/news/press_kits/wise-launchOld.pdf). This correction is negligible relative to WISE's several-arcsecond point-spread function, but the calculation and wording should be fixed.

**Recommendations**

- Propagate each endpoint's full astrometric covariance, including parameter correlations and reference epoch, with a recorded seed and sample count.
- Merge the propagated confidence envelope into the positions actually evaluated by the photometry and injection pipelines.
- Measure coverage empirically and report the largest cross-track uncertainty relative to the WISE response width.
- Use the topocentric observer position or explicitly demonstrate that the geocentric approximation is negligible.

### 4. The injection model is not an end-to-end WISE instrument model

**Severity: High**

Injections are applied analytically to sampled tensors after image acquisition, masking, astrometric projection, interpolation, background estimation, and photometric calibration. They therefore cannot measure failures or attenuation introduced by those stages. They also assume a fixed circular Gaussian response for each band.

The official WISE pipeline used empirical point-spread functions that vary significantly across the focal plane, together with local background estimation and deblending. See the [WISE source-extraction documentation](https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4c.html). WISE photometric interpretation also requires bandpass and color corrections for sufficiently non-stellar spectra, with notable effects in W3 and W4; see the [WISE photometric-calibration documentation](https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html).

The flux-normalization correction is sound for the chosen Gaussian model, but agreement with that model is not evidence that a real WISE source would have the same throughput. The physical-radius table should therefore be treated as an approximate monochromatic interpretation until empirical PRF throughput and source-spectrum corrections are included.

The 32 repetitions per calibration point are also too few to establish a stable 90th-percentile recovery boundary without binomial or bootstrap uncertainty intervals.

**Recommendations**

- Inject sources into individual calibrated images before background subtraction and sampling.
- Draw image-appropriate empirical PRFs, subpixel phases, scan orientations, focal-plane positions, and artifact environments.
- Propagate an explicit source spectrum through the WISE relative spectral response and calibration convention.
- Run enough injections to quote confidence intervals on completeness, not only a point estimate.
- Demonstrate recovery of at least one known moving source processed blindly through the full pipeline.

### 5. The distance, motion, and duty-cycle implementations do not exactly match the hypothesis

**Severity: Moderate to High**

The sampler places 64 distance nodes uniformly in reciprocal distance ([grid construction](scripts/sample_tensor.py#L46-L58)), although the hypothesis states that a log-uniform distance prior governs injection placement. Pointwise sensitivity does not require a prior, but the implementation should not be described as sampling the frozen prior unless injections are actually weighted or drawn from it.

Injections are placed on distance and residual-motion grid nodes. Off-grid sensitivity is represented using an average analytic Gaussian mismatch factor rather than by injecting trajectories continuously between nodes. This cannot reproduce spatially variable confusion or artifacts between neighboring tracks.

Constraint records are created by grouping eight nodes and assigning the first and last node as the interval boundaries ([constraint construction](scripts/injection_calibrate.py#L252-L276)). The resulting ledger leaves gaps between blocks. Those gaps total approximately 9.2% of the declared log-uniform prior mass even though the report describes continuous coverage from 550 to 10,000 AU.

The hypothesis specifies a vector bound of \(|\boldsymbol{\mu}| \leq 1\) arcsec yr\(^{-1}\), while the implementation samples each component from \([-1,1]\). Corner grid points therefore reach \(\sqrt{2}\) arcsec yr\(^{-1}\). This is harmless if an \(L_\infty\) component-wise bound was intended, but the hypothesis and ledger currently describe an \(L_2\) vector bound.

For 50% duty cycle, each exposure is independently activated with Bernoulli probability 0.5. The resulting completeness applies to that independent-flicker model, not to every source with time-averaged duty cycle at least 0.5. Long on/off intervals, phase locking, and structured intermittency can behave differently under the phase and visit vetoes.

**Recommendations**

- Define bin boundaries at appropriate reciprocal-distance midpoints so the constraint intervals tile the full domain.
- Perform continuous off-grid injections in distance and both motion components.
- Either enforce a circular vector-motion bound or explicitly redefine it as a component-wise bound.
- Distinguish an assumed distance prior from the numerical grid used to tabulate pointwise sensitivity.
- Calibrate multiple temporal models, including persistent, independent-flicker, visit-scale, and long-block intermittency.

### 6. Candidate vetoes are useful flags but are not calibrated rejection tests

**Severity: High**

Static-source proximity, bright-star halo, single-phase concentration, and inter-band significance ratios are reasonable reasons to inspect a candidate. They are not individually decisive evidence that the signal cannot be an SGL object. For example:

- A real projected trajectory may cross a static source or bright-star halo by chance.
- An intermittent source may legitimately appear in one observing phase.
- A ratio of detection significances is not an astrophysical color because the bands have different noise, cadence, and backgrounds.
- A source appearing "star-like" does not by itself distinguish natural emission from the technosignature being tested.

The rules should therefore be treated as candidate annotations until their true-positive loss and false-positive rejection rates are measured using complete injections and matched null trials.

**Recommendations**

- Replace hard vetoes with pre-specified likelihood components or candidate-quality flags where feasible.
- Measure the selection function of every veto using injected signals.
- Preserve and report candidates that are ambiguous because of confusion rather than converting ambiguity into a null constraint.
- Require independent evidence for final rejection when a rule can also reject an in-scope signal.

### 7. Exposure-quality filtering needs a sensitivity analysis

**Severity: Moderate**

The precise pass rejects `qual_frame == 0` but does not comprehensively exclude or model all observations with poor scan-quality flags, small South Atlantic Anomaly separation, or close Moon geometry. The production observation set contains nontrivial subsets with these conditions. Official metadata specifically provides `qual_frame`, `qual_scan`, SAA separation, and Moon masking information for quality assessment; see the [NEOWISE image-metadata documentation](https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec3_1c.html).

Local controls reduce some shared contamination, but they do not guarantee that localized trails, halos, latent images, or moving artifacts affect the nominal and control tracks equally.

**Recommendations**

- Define a primary, pre-specified quality mask using official WISE recommendations.
- Repeat the analysis with stricter and looser quality selections and report changes in thresholds, candidates, and recovery limits.
- Include quality variables in null-trial matching so target and control trajectories have comparable exposure conditions.

### 8. Target selection supports a targeted survey, not a population test

**Severity: Moderate**

The Picky Network configuration uses explicit but judgment-based weights and thresholds. It is a prioritization heuristic, not a generative model with validated probabilities or inclusion weights. The universal-target document also identifies seven tier-2 science targets deferred by the track budget. The completed survey is therefore complete for its frozen 88-endpoint portfolio, not for every potentially relevant target within 10 pc.

This does not invalidate the search, but it limits interpretation. A null result cannot be converted into the prevalence of relays or evidence against a network-selection model without a population likelihood and known target inclusion probabilities.

**Recommendations**

- Describe the survey consistently as targeted coverage of a frozen 88-endpoint sample.
- Separate physical search completeness from target-population completeness.
- If population inference is desired later, define a generative target-selection model, inclusion probabilities, and an occurrence-rate likelihood before examining results.

### 9. Reproducibility and bookkeeping claims exceed the implementation

**Severity: Moderate**

The analysis-run identifier hashes configuration and registry information but includes sampled-tensor filenames rather than the tensor contents ([run construction](scripts/injection_calibrate.py#L204-L220)). Observation and intersection hashes are represented by a reference to a tensor manifest, but a complete content-addressed tensor manifest is not part of the present production chain. The `--only-missing` path can also retain old tensors after upstream changes.

Consequently, the report's statements that every stage is immutable, content-addressed, and fully checksummed are too strong. Candidate-record counts and some hypothesis-version endpoint counts also require reconciliation. Automated tests are notably absent despite the earlier flux-normalization error.

**Recommendations**

- Hash the byte content and schema/version of every tensor, cutout, mask, configuration, and upstream record used by a run.
- Make stale-product detection mandatory rather than relying only on filenames and `--only-missing` behavior.
- Add invariant tests for Gaussian normalization, magnitude conversion, interval tiling, search-cell boundaries, record counts, and run hashes.
- Generate report tables directly from validated records and fail the build on mismatched counts or unresolved supersession links.

## Interpretation supported by the present analysis

The current survey establishes that:

- The frozen target trajectories were searched with the implemented nominal-locus, residual-motion, and Gaussian matched-filter model.
- Seventy target-band cells exceeded thresholds defined by eight local controls.
- Those crossings were consistent in number with the null expectation of the threshold construction.
- No crossing remained compelling after the survey's current artifact, cadence, phase, and source-proximity review.
- The corrected tables describe threshold sensitivity for analytic Gaussian injections into sampled tensors.

It does not yet establish that an in-scope source had at least a 90% probability of surviving the complete pipeline or that the probability of at least one false survey detection was controlled at a stated level.

## Priority remediation plan

### Priority 1: Repair the statistical experiment

1. Freeze the complete decision rule.
2. Define a global detection statistic and error-rate target.
3. Construct a large, exchangeable null ensemble matched in sky position, cadence, coverage, and contamination.
4. Reserve independent data for confirmation.

### Priority 2: Measure end-to-end completeness

1. Inject realistic sources into calibrated WISE images.
2. Run acquisition, masking, photometry, scoring, and every veto unchanged.
3. Sample continuous distance, motion, phase, flux, spectrum, and temporal behavior.
4. Quote confidence intervals and distinguish threshold completeness from final-candidate completeness.

### Priority 3: Correct the W3/W4 interpretation

1. Remove expected low cadence as an automatic rejection criterion.
2. Define an external confirmation procedure for single-visit cryogenic detections.
3. Until that is validated, label W3/W4 results as raw sensitivity estimates and withdraw end-to-end warm-structure exclusions.

### Priority 4: Validate geometry and implementation

1. Propagate target covariance into the searched locus.
2. Correct observer-position geometry and motion-bound definitions.
3. Close distance-interval gaps and test continuous off-grid recovery.
4. Recover a known moving object as a positive control.

### Priority 5: Strengthen reproducibility

1. Content-hash all derived products and inputs.
2. Add automated scientific and schema invariants.
3. Rebuild report claims and tables from the validated ledger.

## Final recommendation

The project should be presented as a strong pilot and exploratory non-detection pending these revisions. W1/W2 are the most promising bands for a defensible multi-epoch exclusion analysis. W3/W4 remain scientifically useful for discovery sensitivity, but the present cadence veto prevents them from supporting the report's strongest physical constraints.

After the false-alarm model, end-to-end injections, covariance propagation, and W3/W4 adjudication are repaired, the same transparent infrastructure could support a credible quantitative null result.
