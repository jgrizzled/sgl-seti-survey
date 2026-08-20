# Erratum — WISE depth labels (found 2026-08-20 by the ZTF asteroid control)

`sglsurvey/photometry.py` built its Gaussian matched-filter kernel with
unit *peak*, so the least-squares "flux" was the PSF peak amplitude,
not the total flux. S/N, control thresholds, injection-recovery
fractions, vetoes and candidate decisions are invariant to this scale.
Every **magnitude** derived from it (the `m90` depths on all WISE
Constraint records through AnalysisRun `run-76251a7d4ad3`, and the
report figures) is too faint — i.e. too optimistic — by the constant

    Δm = 2.5 log10(2π σ_pix²),  σ_pix = FWHM / 2.3548 / 2.75″

| band | FWHM | Δm |
|---|---|---|
| W1 | 6.1″ | 1.87 mag |
| W2 | 6.4″ | 1.97 mag |
| W3 | 6.5″ | 2.01 mag |
| W4 | 12.0″ | 3.34 mag |

Corrected depth = quoted depth − Δm (e.g. W1 "17.5" → 15.6 Vega). The
physical size limits in `report/wise_survey_v2.md` (e.g. ≥ 39 km warm
structures via W3) scale accordingly (flux ×6.4 → size ×2.5 at fixed
temperature). The kernel is now unit-sum (total-flux convention); any
regenerated WISE tensors/constraints will carry the corrected scale
automatically. The existing records have NOT been rewritten — decide
whether to rerun `sample_tensor.py` + `injection_calibrate.py` (tensors
only; cutouts cached) or apply the constant per-band correction in the
report.
