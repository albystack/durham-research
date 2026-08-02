# Retained numerical results

All tracked intervals use 10,000 deterministic within-size bootstrap resamples
with seed `20260802`. A positive `delta_bic_log_minus_log2` favours the affine
squared-log curve; a negative value favours the affine ordinary-log curve.

- [`analysis_report.md`](analysis_report.md) gives the scientific summary and
  finite-size caveats.
- `height/` contains the 35,536-sample single-height summary, fit report, fitted
  curves, and figure.
- `double_dimer/` contains the 28,304-pair summary, difference and component
  fits, fitted curves, and three figures.
- [`model_comparison_by_cutoff.csv`](model_comparison_by_cutoff.csv) collects
  cutoff sensitivity for the single height, double difference, paired
  covariance, and independent-data remainder.

The fit reports are plain `key=value` files so that plotting and downstream
analysis can read them without external packages. CSV tables retain more
digits than the prose report.

The affine BIC comparison is based on unweighted least squares across
size-level variance or covariance estimates. It should be read as a diagnostic,
not as a formal raw-sample likelihood comparison.
