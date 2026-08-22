# Retained numerical results

The retained centre-height intervals use 10,000 deterministic within-size
bootstrap resamples with seed `20260802`; the spatial experiment uses seed
`20260805` and preserves whole environments when pooling separations. A
positive `delta_bic_log_minus_log2` favours the log-squared extension; a
negative value favours the ordinary-log curve.

- [`analysis_report.md`](analysis_report.md) gives the scientific summary and
  finite-size caveats.
- `height/` contains the 35,536-sample single-height summary, fit report, fitted
  curves, and figure.
- `double_dimer/` contains the 28,304-pair summary, difference and component
  fits, fitted curves, and three figures.
- `spatial/` contains the 6,900-environment Gamma spatial experiment, the
  5,500-pair uniform control, pooled and weighted comparisons, report, and
  attachment-ready figures.
- [`model_comparison_by_cutoff.csv`](model_comparison_by_cutoff.csv) collects
  cutoff sensitivity for the single height, double difference, paired
  covariance, and independent-data remainder.
- `hamilton_20260809/` contains the cross-parameter Aztec spatial comparison.
- `hamilton_square_grid_20260811/` contains the five-law structured
  square-grid comparison, its sensitivity windows, vector figures, and audited
  scheduler provenance.
- [`glauber_square_grid_20260822/`](glauber_square_grid_20260822/) contains the
  direct frozen-edge production components, diagnostics, scaling comparisons,
  and interpretation.

The newer direct Glauber production traces remain in external project storage
because they are substantially larger than the retained Git datasets. Their
environment-blocked estimates and scientific conclusions are recorded in
[`../../docs/RESULTS.md`](../../docs/RESULTS.md), with reproduction commands in
[`../../docs/REPRODUCIBILITY.md`](../../docs/REPRODUCIBILITY.md).

The fit reports are plain `key=value` files so that plotting and downstream
analysis can read them without external packages. CSV tables retain more
digits than the prose report.

Unweighted and inverse-bootstrap-variance-weighted BIC comparisons are based
on size-level variance or covariance estimates. They should be read as
diagnostics, not as formal raw-sample likelihood comparisons.
