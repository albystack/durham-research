# Spatial Gamma-disorder experiment

This folder retains the publication-focused spatial height-increment analysis:

- 6,900 independent Gamma environments, with two conditional tilings per
  environment;
- 5,500 independent uniform-control replica pairs;
- orders 128 through 1,300;
- symmetric central-row separations `L/32`, `L/16`, `L/8`, and `L/4`;
- 10,000 environment-level bootstrap repetitions.

The primary pooled comparison gives each separation its own intercept and
ordinary-log slope while fitting one common coefficient `c` in

```text
a_fraction + b_fraction log(r) + c (log(r))^2.
```

Environment indices are resampled jointly across the four separations. The
unweighted Gamma disorder estimate is `c = 0.5843`, with 95% bootstrap interval
`[0.1957, 0.9520]`; the inverse-variance-weighted sensitivity estimate is
`c = 0.4948`, with interval `[0.1932, 0.7977]`. The Gamma conditional and
uniform marginal-control intervals both contain zero.

Recommended email attachments:

1. `pooled_log2_coefficients.png` — primary coefficient comparison;
2. `gamma_disorder_spatial_fits.png` — results at all four separations;
3. `spatial_variance_decomposition.png` — physical decomposition;
4. `uniform_control_spatial_fits.png` — optional control detail.

`spatial_analysis_report.md` contains the complete concise report. CSV files
retain the summaries, per-separation comparisons, pooled comparison, weighted
sensitivity, cutoff sensitivity, and fitted curves.
