# Retained LERW results

These are compact, authoritative outputs from the completed campaigns. Raw
per-batch walk records are not retained in this repository.

## Common files

- `summary.csv`: estimates by distribution and lattice size.
- `loglog_fits.csv`: fitted exponent \(p\) in
  \(\operatorname{Var}(W_L)=C(\log L)^p\), including bootstrap intervals.
- `scaling_model_comparison.csv`: affine \(\log L\) versus
  \((\log L)^2\) fits with AIC and BIC.

## Result groups

- `site_iid_single/`: 379,000 walks in independent fixed environments.
- `site_iid_paired/`: 758,000 walks in 379,000 shared environments; the
  reported winding observable is the within-environment difference.
- `temporal_gamma_length/`: refreshed Gamma weights, including the
  \(L^{5/4}\) path-length validation.

The paired directory also contains the three presentation figures used in the
top-level README. The temporal-Gamma directory adds:

- `path_length_fits.csv`: global LERW-length exponent and confidence interval;
- `path_length_pointwise.csv`: size-by-size mean lengths and ratios; and
- `temporal_direction_diagnostics.csv`: empirical direction-frequency checks.
