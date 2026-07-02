# Final aggregate results

These tables summarize the combined 379,000-walk dataset used for the final
finite-size analysis. The `L=16` observations were produced with the optimized
Julia implementation, `L=32` through `L=1024` with the archived Python
implementation, and `L=2048` and `L=4096` with the first Julia implementation.
A targeted optimized Julia run adds 500 baseline and 500 `Gamma(shape=1)`
observations at `L=8192`.

The published confidence intervals use 1,000 environment-clustered bootstrap
replicates with annealed seed `20260702` and quenched seed `20268621`.

## Files

- `summary.csv`: per-distribution, per-size sample counts, annealed and quenched
  variance estimates, cluster-jackknife standard errors, path statistics, exit
  diagnostics, and runtime.
- `loglog_fits.csv`: power fits for `Var(W_L) = C(log L)^p`, including
  environment-clustered bootstrap intervals.
- `scaling_model_comparison.csv`: additive `a + b log L` versus
  `a + b(log L)^2` comparisons using SSE, R-squared, AIC, and BIC.
- `local_effective_exponents.csv`: adjacent-size effective powers.
- `pointwise_ratios.csv`: pointwise log-variance diagnostics.
- `figures/professor_scaling_overview.pdf`: four-page overview of the scaling
  curves, exponent intervals, and BIC comparisons.
- `figures/all_distributions_detailed.pdf`: one detailed page per distribution.

The raw walk-level CSVs are excluded from Git because of their size. These
aggregate files are sufficient to inspect the reported estimates, but exact
bootstrap reproduction requires the raw environment-level observations. A
separate versioned data archive should accompany any formal research release.

## Interpretation

The ten-size annealed estimates are `p=1.045` for baseline and `p=1.077` for
`Gamma(shape=1)`, with clustered-bootstrap intervals `[0.968,1.112]` and
`[0.999,1.149]`, respectively. All fifteen annealed and quenched additive model
comparisons now favor ordinary logarithmic growth by BIC over their available
size ranges. The evidence therefore does not support a robust super-rough
regime over the simulated sizes.
