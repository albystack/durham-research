# Final aggregate results

These tables summarize the combined 303,000-walk dataset used for the final
finite-size analysis. The smaller sizes (`L=32` through `L=1024`) were produced
with the archived Python implementation; `L=2048` and `L=4096` were produced
with the Julia implementation.

The published confidence intervals use 1,000 environment-clustered bootstrap
replicates with annealed seed `20260628` and quenched seed `20268547`.

## Files

- `summary.csv`: per-distribution, per-size sample counts, annealed and quenched
  variance estimates, cluster-jackknife standard errors, path statistics, exit
  diagnostics, and runtime.
- `loglog_fits.csv`: power fits for `Var(W_L) = C(log L)^p`, including
  environment-clustered bootstrap intervals.
- `scaling_model_comparison.csv`: additive `a + b log L` versus
  `a + b(log L)^2` comparisons using SSE, R-squared, AIC, and BIC.
- `local_effective_exponents.csv`: adjacent-size effective powers.

The raw walk-level CSVs are excluded from Git because of their size. These
aggregate files are sufficient to inspect the reported estimates, but exact
bootstrap reproduction requires the raw environment-level observations. A
separate versioned data archive should accompany any formal research release.

## Interpretation

The baseline power estimate is close to one, and fourteen of fifteen additive
model comparisons favor ordinary logarithmic growth by BIC. `Gamma(shape=1)`
is the only additive squared-log BIC preference, but its direct power estimate
and largest-scale local exponent remain close to one. The evidence therefore
does not support a robust super-rough regime over the simulated sizes.
