# Kasteleyn replay of the direct square-grid campaign

The finite-volume determinantal calculation confirms the previous central-
height conclusion without a Markov-chain mixing assumption. Across the same
960 Gamma environments and sizes `L=2--20`, there is no robust positive
quadratic-log coefficient in the conditional, disorder, or total component.

For the primary unweighted full-range fit:

| Component | c | 95% environment-bootstrap interval | Delta BIC | LOOCV RMSE log / quadratic |
|---|---:|---:|---:|---:|
| Conditional | -0.145 | [-0.379, 0.090] | 2.454 | 0.179 / 0.145 |
| Disorder | 0.130 | [-0.552, 0.785] | -0.625 | 0.241 / 0.283 |
| Total | -0.015 | [-0.801, 0.714] | -2.046 | 0.145 / 0.251 |

Positive Delta BIC favours the quadratic extension. The conditional quadratic
fit has modest model-selection support, but its coefficient is negative and
its interval includes zero; it is not evidence for super-rough growth.

## Frontier comparison

| L | Component | Kasteleyn | MCMC |
|---:|---|---:|---:|
| 16 | Conditional | 2.241 | 2.301 |
| 16 | Disorder | 2.421 | 2.199 |
| 16 | Total | 4.662 | 4.500 |
| 20 | Conditional | 2.057 | 2.302 |
| 20 | Disorder | 3.066 | 3.137 |
| 20 | Total | 5.124 | 5.439 |

Individual MCMC environments can differ materially from their determinantal
moments, especially in rare-event environments and at the frontier. Those
errors largely cancel in the current aggregate estimates. The exact replay
therefore strengthens the finite-size null conclusion but does not establish
asymptotic ordinary-log behaviour.

## Numerical validation

- weighted partition functions and central-height moments agree with complete
  `L=1,2` enumeration to approximately `1e-12`;
- the uniform `L=2` determinant gives all 36 tilings;
- all production environment seeds match exactly;
- one of 960 environments used the automatic 256-bit fallback;
- the maximum retained relative solve residual is `7.39e-12`;
- probability and covariance imaginary residuals are zero at output precision.

## Files

- `kasteleyn_component_summary.csv`: size-level component estimates and CIs.
- `kasteleyn_scaling_comparison.csv`: cutoff, weighting, BIC, prediction, and
  coefficient comparisons.
- `kasteleyn_mcmc_comparison.csv`: environment-by-environment MCMC comparison.
- `kasteleyn_mcmc_size_summary.csv`: size-level MCMC discrepancy diagnostics.
- `NUMERICAL_DIAGNOSTICS.txt`: factorization diagnostics.
- `ANALYSIS_METHOD.txt`: estimator and bootstrap declaration.
- `CHECKSUMS.sha256`: SHA-256 hashes for the retained analysis artifacts.

The environment-level deterministic moments are retained in
[`../../data/glauber_square_grid_kasteleyn_20260822/`](../../data/glauber_square_grid_kasteleyn_20260822/).
