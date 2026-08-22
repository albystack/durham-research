# Direct square-grid Glauber production result

The completed frozen-edge Gamma campaign does not provide robust evidence for
a positive quadratic-log contribution to central-height variance over
`L=2--20`.

For the primary unweighted all-environment fit:

| Component | c | 95% environment-bootstrap interval | Delta BIC | LOOCV RMSE log / quadratic |
|---|---:|---:|---:|---:|
| Conditional | -0.048 | [-0.320, 0.227] | -1.283 | 0.114 / 0.142 |
| Disorder | 0.108 | [-0.545, 0.794] | -1.632 | 0.329 / 0.386 |
| Total | 0.060 | [-0.721, 0.822] | -1.845 | 0.242 / 0.279 |

Positive Delta BIC favours the quadratic extension. All three point estimates
and leave-one-size-out prediction errors favour the simpler ordinary-log
description. Weighting, lower-size cutoffs, and the prespecified start-gap
sensitivity filter do not produce a robust positive coefficient. The all-one
control's disorder covariance remains numerically near zero.

## Mixing qualification

The largest standardized extremal-start gaps are 6.28 at `L=16` and 7.76 at
`L=20`. The lowest adjacent-beta pair has pooled exchange acceptance 0.062 at
`L=20`, while target-adjacent exchange remains healthy. Full replica round
trips were not recorded, so the larger-size result is not a mixing certificate.

The defensible interpretation is a **finite-size null result with a frontier
mixing caveat**, not a proof of asymptotic ordinary-log behaviour.

## Files

- `production_component_summary.csv`: size-level estimates and bootstrap
  intervals for all and diagnostic-filtered environments.
- `production_diagnostic_summary.csv`: start gaps, ESS, and exchange summaries.
- `production_scaling_comparison.csv`: all cutoff, weighting, component, and
  filter comparisons.
- `ANALYSIS_METHOD.txt`: estimator and resampling declaration.

The retained environment blocks and exact reproduction command are in
[`../../data/glauber_square_grid_20260822/`](../../data/glauber_square_grid_20260822/).
