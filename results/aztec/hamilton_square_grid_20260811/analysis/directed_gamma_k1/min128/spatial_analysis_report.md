# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **128**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.004675 | [-0.01157, 0.02028] | 0.693 | -3.227 | 0.108 / 0.107 |
| Disorder covariance | weighted | 0.001274 | [-0.003991, 0.006204] | 0.675 | -3.589 | 0.096 / 0.096 |
| Disordered conditional | unweighted | -0.007352 | [-0.02895, 0.01073] | 0.163 | -3.185 | 0.201 / 0.199 |
| Disordered conditional | weighted | -0.004126 | [-0.01073, 0.00124] | 0.063 | -2.194 | 0.187 / 0.185 |
| No-disorder marginal control | unweighted | 0.006992 | [-0.01264, 0.02043] | 0.676 | -2.695 | 0.151 / 0.151 |
| No-disorder marginal control | weighted | 0.002269 | [-0.003396, 0.006564] | 0.730 | -2.657 | 0.158 / 0.154 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.0001029 | [-0.03494, 0.03399] | 0.492 | -2.485 | 0.069 / 0.086 |
| 1/16 | 0.009665 | [-0.01422, 0.03277] | 0.766 | -1.929 | 0.137 / 0.176 |
| 1/8 | 0.009982 | [-0.01828, 0.03867] | 0.743 | -1.800 | 0.123 / 0.117 |
| 1/4 | -0.001049 | [-0.03178, 0.03062] | 0.452 | -2.478 | 0.090 / 0.140 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | -0.0007453 | [-0.02901, 0.02313] | -2.450 | 0.014 / 0.024 |
| 1/16 | 0.03913 | [-0.009104, 0.0833] | 1.816 | 0.260 / 0.243 |
| 1/8 | 0.00119 | [-0.02846, 0.02548] | -2.402 | 0.029 / 0.029 |
| 1/4 | -0.0116 | [-0.04264, 0.01411] | -1.848 | 0.150 / 0.147 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.00444 | [-0.01085, 0.0009346] | 0.050 | -1.678 | 0.187 / 0.185 |
| Disorder covariance | 0.001144 | [-0.003937, 0.005915] | 0.659 | -3.670 | 0.096 / 0.096 |
| No-disorder marginal | 0.002409 | [-0.003329, 0.006846] | 0.754 | -2.984 | 0.158 / 0.155 |
| No-disorder replica covariance | -0.0006665 | [-0.005831, 0.004416] | 0.387 | -3.804 | 0.063 / 0.067 |
