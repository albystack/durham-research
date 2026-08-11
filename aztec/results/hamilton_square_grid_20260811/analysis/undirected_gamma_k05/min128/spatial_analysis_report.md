# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **128**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.006922 | [-0.007847, 0.02068] | 0.818 | -2.288 | 0.111 / 0.116 |
| Disorder covariance | weighted | 0.003539 | [-0.001453, 0.008401] | 0.914 | -1.002 | 0.117 / 0.115 |
| Disordered conditional | unweighted | -0.01689 | [-0.03545, -0.002268] | 0.014 | -1.098 | 0.227 / 0.222 |
| Disordered conditional | weighted | -0.005818 | [-0.0129, -6.47e-06] | 0.025 | -0.933 | 0.232 / 0.229 |
| No-disorder marginal control | unweighted | 0.006992 | [-0.01271, 0.02055] | 0.671 | -2.695 | 0.151 / 0.151 |
| No-disorder marginal control | weighted | 0.002253 | [-0.003362, 0.006699] | 0.741 | -2.661 | 0.158 / 0.154 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.01145 | [-0.01352, 0.03725] | 0.807 | 0.905 | 0.045 / 0.064 |
| 1/16 | 0.002451 | [-0.02801, 0.03479] | 0.557 | -2.354 | 0.046 / 0.074 |
| 1/8 | -0.008823 | [-0.03231, 0.01452] | 0.228 | -2.121 | 0.165 / 0.168 |
| 1/4 | 0.02261 | [-0.00698, 0.05102] | 0.931 | 0.656 | 0.135 / 0.115 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | -0.0007453 | [-0.0293, 0.02315] | -2.450 | 0.014 / 0.024 |
| 1/16 | 0.03913 | [-0.008654, 0.08338] | 1.816 | 0.260 / 0.243 |
| 1/8 | 0.00119 | [-0.02799, 0.02573] | -2.402 | 0.029 / 0.029 |
| 1/4 | -0.0116 | [-0.04251, 0.01422] | -1.848 | 0.150 / 0.147 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.00525 | [-0.01233, 0.0005029] | 0.036 | -1.341 | 0.232 / 0.229 |
| Disorder covariance | 0.003652 | [-0.001268, 0.008608] | 0.925 | -1.778 | 0.117 / 0.115 |
| No-disorder marginal | 0.002471 | [-0.003143, 0.006937] | 0.775 | -2.936 | 0.158 / 0.155 |
| No-disorder replica covariance | -0.0006318 | [-0.005671, 0.004434] | 0.412 | -3.810 | 0.062 / 0.067 |
