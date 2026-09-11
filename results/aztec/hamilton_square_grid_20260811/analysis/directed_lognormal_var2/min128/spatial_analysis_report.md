# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **128**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.001476 | [-0.01431, 0.01705] | 0.550 | -3.832 | 0.174 / 0.178 |
| Disorder covariance | weighted | 0.001364 | [-0.003765, 0.006342] | 0.703 | -3.594 | 0.147 / 0.147 |
| Disordered conditional | unweighted | -0.00324 | [-0.02331, 0.01248] | 0.255 | -3.726 | 0.173 / 0.182 |
| Disordered conditional | weighted | -0.004057 | [-0.01087, 0.001674] | 0.075 | -1.354 | 0.159 / 0.157 |
| No-disorder marginal control | unweighted | 0.006992 | [-0.01261, 0.02079] | 0.679 | -2.695 | 0.151 / 0.151 |
| No-disorder marginal control | weighted | 0.002261 | [-0.003237, 0.006695] | 0.732 | -2.651 | 0.158 / 0.154 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.009132 | [-0.04324, 0.02489] | 0.297 | -1.303 | 0.106 / 0.128 |
| 1/16 | -0.01169 | [-0.03609, 0.01263] | 0.172 | -0.915 | 0.115 / 0.129 |
| 1/8 | 0.02473 | [-0.007955, 0.05768] | 0.925 | -0.762 | 0.249 / 0.325 |
| 1/4 | 0.001997 | [-0.03247, 0.03762] | 0.518 | -2.473 | 0.186 / 0.199 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | -0.0007453 | [-0.02878, 0.02366] | -2.450 | 0.014 / 0.024 |
| 1/16 | 0.03913 | [-0.009534, 0.08334] | 1.816 | 0.260 / 0.243 |
| 1/8 | 0.00119 | [-0.02795, 0.02566] | -2.402 | 0.029 / 0.029 |
| 1/4 | -0.0116 | [-0.04294, 0.01521] | -1.848 | 0.150 / 0.147 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.004474 | [-0.0111, 0.001049] | 0.054 | -1.823 | 0.159 / 0.156 |
| Disorder covariance | 0.001218 | [-0.003915, 0.006297] | 0.665 | -3.649 | 0.148 / 0.147 |
| No-disorder marginal | 0.002489 | [-0.003023, 0.006824] | 0.772 | -2.921 | 0.158 / 0.154 |
| No-disorder replica covariance | -0.000662 | [-0.005695, 0.004294] | 0.394 | -3.805 | 0.063 / 0.067 |
