# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **1024**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | -0.01364 | [-0.09847, 0.06776] | 0.356 | -3.199 | 0.097 / 0.130 |
| Disorder covariance | weighted | 0.002352 | [-0.05036, 0.05271] | 0.520 | -3.444 | 0.079 / 0.087 |
| Disordered conditional | unweighted | 0.03852 | [-0.07631, 0.1336] | 0.694 | -1.885 | 0.121 / 0.134 |
| Disordered conditional | weighted | 0.02659 | [-0.03955, 0.07917] | 0.755 | -1.655 | 0.119 / 0.106 |
| No-disorder marginal control | unweighted | 0.0337 | [-0.0782, 0.1207] | 0.654 | -2.647 | 0.145 / 0.151 |
| No-disorder marginal control | weighted | -0.0126 | [-0.07115, 0.03209] | 0.222 | -2.917 | 0.154 / 0.169 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.02919 | [-0.1853, 0.1211] | 0.356 | -1.756 | 0.089 / 0.102 |
| 1/16 | 0.0608 | [-0.1042, 0.2411] | 0.736 | 0.839 | 0.066 / 0.051 |
| 1/8 | -0.01633 | [-0.2075, 0.185] | 0.411 | -2.010 | 0.119 / 0.296 |
| 1/4 | -0.06981 | [-0.2568, 0.1293] | 0.237 | -0.354 | 0.106 / 0.144 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.008183 | [-0.15, 0.1507] | -1.964 | 0.018 / 0.021 |
| 1/16 | 0.2094 | [-0.06675, 0.4834] | 2.278 | 0.249 / 0.222 |
| 1/8 | 0.01007 | [-0.1512, 0.1511] | -1.899 | 0.033 / 0.030 |
| 1/4 | -0.09291 | [-0.2606, 0.05714] | -0.961 | 0.146 / 0.159 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | 0.03241 | [-0.03155, 0.08645] | 0.811 | -2.311 | 0.118 / 0.105 |
| Disorder covariance | -0.004031 | [-0.0553, 0.04446] | 0.416 | -3.441 | 0.079 / 0.088 |
| No-disorder marginal | -0.006673 | [-0.0627, 0.03747] | 0.295 | -3.399 | 0.155 / 0.168 |
| No-disorder replica covariance | 0.03606 | [-0.0151, 0.08501] | 0.912 | -1.482 | 0.072 / 0.055 |
