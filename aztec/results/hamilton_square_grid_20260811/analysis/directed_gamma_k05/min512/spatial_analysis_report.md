# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **512**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | -0.01342 | [-0.05422, 0.02682] | 0.240 | -2.477 | 0.091 / 0.097 |
| Disorder covariance | weighted | -0.02267 | [-0.04238, -0.004158] | 0.009 | 8.301 | 0.095 / 0.085 |
| Disordered conditional | unweighted | 0.01802 | [-0.0383, 0.06387] | 0.666 | -2.042 | 0.112 / 0.116 |
| Disordered conditional | weighted | 0.009007 | [-0.01524, 0.02982] | 0.721 | -2.586 | 0.110 / 0.102 |
| No-disorder marginal control | unweighted | 0.01396 | [-0.04044, 0.05488] | 0.603 | -3.055 | 0.148 / 0.154 |
| No-disorder marginal control | weighted | -0.003984 | [-0.025, 0.01228] | 0.244 | -3.375 | 0.159 / 0.162 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.02273 | [-0.09962, 0.05121] | 0.273 | -1.344 | 0.087 / 0.097 |
| 1/16 | 0.02094 | [-0.05891, 0.1038] | 0.667 | -0.773 | 0.065 / 0.061 |
| 1/8 | -0.0006216 | [-0.09132, 0.09701] | 0.472 | -2.302 | 0.089 / 0.196 |
| 1/4 | -0.05125 | [-0.1387, 0.04169] | 0.137 | 1.852 | 0.116 / 0.102 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.005163 | [-0.07279, 0.07375] | -2.086 | 0.019 / 0.023 |
| 1/16 | 0.1065 | [-0.03237, 0.2323] | 2.504 | 0.254 / 0.230 |
| 1/8 | -0.0096 | [-0.08915, 0.06126] | -1.578 | 0.028 / 0.038 |
| 1/4 | -0.04624 | [-0.1319, 0.02857] | -0.980 | 0.150 / 0.144 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | 0.01058 | [-0.01342, 0.03072] | 0.777 | -2.807 | 0.110 / 0.101 |
| Disorder covariance | -0.0248 | [-0.04424, -0.006096] | 0.005 | 2.851 | 0.095 / 0.086 |
| No-disorder marginal | -0.002757 | [-0.02402, 0.0131] | 0.282 | -3.604 | 0.160 / 0.163 |
| No-disorder replica covariance | 0.01313 | [-0.005344, 0.03151] | 0.919 | -1.764 | 0.071 / 0.056 |
