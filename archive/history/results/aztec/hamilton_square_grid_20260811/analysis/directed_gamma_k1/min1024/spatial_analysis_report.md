# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **1024**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.007377 | [-0.0853, 0.101] | 0.553 | -3.425 | 0.130 / 0.139 |
| Disorder covariance | weighted | 0.03454 | [-0.01529, 0.08328] | 0.905 | -1.494 | 0.104 / 0.102 |
| Disordered conditional | unweighted | 0.007824 | [-0.1162, 0.1165] | 0.466 | -3.445 | 0.231 / 0.232 |
| Disordered conditional | weighted | -0.06097 | [-0.1237, -0.008268] | 0.011 | 2.219 | 0.209 / 0.216 |
| No-disorder marginal control | unweighted | 0.0337 | [-0.07685, 0.1188] | 0.653 | -2.647 | 0.145 / 0.151 |
| No-disorder marginal control | weighted | -0.01264 | [-0.07106, 0.03064] | 0.219 | -2.918 | 0.154 / 0.169 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.04923 | [-0.2605, 0.1522] | 0.316 | -1.198 | 0.088 / 0.130 |
| 1/16 | 0.09265 | [-0.04256, 0.223] | 0.909 | -0.707 | 0.172 / 0.307 |
| 1/8 | 0.03535 | [-0.1309, 0.2089] | 0.649 | -1.854 | 0.117 / 0.154 |
| 1/4 | -0.04926 | [-0.2303, 0.1393] | 0.290 | -1.674 | 0.127 / 0.283 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.008183 | [-0.1453, 0.1555] | -1.964 | 0.018 / 0.021 |
| 1/16 | 0.2094 | [-0.07376, 0.4783] | 2.278 | 0.249 / 0.222 |
| 1/8 | 0.01007 | [-0.1513, 0.1576] | -1.899 | 0.033 / 0.030 |
| 1/4 | -0.09291 | [-0.2639, 0.0618] | -0.961 | 0.146 / 0.159 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.05869 | [-0.1202, -0.007887] | 0.014 | 0.555 | 0.208 / 0.216 |
| Disorder covariance | 0.02968 | [-0.02061, 0.07804] | 0.874 | -2.118 | 0.104 / 0.102 |
| No-disorder marginal | -0.007771 | [-0.06597, 0.03577] | 0.276 | -3.376 | 0.156 / 0.168 |
| No-disorder replica covariance | 0.03603 | [-0.0144, 0.0869] | 0.919 | -1.478 | 0.072 / 0.055 |
