# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **1024**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.05481 | [-0.02202, 0.1244] | 0.921 | -2.133 | 0.200 / 0.230 |
| Disorder covariance | weighted | -0.003921 | [-0.05009, 0.0411] | 0.419 | -3.441 | 0.163 / 0.203 |
| Disordered conditional | unweighted | -0.04775 | [-0.167, 0.0548] | 0.155 | -2.891 | 0.266 / 0.275 |
| Disordered conditional | weighted | -0.04552 | [-0.1106, 0.007798] | 0.042 | -1.787 | 0.215 / 0.225 |
| No-disorder marginal control | unweighted | 0.0337 | [-0.07597, 0.1191] | 0.656 | -2.647 | 0.145 / 0.151 |
| No-disorder marginal control | weighted | -0.01244 | [-0.0694, 0.03248] | 0.227 | -2.929 | 0.154 / 0.168 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.02474 | [-0.1227, 0.1752] | 0.628 | -1.991 | 0.147 / 0.140 |
| 1/16 | -0.04019 | [-0.2192, 0.1392] | 0.328 | -0.504 | 0.062 / 0.085 |
| 1/8 | 0.01315 | [-0.1337, 0.1565] | 0.563 | -1.972 | 0.055 / 0.119 |
| 1/4 | 0.2216 | [0.04918, 0.3867] | 0.996 | -0.030 | 0.362 / 0.573 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.008183 | [-0.1483, 0.1489] | -1.964 | 0.018 / 0.021 |
| 1/16 | 0.2094 | [-0.07304, 0.4845] | 2.278 | 0.249 / 0.222 |
| 1/8 | 0.01007 | [-0.1526, 0.153] | -1.899 | 0.033 / 0.030 |
| 1/4 | -0.09291 | [-0.2631, 0.057] | -0.961 | 0.146 / 0.159 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.05244 | [-0.1168, 0.0001403] | 0.025 | -0.300 | 0.218 / 0.230 |
| Disorder covariance | -0.001267 | [-0.04652, 0.04222] | 0.453 | -3.463 | 0.163 / 0.200 |
| No-disorder marginal | -0.007498 | [-0.06467, 0.03558] | 0.283 | -3.382 | 0.156 / 0.168 |
| No-disorder replica covariance | 0.03564 | [-0.01523, 0.08536] | 0.919 | -1.512 | 0.072 / 0.054 |
