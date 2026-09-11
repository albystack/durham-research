# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **1024**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.02743 | [-0.06373, 0.1168] | 0.705 | -3.117 | 0.214 / 0.240 |
| Disorder covariance | weighted | -0.03614 | [-0.08594, 0.01371] | 0.077 | -1.071 | 0.198 / 0.211 |
| Disordered conditional | unweighted | -0.02841 | [-0.1376, 0.06405] | 0.223 | -3.186 | 0.205 / 0.272 |
| Disordered conditional | weighted | -0.005798 | [-0.07132, 0.04942] | 0.343 | -3.403 | 0.168 / 0.180 |
| No-disorder marginal control | unweighted | 0.0337 | [-0.07442, 0.1197] | 0.661 | -2.647 | 0.145 / 0.151 |
| No-disorder marginal control | weighted | -0.01214 | [-0.06951, 0.03238] | 0.227 | -2.955 | 0.154 / 0.168 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.05772 | [-0.2636, 0.1484] | 0.285 | -0.799 | 0.129 / 0.180 |
| 1/16 | -0.1211 | [-0.2663, 0.02312] | 0.048 | 4.523 | 0.124 / 0.080 |
| 1/8 | 0.2494 | [0.06434, 0.4374] | 0.997 | 4.052 | 0.336 / 0.404 |
| 1/4 | 0.03912 | [-0.1527, 0.2461] | 0.619 | -1.970 | 0.196 / 0.346 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.008183 | [-0.1523, 0.1464] | -1.964 | 0.018 / 0.021 |
| 1/16 | 0.2094 | [-0.07368, 0.4824] | 2.278 | 0.249 / 0.222 |
| 1/8 | 0.01007 | [-0.15, 0.1534] | -1.899 | 0.033 / 0.030 |
| 1/4 | -0.09291 | [-0.2634, 0.059] | -0.961 | 0.146 / 0.159 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.003985 | [-0.06937, 0.0479] | 0.362 | -3.448 | 0.167 / 0.175 |
| Disorder covariance | -0.03973 | [-0.09073, 0.01054] | 0.058 | -1.074 | 0.197 / 0.208 |
| No-disorder marginal | -0.007315 | [-0.06477, 0.03704] | 0.292 | -3.386 | 0.156 / 0.169 |
| No-disorder replica covariance | 0.03606 | [-0.01554, 0.086] | 0.917 | -1.488 | 0.072 / 0.055 |
