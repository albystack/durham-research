# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **512**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.01865 | [-0.01678, 0.05188] | 0.836 | -2.935 | 0.180 / 0.207 |
| Disorder covariance | weighted | -0.005595 | [-0.02355, 0.01184] | 0.251 | -3.286 | 0.151 / 0.167 |
| Disordered conditional | unweighted | -0.02327 | [-0.08022, 0.02696] | 0.141 | -3.028 | 0.240 / 0.248 |
| Disordered conditional | weighted | -0.02564 | [-0.04959, -0.005965] | 0.007 | 0.213 | 0.197 / 0.194 |
| No-disorder marginal control | unweighted | 0.01396 | [-0.04007, 0.05421] | 0.610 | -3.055 | 0.148 / 0.154 |
| No-disorder marginal control | weighted | -0.004226 | [-0.02547, 0.01171] | 0.240 | -3.334 | 0.159 / 0.163 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.01438 | [-0.05649, 0.08625] | 0.651 | -2.151 | 0.146 / 0.144 |
| 1/16 | -0.02369 | [-0.1081, 0.0655] | 0.296 | 0.145 | 0.065 / 0.070 |
| 1/8 | -0.001791 | [-0.07365, 0.06712] | 0.469 | -2.293 | 0.041 / 0.094 |
| 1/4 | 0.0857 | [0.004725, 0.1641] | 0.981 | -0.864 | 0.321 / 0.463 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.005163 | [-0.07267, 0.07464] | -2.086 | 0.019 / 0.023 |
| 1/16 | 0.1065 | [-0.03095, 0.2343] | 2.504 | 0.254 / 0.230 |
| 1/8 | -0.0096 | [-0.09037, 0.05672] | -1.578 | 0.028 / 0.038 |
| 1/4 | -0.04624 | [-0.1292, 0.02826] | -0.980 | 0.150 / 0.144 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.02667 | [-0.05064, -0.007043] | 0.005 | 1.971 | 0.198 / 0.196 |
| Disorder covariance | -0.004514 | [-0.02134, 0.01212] | 0.290 | -3.413 | 0.152 / 0.167 |
| No-disorder marginal | -0.002903 | [-0.02426, 0.01316] | 0.289 | -3.594 | 0.160 / 0.164 |
| No-disorder replica covariance | 0.01288 | [-0.005758, 0.03125] | 0.908 | -1.820 | 0.071 / 0.057 |
