# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **512**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.003197 | [-0.04054, 0.04754] | 0.533 | -3.666 | 0.189 / 0.204 |
| Disorder covariance | weighted | -0.006533 | [-0.02526, 0.01125] | 0.231 | -3.184 | 0.156 / 0.156 |
| Disordered conditional | unweighted | -0.008397 | [-0.06191, 0.0361] | 0.272 | -3.566 | 0.184 / 0.215 |
| Disordered conditional | weighted | -0.01249 | [-0.03765, 0.007596] | 0.097 | -1.865 | 0.158 / 0.157 |
| No-disorder marginal control | unweighted | 0.01396 | [-0.04045, 0.05348] | 0.604 | -3.055 | 0.148 / 0.154 |
| No-disorder marginal control | weighted | -0.004132 | [-0.02528, 0.01194] | 0.232 | -3.348 | 0.159 / 0.163 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.02526 | [-0.1278, 0.0713] | 0.303 | -1.121 | 0.115 / 0.158 |
| 1/16 | -0.05635 | [-0.1253, 0.01348] | 0.056 | 3.327 | 0.126 / 0.107 |
| 1/8 | 0.08752 | [-0.00292, 0.1772] | 0.970 | 0.703 | 0.281 / 0.398 |
| 1/4 | 0.00688 | [-0.0872, 0.1092] | 0.535 | -2.286 | 0.188 / 0.246 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.005163 | [-0.0739, 0.07213] | -2.086 | 0.019 / 0.023 |
| 1/16 | 0.1065 | [-0.03305, 0.2354] | 2.504 | 0.254 / 0.230 |
| 1/8 | -0.0096 | [-0.08883, 0.06218] | -1.578 | 0.028 / 0.038 |
| 1/4 | -0.04624 | [-0.1308, 0.02986] | -0.980 | 0.150 / 0.144 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.013 | [-0.03712, 0.006385] | 0.084 | -2.362 | 0.158 / 0.156 |
| Disorder covariance | -0.007155 | [-0.02596, 0.01052] | 0.205 | -3.109 | 0.155 / 0.155 |
| No-disorder marginal | -0.003046 | [-0.02395, 0.01313] | 0.276 | -3.585 | 0.160 / 0.163 |
| No-disorder replica covariance | 0.01319 | [-0.005244, 0.03166] | 0.916 | -1.735 | 0.071 / 0.056 |
