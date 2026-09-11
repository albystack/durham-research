# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **512**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.01073 | [-0.0343, 0.056] | 0.665 | -3.263 | 0.117 / 0.120 |
| Disorder covariance | weighted | 0.01259 | [-0.006016, 0.03073] | 0.907 | -1.714 | 0.103 / 0.098 |
| Disordered conditional | unweighted | -0.01086 | [-0.07247, 0.04175] | 0.269 | -3.497 | 0.211 / 0.211 |
| Disordered conditional | weighted | -0.0309 | [-0.05379, -0.01197] | 0.001 | 4.346 | 0.184 / 0.182 |
| No-disorder marginal control | unweighted | 0.01396 | [-0.0402, 0.05343] | 0.607 | -3.055 | 0.148 / 0.154 |
| No-disorder marginal control | weighted | -0.004298 | [-0.02548, 0.01155] | 0.231 | -3.323 | 0.159 / 0.163 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.005908 | [-0.1058, 0.09086] | 0.453 | -2.242 | 0.073 / 0.116 |
| 1/16 | 0.0347 | [-0.03169, 0.09847] | 0.840 | -1.380 | 0.151 / 0.235 |
| 1/8 | 0.0327 | [-0.04604, 0.1147] | 0.783 | -1.353 | 0.123 / 0.120 |
| 1/4 | -0.01858 | [-0.1032, 0.07232] | 0.326 | -2.016 | 0.109 / 0.196 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.005163 | [-0.07341, 0.07274] | -2.086 | 0.019 / 0.023 |
| 1/16 | 0.1065 | [-0.03203, 0.2361] | 2.504 | 0.254 / 0.230 |
| 1/8 | -0.0096 | [-0.08964, 0.0587] | -1.578 | 0.028 / 0.038 |
| 1/4 | -0.04624 | [-0.1319, 0.02703] | -0.980 | 0.150 / 0.144 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.03128 | [-0.0548, -0.01238] | 0.001 | 4.707 | 0.183 / 0.182 |
| Disorder covariance | 0.01189 | [-0.00676, 0.03014] | 0.885 | -2.091 | 0.103 / 0.098 |
| No-disorder marginal | -0.0031 | [-0.02414, 0.01288] | 0.269 | -3.581 | 0.160 / 0.163 |
| No-disorder replica covariance | 0.01299 | [-0.005765, 0.03141] | 0.911 | -1.785 | 0.071 / 0.056 |
