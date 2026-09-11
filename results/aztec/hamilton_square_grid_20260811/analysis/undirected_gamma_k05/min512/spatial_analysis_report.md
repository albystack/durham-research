# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **512**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.01311 | [-0.02814, 0.05151] | 0.728 | -2.933 | 0.111 / 0.123 |
| Disorder covariance | weighted | 0.00636 | [-0.01204, 0.02394] | 0.741 | -2.936 | 0.109 / 0.110 |
| Disordered conditional | unweighted | -0.0445 | [-0.09436, -0.001927] | 0.021 | -1.152 | 0.224 / 0.220 |
| Disordered conditional | weighted | -0.03482 | [-0.05985, -0.01377] | 0.001 | 5.671 | 0.225 / 0.212 |
| No-disorder marginal control | unweighted | 0.01396 | [-0.04087, 0.0539] | 0.606 | -3.055 | 0.148 / 0.154 |
| No-disorder marginal control | weighted | -0.004124 | [-0.02526, 0.01219] | 0.247 | -3.351 | 0.159 / 0.162 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.01036 | [-0.05907, 0.08369] | 0.596 | -1.909 | 0.050 / 0.082 |
| 1/16 | -0.008341 | [-0.09284, 0.08451] | 0.411 | -2.106 | 0.059 / 0.087 |
| 1/8 | -0.01076 | [-0.07359, 0.05337] | 0.370 | -2.234 | 0.165 / 0.182 |
| 1/4 | 0.06121 | [-0.02028, 0.142] | 0.926 | 0.985 | 0.128 / 0.095 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.005163 | [-0.07328, 0.07252] | -2.086 | 0.019 / 0.023 |
| 1/16 | 0.1065 | [-0.03008, 0.2355] | 2.504 | 0.254 / 0.230 |
| 1/8 | -0.0096 | [-0.09052, 0.05901] | -1.578 | 0.028 / 0.038 |
| 1/4 | -0.04624 | [-0.1308, 0.02859] | -0.980 | 0.150 / 0.144 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.03112 | [-0.05715, -0.0105] | 0.002 | 3.307 | 0.224 / 0.212 |
| Disorder covariance | 0.006979 | [-0.01159, 0.02509] | 0.770 | -3.114 | 0.109 / 0.110 |
| No-disorder marginal | -0.002729 | [-0.02427, 0.01394] | 0.296 | -3.607 | 0.160 / 0.163 |
| No-disorder replica covariance | 0.01277 | [-0.005669, 0.03143] | 0.913 | -1.849 | 0.071 / 0.057 |
