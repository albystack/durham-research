# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **1024**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.02611 | [-0.05839, 0.1063] | 0.722 | -2.842 | 0.112 / 0.140 |
| Disorder covariance | weighted | 0.02153 | [-0.02812, 0.06811] | 0.800 | -2.285 | 0.107 / 0.121 |
| Disordered conditional | unweighted | -0.08921 | [-0.1882, -0.003551] | 0.021 | -1.354 | 0.223 / 0.224 |
| Disordered conditional | weighted | -0.06451 | [-0.1307, -0.009629] | 0.013 | 1.229 | 0.218 / 0.213 |
| No-disorder marginal control | unweighted | 0.0337 | [-0.07586, 0.1188] | 0.661 | -2.647 | 0.145 / 0.151 |
| No-disorder marginal control | weighted | -0.01208 | [-0.07017, 0.03169] | 0.228 | -2.962 | 0.154 / 0.168 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.01476 | [-0.1292, 0.1669] | 0.561 | -1.919 | 0.056 / 0.114 |
| 1/16 | -0.04404 | [-0.2234, 0.1537] | 0.312 | -0.926 | 0.074 / 0.071 |
| 1/8 | 0.01061 | [-0.1156, 0.1402] | 0.563 | -2.066 | 0.166 / 0.217 |
| 1/4 | 0.1231 | [-0.04462, 0.289] | 0.921 | 0.979 | 0.118 / 0.092 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.008183 | [-0.1499, 0.1471] | -1.964 | 0.018 / 0.021 |
| 1/16 | 0.2094 | [-0.07142, 0.4747] | 2.278 | 0.249 / 0.222 |
| 1/8 | 0.01007 | [-0.1504, 0.1517] | -1.899 | 0.033 / 0.030 |
| 1/4 | -0.09291 | [-0.2599, 0.06143] | -0.961 | 0.146 / 0.159 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.05811 | [-0.1246, -0.003556] | 0.020 | 0.140 | 0.217 / 0.213 |
| Disorder covariance | 0.02425 | [-0.02417, 0.07264] | 0.834 | -2.502 | 0.108 / 0.122 |
| No-disorder marginal | -0.007372 | [-0.06519, 0.03746] | 0.289 | -3.385 | 0.155 / 0.168 |
| No-disorder replica covariance | 0.03564 | [-0.01429, 0.08345] | 0.916 | -1.512 | 0.072 / 0.054 |
