# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **128**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | 0.008425 | [-0.004333, 0.02014] | 0.900 | -2.689 | 0.175 / 0.183 |
| Disorder covariance | weighted | 0.003129 | [-0.001898, 0.007922] | 0.889 | -2.294 | 0.165 / 0.161 |
| Disordered conditional | unweighted | -0.01412 | [-0.0347, 0.002797] | 0.044 | -2.046 | 0.233 / 0.231 |
| Disordered conditional | weighted | -0.01405 | [-0.02087, -0.008415] | 0.000 | 8.858 | 0.222 / 0.205 |
| No-disorder marginal control | unweighted | 0.006992 | [-0.01273, 0.02095] | 0.674 | -2.695 | 0.151 / 0.151 |
| No-disorder marginal control | weighted | 0.002274 | [-0.003305, 0.006783] | 0.740 | -2.636 | 0.158 / 0.154 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.007791 | [-0.01788, 0.03288] | 0.715 | -2.134 | 0.146 / 0.145 |
| 1/16 | -0.005185 | [-0.03628, 0.0254] | 0.371 | -1.631 | 0.060 / 0.070 |
| 1/8 | -0.001544 | [-0.02703, 0.02346] | 0.442 | -2.427 | 0.032 / 0.059 |
| 1/4 | 0.03264 | [0.00397, 0.06045] | 0.988 | -0.913 | 0.311 / 0.358 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | -0.0007453 | [-0.02935, 0.02405] | -2.450 | 0.014 / 0.024 |
| 1/16 | 0.03913 | [-0.009468, 0.0835] | 1.816 | 0.260 / 0.243 |
| 1/8 | 0.00119 | [-0.02849, 0.02498] | -2.402 | 0.029 / 0.029 |
| 1/4 | -0.0116 | [-0.04328, 0.01448] | -1.848 | 0.150 / 0.147 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | -0.01448 | [-0.02135, -0.008921] | 0.000 | 16.888 | 0.222 / 0.206 |
| Disorder covariance | 0.00337 | [-0.001555, 0.008088] | 0.909 | -1.992 | 0.165 / 0.161 |
| No-disorder marginal | 0.00243 | [-0.003246, 0.006685] | 0.758 | -2.972 | 0.158 / 0.155 |
| No-disorder replica covariance | -0.000663 | [-0.005615, 0.004347] | 0.391 | -3.804 | 0.062 / 0.067 |
