# Spatial height-increment analysis

Independent disordered environments: **39570**.
Independent no-disorder replica pairs: **39570**.
Bootstrap repetitions: **10000**; fitted orders start at **128**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Disorder covariance | unweighted | -0.005109 | [-0.01986, 0.008955] | 0.224 | -2.590 | 0.088 / 0.090 |
| Disorder covariance | weighted | -0.001901 | [-0.007354, 0.00333] | 0.227 | -3.117 | 0.096 / 0.093 |
| Disordered conditional | unweighted | 0.007199 | [-0.01284, 0.02318] | 0.683 | -1.873 | 0.107 / 0.107 |
| Disordered conditional | weighted | 0.001816 | [-0.004504, 0.007236] | 0.663 | -3.249 | 0.106 / 0.102 |
| No-disorder marginal control | unweighted | 0.006992 | [-0.01256, 0.021] | 0.674 | -2.695 | 0.151 / 0.151 |
| No-disorder marginal control | weighted | 0.002241 | [-0.003509, 0.00673] | 0.730 | -2.677 | 0.158 / 0.154 |

## Primary disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | -0.01245 | [-0.04109, 0.01393] | 0.180 | -0.324 | 0.089 / 0.089 |
| 1/16 | 0.008882 | [-0.01928, 0.03802] | 0.707 | -0.450 | 0.067 / 0.061 |
| 1/8 | 0.003398 | [-0.0292, 0.03735] | 0.550 | -2.367 | 0.065 / 0.123 |
| 1/4 | -0.02027 | [-0.05149, 0.01218] | 0.103 | 1.906 | 0.121 / 0.108 |

## No-disorder control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | -0.0007453 | [-0.02907, 0.02327] | -2.450 | 0.014 / 0.024 |
| 1/16 | 0.03913 | [-0.008367, 0.08399] | 1.816 | 0.260 / 0.243 |
| 1/8 | 0.00119 | [-0.02779, 0.02516] | -2.402 | 0.029 / 0.029 |
| 1/4 | -0.0116 | [-0.04251, 0.0152] | -1.848 | 0.150 / 0.147 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.

## Covariance-aware pooled sensitivity

Block GLS uses one bootstrap-estimated cross-fraction covariance matrix per order. Independent orders form separate covariance blocks; eigenvalues below `1e-8` of the largest block eigenvalue are regularized.

| observable | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---:|---:|---:|---:|---:|
| Disordered conditional | 0.002163 | [-0.004233, 0.007416] | 0.706 | -3.346 | 0.106 / 0.102 |
| Disorder covariance | -0.002513 | [-0.007742, 0.002546] | 0.166 | -2.977 | 0.096 / 0.092 |
| No-disorder marginal | 0.002611 | [-0.003061, 0.006926] | 0.780 | -2.826 | 0.158 / 0.154 |
| No-disorder replica covariance | -0.0006546 | [-0.005699, 0.004442] | 0.397 | -3.806 | 0.062 / 0.067 |
