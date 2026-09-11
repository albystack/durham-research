# Spatial height-increment analysis

Independent Gamma environments: **6900**.
Independent uniform replica pairs: **5500**.
Bootstrap repetitions: **10000**; fitted orders start at **128**; the largest **2** orders test prediction.

## Joint environment-clustered result

The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.

| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---|---|---:|---:|---:|---:|---:|
| Gamma disorder covariance | unweighted | 0.5843 | [0.1957, 0.952] | 0.999 | 10.848 | 0.886 / 0.588 |
| Gamma disorder covariance | weighted | 0.4948 | [0.1932, 0.7977] | 0.999 | 14.569 | 1.026 / 0.608 |
| Gamma conditional | unweighted | -0.2013 | [-0.541, 0.1411] | 0.121 | -0.814 | 0.370 / 0.476 |
| Gamma conditional | weighted | -0.1873 | [-0.4866, 0.1116] | 0.106 | -0.572 | 0.475 / 0.492 |
| Uniform marginal control | unweighted | 0.00882 | [-0.2623, 0.2835] | 0.520 | -3.578 | 0.374 / 0.431 |
| Uniform marginal control | weighted | 0.06353 | [-0.1847, 0.316] | 0.694 | -3.172 | 0.415 / 0.406 |

## Primary Gamma disorder result

The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.

| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|---:|
| 1/32 | 0.5655 | [0.02696, 1.111] | 0.980 | 8.102 | 0.673 / 0.156 |
| 1/16 | 0.338 | [-0.3087, 0.9744] | 0.841 | 0.310 | 0.790 / 0.925 |
| 1/8 | 1.045 | [0.2372, 1.874] | 0.994 | 5.382 | 1.260 / 0.846 |
| 1/4 | 0.3891 | [-0.4762, 1.245] | 0.817 | -1.064 | 0.688 / 0.765 |

## Uniform control

| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |
|---:|---:|---:|---:|---:|
| 1/32 | 0.158 | [-0.2758, 0.5886] | -1.232 | 0.321 / 0.411 |
| 1/16 | -0.05249 | [-0.5536, 0.454] | -2.159 | 0.278 / 0.307 |
| 1/8 | 0.2162 | [-0.3546, 0.7886] | -1.333 | 0.538 / 0.512 |
| 1/4 | -0.2867 | [-0.8739, 0.3019] | -1.090 | 0.302 / 1.094 |

These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.

The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.
