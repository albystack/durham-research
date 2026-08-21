# Spatial block-GLS campaign comparison

The table below selects the disorder covariance from each analysis.
Intervals and model comparisons come from joint environment bootstraps.

| Campaign | min-order tag | c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / log2 |
|---|---:|---:|---:|---:|---:|---:|
| original_min128 | 128 | 0.5019 | [0.2022, 0.807] | 0.9994 | 6.928 | 1.034 / 0.614 |
| strong_min128 | 128 | 0.6077 | [0.2723, 0.9425] | 0.9999 | 9.093 | 1.106 / 0.638 |
| weak_min128 | 128 | 0.08045 | [-0.1837, 0.3461] | 0.7224 | -3.227 | 0.785 / 1.121 |
| symmetric_min128 | 128 | 0.2559 | [-0.04009, 0.5547] | 0.9546 | -0.774 | 1.155 / 1.349 |
| original_min384 | 384 | 0.9503 | [-0.6456, 2.519] | 0.8730 | -1.814 | 0.774 / 0.718 |
| strong_min384 | 384 | 0.9179 | [-0.8368, 2.698] | 0.8458 | -2.144 | 0.539 / 0.408 |
| weak_min384 | 384 | -1.111 | [-2.442, 0.2395] | 0.0537 | -0.509 | 0.990 / 1.205 |
| symmetric_min384 | 384 | -1.608 | [-3.143, -0.04565] | 0.0220 | 0.985 | 1.544 / 1.501 |
| original_min512 | 512 | 1.697 | [-1.327, 4.691] | 0.8633 | -1.759 | 0.937 / 1.161 |
| strong_min512 | 512 | 0.3263 | [-3.003, 3.573] | 0.5724 | -2.958 | 1.073 / 1.205 |
| weak_min512 | 512 | -1.464 | [-3.929, 1.017] | 0.1272 | -1.636 | 1.029 / 2.036 |
| symmetric_min512 | 512 | -0.9636 | [-3.735, 1.864] | 0.2480 | -2.544 | 1.309 / 1.285 |

These are finite-size numerical comparisons, not asymptotic proofs.
