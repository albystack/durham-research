# Numerical Analysis Dossier

Source data: `results/research_batch_main.csv`

## Dataset

- Models: gamma, gamma4, symmetric
- Box sizes: 16, 32, 64, 128, 256
- Total grouped cells: 65

## Main Diagnostics

- `symmetric` k=None at L=256: Var(W)=4.171, mean W=0.053, mean hit/L=(-0.055, 0.010), mean raw steps=78742.3.
- `gamma` k=1.0 at L=256: Var(W)=1.889, mean W=0.017, mean hit/L=(0.884, 0.888), mean raw steps=4665.9.
- `gamma4` k=1.0 at L=256: Var(W)=4.270, mean W=0.023, mean hit/L=(-0.020, 0.015), mean raw steps=93050.7.

The two-random-weight `gamma` model should be treated carefully: the boundary-hit diagnostics measure whether the normalised probabilities create an effective drift. The `gamma4` model is a balanced comparison, not a replacement for the supervisor-suggested model.

## Scaling Fits

The fits below are linear regressions of winding variance against either `log L` or `(log L)^2`. They are descriptive diagnostics, not proof of asymptotic scaling.

| model | k | fit | slope | intercept | R^2 | RSS |
|---|---:|---|---:|---:|---:|---:|
| gamma | 0.5 | logL | 0.05961 | 1.225 | 0.1951 | 0.07044 |
| gamma | 0.5 | logL_squared | 0.007688 | 1.333 | 0.2267 | 0.06767 |
| gamma | 1.0 | logL | 0.07363 | 1.695 | 0.1795 | 0.1191 |
| gamma | 1.0 | logL_squared | 0.007062 | 1.872 | 0.1153 | 0.1284 |
| gamma | 10.0 | logL | 0.537 | 0.7595 | 0.8110 | 0.3227 |
| gamma | 10.0 | logL_squared | 0.06154 | 1.869 | 0.7443 | 0.4368 |
| gamma | 2.0 | logL | 0.1578 | 1.457 | 0.3970 | 0.1818 |
| gamma | 2.0 | logL_squared | 0.01827 | 1.78 | 0.3715 | 0.1894 |
| gamma | 20.0 | logL | 0.6285 | 0.1219 | 0.9542 | 0.09111 |
| gamma | 20.0 | logL_squared | 0.07419 | 1.381 | 0.9286 | 0.142 |
| gamma | 5.0 | logL | 0.4699 | 0.7657 | 0.7729 | 0.3117 |
| gamma | 5.0 | logL_squared | 0.05429 | 1.729 | 0.7207 | 0.3834 |
| gamma4 | 0.5 | logL | 0.9047 | -0.6575 | 0.9182 | 0.3502 |
| gamma4 | 0.5 | logL_squared | 0.1082 | 1.13 | 0.9168 | 0.3562 |
| gamma4 | 1.0 | logL | 0.6296 | 0.5611 | 0.9327 | 0.1374 |
| gamma4 | 1.0 | logL_squared | 0.07691 | 1.775 | 0.9723 | 0.05648 |
| gamma4 | 10.0 | logL | 0.8589 | -0.6653 | 0.8733 | 0.5142 |
| gamma4 | 10.0 | logL_squared | 0.1044 | 1 | 0.9016 | 0.3994 |
| gamma4 | 2.0 | logL | 0.9397 | -0.797 | 0.9319 | 0.3099 |
| gamma4 | 2.0 | logL_squared | 0.1134 | 1.041 | 0.9479 | 0.2374 |
| gamma4 | 20.0 | logL | 0.8259 | -0.4162 | 0.9953 | 0.01531 |
| gamma4 | 20.0 | logL_squared | 0.09877 | 1.215 | 0.9944 | 0.01841 |
| gamma4 | 5.0 | logL | 0.7354 | -0.06997 | 0.9742 | 0.06873 |
| gamma4 | 5.0 | logL_squared | 0.08857 | 1.371 | 0.9873 | 0.03395 |
| symmetric | None | logL | 0.9063 | -0.7222 | 0.9303 | 0.2956 |
| symmetric | None | logL_squared | 0.1073 | 1.088 | 0.9111 | 0.3772 |

## Figures

![variance_vs_logL_selected.svg](figures/variance_vs_logL_selected.svg)

![variance_vs_logL2_selected.svg](figures/variance_vs_logL2_selected.svg)

![gamma_variance_by_k_logL.svg](figures/gamma_variance_by_k_logL.svg)

![gamma4_variance_by_k_logL.svg](figures/gamma4_variance_by_k_logL.svg)

![mean_hit_x_over_L.svg](figures/mean_hit_x_over_L.svg)

![mean_hit_y_over_L.svg](figures/mean_hit_y_over_L.svg)

![raw_walk_length.svg](figures/raw_walk_length.svg)

![lerw_path_length.svg](figures/lerw_path_length.svg)

![disorder_sweep_L128.svg](figures/disorder_sweep_L128.svg)

![disorder_sweep_L256.svg](figures/disorder_sweep_L256.svg)

![boundary_hits_L128.svg](figures/boundary_hits_L128.svg)

![winding_histograms_L128.svg](figures/winding_histograms_L128.svg)

![sample_lerw_paths.svg](figures/sample_lerw_paths.svg)

![environment_vector_fields.svg](figures/environment_vector_fields.svg)

## Preliminary Interpretation

- The symmetric baseline is the calibration case and should show roughly logarithmic growth, subject to finite-size noise.
- A strong mean boundary-hit displacement indicates effective drift, in which case winding variance may flatten rather than grow like `log L` or `(log L)^2`.
- The balanced `gamma4` diagnostic helps distinguish disorder effects from drift effects.
- The current run should be used to decide which models and size ranges deserve larger sample counts before paper drafting.
