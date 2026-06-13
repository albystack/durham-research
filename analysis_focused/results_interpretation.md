# Results Interpretation: Focused Batch

Source data: `results/research_batch_focused.csv`

This focused batch used:

- box sizes `L = 16, 32, 64, 128, 256`
- `700` samples for each model/size/disorder configuration
- models: `symmetric`, two-weight `gamma`, balanced `gamma4`
- disorder levels: `k = 20, 1, 0.5`

The purpose of this run was not to prove an asymptotic law, but to get a cleaner numerical picture after the exploratory all-k sweep.

## 1. Symmetric Baseline

The symmetric walk behaves as expected and is a useful calibration.

| L | Var(W) | mean raw steps | mean LERW steps |
|---:|---:|---:|---:|
| 16 | 1.894 | 296.9 | 40.3 |
| 32 | 2.353 | 1239.1 | 96.2 |
| 64 | 3.148 | 4861.2 | 231.8 |
| 128 | 3.572 | 20045.7 | 559.1 |
| 256 | 4.300 | 77271.2 | 1329.5 |

The linear fit against `log L` gives:

```text
Var(W) = 0.870 log L - 0.565,   R^2 = 0.992
```

The fit against `(log L)^2` is also visually reasonable over only five sizes, but slightly worse:

```text
Var(W) = 0.104 (log L)^2 + 1.158,   R^2 = 0.987
```

Interpretation: the baseline is consistent with the expected `C log L` law, and the simulation/winding code is behaving sensibly.

## 2. Two-Weight Gamma Model

This is the model closest to Prof. Chhita's latest email:

```text
w_N = 1
w_E = 1
w_S = u
w_W = v
```

with `u, v ~ Gamma(k, 1/k)`.

The key finding is that this model develops a strong north-east effective drift after normalisation, especially for smaller `k`.

### Drift Diagnostic

At `L = 256`, the mean boundary hit locations are:

| model | k | mean hit x/L | mean hit y/L |
|---|---:|---:|---:|
| gamma | 20 | 0.589 | 0.591 |
| gamma | 1 | 0.894 | 0.889 |
| gamma | 0.5 | 0.916 | 0.918 |

This is not a small finite-sample effect. The walks are overwhelmingly exiting near the north-east corner for `k = 1` and `k = 0.5`.

### Winding Variance

| k | Var(W), L=16 | Var(W), L=256 | logL fit R^2 | logL^2 fit R^2 |
|---:|---:|---:|---:|---:|
| 20 | 1.973 | 3.735 | 0.952 | 0.908 |
| 1 | 1.755 | 2.099 | 0.524 | 0.587 |
| 0.5 | 1.358 | 1.548 | 0.950 | 0.916 |

Interpretation:

- `k = 20` is weak disorder, so the drift is less extreme and the variance still grows with `L`.
- `k = 1` and `k = 0.5` are dominated by effective drift, giving much shorter raw walks and nearly flat winding variance.
- These drifted cases are not good evidence for either `C log L` or `C(log L)^2`; they are closer to the finite-variance drift scenario Prof. Chhita warned about.

## 3. Balanced Gamma4 Diagnostic

The balanced comparison model uses four independent mean-one Gamma weights:

```text
w_N, w_E, w_S, w_W iid Gamma(k, 1/k)
```

This is not the exact supervisor-suggested model, but it is useful because it removes the structural north-east preference.

### Drift Diagnostic

At `L = 256`, the mean boundary hit locations are:

| model | k | mean hit x/L | mean hit y/L |
|---|---:|---:|---:|
| gamma4 | 20 | 0.017 | -0.027 |
| gamma4 | 1 | 0.028 | 0.026 |
| gamma4 | 0.5 | 0.017 | -0.003 |

This is well balanced at the sample sizes tested.

### Winding Variance

| k | Var(W), L=16 | Var(W), L=256 | logL fit R^2 | logL^2 fit R^2 |
|---:|---:|---:|---:|---:|
| 20 | 2.034 | 4.104 | 0.963 | 0.938 |
| 1 | 1.987 | 4.001 | 0.978 | 0.943 |
| 0.5 | 1.965 | 4.184 | 0.976 | 0.946 |

Interpretation:

- The balanced random environment shows clear growth of winding variance with `L`.
- Over `L <= 256`, the data are still closer to `C log L` than `C(log L)^2`.
- There is no convincing evidence yet for anomalous `(log L)^2` growth in this balanced diagnostic model.
- The disorder strength does not dramatically increase the variance at this size range; all three `k` values end near variance `4.0-4.2` at `L=256`.

## 4. Main Scientific Takeaway

The most important result is methodological:

> The two-weight edge model from the email, although it uses mean-one random weights, appears to create strong effective drift after probability normalisation.

This means the exact model as currently implemented may not be testing the intended "mean-zero random environment" winding question for moderate or strong disorder. Instead, for `k = 1` and `k = 0.5`, it behaves like a drifted walk with flattened winding variance.

The balanced `gamma4` model is a useful control. It stays centred and gives winding variance growth similar to the symmetric baseline over the available sizes.

## 5. What We Can Safely Claim

Safe claims:

- The simulation pipeline reproduces the expected qualitative behaviour of the symmetric case.
- The two-weight model has strong boundary-hit drift for moderate/strong disorder.
- Drift correlates with shorter raw walks and suppressed winding variance.
- The balanced four-weight diagnostic does not show clear `(log L)^2` growth up to `L=256`.
- Larger sizes and more samples are required before making asymptotic claims.

Claims to avoid:

- Do not claim that the open problem is resolved.
- Do not claim that random environments never produce `(log L)^2` growth.
- Do not describe the two-weight model as drift-free merely because the random weights have mean one.

## 6. Recommended Next Scientific Step

Before spending compute on much larger runs, ask Prof. Chhita to confirm the intended no-drift model.

Suggested question:

> In the edge-weight model with `w_N = w_E = 1`, `w_S = u`, `w_W = v`, I observe strong north-east boundary drift after normalisation, especially for `u, v` drawn from mean-one Gamma distributions. Should the intended random-environment model instead randomise all four outgoing weights symmetrically, or impose a local balance condition to remove this effective drift?

If the two-weight model is still desired, the paper should frame the main result as a drift diagnostic and finite-variance observation, not as evidence for the conjectured `(log L)^2` law.

