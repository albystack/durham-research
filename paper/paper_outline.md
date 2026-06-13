# Paper Outline

## Working Title

**Numerical Experiments on Winding of Loop-Erased Random Walks in Fixed Random Environments**

## Abstract Draft

We study the winding of loop-erased random walks from the origin to the boundary of a square box in two-dimensional lattice random environments. In the symmetric case, the winding variance is expected to grow logarithmically with the box size. Motivated by a proposed random-environment extension of the uniform spanning tree model, we implement fixed site-dependent transition probabilities generated from positive edge weights and compare the observed winding variance against `C log L` and `C(log L)^2` scaling forms. Our simulations reproduce the expected qualitative logarithmic growth in the symmetric baseline. For the two-random-weight model `w_N = w_E = 1`, `w_S = u`, `w_W = v`, we find substantial effective north-east drift after normalisation, especially for moderate and strong disorder. This drift shortens walks and suppresses winding variance growth. As a balanced diagnostic, we also simulate a four-random-weight model with independent mean-one Gamma weights in all directions; this model remains centred and shows variance growth broadly comparable to the symmetric baseline for `L <= 256`. These results suggest that drift diagnostics are essential before interpreting random-environment winding simulations as evidence for anomalous scaling.

## Core Thesis

The numerical project currently supports a careful methodological conclusion:

> The proposed random-environment simulation is highly sensitive to how the local weights are normalised. A mean-one edge-weight rule does not automatically produce a drift-free transition kernel.

The paper should present this as a useful finding, not as a proof or disproof of the conjectured `(log L)^2` winding variance.

## Proposed Structure

### 1. Introduction

Goals:

- Introduce USTs and LERW.
- Explain the winding observable.
- State the known symmetric expectation `Var(W_L) ~ C log L`.
- Motivate the fixed random-environment extension.
- State the numerical question: can random environments produce stronger, possibly `(log L)^2`, winding variance?

### 2. Background

Include:

- square lattice boxes `B_L`
- random walks and boundary hitting
- loop erasure
- Wilson's algorithm and the UST/LERW connection
- winding as left turns minus right turns

### 3. Random-Environment Models

#### 3.1 Two-Weight Model

Define:

```text
w_N = 1
w_E = 1
w_S = u_x
w_W = v_x
```

with:

```text
u_x, v_x ~ Gamma(k, 1/k)
```

and transition probabilities given by normalisation.

#### 3.2 Balanced Four-Weight Diagnostic

Define:

```text
w_N, w_E, w_S, w_W iid Gamma(k, 1/k)
```

Explain that this is a diagnostic comparison model, not the primary supervisor-suggested model.

### 4. Simulation Method

Include:

- box sizes `L = 16, 32, 64, 128, 256`
- fixed quenched environment
- one walk from origin to boundary
- chronological loop erasure
- winding computation
- sample counts
- random seeds
- data recorded per trial

Mention that the focused batch used `700` samples per model/size/disorder cell.

### 5. Correctness And Baseline Checks

Use:

- symmetric baseline variance table
- symmetric variance vs `log L` figure
- path-length sanity checks

Main point:

```text
The symmetric baseline is consistent with logarithmic growth and validates the pipeline qualitatively.
```

### 6. Results: Two-Weight Gamma Model

Use figures:

- `mean_hit_x_over_L.svg`
- `mean_hit_y_over_L.svg`
- `boundary_hits_L128.svg`
- `gamma_variance_by_k_logL.svg`
- `raw_walk_length.svg`

Main claims:

- the model has strong north-east drift for `k = 1` and `k = 0.5`
- drift grows with box size
- raw walk lengths become much shorter than symmetric/balanced cases
- winding variance flattens for stronger disorder

### 7. Results: Balanced Gamma4 Diagnostic

Use figures:

- `gamma4_variance_by_k_logL.svg`
- `variance_vs_logL_selected.svg`
- `variance_vs_logL2_selected.svg`
- `disorder_sweep_L256.svg`

Main claims:

- boundary-hit means stay close to zero
- winding variance grows with `L`
- no strong evidence for `(log L)^2` growth up to `L = 256`
- disorder strength has weaker effect than expected over this size range

### 8. Discussion

Discuss:

- why mean-one weights are not sufficient to remove drift
- why the two-weight model may not address the intended no-drift open problem
- why the balanced model is useful but should be confirmed with the supervisor
- finite-size limitations
- sample-size limitations
- need for variance decomposition

### 9. Conclusion

Possible conclusion:

> The simulations identify effective drift as a central modelling issue in fixed random-environment LERW winding experiments. The exact two-weight model produces strong directional exit bias under moderate and strong Gamma disorder, suppressing winding variance. A balanced four-weight diagnostic remains centred and displays variance growth close to the symmetric baseline over the tested size range. Further work should clarify the intended drift-free random-environment model before larger-scale simulations are used to test the conjectured `(log L)^2` law.

## Figures To Use First

From `analysis_focused/figures`:

1. `variance_vs_logL_selected.svg`
2. `variance_vs_logL2_selected.svg`
3. `mean_hit_x_over_L.svg`
4. `mean_hit_y_over_L.svg`
5. `boundary_hits_L128.svg`
6. `gamma_variance_by_k_logL.svg`
7. `gamma4_variance_by_k_logL.svg`
8. `raw_walk_length.svg`
9. `lerw_path_length.svg`
10. `sample_lerw_paths.svg`
11. `environment_vector_fields.svg`

## Tables To Use First

From `analysis_focused/tables`:

- `group_stats.csv`
- `scaling_fits.csv`

## Claims We Can Defend

- The code reproduces the expected qualitative symmetric baseline.
- The two-weight model is strongly drifted for moderate and strong disorder.
- Strong drift suppresses winding variance growth.
- The balanced diagnostic model remains centred.
- The balanced diagnostic does not provide compelling `(log L)^2` evidence up to `L = 256`.

## Claims To Avoid

- Avoid saying the conjecture is false.
- Avoid saying the random-environment problem is solved.
- Avoid saying the two-weight model is drift-free.
- Avoid treating `gamma4` as the official model without asking Prof. Chhita.

