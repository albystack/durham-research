# Background And Write-Up Outline

## Working Title

**Winding Variance of Loop-Erased Random Walks in Fixed Random Environments**

## Abstract Draft Structure

- Introduce loop-erased random walks and their connection to uniform spanning trees.
- State the known symmetric behaviour: winding variance grows like `C log L`.
- Introduce the fixed site-dependent random environment.
- State the numerical question: whether disorder changes the growth to `C(log L)^2`.
- Summarise the simulation method and the evidence obtained.

## 1. Introduction

Purpose:

- Explain why USTs and LERW are natural objects to study.
- Explain why winding is an interesting observable.
- Motivate the random-environment extension from Prof. Chhita's suggested project.

Key points to include:

- Wilson's algorithm samples USTs using loop-erased random walks.
- A UST branch has the same law as a LERW path.
- Therefore we can study the branch from the centre to the boundary by simulating one walk and loop-erasing it.
- The ordinary symmetric case has known logarithmic winding variance.
- The random-environment case appears to be open and suitable for numerical exploration.

## 2. Mathematical Background

### 2.1 Square Lattice And Finite Boxes

Define:

- `Z^2`
- nearest-neighbour edges
- finite box `B_L`
- boundary `partial B_L`
- origin/centre starting point

### 2.2 Random Walks

Define:

- simple symmetric random walk
- transition probabilities
- hitting the boundary
- raw path length

### 2.3 Loop Erasure

Define:

- chronological loop erasure
- self-avoiding path
- loop-erased path length

Implementation note:

- The loop-erasure algorithm should keep the most recent index of each visited vertex and delete loops when a vertex is revisited.

### 2.4 Uniform Spanning Trees And Wilson's Algorithm

Explain:

- spanning tree
- uniform spanning tree
- Wilson's algorithm
- why a UST branch can be sampled as a LERW

Keep this section focused; the simulations do not need to generate the whole tree.

### 2.5 Winding

Define:

```text
W_L = number of left turns - number of right turns
```

Implementation convention:

- Encode directions cyclically, for example `E = 0`, `N = 1`, `W = 2`, `S = 3`.
- Difference `+1 mod 4` means left turn.
- Difference `-1 mod 4` or `3 mod 4` means right turn.
- Difference `0` means straight.

## 3. Known Symmetric Case

State the known or expected behaviour:

```text
Var(W_L) ~ C log L
```

Things to research/cite:

- Schramm's scaling limit result for LERW/SLE_2.
- Winding variance results for SLE or LERW.
- Numerical or physics literature on winding angle variance.

This section should be used to justify the symmetric baseline simulation.

## 4. Random-Environment Model

### 4.1 Quenched Environment

Define the environment as fixed in time:

- Sample random variables at each site once.
- During a walk, every revisit to the same site uses the same local transition probabilities.

Contrast with annealed randomness:

- Annealed would resample at each step.
- That is not the model for this project.

### 4.2 Edge-Weight Definition

At each site:

```text
w_N = 1
w_E = 1
w_S = u_x
w_W = v_x
```

Normalise:

```text
Z_x = 2 + u_x + v_x
p_N = 1 / Z_x
p_E = 1 / Z_x
p_S = u_x / Z_x
p_W = v_x / Z_x
```

### 4.3 Disorder Distributions

Start with:

```text
u_x, v_x ~ Gamma(shape = k, scale = 1/k)
```

Then:

```text
E[u_x] = E[v_x] = 1
Var(u_x) = Var(v_x) = 1/k
```

Interpretation:

- Large `k`: weak disorder.
- Small `k`: strong disorder.
- `k = 1`: exponential distribution with mean `1`.

## 5. Numerical Method

For each pair `(L, k)`:

1. Generate a fixed random environment.
2. Run a random walk from the origin to the boundary.
3. Loop-erase the path.
4. Compute the winding.
5. Repeat for many samples.
6. Estimate the mean and variance of winding.

Record:

- seed
- box size
- distribution and parameter
- environment id
- walk id
- raw walk length
- loop-erased path length
- boundary hit point
- winding

## 6. Scaling Analysis

Fit and compare:

```text
Var(W_L) = a log L + b
Var(W_L) = a (log L)^2 + b
```

Use:

- scatter plots
- linear fits
- residual plots
- comparison after dropping the smallest `L`
- confidence intervals or bootstrap error bars if time permits

## 7. Optional Variance Decomposition

If using multiple walks per environment, apply:

```text
Var(W) = E[Var(W | environment)] + Var(E[W | environment])
```

This can help determine whether anomalous scaling comes from:

- random walk noise within one fixed environment
- variation between different fixed environments

## 8. Results Section Template

Suggested order:

1. Correctness checks.
2. Symmetric baseline.
3. Random environment with `k = 1`.
4. Disorder sweep over several `k`.
5. Scaling comparison.
6. Variance decomposition if available.

## 9. Discussion

Address:

- Whether the numerical evidence favours `log L` or `(log L)^2`.
- How strong the evidence is.
- Finite-size limitations.
- Runtime limitations.
- Whether stronger disorder changes the apparent scaling.
- How the results compare with Prof. Chhita's expectation.

## 10. Further Work

Possible extensions:

- Larger box sizes.
- More samples.
- Better variance decomposition.
- Triangular lattice extension.
- Drifted walks as a contrast case.
- Other mean-one positive distributions.

