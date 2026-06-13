# Project Brief: Fixed Random Environments for LERW Winding

## One-Sentence Aim

Simulate loop-erased random walks from the centre of a finite square box to the boundary in a fixed site-dependent random environment, then estimate how the winding variance scales with box size.

## Supervisor Guidance From `outlook.txt`

The project should focus on **Extension 3**, especially the case where random transition probabilities are sampled once at each site and then kept fixed throughout the simulation.

Prof. Chhita's latest modelling advice is to avoid sampling transition probabilities directly. Instead, assign positive edge weights at each site and normalise them into probabilities.

For a site with compass directions `N, E, S, W`, use:

```text
weight(N) = 1
weight(E) = 1
weight(S) = u
weight(W) = v
```

where `u` and `v` are sampled from a positive distribution with mean `1`, such as an exponential or Gamma distribution.

The transition probabilities are:

```text
p_N = 1 / (2 + u + v)
p_E = 1 / (2 + u + v)
p_S = u / (2 + u + v)
p_W = v / (2 + u + v)
```

This is the model we should implement first.

## Main Observable

The observable is the winding of the loop-erased path from the origin to the boundary:

```text
winding = number of left turns - number of right turns
```

We should focus on the path from the origin rather than generating the whole spanning tree, because this is computationally cheaper and still targets the branch of interest.

## Main Scaling Question

For the ordinary symmetric case, where each step has probability `1/4`, the winding variance is expected to scale like:

```text
Var(W_L) ~ C log L
```

The open/numerical question is whether the fixed random-environment case instead scales like:

```text
Var(W_L) ~ C (log L)^2
```

The simulation should compare these two candidate behaviours over increasing box sizes.

## Important Distinctions

- The environment is **quenched**: sampled once, then fixed during the walk.
- The model should use **positive edge weights**, not direct normal samples for probabilities.
- The random weights should have mean `1`, matching the latest email guidance. This balances the raw edge-weight scale, but the normalised transition probabilities should still be checked empirically for effective drift.
- If the environment has a genuine non-zero drift, the winding variance may become finite rather than growing with `L`.
- Once the walk is no longer symmetric, it is better to describe the simulation as **LERW in a random environment** or a **weighted random-walk spanning tree model**, rather than simply as the ordinary UST.

## Early Modelling Caveat

The exact two-random-weight model suggested in the email,

```text
w_N = 1
w_E = 1
w_S = u
w_W = v
```

can show an effective drift after normalisation, even when `E[u] = E[v] = 1`. For that reason, the simulator also includes an optional diagnostic model called `gamma4`, where all four directions receive independent mean-one Gamma weights. This is not a replacement for Prof. Chhita's suggested model, but it is useful as a rotation-symmetric comparison case.

## Initial Implementation Target

Start with the smallest useful experiment:

```text
L = 16, 32, 64
distribution = Gamma(shape = 1, scale = 1)
samples = 50 to 100 per L
```

This is not intended to answer the research question conclusively. Its purpose is to check correctness, runtime, output format, and plotting.

## First Evidence To Produce

1. A symmetric baseline plot of winding variance versus `log L`.
2. A random-environment plot of winding variance versus `log L`.
3. A random-environment plot of winding variance versus `(log L)^2`.
4. A small table containing `L`, sample count, mean winding, variance, mean raw walk length, and mean loop-erased path length.
