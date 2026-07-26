# Winding fluctuations of loop-erased random walks

A computational study of loop-erased random walks (LERWs) in random
environments. The central question is whether local disorder changes the
winding variance from logarithmic growth,

```math
\operatorname{Var}(W_L) \sim C\log L,
```

to a squared-logarithmic regime,

```math
\operatorname{Var}(W_L) \sim C(\log L)^2.
```

## Model

A walk starts at the origin and stops when it reaches the boundary of
the square \([-L,L]^2\). At each lattice site \(x\), four independent positive
weights \(w_{x,d}\) determine the next step:

```math
\mathbb P(X_{n+1}=x+e_d\mid X_n=x,\omega)
=\frac{w_{x,d}}{\sum_{d'}w_{x,d'}}.
```

Loops are erased chronologically. Along the resulting self-avoiding path,
\(W_L\) is the number of left quarter-turns minus right quarter-turns; the
angular winding is \(\Theta_L=(\pi/2)W_L\).

The simulations cover 15 bounded, light-tailed, and heavy-tailed weight
specifications. Two observables are studied:

- one LERW in each independent environment;
- a pair of conditionally independent LERWs in the same environment, with
  \(\Delta W_L=W_L^{(1)}-W_L^{(2)}\).

The implemented environment is directed and site-i.i.d. This repository
currently studies LERWs only; it does not yet model domino tilings or Aztec
diamonds.

## Results

The effective exponent is fitted through

```math
\operatorname{Var}(W_L)=C(\log L)^p,
\qquad
\log\operatorname{Var}(W_L)=\log C+p\log\log L.
```

Thus \(p=1\) represents logarithmic growth and \(p=2\) squared-logarithmic
growth.

| Experiment | Sample | Fitted \(p\) across 15 specifications | Log model preferred by BIC |
|---|---:|---:|---:|
| Single walk | 379,000 environments and walks | 1.011–1.160 | 14/15 |
| Paired walks | 379,000 environments, 758,000 walks | 1.036–1.187 | 12/15 |

Every bootstrap confidence interval remains far below \(p=2\). Over the
simulated range \(L=16,\ldots,8192\), both experiments therefore show no robust
evidence for squared-logarithmic growth. These are finite-size Monte Carlo
results, not an asymptotic proof.

The tables and figures in [`reports/`](reports/) show the completed paired-walk
experiment. Raw batch files and the earlier single-walk results are preserved
in [`legacy/`](legacy/).

### Paired winding variance

Points are variance estimates with 95% error bars. The blue and orange curves
fit \(a+b\log L\) and \(a+b(\log L)^2\), respectively.

<p align="center">
  <img src="reports/figures/annealed_scaling_all_distributions.png" alt="Paired LERW winding-difference variance across all weight distributions" width="100%">
</p>

### Effective exponents

The green line marks \(p=1\); the dashed red line marks \(p=2\).

<p align="center">
  <img src="reports/figures/scaling_exponent_forest.png" alt="Effective exponent estimates for paired LERW winding differences" width="82%">
</p>

### Direct model comparison

Positive \(\Delta\mathrm{BIC}=\mathrm{BIC}(\log^2)-\mathrm{BIC}(\log)\)
favors ordinary logarithmic growth.

<p align="center">
  <img src="reports/figures/bic_model_comparison.png" alt="BIC comparison of logarithmic and squared-logarithmic fits" width="82%">
</p>

## Setup

Julia 1.10 or newer is required.

```bash
git clone https://github.com/albystack/lerw-random-environment-research.git
cd lerw-random-environment-research/julia

julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

The simulation, analysis, and plotting entry points are in [`julia/scripts/`](julia/scripts/).
