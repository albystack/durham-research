# Random walks and dimers in random environments

This repository contains two computational probability projects supervised by
Prof. Sunil Chhita at Durham University:

1. **Loop-erased random walks (LERWs):** does spatial disorder change winding
   fluctuations from \(C\log L\) to \(C(\log L)^2\)?
2. **Random-weight Aztec diamonds:** how does the height at the central face
   fluctuate under Gamma-disordered domino weights?

Both studies use reproducible Julia simulations, deterministic seeds,
bootstrap inference and direct model comparison. The numerical conclusions
are finite-size evidence, not asymptotic proofs.

## Main results

| Study | Scale | Main numerical result |
|---|---:|---|
| LERW, one walk per fixed environment | 379,000 walks; 15 weight laws; \(L\le8192\) | Fitted winding exponents \(p=1.01\)--\(1.16\); logarithmic model preferred by BIC for 14/15 specifications |
| LERW, two walks sharing each environment | 758,000 walks in 379,000 pairs | Fitted exponents \(p=1.04\)--\(1.19\); logarithmic model preferred for 12/15 specifications |
| Temporally refreshed Gamma LERW | \(L\le8192\) | Recovered the length law \(E[|\mathrm{LERW}_L|]\propto L^{1.2518}\), with 95% bootstrap interval \(1.2456\)--\(1.2573\) |
| Gamma-disordered Aztec diamond | 26,050 tilings; \(L\le600\) | The affine \((\log L)^2\) height-variance fit beats the affine \(\log L\) fit by 5.260 BIC units and wins 94.2% of bootstrap resamples |

The Aztec result is promising but not decisive: the 95% bootstrap interval for
the BIC difference is \([-1.273,6.691]\), and the two fitted curves remain
close over the simulated range.

## Repository layout

```text
.
├── random_walk/   LERW package, experiment configs and compact results
├── aztec/         weighted domino sampler, datasets and results
└── .github/       CI for both Julia test suites
```

Each study is self-contained:

- [`random_walk/README.md`](random_walk/README.md) explains the LERW models,
  retained results and campaign commands.
- [`aztec/README.md`](aztec/README.md) explains the domino sampler, height
  convention, retained datasets and reproduction commands.

## Quick start

Julia 1.10 or newer is recommended.

### Random-walk tests

```bash
cd random_walk
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

### Aztec-diamond tests and saved-data analysis

From the repository root:

```bash
julia --startup-file=no aztec/test/runtests.jl
julia --startup-file=no aztec/scripts/analyze_height_campaign.jl
julia --startup-file=no aztec/scripts/plot_height_campaign.jl
```

Generated batches and scratch analyses are written below `output/` and are
ignored by Git. The repository tracks source code, tests, frozen
configurations, compact datasets and the result files needed to inspect the
reported conclusions.

## Selected figures

### Paired LERW winding variance

<p align="center">
  <img src="random_walk/results/site_iid_paired/figures/annealed_scaling_all_distributions.png" alt="Paired LERW winding-difference variance across random environments" width="92%">
</p>

### Gamma-disordered Aztec center-height variance

<p align="center">
  <img src="aztec/results/height/center_height_variance.png" alt="Central-height variance for Gamma-disordered Aztec diamonds" width="82%">
</p>
