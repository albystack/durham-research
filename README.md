# Gamma-disordered Aztec diamonds

Reproducible Julia simulations for central-height fluctuations in the
Gamma-disordered Aztec diamond. The repository compares ordinary logarithmic
growth with the super-rough squared-logarithmic prediction and separates
conditional tiling noise from disorder-induced fluctuations using paired
double-dimer samples.

The active project is [`aztec/`](aztec/). Earlier loop-erased-random-walk work
is preserved under [`old/`](old/) and is not required by the Aztec code.

## Quick start

Requirements: Julia 1.10 or newer. Only Julia standard libraries are used.

```bash
# Run all unit tests.
julia --project=aztec -e 'using Pkg; Pkg.test()'

# Reanalyse the retained single-height data.
julia --project=aztec aztec/scripts/analyze_height_campaign.jl

# Reanalyse the retained shared-environment pairs.
julia --project=aztec aztec/scripts/analyze_double_dimer_campaign.jl
```

Generated batches and scratch analyses go to the ignored `aztec/output/`
directory. Retained data, checksums, tables, and figures live in
`aztec/data/` and `aztec/results/`.

## Retained experiment

- 35,536 independent single-height samples through order 1,300.
- 28,304 independent environments with two conditional tilings per
  environment, also through order 1,300.
- Deterministic per-sample seeds, atomic resumable batches, and 10,000-replicate
  within-size bootstrap analyses.

The total single-height data do not decisively distinguish the two affine
curves. The double-dimer difference favours ordinary log over the full fitted
range. The paired covariance, which directly estimates the disorder-induced
variance, has suggestive upward curvature but remains noisy at the largest
orders.

See the [project guide](aztec/README.md) for the model and commands, and the
[analysis report](aztec/results/analysis_report.md) for numerical results and
finite-size caveats.
