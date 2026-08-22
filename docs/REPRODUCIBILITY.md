# Reproducibility

## Environment

The active package supports Julia 1.10 or newer and uses only Julia standard
libraries. A clean checkout can run the full suite with:

```bash
julia --project=aztec -e 'using Pkg; Pkg.test()'
sh aztec/test/smoke_workflows.sh
```

The CI workflow runs the package tests on Julia 1.10 and the latest stable
Julia release.

## Randomness and pairing

Campaign seeds are explicit CLI arguments. Per-observation seeds are derived
deterministically from the campaign seed, model parameters, lattice size,
environment identifier, and replica identifier. A restart validates existing
metadata and seeds before skipping a complete batch.

For shared-environment experiments, the frozen environment is the independent
statistical block. Both conditional replicas and every observable recorded
from that environment remain paired during bootstrap resampling.

## Output integrity

Runners write one atomic batch at a time. A batch is considered complete only
after its header, row count, identifiers, and deterministic seeds validate.
Malformed or incompatible output stops the run instead of being overwritten.

Production directories contain three artifact classes:

- `batch_*.csv`: retained observables or traces;
- `diagnostic_*.csv`: mixing, ESS, and exchange diagnostics;
- `execution_*.txt`: runtime and scheduler provenance.

Temporary files are not valid results. A completed campaign should have no
residual `.tmp` files and no non-empty error logs.

## Retained analyses

Single-height and paired Aztec analyses can be regenerated with the commands
in [`aztec/README.md`](../aztec/README.md). The direct Glauber pipeline uses:

```bash
julia --project=aztec aztec/scripts/analyze_glauber_square_grid_production.jl \
  --production-dir /path/to/production \
  --analysis-dir aztec/output/glauber_environment_blocks

julia --project=aztec aztec/scripts/analyze_glauber_square_grid_scaling.jl \
  --gamma-blocks /path/to/core_blocks.csv,/path/to/frontier_blocks.csv \
  --control-blocks /path/to/control_blocks.csv \
  --output-dir aztec/output/glauber_scaling \
  --bootstrap-reps 5000 \
  --bootstrap-seed 20260822
```

The first command converts raw traces into one row per independent frozen
environment. The second resamples those rows within lattice size and compares
`a+b*log(L)` against `a+b*log(L)+c*log(L)^2` over prespecified cutoffs.

The mixing-independent replay of the 960 Gamma environments uses the same
environment seeds:

```bash
julia --project=aztec aztec/scripts/run_glauber_kasteleyn_campaign.jl \
  --config aztec/configs/glauber_square_grid_kasteleyn_replay.csv \
  --output-dir aztec/output/glauber_square_grid_kasteleyn_replay_20260822 \
  --distribution gamma --parameter 0.5 --base-seed 2026082201

julia --project=aztec aztec/scripts/analyze_glauber_kasteleyn_campaign.jl \
  --input-root aztec/output/glauber_square_grid_kasteleyn_replay_20260822 \
  --output-dir aztec/output/glauber_square_grid_kasteleyn_analysis_20260822 \
  --mcmc-blocks aztec/data/glauber_square_grid_20260822/gamma_core_environment_blocks.csv,aztec/data/glauber_square_grid_20260822/gamma_frontier_environment_blocks.csv \
  --bootstrap-reps 5000 --bootstrap-seed 20260823 --cutoffs 2,4,6,8
```

The Kasteleyn runner records solve and real-projection residuals for every
environment and automatically retries at 256-bit precision when the Float64
solve residual exceeds its acceptance threshold.

## Reproduction checklist

Before citing a result, verify:

- the model label and disorder parameterization;
- environment and replica seed uniqueness;
- sample counts by lattice size;
- absence of incomplete batches;
- estimator definitions and bootstrap block;
- BIC sign convention and prediction metric;
- fit-window and weighting sensitivity;
- negative-control behaviour;
- mixing diagnostics for every included size;
- checksums for retained input and derived tables.

## Large data

Raw HPC traces are too large for the Git repository. Keep them in managed
project storage and preserve a compact manifest containing the campaign config,
base seed, code revision, file counts, checksums, scheduler outcome, and the
exact command used to derive retained summaries.
