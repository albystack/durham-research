# Julia LERW simulation and analysis pipeline

This directory is the active implementation of the quenched random-environment
experiment described in the [repository README](../README.md).

## What the code does

- lazily samples four positive outgoing weights at each exposed lattice site
  from a reproducible point-addressed random stream;
- reuses each site's weights across walks in the same quenched environment;
- bounds the site cache and regenerates evicted weights exactly, avoiding
  unbounded memory growth at large `L`;
- maintains chronological loop erasure online;
- records winding, path lengths, exit location, runtime, and reproducibility
  metadata;
- writes restart-safe, atomic batch CSVs;
- validates task completeness before analysis;
- computes annealed and quenched variance summaries;
- fits `Var(W_L) = C(log L)^p` and compares additive `log L` and `(log L)^2`
  finite-size models;
- uses environment-clustered bootstrap intervals;
- runs independent quenched environments concurrently when Julia has multiple
  threads.

The implementation supports baseline, exponential, gamma, lognormal, Pareto,
uniform, beta, Weibull, inverse-gamma, Bernoulli, and triangular site-weight
distributions. Parameterized specifications use the `name:value` form, such as
`gamma:0.5` or `lognormal:1.0`.

## Requirements

- Julia 1.10 or newer
- packages pinned by `Project.toml` and `Manifest.toml`

Install dependencies and run the test suite:

```bash
cd research-julia
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Local smoke test

Run all four tasks in the bundled smoke configuration:

```bash
for task in 0 1 2 3; do
  julia --project=. scripts/run_batch.jl \
    --task-id "$task" \
    --config configs/hpc_smoke.csv \
    --output-dir results_hpc_smoke
done
```

Analyse the completed batches:

```bash
julia --project=. scripts/analyze_results.jl \
  --config configs/hpc_smoke.csv \
  --results-dir results_hpc_smoke \
  --output-dir analysis_hpc_smoke \
  --bootstrap-reps 100
```

The first Julia invocation includes compilation time. Benchmark subsequent
tasks rather than startup.

For local multicore execution, start Julia with all available threads:

```bash
julia --threads=auto --project=. scripts/run_batch.jl ...
```

Parallelism is across environments; walks within one environment remain
sequential because they share its lazily generated site weights.

## Generate an experiment grid

For example, the following command writes 180 tasks: three distributions,
three sizes, and twenty batches per distribution-size pair.

```bash
julia --project=. scripts/generate_config.jl \
  --output configs/hpc_large.csv \
  --sizes 1024,2048,4096 \
  --distributions baseline,gamma:0.5,lognormal:1.0 \
  --batches 20 \
  --num-environments 10 \
  --walks-per-environment 5
```

Time one representative task at every size before launching a full grid.
Memory use grows with the number of distinct sites exposed in an environment.

## Batch semantics

Each configuration row defines one task and includes its distribution, lattice
size, batch ID, number of environments, walks per environment, and base seed.
Output paths have the form

```text
OUTPUT/site_iid/DISTRIBUTION/L_NNNN/batch_NNNN.csv
```

Completed files are skipped. An existing partial or failed file is never
silently overwritten:

```bash
# Rerun an incomplete task intentionally
julia --project=. scripts/run_batch.jl ... --rerun-failed

# Overwrite any existing task output intentionally
julia --project=. scripts/run_batch.jl ... --force
```

## Analysis outputs

`scripts/analyze_results.jl` writes:

```text
validation.csv
combined_raw.csv
summary.csv
loglog_fits.csv
scaling_model_comparison.csv
pointwise_ratios.csv
local_effective_exponents.csv
```

`loglog_fits.csv` fits

```text
log Var(W_L) = log C + p log(log L).
```

`scaling_model_comparison.csv` separately compares `a + b log L` and
`a + b(log L)^2`. Annealed bootstrap intervals resample entire environments,
because walks sharing one environment are correlated.

Analysis fails on missing, partial, or failed tasks by default. Use
`--allow-incomplete` only for an explicitly provisional analysis.

To combine compatible result trees with different batch sizes (for example,
the Python and Julia runs used in the final report), supply the expected row
count for each tree:

```bash
julia --project=. scripts/analyze_combined_results.jl \
  --results-dirs path/to/python/site_iid,path/to/julia/site_iid \
  --expected-rows 100,50 \
  --output-dir analysis_combined \
  --bootstrap-reps 1000 \
  --bootstrap-seed 20260628
```

## Slurm

The scripts in `slurm/` submit simulation arrays and analysis jobs. See
[`HPC_README.md`](HPC_README.md) for a cluster workflow and resource-checking
guidance.
