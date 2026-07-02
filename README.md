# Quenched LERW winding in random environments

Monte Carlo research on loop-erased random walks (LERW), winding fluctuations,
and random dimer height functions. The numerical question is whether quenched
disorder changes winding variance from ordinary logarithmic growth,

$$
\operatorname{Var}(W_L) \asymp C\log L,
$$

to a super-rough regime,

$$
\operatorname{Var}(W_L) \asymp C(\log L)^2.
$$

This repository contains a restart-safe Julia simulation and analysis pipeline,
plus the archived Python implementation used for the smaller-size runs.

## Main result

The combined study contains **303,000 simulated walks**, **15 weight
distributions**, and **8 lattice sizes** from `L=32` to `L=4096`.

- The simple-random-walk baseline gives a fitted power `p = 1.009` in
  `Var(W_L) = C(log L)^p`.
- Fourteen of fifteen disorder specifications prefer the additive
  `a + b log L` model to `a + b(log L)^2` by BIC.
- The sole additive-fit exception, `Gamma(shape=1)`, still has fitted power
  close to one and a largest-scale local exponent of `0.979`.
- The experiments therefore show no robust evidence for squared-log growth in
  the tested regime.

These are finite-size numerical results, not a proof of asymptotic scaling.
Compact aggregate outputs are stored in [`reports/`](reports/).

## Model

The walk starts at the origin in `[-L,L]^2`. Each lattice site receives four
independent positive outgoing weights in north, east, south, and west order.
Weights are sampled when a site is first exposed and then reused by every walk
in that quenched environment. Steps are selected proportionally to those
weights, the walk stops on first hitting the square boundary, and chronological
loop erasure is maintained online. The observable is left quarter-turns minus
right quarter-turns along the loop-erased path.

## Repository layout

```text
research-julia/   Active Julia package, tests, configs, and Slurm scripts
research-python/  Archived Python implementation and earlier experiment code
reports/          Small aggregate tables from the final combined analysis
```

Generated raw walks, analysis directories, logs, and local environments are
excluded from version control.

## Quick start

Julia 1.10 or newer is required.

```bash
cd research-julia
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

julia --project=. scripts/run_batch.jl \
  --task-id 0 \
  --config configs/hpc_smoke.csv \
  --output-dir results_hpc_smoke

julia --project=. scripts/analyze_results.jl \
  --config configs/hpc_smoke.csv \
  --results-dir results_hpc_smoke \
  --output-dir analysis_hpc_smoke \
  --allow-incomplete
```

See [`research-julia/README.md`](research-julia/README.md) for the complete
workflow and [`research-julia/HPC_README.md`](research-julia/HPC_README.md) for
Slurm execution.

## Reproducibility

- Seeds are derived with SHA-256 and are stable across platforms.
- Result files record code, Julia, task, environment, and walk metadata.
- Batch CSVs are written atomically and completed tasks are not overwritten.
- Annealed uncertainty resamples whole environments rather than treating walks
  from one environment as independent.

The full raw dataset is intentionally not stored in Git because of its size. A
versioned data release can be attached separately when the repository is made
public.
