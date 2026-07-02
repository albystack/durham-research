# Slurm workflow

This guide runs the Julia simulation pipeline on a generic Slurm cluster. Run
all commands from `research-julia/` and adapt module, partition, account, time,
and memory settings to the target system.

## 1. Prepare the environment

Load Julia 1.10 or newer using the cluster's module system, for example:

```bash
module spider julia
module load julia/1.11
julia --version
```

Install and precompile dependencies once from a login or development node:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Do not compile Julia from source or run large simulations on a login node.

## 2. Run a smoke test

Run one small task interactively:

```bash
julia --project=. scripts/run_batch.jl \
  --task-id 0 \
  --config configs/hpc_smoke.csv \
  --output-dir results_hpc_smoke
```

Then submit the four-task smoke array:

```bash
mkdir -p logs
sbatch --array=0-3 slurm/run_lerw_array.slurm
```

Monitor jobs and inspect resource use:

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

## 3. Benchmark large sizes

Generate an experiment grid as described in `README.md`. Before launching the
full array, submit one representative task at each lattice size:

```bash
CONFIG=configs/hpc_large.csv OUTPUT_DIR=results_hpc_large \
  sbatch --array=0,20,40 slurm/run_lerw_array.slurm
```

Use the measured elapsed time and peak RSS to set realistic Slurm limits. Each
active disordered environment can use roughly 320 MiB for its bounded site
cache, plus the loop-erased path map. Memory therefore grows with Julia thread
count. Baseline runs do not allocate the site cache.

## 4. Submit the experiment

After resource validation, submit the complete task range:

```bash
CONFIG=configs/hpc_large.csv OUTPUT_DIR=results_hpc_large \
  sbatch --array=0-179 slurm/run_lerw_array.slurm
```

The runner is restart-safe: complete files are skipped. To intentionally
replace failed or partial output, pass `--rerun-failed` to the batch command.

## 5. Analyse results

```bash
CONFIG=configs/hpc_large.csv \
RESULTS_DIR=results_hpc_large \
OUTPUT_DIR=analysis_hpc_large \
FIT_MIN_L=1024 \
BOOTSTRAP_REPS=1000 \
sbatch slurm/run_analysis.slurm
```

Analysis rejects missing, incomplete, failed, or duplicate tasks by default.
Use `--allow-incomplete` only for explicitly provisional work.

## Site-specific configuration

The templates intentionally omit cluster-specific partition and account
directives. Add them at submission time or copy the templates into a local,
untracked configuration:

```bash
sbatch --partition=PARTITION --account=ACCOUNT \
  --array=0-179 slurm/run_lerw_array.slurm
```

Keep credentials, allocation identifiers, and private filesystem paths out of
version control.
