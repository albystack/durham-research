# Slurm/Hamilton workflow for the Julia implementation

Run every command below from `research-julia/` on Hamilton.

## 1. Copy the project

From a local checkout, replace `<username>` and the destination as needed:

```bash
scp -r research-julia \
  <username>@hamilton8.dur.ac.uk:~/research-julia
```

Then connect:

```bash
ssh <username>@hamilton8.dur.ac.uk
cd ~/research-julia
```

## 2. Load Julia

First inspect Hamilton's installed modules:

```bash
module spider julia
```

Load the newest available Julia version that is at least 1.10. The exact module
name is determined by the previous command, for example:

```bash
module load julia/1.11
julia --version
```

If Hamilton has no suitable Julia module, use the official instructions at
<https://julialang.org/downloads/> or ask Durham ARC support which Julia module
they support. Do not compile Julia from source on a login node.

## 3. Install and precompile packages once

This step downloads packages, so run it interactively on the login node rather
than inside every array task:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

After the manifest is created on Hamilton, keep it with the copied project for
reproducible subsequent runs.

## 4. Smoke test

Run one tiny task on the login node only:

```bash
julia --project=. scripts/run_batch.jl \
  --task-id 0 \
  --config configs/hpc_smoke.csv \
  --output-dir results_hpc_smoke
```

Then submit the four-task smoke array:

```bash
mkdir -p logs
sbatch slurm/run_lerw_array.slurm
```

Monitor and inspect resource use:

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

## 5. Time the large sizes before full submission

Generate the proposed config as described in `README.md`. Submit one task for
each size first. Find their task IDs in the CSV, then run a small array such as:

```bash
CONFIG=configs/hpc_large.csv OUTPUT_DIR=results_hpc_large \
  sbatch --array=0,20,40 slurm/run_lerw_array.slurm
```

Use `sacct` to set realistic `--time` and `--mem` values. The disordered
environment cache is bounded at roughly 320 MiB per concurrently active
environment, plus the loop-erased path map. Multiple Julia threads therefore
need proportionally more memory; baseline runs do not allocate this cache.

Once those checks pass:

```bash
CONFIG=configs/hpc_large.csv OUTPUT_DIR=results_hpc_large \
  sbatch --array=0-179 slurm/run_lerw_array.slurm
```

The batch runner is restart-safe: completed batch files are detected and left
unchanged. Failed or partial tasks can be intentionally rerun with the command
line option `--rerun-failed`.

## 6. Analyse

```bash
CONFIG=configs/hpc_large.csv \
RESULTS_DIR=results_hpc_large \
OUTPUT_DIR=analysis_hpc_large \
FIT_MIN_L=1024 \
BOOTSTRAP_REPS=1000 \
sbatch slurm/run_analysis.slurm
```

## 7. Copy results back

From the local machine:

```bash
scp -r <username>@hamilton8.dur.ac.uk:~/research-julia/results_hpc_large \
  research-julia/
scp -r <username>@hamilton8.dur.ac.uk:~/research-julia/analysis_hpc_large \
  research-julia/
```
