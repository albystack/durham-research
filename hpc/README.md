# Durham Hamilton 8 workflow

Verified account details for this project:

```text
host: hamilton8.dur.ac.uk
username: fvkl37
home: /home/fvkl37
Julia module: julia/1.10.4
partitions: test, shared, long, multi, bigmem, cuda
home quota: 10 GB
nobackup quota: 600 GB
```

Keep source and small logs under `/home/fvkl37`. Put large batch output under
`/nobackup/fvkl37` and copy retained results back to the Mac because
`/nobackup` is not the backed-up source repository.

## Direct weighted-dimer Glauber pilot

This is the independent square-grid dimer sampler requested in Sunil's latest
email. It uses the supplied fixed height boundary, frozen i.i.d. mean-one edge
weights, and exact random-face heat-bath dynamics. Its `accelerated` option
skips self-loops while preserving the distribution at fixed attempted-update
times; it is not the biased active-site jump chain.

Prepare the Julia environment, then begin with the (L=2) calibration.  This
is the only initial size where every frozen weighted environment can be checked
against exact enumeration.  Do **not** submit the broader size pilot until this
calibration has passed review.

```bash
bash hpc/setup_hamilton_environment.sh
mkdir -p logs /nobackup/fvkl37/glauber_dimer

sbatch --array=1-2%2 --export=ALL,CONFIG=aztec/configs/glauber_square_grid_L2_calibration.csv,OUTPUT_DIR=/nobackup/fvkl37/glauber_dimer/L2_constant_calibration,BASE_SEED=2026082101,DISTRIBUTION=constant \
  hpc/glauber_square_grid_pilot.slurm

sbatch --array=1-2%2 --export=ALL,CONFIG=aztec/configs/glauber_square_grid_L2_calibration.csv,OUTPUT_DIR=/nobackup/fvkl37/glauber_dimer/L2_gamma_k05_calibration,BASE_SEED=2026082102,DISTRIBUTION=gamma,PARAMETER=0.5 \
  hpc/glauber_square_grid_pilot.slurm
```

After all eight array elements in one arm finish, generate the diagnostic
summary (including the exact (L=2) check):

```bash
sbatch --export=ALL,PILOT_DIR=/nobackup/fvkl37/glauber_dimer/L2_gamma_k05_calibration,ANALYSIS_DIR=/nobackup/fvkl37/glauber_dimer/L2_gamma_k05_calibration_analysis \
  hpc/analyze_glauber_square_grid_pilot.slurm
```

The analysis deliberately creates `REVIEW_REQUIRED.txt`, rather than a
production approval. Inspect start gaps, ESS, exact (L=2) agreement, raw
traces, and Slurm accounting before choosing the production burn-in and
sampling schedule.

## Prepare Julia 1.10.4

The retained `aztec/Manifest.toml` was generated with Julia 1.12.6. Preserve it and create a separate Hamilton environment:

```bash
bash hpc/setup_hamilton_environment.sh
```

This runs the full package test suite using Hamilton Julia 1.10.4.

## First Hamilton smoke test

From the repository root on Hamilton:

```bash
mkdir -p logs
sbatch hpc/square_grid_smoke.slurm
squeue -u "$USER"
```

After completion:

```bash
cat logs/sqgrid-smoke-*.out
cat logs/sqgrid-smoke-*.err
```

## Pilot task table

```bash
module load julia/1.10.4
julia --project=aztec aztec/scripts/run_square_grid_campaign.jl \
  --config aztec/configs/square_grid_pilot.csv --list-tasks
```

The frozen pilot expands to 50 tasks, matching the Slurm array `1-50`.

## Submit baseline pilot

```bash
mkdir -p /nobackup/fvkl37/square_grid/baseline
sbatch \
  --export=ALL,ENVIRONMENT_MODEL=baseline,PARAMETER=0.5,BASE_SEED=20260840,OUTPUT_DIR=/nobackup/fvkl37/square_grid/baseline \
  hpc/square_grid_pilot.slurm
```

## Submit directed Gamma pilot

```bash
mkdir -p /nobackup/fvkl37/square_grid/directed_gamma_k05
sbatch \
  --export=ALL,ENVIRONMENT_MODEL=directed_site_iid,PARAMETER=0.5,BASE_SEED=20260841,OUTPUT_DIR=/nobackup/fvkl37/square_grid/directed_gamma_k05 \
  hpc/square_grid_pilot.slurm
```

Do not submit the undirected pilot until the directed pilot, tests, and
baseline analysis have been checked.

## Monitor and inspect

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
scancel JOB_ID
```

## Analyse after both arrays finish

```bash
sbatch \
  --export=ALL,DISORDER_RESULTS=/nobackup/fvkl37/square_grid/directed_gamma_k05,BASELINE_RESULTS=/nobackup/fvkl37/square_grid/baseline,ANALYSIS_DIR=/nobackup/fvkl37/square_grid/analysis_directed_k05 \
  hpc/analyze_square_grid.slurm
```

## Robustness campaigns

The hardened runner gives every law/model pair a unique label, includes the
law in deterministic seed derivation, materializes one transition table shared
by both replicas, and writes a hashed TOML campaign manifest plus per-task
Slurm provenance. List the 899 production tasks with:

```bash
julia --project=.hamilton_env aztec/scripts/run_square_grid_campaign.jl \
  --config aztec/configs/square_grid_robustness.csv --list-tasks
```

Submit with an explicit array and concurrency cap. For example, directed
Gamma `k=1` is:

```bash
sbatch --array=1-899%20 \
  --export=ALL,ENVIRONMENT_MODEL=directed_site_iid,DISTRIBUTION=gamma,PARAMETER=1.0,BASE_SEED=20260881,CONFIG=aztec/configs/square_grid_robustness.csv,OUTPUT_DIR=/nobackup/fvkl37/square_grid_robustness/directed_gamma_k1 \
  hpc/square_grid_robustness.slurm
```

Use distinct output directories and base seeds for directed matched-variance
lognormal (`sigma=1.048147073968205`), directed `Uniform(0,2)`, and
undirected-conductance campaigns. Analyse each separately against the frozen
legacy baseline; pass its actual campaign labels through `DISORDER_LABEL` and
`BASELINE_LABEL`. The analysis now writes a block-GLS sensitivity table using
the bootstrap-estimated `4 x 4` cross-fraction covariance at each order.

For a node-QOS-limited recovery of already submitted one-sample arrays, see
[`PACKED_CAMPAIGNS.md`](PACKED_CAMPAIGNS.md). The packed workflow preserves the
original task IDs and dependency chain; it must not be used without first
holding the corresponding pending original tasks.
