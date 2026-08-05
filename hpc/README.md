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
