# Hamilton and Slurm workflows

This directory contains restart-safe Slurm wrappers for the repository's CPU
campaigns. The scripts target Durham's Hamilton cluster, but their separation
of configs, source, and large output is applicable to other Slurm systems.

## Principles

- Use login nodes only for setup, inspection, and short validations.
- Put source and compact metadata in backed-up project storage.
- Put large batches below an explicit project-data root such as
  `/nobackup/$USER/...`.
- Submit a smoke test, then a pilot, then production.
- Give every campaign a unique output directory and base seed.
- Never overwrite an existing production root with a different config or seed.
- Inspect error logs, exit codes, temporary files, and scientific diagnostics
  before accepting a campaign.

## Julia environment

Hamilton provides a fixed Julia executable that may differ from a developer's
local version. The setup helper creates an isolated `.hamilton_env`, develops
the local package into it, and runs the package tests:

```bash
bash hpc/setup_hamilton_environment.sh
```

The resulting directory is machine-local and ignored by Git.

## Wrapper inventory

| Workflow | Wrappers |
|---|---|
| Aztec spatial increments | `aztec_spatial_smoke.slurm`, `aztec_spatial_parameter.slurm`, `analyze_aztec_parameter.slurm` |
| Temperley square-grid experiments | `square_grid_smoke.slurm`, `square_grid_pilot.slurm`, `square_grid_robustness.slurm`, `square_grid_high_l.slurm`, `analyze_square_grid.slurm` |
| Direct square-grid Glauber dynamics | `glauber_square_grid_pilot.slurm`, `analyze_glauber_square_grid_pilot.slurm`, `glauber_square_grid_production.slurm` |
| Refreshed-transition checks | `temporal_square_grid_pilot.slurm`, `temporal_ust_confirmation.slurm`, `temporal_ust_thread_benchmark.slurm` |
| Historical packed recovery | `../archive/hpc_recovery_202608/` |

Packed recovery scripts encode a specific audited campaign history. They are
kept outside the reusable wrapper directory and should not be treated as
generic launchers.

## First smoke test

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
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

Proceed only when the job exits successfully and the error log is empty.

## Direct weighted-dimer workflow

The direct sampler freezes independent positive mean-one weights on square-grid
dimer edges and uses exact random-face heat-bath dynamics. The `accelerated`
kernel skips self-loops while preserving literal attempted-update time. The
`tempered` kernel runs exact replica exchange in the same frozen environment
and records the physical observable only at `beta=1`.

List a config's tasks before submission:

```bash
/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia \
  --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_glauber_square_grid_campaign.jl \
  --config aztec/configs/glauber_square_grid_production_smoke.csv \
  --output-dir /nobackup/$USER/glauber_smoke \
  --phase production \
  --base-seed 2026082203 \
  --list-tasks
```

Submit a small production-path smoke array with explicit parameters:

```bash
sbatch --array=1-2%2 \
  --export=ALL,CONFIG=aztec/configs/glauber_square_grid_production_smoke.csv,OUTPUT_DIR=/nobackup/$USER/glauber_smoke,BASE_SEED=2026082203,DISTRIBUTION=gamma,PARAMETER=0.5,ALGORITHM=tempered \
  hpc/glauber_square_grid_production.slurm
```

Production configs are split by scientific role:

- `glauber_square_grid_production_core.csv`: calibrated size range;
- `glauber_square_grid_production_frontier.csv`: expensive larger sizes;
- `glauber_square_grid_production_control.csv`: all-one negative control.

Keeping the frontier separate prevents a timeout or failed mixing diagnostic
from invalidating the core campaign.

## Analysis

Convert one completed raw root to independent environment blocks:

```bash
/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia \
  --project=.hamilton_env --startup-file=no \
  aztec/scripts/analyze_glauber_square_grid_production.jl \
  --production-dir /nobackup/$USER/glauber_production \
  --analysis-dir /nobackup/$USER/glauber_production_analysis
```

The size-scaling analysis can then consume separate core/frontier block tables
while preserving the environment as the bootstrap unit. See
[`../docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md).

## Completion audit

For every array, record:

```bash
sacct -j JOB_ID -X --format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS
find /path/to/output -type f -name '*.tmp' -print
find /path/to/output -type f -name '*.err' -size +0c -print
```

Then verify expected counts for batches, diagnostics, execution provenance,
and campaign metadata. A successful scheduler exit is necessary but not
sufficient: start-state agreement, ESS, exchange flow, and controls determine
whether the scientific output is usable.

## Transfer and retention

Copy compact derived tables, configs, metadata, and checksums back to the main
repository or managed archival storage. Large `/nobackup` roots are not a
substitute for a retained manifest. Do not commit credentials, account IDs,
private absolute paths, or raw scheduler logs containing personal metadata.
