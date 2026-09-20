#!/bin/bash
#SBATCH --job-name=sqgrid-smoke
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}"
SMOKE_DIR="output/hpc_square_grid_smoke_${SLURM_JOB_ID}"
mkdir -p logs "${SMOKE_DIR}"

srun bash hpc/setup_hamilton_environment.sh

srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model baseline \
  --config aztec/configs/square_grid_smoke.csv \
  --output-dir "${SMOKE_DIR}/baseline" \
  --base-seed 20260820

srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model directed_site_iid \
  --distribution gamma --parameter 0.5 \
  --config aztec/configs/square_grid_smoke.csv \
  --output-dir "${SMOKE_DIR}/directed_gamma_k05" \
  --base-seed 20260821
