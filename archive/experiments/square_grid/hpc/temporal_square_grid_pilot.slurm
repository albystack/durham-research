#!/bin/bash
#SBATCH --job-name=sqgrid-temporal-pilot
#SBATCH --partition=shared
#SBATCH --array=1-6%3
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

# Prepared only: do not submit until the local pilot and its scientific review
# have been accepted.  Each array element owns one deterministic output batch.
set -euo pipefail

JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

: "${OUTPUT_DIR:?Set OUTPUT_DIR under /nobackup for the campaign result directory}"
: "${BASE_SEED:?Set BASE_SEED}"
CONFIG="${CONFIG:-aztec/configs/temporal_square_grid_local_pilot.csv}"
DISTRIBUTION="${DISTRIBUTION:-gamma}"
PARAMETER="${PARAMETER:-0.5}"

[[ "${OUTPUT_DIR}" == /nobackup/* ]] || {
  echo "OUTPUT_DIR must be under /nobackup on Hamilton8" >&2
  exit 2
}

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || {
  echo "Run bash hpc/setup_hamilton_environment.sh first" >&2
  exit 1
}
mkdir -p logs "${OUTPUT_DIR}"

srun --cpu-bind=cores "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model temporal_iid \
  --distribution "${DISTRIBUTION}" \
  --parameter "${PARAMETER}" \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}" \
  --base-seed "${BASE_SEED}" \
  --task-id "${SLURM_ARRAY_TASK_ID}"
