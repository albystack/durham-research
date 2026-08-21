#!/bin/bash
#SBATCH --job-name=glauber-dimer-pilot
#SBATCH --partition=shared
#SBATCH --array=1-8%4
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL

# This is an equilibration/validation pilot, not a production size-scaling run.
set -euo pipefail

JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

: "${OUTPUT_DIR:?Set OUTPUT_DIR under /nobackup for the pilot output}"
: "${BASE_SEED:?Set BASE_SEED}"
CONFIG="${CONFIG:-aztec/configs/glauber_square_grid_pilot.csv}"
DISTRIBUTION="${DISTRIBUTION:-gamma}"
PARAMETER="${PARAMETER:-0.5}"
ALGORITHM="${ALGORITHM:-accelerated}"
TEMPERING_BETAS="${TEMPERING_BETAS:-0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}"
SWAP_INTERVAL_ATTEMPTS="${SWAP_INTERVAL_ATTEMPTS:-10000}"

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
  aztec/scripts/run_glauber_square_grid_campaign.jl \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}" \
  --distribution "${DISTRIBUTION}" \
  --parameter "${PARAMETER}" \
  --base-seed "${BASE_SEED}" \
  --algorithm "${ALGORITHM}" \
  --tempering-betas "${TEMPERING_BETAS}" \
  --swap-interval-attempts "${SWAP_INTERVAL_ATTEMPTS}" \
  --task-id "${SLURM_ARRAY_TASK_ID}"
