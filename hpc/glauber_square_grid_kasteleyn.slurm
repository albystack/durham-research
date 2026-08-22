#!/bin/bash
#SBATCH --job-name=glauber-kasteleyn
#SBATCH --partition=shared
#SBATCH --array=1-18%4
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL

# Finite-volume determinantal validation, not an MCMC production campaign.
set -euo pipefail

JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

: "${OUTPUT_DIR:?Set OUTPUT_DIR under /nobackup for pilot output}"
: "${BASE_SEED:?Set BASE_SEED; use 2026082201 to replay the Gamma production environments}"
CONFIG="${CONFIG:-aztec/configs/glauber_square_grid_kasteleyn_pilot.csv}"
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
  aztec/scripts/run_glauber_kasteleyn_campaign.jl \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}" \
  --distribution "${DISTRIBUTION}" \
  --parameter "${PARAMETER}" \
  --base-seed "${BASE_SEED}" \
  --task-id "${SLURM_ARRAY_TASK_ID}"
