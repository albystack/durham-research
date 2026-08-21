#!/bin/bash
#SBATCH --job-name=sqgrid-pilot
#SBATCH --partition=shared
#SBATCH --array=1-50%10
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1

: "${ENVIRONMENT_MODEL:?Set ENVIRONMENT_MODEL to baseline, directed_site_iid, or undirected_conductance}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to the campaign result directory}"
PARAMETER="${PARAMETER:-0.5}"
DISTRIBUTION="${DISTRIBUTION:-gamma}"
BASE_SEED="${BASE_SEED:-20260830}"

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || { echo "Run bash hpc/setup_hamilton_environment.sh first" >&2; exit 1; }
mkdir -p logs "${OUTPUT_DIR}"

srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model "${ENVIRONMENT_MODEL}" \
  --distribution "${DISTRIBUTION}" \
  --parameter "${PARAMETER}" \
  --config aztec/configs/square_grid_pilot.csv \
  --output-dir "${OUTPUT_DIR}" \
  --base-seed "${BASE_SEED}" \
  --task-id "${SLURM_ARRAY_TASK_ID}"
