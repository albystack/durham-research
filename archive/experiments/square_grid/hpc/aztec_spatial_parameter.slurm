#!/bin/bash
#SBATCH --job-name=aztec-spatial
#SBATCH --partition=long
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

: "${CONFIG:?Set CONFIG}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR}"
: "${BASE_SEED:?Set BASE_SEED}"
: "${ALPHA:?Set ALPHA}"
: "${BETA:?Set BETA}"
: "${CAMPAIGN_LABEL:?Set CAMPAIGN_LABEL}"

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || {
  echo "Run bash hpc/setup_hamilton_environment.sh first" >&2
  exit 1
}
mkdir -p logs "${OUTPUT_DIR}"

/usr/bin/time -v srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model gamma \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}" \
  --fractions 1/32,1/16,1/8,1/4 \
  --base-seed "${BASE_SEED}" \
  --alpha "${ALPHA}" \
  --beta "${BETA}" \
  --campaign-label "${CAMPAIGN_LABEL}" \
  --task-id "${SLURM_ARRAY_TASK_ID}"
