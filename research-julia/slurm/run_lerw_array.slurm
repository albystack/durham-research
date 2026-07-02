#!/usr/bin/env bash
#SBATCH --job-name=lerw-julia
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-3
#SBATCH --output=logs/lerw_%A_%a.out
#SBATCH --error=logs/lerw_%A_%a.err
## Add site-specific --partition and --account directives when required.

set -euo pipefail

CONFIG="${CONFIG:-configs/hpc_smoke.csv}"
OUTPUT_DIR="${OUTPUT_DIR:-results_hpc}"
JULIA="${JULIA:-julia}"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

"${JULIA}" --project=. scripts/run_batch.jl \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}"
