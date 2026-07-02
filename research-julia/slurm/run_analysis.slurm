#!/usr/bin/env bash
#SBATCH --job-name=lerw-analysis-julia
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err
##SBATCH --partition=CHANGE_ME
##SBATCH --account=CHANGE_ME

set -euo pipefail

CONFIG="${CONFIG:-configs/hpc_smoke.csv}"
RESULTS_DIR="${RESULTS_DIR:-results_hpc}"
OUTPUT_DIR="${OUTPUT_DIR:-analysis_hpc}"
JULIA="${JULIA:-julia}"
FIT_MIN_L="${FIT_MIN_L:-}"
BOOTSTRAP_REPS="${BOOTSTRAP_REPS:-500}"

ARGS=(
  --config "${CONFIG}"
  --results-dir "${RESULTS_DIR}"
  --output-dir "${OUTPUT_DIR}"
  --bootstrap-reps "${BOOTSTRAP_REPS}"
)

if [[ -n "${FIT_MIN_L}" ]]; then
  ARGS+=(--fit-min-L "${FIT_MIN_L}")
fi

"${JULIA}" --project=. scripts/analyze_results.jl "${ARGS[@]}"
