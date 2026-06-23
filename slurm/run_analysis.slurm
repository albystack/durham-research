#!/usr/bin/env bash
# Hamilton 8 Slurm job for combining completed LERW array results.
#
# Edit resources after checking result size:
#   --time
#   --mem
#   --cpus-per-task
#   --partition
#   --account
#
#SBATCH --job-name=lerw-analysis
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err
##SBATCH --partition=CHANGE_ME
##SBATCH --account=CHANGE_ME

set -euo pipefail

CONFIG="${CONFIG:-configs/hpc_test.csv}"
RESULTS_DIR="${RESULTS_DIR:-results_hpc}"
OUTPUT_DIR="${OUTPUT_DIR:-analysis_hpc}"
PYTHON="${PYTHON:-.venv/bin/python}"
FIT_MIN_L="${FIT_MIN_L:-}"
BOOTSTRAP_REPS="${BOOTSTRAP_REPS:-0}"

ARGS=(
  --config "${CONFIG}"
  --results-dir "${RESULTS_DIR}"
  --output-dir "${OUTPUT_DIR}"
  --bootstrap-reps "${BOOTSTRAP_REPS}"
)

if [[ -n "${FIT_MIN_L}" ]]; then
  ARGS+=(--fit-min-L "${FIT_MIN_L}")
fi

"${PYTHON}" analyze_hpc_results.py "${ARGS[@]}"
