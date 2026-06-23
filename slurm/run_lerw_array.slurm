#!/usr/bin/env bash
# Hamilton 8 Slurm array for one LERW batch per task.
#
# Before submitting, create the log directory from the repository root:
#   mkdir -p logs
#
# Edit these resource lines after a small Hamilton smoke test:
#   --time              wall clock limit
#   --mem               memory per task
#   --cpus-per-task     keep at 1 unless profiling supports internal parallelism
#   --array             task-id range from generate_hpc_config.py
#   --partition         Hamilton partition, if required
#   --account           Hamilton account, if required
#
#SBATCH --job-name=lerw-hpc
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-3
#SBATCH --output=logs/lerw_%A_%a.out
#SBATCH --error=logs/lerw_%A_%a.err
##SBATCH --partition=CHANGE_ME
##SBATCH --account=CHANGE_ME

set -euo pipefail

CONFIG="${CONFIG:-configs/hpc_test.csv}"
OUTPUT_DIR="${OUTPUT_DIR:-results_hpc}"
PYTHON="${PYTHON:-.venv/bin/python}"

"${PYTHON}" run_hpc_batch.py \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}"
