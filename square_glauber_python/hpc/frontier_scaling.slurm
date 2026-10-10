#!/bin/bash
#SBATCH --job-name=square-frontier-scale
#SBATCH --partition=shared
#SBATCH --array=0-17%2
#SBATCH --time=00:04:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit the verified source snapshot}"
: "${FRONTIER_OUTPUT:?Set a new /nobackup output root}"
: "${ORIGINAL_CASES:?Set retained original pilot cases root}"
: "${RESEARCH_PYTHON:?Set validated Hamilton Python interpreter}"
[[ "$FRONTIER_OUTPUT" == /nobackup/* ]] || exit 2
module load python/3.12.6
export PYTHONPATH="$PWD/square_glauber_python/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
srun "$RESEARCH_PYTHON" square_glauber_python/scripts/run_frontier_scaling.py \
  --plan square_glauber_python/configs/frontier_scaling_20261003.json \
  --case "$SLURM_ARRAY_TASK_ID" --output "$FRONTIER_OUTPUT" \
  --original-cases "$ORIGINAL_CASES"
