#!/bin/bash
#SBATCH --job-name=frontier-validation
#SBATCH --partition=test
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit the verified source snapshot}"
: "${VALIDATION_OUTPUT:?Set a fresh validation directory}"
: "${RESEARCH_PYTHON:?Set Hamilton Python interpreter}"
module load python/3.12.6
export PYTHONPATH="$PWD/square_glauber_python/src"
export NUMBA_CACHE_DIR="$VALIDATION_OUTPUT/numba-cache"
export MPLCONFIGDIR="$VALIDATION_OUTPUT/matplotlib-cache"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"
"$RESEARCH_PYTHON" -m pytest square_glauber_python/tests -q
# The runner's small complete/censored/original-input paths are covered by
# test_frontier_scaling_runner.py. The larger cohort starts only after this job.
