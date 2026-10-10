#!/bin/bash
#SBATCH --job-name=september-remeasure
#SBATCH --partition=shared
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Read-only use of historical weights/states; new observables go elsewhere.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the isolated code snapshot}"
: "${RESEARCH_PYTHON:?Set the existing Hamilton Python executable}"
: "${HISTORICAL_CAMPAIGN:?Set the September campaign directory}"
: "${REMEASUREMENT_OUTPUT:?Set a new output directory}"
[[ "$REMEASUREMENT_OUTPUT" != "$HISTORICAL_CAMPAIGN"* ]] || {
  echo "Keep derived output outside the historical campaign" >&2; exit 2;
}
module load "${RESEARCH_PYTHON_MODULE:-python/3.12.6}"
export PYTHONPATH="$PWD/square_glauber_python/src"
export MPLCONFIGDIR="$REMEASUREMENT_OUTPUT/matplotlib-cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR"
srun "$RESEARCH_PYTHON" -m square_glauber.remeasure \
  --snapshot "$HISTORICAL_CAMPAIGN/center_measurements.csv" \
  --states "$HISTORICAL_CAMPAIGN/cases" \
  --output "$REMEASUREMENT_OUTPUT"
