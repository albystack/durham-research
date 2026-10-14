#!/bin/bash
#SBATCH --job-name=sq-perfect-finish
#SBATCH --partition=shared
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
: "${MANIFEST:?Set the frozen campaign manifest}"
: "${OUTPUT_DIR:?Set the existing campaign root}"
: "${ANALYSIS_DIR:?Set a new analysis directory}"
cd "$SLURM_SUBMIT_DIR"
export PYTHONPATH="$SLURM_SUBMIT_DIR/src"
export MPLCONFIGDIR="$OUTPUT_DIR/matplotlib_cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test ! -e "$ANALYSIS_DIR"

# Strict aggregation writes the status/retry tables, then fails if any pair
# remains missing or censored. With set -e, analysis cannot run in that case.
srun .hamilton_python/bin/python -m square_glauber.cli aggregate-perfect-campaign \
  --manifest "$MANIFEST" --output "$OUTPUT_DIR"
srun .hamilton_python/bin/python -m square_glauber.cli analyze \
  --input "$OUTPUT_DIR" --output "$ANALYSIS_DIR" \
  --bootstrap-replicates 10000 --bootstrap-seed 2026090902
