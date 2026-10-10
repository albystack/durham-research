#!/bin/bash
#SBATCH --job-name=larger-matched-analysis
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from frozen code directory}"
: "${LARGER_OUTPUT:?}" "${RESEARCH_PYTHON:?}"
module load python/3.12.6
export PYTHONPATH="$PWD/square_glauber_python/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$RESEARCH_PYTHON" square_glauber_python/scripts/larger_campaign.py analyse \
  --plan square_glauber_python/configs/larger_matched_20261007.json \
  --root "$LARGER_OUTPUT/batches" --output "$LARGER_OUTPUT/analysis"
