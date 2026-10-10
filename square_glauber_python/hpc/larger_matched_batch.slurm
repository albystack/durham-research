#!/bin/bash
#SBATCH --job-name=larger-matched
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
# Array range and fixed wall limit are supplied by the launch script.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from frozen code directory}"
: "${LARGER_OUTPUT:?}" "${RESEARCH_PYTHON:?}"
module load python/3.12.6
export PYTHONPATH="$PWD/square_glauber_python/src"
export RESEARCH_JULIA=/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia
export JULIA_DEPOT_PATH="$LARGER_OUTPUT/julia-depot-$SLURM_ARRAY_TASK_ID"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$JULIA_DEPOT_PATH"
srun "$RESEARCH_PYTHON" square_glauber_python/scripts/larger_campaign.py run \
  --plan square_glauber_python/configs/larger_matched_20261007.json \
  --root "$LARGER_OUTPUT/batches" --batch "$SLURM_ARRAY_TASK_ID"
