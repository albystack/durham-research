#!/bin/bash
#SBATCH --job-name=matched-shuffle
#SBATCH --partition=shared
#SBATCH --array=0-35%2
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from verified source snapshot}"
: "${MATCHED_OUTPUT:?Set isolated output directory}"
: "${RESEARCH_PYTHON:?Set Hamilton Python interpreter}"
[[ "$MATCHED_OUTPUT" == /nobackup/* ]] || exit 2
module load python/3.12.6
export PYTHONPATH="$PWD/square_glauber_python/src"
export RESEARCH_JULIA=/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia
export JULIA_DEPOT_PATH="$MATCHED_OUTPUT/julia-depot-$SLURM_ARRAY_TASK_ID"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$JULIA_DEPOT_PATH"
srun "$RESEARCH_PYTHON" square_glauber_python/scripts/run_matched_shuffle.py \
  --plan square_glauber_python/configs/matched_shuffle_20261004.json \
  --cell "$SLURM_ARRAY_TASK_ID" --output "$MATCHED_OUTPUT/cells"
