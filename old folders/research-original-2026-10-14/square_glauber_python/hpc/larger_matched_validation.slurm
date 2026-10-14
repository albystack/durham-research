#!/bin/bash
#SBATCH --job-name=larger-matched-check
#SBATCH --partition=shared
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from frozen code directory}"
: "${VALIDATION_OUTPUT:?}" "${RESEARCH_PYTHON:?}"
module load python/3.12.6
export PYTHONPATH="$PWD/square_glauber_python/src"
export RESEARCH_JULIA=/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia
export JULIA_DEPOT_PATH="$VALIDATION_OUTPUT/julia-depot"
export NUMBA_CACHE_DIR="$VALIDATION_OUTPUT/numba-cache"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$JULIA_DEPOT_PATH" "$NUMBA_CACHE_DIR"
"$RESEARCH_PYTHON" -m pytest square_glauber_python/tests -q
"$RESEARCH_JULIA" --startup-file=no --project=aztec aztec/test/runtests.jl
for batch in $(seq 0 41); do
  "$RESEARCH_PYTHON" square_glauber_python/scripts/larger_campaign.py run \
    --plan square_glauber_python/configs/larger_matched_20261007_validation.json \
    --root "$VALIDATION_OUTPUT/batches" --batch "$batch"
done
# Exercise reuse of a completed batch without regenerating it.
"$RESEARCH_PYTHON" square_glauber_python/scripts/larger_campaign.py run \
  --plan square_glauber_python/configs/larger_matched_20261007_validation.json \
  --root "$VALIDATION_OUTPUT/batches" --batch 0
"$RESEARCH_PYTHON" square_glauber_python/scripts/larger_campaign.py analyse \
  --plan square_glauber_python/configs/larger_matched_20261007_validation.json \
  --root "$VALIDATION_OUTPUT/batches" --output "$VALIDATION_OUTPUT/analysis"
