#!/bin/bash
#SBATCH --job-name=matched-shuffle-check
#SBATCH --partition=test
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit verified source snapshot}"
: "${VALIDATION_OUTPUT:?Set fresh validation directory}"
: "${RESEARCH_PYTHON:?Set Hamilton Python interpreter}"
module load python/3.12.6
export PYTHONPATH="$PWD/square_glauber_python/src"
export RESEARCH_JULIA=/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia
export JULIA_DEPOT_PATH="$VALIDATION_OUTPUT/julia-depot"
export NUMBA_CACHE_DIR="$VALIDATION_OUTPUT/numba-cache"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$JULIA_DEPOT_PATH" "$NUMBA_CACHE_DIR"
"$RESEARCH_PYTHON" -m pytest square_glauber_python/tests -q
"$RESEARCH_JULIA" --startup-file=no --project=aztec aztec/test/runtests.jl
"$RESEARCH_JULIA" --startup-file=no --project=aztec aztec/scripts/audit_robust_gamma.jl "$VALIDATION_OUTPUT/aztec_numeric_audit.json"
"$RESEARCH_PYTHON" square_glauber_python/scripts/audit_embedded_shuffle.py --output "$VALIDATION_OUTPUT/numeric_audit.json"
for cell in $(seq 0 35); do
  "$RESEARCH_PYTHON" square_glauber_python/scripts/run_matched_shuffle.py \
    --plan square_glauber_python/configs/matched_shuffle_20261004.json \
    --cell "$cell" --output "$VALIDATION_OUTPUT/cells" --validation-only
done
"$RESEARCH_PYTHON" square_glauber_python/scripts/analyse_matched_shuffle.py \
  --plan square_glauber_python/configs/matched_shuffle_20261004.json \
  --input "$VALIDATION_OUTPUT/cells" --output "$VALIDATION_OUTPUT/verification" --validation-only
"$RESEARCH_PYTHON" square_glauber_python/scripts/verify_matched_shuffle_raw.py \
  --plan square_glauber_python/configs/matched_shuffle_20261004.json \
  --input "$VALIDATION_OUTPUT/cells" --output "$VALIDATION_OUTPUT/raw_verification.json" --validation-only
