#!/bin/bash
#SBATCH --job-name=robust-validation
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit isolated verified source snapshot}"
: "${VALIDATION_OUTPUT:?Set new output location}"
: "${RESEARCH_PYTHON:?Set Hamilton Python executable}"
module load python/3.12.6
comparison_julia="${RESEARCH_JULIA:-/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia}"
export PYTHONPATH="$PWD/square_glauber_python/src"
export JULIA_DEPOT_PATH="$VALIDATION_OUTPUT/julia-depot"
export NUMBA_CACHE_DIR="$VALIDATION_OUTPUT/numba-cache"
export MPLCONFIGDIR="$VALIDATION_OUTPUT/matplotlib-cache"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$JULIA_DEPOT_PATH" "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"
"$RESEARCH_PYTHON" -m pytest square_glauber_python/tests -q
"$comparison_julia" --startup-file=no --project=aztec aztec/test/runtests.jl
"$comparison_julia" --startup-file=no --project=aztec aztec/scripts/audit_robust_gamma.jl "$VALIDATION_OUTPUT/aztec_numeric_audit.json"
"$comparison_julia" --startup-file=no --project=aztec aztec/scripts/run_robust_refinement.jl --output "$VALIDATION_OUTPUT/benchmark/n128" --size 128 --environments 4 --seed 2026100102
