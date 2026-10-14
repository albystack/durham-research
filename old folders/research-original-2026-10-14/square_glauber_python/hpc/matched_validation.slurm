#!/bin/bash
#SBATCH --job-name=matched-validation
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Small mathematical tests on a compute node; no production sampling.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the isolated code snapshot}"
: "${RESEARCH_PYTHON:?Set the existing Hamilton Python executable}"
: "${VALIDATION_OUTPUT:?Set the isolated validation output directory}"
module load "${RESEARCH_PYTHON_MODULE:-python/3.12.6}"
comparison_julia="${RESEARCH_JULIA:-/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia}"
export PYTHONPATH="$PWD/square_glauber_python/src"
export NUMBA_CACHE_DIR="$VALIDATION_OUTPUT/numba-cache"
export MPLCONFIGDIR="$VALIDATION_OUTPUT/matplotlib-cache"
export JULIA_DEPOT_PATH="$VALIDATION_OUTPUT/julia-depot"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR" "$JULIA_DEPOT_PATH"
"$RESEARCH_PYTHON" -m pytest square_glauber_python/tests -q
"$comparison_julia" --startup-file=no --project=aztec aztec/test/runtests.jl
"$comparison_julia" --startup-file=no --project=aztec \
  aztec/scripts/audit_comparison_mapping.jl "$VALIDATION_OUTPUT/aztec_mapping_audit.json"
"$RESEARCH_PYTHON" -c 'import json,sys,numpy,pandas,numba; from pathlib import Path; Path(sys.argv[1]).write_text(json.dumps(dict(python=sys.version,numpy=numpy.__version__,pandas=pandas.__version__,numba=numba.__version__,validation_passed=True),indent=2)+"\n")' \
  "$VALIDATION_OUTPUT/runtime.json"
