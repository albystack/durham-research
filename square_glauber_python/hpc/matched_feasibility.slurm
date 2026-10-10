#!/bin/bash
#SBATCH --job-name=matched-feasibility
#SBATCH --partition=shared
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Submit from an isolated repository snapshot after authenticated Hamilton login.
# Confirm partition/account choices against the live Hamilton setup first.
# This runs the capped feasibility design, never a production scaling campaign.
set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the research repository root}"
: "${MATCHED_OUTPUT:?Set an explicit output directory under /nobackup}"
[[ "$MATCHED_OUTPUT" == /nobackup/* ]] || { echo "Use /nobackup for MATCHED_OUTPUT" >&2; exit 2; }
cd "$SLURM_SUBMIT_DIR"
module load "${RESEARCH_PYTHON_MODULE:-python/3.12.6}"
comparison_python="${RESEARCH_PYTHON:-square_glauber_python/.hamilton_python/bin/python}"
test -x "$comparison_python"
export PYTHONPATH="$SLURM_SUBMIT_DIR/square_glauber_python/src"
export NUMBA_CACHE_DIR="$MATCHED_OUTPUT/numba-cache"
export MPLCONFIGDIR="$MATCHED_OUTPUT/matplotlib-cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"
srun "$comparison_python" -m square_glauber.comparison_pilot \
  --config square_glauber_python/configs/matched_feasibility_hamilton.json \
  --output "$MATCHED_OUTPUT"
