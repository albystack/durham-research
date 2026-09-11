#!/bin/bash
#SBATCH --job-name=sq-perfect-pairs
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --mail-type=FAIL

set -euo pipefail

: "${MANIFEST:?Set MANIFEST to the frozen production manifest}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR under /nobackup}"
: "${SLURM_ARRAY_TASK_ID:?Submit this wrapper as an array}"

[[ "$OUTPUT_DIR" == /nobackup/* ]] || {
  echo "OUTPUT_DIR must be under /nobackup on Hamilton" >&2
  exit 2
}

cd "$SLURM_SUBMIT_DIR"
test -x .hamilton_python/bin/python || {
  echo "Run hpc/setup_hamilton_python.sh before submission" >&2
  exit 1
}
test -f "$MANIFEST"
mkdir -p logs "$OUTPUT_DIR" "$OUTPUT_DIR/numba_cache"

export PYTHONPATH="$SLURM_SUBMIT_DIR/src"
export NUMBA_CACHE_DIR="$OUTPUT_DIR/numba_cache"
export MPLCONFIGDIR="$OUTPUT_DIR/matplotlib_cache"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR"

arguments=(
  -m square_glauber.cli perfect-campaign-task
  --manifest "$MANIFEST"
  --task-id "$SLURM_ARRAY_TASK_ID"
  --output "$OUTPUT_DIR"
)
if [[ -n "${MAX_HORIZON_SWEEPS:-}" ]]; then
  arguments+=(--max-horizon-sweeps "$MAX_HORIZON_SWEEPS")
fi

srun --cpu-bind=cores .hamilton_python/bin/python "${arguments[@]}"
