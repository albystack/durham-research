#!/bin/bash
#SBATCH --job-name=sq-perfect-smoke
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-type=FAIL

set -euo pipefail

: "${SMOKE_OUTPUT:?Set a new SMOKE_OUTPUT under /nobackup}"
[[ "$SMOKE_OUTPUT" == /nobackup/* ]] || {
  echo "SMOKE_OUTPUT must be under /nobackup on Hamilton" >&2
  exit 2
}

cd "$SLURM_SUBMIT_DIR"
test -x .hamilton_python/bin/python
test ! -e "$SMOKE_OUTPUT" || {
  echo "refusing to overwrite smoke output $SMOKE_OUTPUT" >&2
  exit 2
}
mkdir -p logs "$SMOKE_OUTPUT"
export PYTHONPATH="$SLURM_SUBMIT_DIR/src"
export NUMBA_CACHE_DIR="$SMOKE_OUTPUT/numba_cache"
export MPLCONFIGDIR="$SMOKE_OUTPUT/matplotlib_cache"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR"

manifest="$SMOKE_OUTPUT/smoke_manifest.csv"
srun .hamilton_python/bin/python -m square_glauber.cli \
  create-perfect-campaign-manifest \
  --schedule configs/perfect_covariance_smoke.csv \
  --campaign-id "hamilton-perfect-smoke-${SLURM_JOB_ID}" \
  --seed 2026090907 --backend numba --output "$manifest"

for task_id in 0 1; do
  srun .hamilton_python/bin/python -m square_glauber.cli perfect-campaign-task \
    --manifest "$manifest" --task-id "$task_id" --output "$SMOKE_OUTPUT"
done

srun .hamilton_python/bin/python -m square_glauber.cli aggregate-perfect-campaign \
  --manifest "$manifest" --output "$SMOKE_OUTPUT"

.hamilton_python/bin/python -c \
  'import json,sys; p=json.load(open(sys.argv[1])); assert p["complete_environments"] == 2; assert p["all_pairs_cftp_certified"]' \
  "$SMOKE_OUTPUT/campaign_summary.json"
