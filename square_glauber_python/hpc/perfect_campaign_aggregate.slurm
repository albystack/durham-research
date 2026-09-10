#!/bin/bash
#SBATCH --job-name=sq-perfect-audit
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail
: "${MANIFEST:?Set MANIFEST}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR}"
cd "$SLURM_SUBMIT_DIR"
export PYTHONPATH="$SLURM_SUBMIT_DIR/src"
export MPLCONFIGDIR="$OUTPUT_DIR/matplotlib_cache"
mkdir -p "$MPLCONFIGDIR"

# This audit deliberately succeeds on an incomplete campaign and records every
# missing/censored task.  Scientific analysis is a separate manual step only
# after all expected pairs certify.
srun .hamilton_python/bin/python -m square_glauber.cli aggregate-perfect-campaign \
  --manifest "$MANIFEST" --output "$OUTPUT_DIR" --allow-incomplete
