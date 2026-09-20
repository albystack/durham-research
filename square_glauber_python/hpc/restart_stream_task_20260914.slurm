#!/bin/bash
#SBATCH --job-name=sq-restart-stream
#SBATCH --partition=shared
#SBATCH --time=1-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONPATH="$SLURM_SUBMIT_DIR/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_CACHE_DIR="$OUTPUT_DIR/numba_cache"
export MPLCONFIGDIR="$OUTPUT_DIR/matplotlib_cache"
srun --cpu-bind=cores .hamilton_python/bin/python hpc/restart_stage_20260914.py worker
