#!/bin/bash
#SBATCH --job-name=sq-restart-finish
#SBATCH --partition=shared
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONPATH="$SLURM_SUBMIT_DIR/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$OUTPUT_DIR/matplotlib_cache"
srun .hamilton_python/bin/python hpc/restart_stage_20260913.py consolidate
# Existing strict aggregation then frozen analysis; no analysis on censoring.
bash hpc/perfect_campaign_finish.slurm
