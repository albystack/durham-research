#!/bin/bash
#SBATCH --job-name=sqgrid-a4
#SBATCH --partition=shared
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -eo pipefail

JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}"

CUTS=(32 128 256 512)
MIN_ORDER="${CUTS[$((SLURM_ARRAY_TASK_ID - 1))]}"

GAMMA_RESULTS="output/square_grid_max15h/directed_gamma_k05"
UNIFORM_RESULTS="output/square_grid_max15h/baseline"
ANALYSIS_DIR="output/square_grid_max15h_analysis_min${MIN_ORDER}"

mkdir -p "$ANALYSIS_DIR"

echo "Starting analysis for minimum order ${MIN_ORDER}"
echo "Host: $(hostname)"
echo "Started: $(date -Is)"

/usr/bin/time -v "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/analyze_spatial_campaign.jl \
  --gamma-results "$GAMMA_RESULTS" \
  --uniform-results "$UNIFORM_RESULTS" \
  --output-dir "$ANALYSIS_DIR" \
  --bootstrap-reps 2000 \
  --bootstrap-seed "$((20261000 + MIN_ORDER))" \
  --min-order "$MIN_ORDER" \
  --holdout-orders 2

"$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/plot_spatial_campaign.jl \
  --analysis-dir "$ANALYSIS_DIR" \
  --output-dir "$ANALYSIS_DIR"

echo "Completed analysis for minimum order ${MIN_ORDER}"
echo "Finished: $(date -Is)"
