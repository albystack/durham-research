#!/bin/bash
#SBATCH --job-name=sqgrid-analysis
#SBATCH --partition=shared
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
source /etc/profile
module purge
module load julia/1.10.4

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1

: "${DISORDER_RESULTS:?Set DISORDER_RESULTS to the directed/undirected Gamma directory}"
: "${BASELINE_RESULTS:?Set BASELINE_RESULTS to the baseline directory}"
ANALYSIS_DIR="${ANALYSIS_DIR:-output/square_grid_analysis}"

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || { echo "Run bash hpc/setup_hamilton_environment.sh first" >&2; exit 1; }
mkdir -p logs "${ANALYSIS_DIR}"

srun julia --project=.hamilton_env --startup-file=no \
  aztec/scripts/analyze_spatial_campaign.jl \
  --gamma-results "${DISORDER_RESULTS}" \
  --uniform-results "${BASELINE_RESULTS}" \
  --output-dir "${ANALYSIS_DIR}" \
  --bootstrap-reps 10000 \
  --bootstrap-seed 20260831 \
  --min-order 16 \
  --holdout-orders 2

srun julia --project=.hamilton_env --startup-file=no \
  aztec/scripts/plot_spatial_campaign.jl \
  --analysis-dir "${ANALYSIS_DIR}" \
  --output-dir "${ANALYSIS_DIR}"
