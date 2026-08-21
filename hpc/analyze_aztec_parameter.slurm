#!/bin/bash
#SBATCH --job-name=aztec-analysis
#SBATCH --partition=shared
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1

: "${GAMMA_RESULTS:?Set GAMMA_RESULTS}"
: "${UNIFORM_RESULTS:?Set UNIFORM_RESULTS}"
: "${GAMMA_LABEL:?Set GAMMA_LABEL}"
: "${ANALYSIS_DIR:?Set ANALYSIS_DIR}"
BOOTSTRAP_REPS="${BOOTSTRAP_REPS:-10000}"
BOOTSTRAP_SEED="${BOOTSTRAP_SEED:-20260901}"
MIN_ORDER="${MIN_ORDER:-128}"
HOLDOUT_ORDERS="${HOLDOUT_ORDERS:-2}"

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || { echo "Run bash hpc/setup_hamilton_environment.sh first" >&2; exit 1; }
mkdir -p logs "${ANALYSIS_DIR}"

srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/analyze_spatial_campaign.jl \
  --gamma-results "${GAMMA_RESULTS}" \
  --uniform-results "${UNIFORM_RESULTS}" \
  --gamma-model "${GAMMA_LABEL}" \
  --uniform-model uniform \
  --output-dir "${ANALYSIS_DIR}" \
  --bootstrap-reps "${BOOTSTRAP_REPS}" \
  --bootstrap-seed "${BOOTSTRAP_SEED}" \
  --min-order "${MIN_ORDER}" \
  --holdout-orders "${HOLDOUT_ORDERS}"

srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/plot_spatial_campaign.jl \
  --analysis-dir "${ANALYSIS_DIR}" \
  --output-dir "${ANALYSIS_DIR}" \
  --geometry-label "Aztec diamond"
