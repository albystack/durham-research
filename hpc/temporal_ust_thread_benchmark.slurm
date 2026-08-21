#!/bin/bash
#SBATCH --job-name=sqgrid-temporal-thread-bench
#SBATCH --partition=shared
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Submit this tiny temporal-only scaling benchmark four times, overriding
# --cpus-per-task with 1, 2, 4, and 8.  It uses the same L=1536 sample IDs
# for every CPU count, but a distinct output directory for each configuration.
set -euo pipefail

JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

: "${BENCHMARK_ROOT:?Set BENCHMARK_ROOT under /nobackup}"
: "${BASE_SEED:?Set the fixed public benchmark seed}"
CONFIG="${CONFIG:-aztec/configs/temporal_ust_thread_benchmark.csv}"
BENCHMARK_LABEL="${BENCHMARK_LABEL:-cpus_${SLURM_CPUS_PER_TASK}}"
OUTPUT_DIR="${BENCHMARK_ROOT}/${BENCHMARK_LABEL}"

[[ "${BENCHMARK_ROOT}" == /nobackup/* ]] || {
  echo "BENCHMARK_ROOT must be under /nobackup on Hamilton8" >&2
  exit 2
}
[[ -f "${CONFIG}" ]] || { echo "missing config: ${CONFIG}" >&2; exit 2; }

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || {
  echo "Run bash hpc/setup_hamilton_environment.sh first" >&2
  exit 1
}
mkdir -p logs "${OUTPUT_DIR}"

echo "benchmark_label=${BENCHMARK_LABEL} cpus=${SLURM_CPUS_PER_TASK} seed=${BASE_SEED}"
echo "config=${CONFIG} output_dir=${OUTPUT_DIR}"
echo "git_commit=$(git rev-parse HEAD)"
"${JULIA}" --version

srun --cpu-bind=cores "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model temporal_iid \
  --distribution gamma \
  --parameter 0.5 \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}" \
  --base-seed "${BASE_SEED}" \
  --task-id 0
