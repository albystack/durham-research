#!/bin/bash
#SBATCH --job-name=sqgrid-temporal-ust
#SBATCH --partition=shared
#SBATCH --array=1-150
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

# Production defaults are benchmarked for temporal_iid.  Submit ordinary_ust
# with `--array=1-79`; it uses fewer deterministic batches.  No %N cap is
# embedded because the current shared QoS has no per-user job cap and broad
# one-core placement maximizes useful work per CPU-second.
set -euo pipefail

JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

: "${OUTPUT_DIR:?Set OUTPUT_DIR under /nobackup for the selected model result directory}"
: "${BASE_SEED:?Set one public master seed for both model submissions}"
: "${CAMPAIGN_MODE:?Set temporal_iid or ordinary_ust}"
CONFIG="${CONFIG:-}"

[[ "${OUTPUT_DIR}" == /nobackup/* ]] || {
  echo "OUTPUT_DIR must be under /nobackup on Hamilton8" >&2
  exit 2
}
case "${CAMPAIGN_MODE}" in
  temporal_iid)
    MODEL_ARGS=(--environment-model temporal_iid --distribution gamma --parameter 0.5)
    [[ -n "${CONFIG}" ]] || CONFIG="aztec/configs/temporal_ust_confirmation_temporal.csv"
    ;;
  ordinary_ust)
    MODEL_ARGS=(--environment-model ordinary_ust --distribution baseline --parameter 0.5)
    [[ -n "${CONFIG}" ]] || CONFIG="aztec/configs/temporal_ust_confirmation_ordinary_ust.csv"
    ;;
  *)
    echo "CAMPAIGN_MODE must be temporal_iid or ordinary_ust" >&2
    exit 2
    ;;
esac

[[ -f "${CONFIG}" ]] || { echo "missing config: ${CONFIG}" >&2; exit 2; }

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || {
  echo "Run bash hpc/setup_hamilton_environment.sh first" >&2
  exit 1
}
mkdir -p logs "${OUTPUT_DIR}"

# A job-specific sidecar supplements the runner's immutable batch and
# execution files.  It never shares a mutable filename with another array
# element or submission.
PROVENANCE_DIR="${OUTPUT_DIR}/hpc_provenance"
mkdir -p "${PROVENANCE_DIR}"
PROVENANCE_PATH="${PROVENANCE_DIR}/job_${SLURM_ARRAY_JOB_ID}_task_${SLURM_ARRAY_TASK_ID}.txt"
PROVENANCE_TMP="${PROVENANCE_PATH}.$$.tmp"
{
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID}"
  echo "slurm_array_job_id=${SLURM_ARRAY_JOB_ID}"
  echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID}"
  echo "hostname=$(hostname)"
  echo "allocated_cpus=${SLURM_CPUS_PER_TASK}"
  echo "slurm_mem_per_node=${SLURM_MEM_PER_NODE:-not_exported}"
  echo "julia_threads=${JULIA_NUM_THREADS}"
  echo "julia_version=$(${JULIA} --version)"
  echo "campaign_mode=${CAMPAIGN_MODE}"
  echo "base_seed=${BASE_SEED}"
  echo "config=${CONFIG}"
  sha256sum "${CONFIG}"
} > "${PROVENANCE_TMP}"
mv "${PROVENANCE_TMP}" "${PROVENANCE_PATH}"

echo "campaign_mode=${CAMPAIGN_MODE} task=${SLURM_ARRAY_TASK_ID} base_seed=${BASE_SEED}"
echo "config=${CONFIG} output_dir=${OUTPUT_DIR} threads=${JULIA_NUM_THREADS}"

srun --cpu-bind=cores "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_square_grid_campaign.jl \
  "${MODEL_ARGS[@]}" \
  --config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}" \
  --base-seed "${BASE_SEED}" \
  --task-id "${SLURM_ARRAY_TASK_ID}"
