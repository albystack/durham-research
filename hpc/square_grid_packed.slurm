#!/bin/bash
#SBATCH --job-name=sqgrid-pack
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

# Run several independent square-grid task IDs inside one Slurm allocation.
# The submission command must set --ntasks to at least GROUP_SIZE and request
# enough memory for all workers.  TASK_FILE contains one original task ID per
# line, preserving the exact task/seed identity of the original campaign.

set -euo pipefail

JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

: "${TASK_FILE:?Set TASK_FILE to an absolute task-manifest path}"
: "${GROUP_SIZE:?Set GROUP_SIZE to the number of concurrent task IDs}"
: "${WORKER_MEMORY:?Set WORKER_MEMORY (for example, 22G)}"
: "${ENVIRONMENT_MODEL:?Set ENVIRONMENT_MODEL}"
: "${DISTRIBUTION:?Set DISTRIBUTION}"
: "${PARAMETER:?Set PARAMETER}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR}"
: "${CONFIG:?Set CONFIG}"
: "${BASE_SEED:?Set BASE_SEED}"

REPOSITORY_ROOT="${REPOSITORY_ROOT:-${SLURM_SUBMIT_DIR}}"
BOOSTER_LOG_DIR="${BOOSTER_LOG_DIR:-${REPOSITORY_ROOT}/logs/packed-square-grid}"

[[ "${GROUP_SIZE}" =~ ^[1-9][0-9]*$ ]] || {
  echo "GROUP_SIZE must be a positive integer" >&2
  exit 2
}
[[ "${SLURM_ARRAY_TASK_ID:-}" =~ ^[1-9][0-9]*$ ]] || {
  echo "This wrapper must run as a one-based Slurm array" >&2
  exit 2
}
test -r "${TASK_FILE}" || { echo "Cannot read TASK_FILE=${TASK_FILE}" >&2; exit 2; }
test -f "${REPOSITORY_ROOT}/.hamilton_env/Project.toml" || {
  echo "Hamilton Julia environment is missing under ${REPOSITORY_ROOT}" >&2
  exit 2
}

first_line=$(( (SLURM_ARRAY_TASK_ID - 1) * GROUP_SIZE + 1 ))
last_line=$(( SLURM_ARRAY_TASK_ID * GROUP_SIZE ))
mapfile -t task_ids < <(sed -n "${first_line},${last_line}p" "${TASK_FILE}")

((${#task_ids[@]} > 0)) || {
  echo "No task IDs for packed group ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
}
((${#task_ids[@]} <= ${SLURM_NTASKS:-1})) || {
  echo "Packed group has ${#task_ids[@]} tasks but allocation has ${SLURM_NTASKS:-1}" >&2
  exit 2
}

mkdir -p "${BOOSTER_LOG_DIR}"
cd "${REPOSITORY_ROOT}"

for task_id in "${task_ids[@]}"; do
  [[ "${task_id}" =~ ^[1-9][0-9]*$ ]] || {
    echo "Invalid task ID in ${TASK_FILE}: ${task_id}" >&2
    exit 2
  }
done

export JULIA REPOSITORY_ROOT BOOSTER_LOG_DIR
export PACKED_FIRST_LINE="${first_line}"

srun --exact --nodes=1 --ntasks="${#task_ids[@]}" --cpus-per-task=1 \
  --mem-per-cpu="${WORKER_MEMORY}" --kill-on-bad-exit=1 \
  /bin/bash -c '
    set -euo pipefail
    line_number=$((PACKED_FIRST_LINE + SLURM_PROCID))
    task_id=$(sed -n "${line_number}p" "${TASK_FILE}")
    [[ "${task_id}" =~ ^[1-9][0-9]*$ ]]
    task_log="${BOOSTER_LOG_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_task_${task_id}.log"
    exec >"${task_log}" 2>&1
    echo "packed_rank=${SLURM_PROCID} original_task_id=${task_id}"
    exec "${JULIA}" --project=.hamilton_env --startup-file=no \
    aztec/scripts/run_square_grid_campaign.jl \
    --environment-model "${ENVIRONMENT_MODEL}" \
    --distribution "${DISTRIBUTION}" \
    --parameter "${PARAMETER}" \
    --config "${CONFIG}" \
    --output-dir "${OUTPUT_DIR}" \
    --base-seed "${BASE_SEED}" \
    --task-id "${task_id}"
  '

echo "Packed square-grid group ${SLURM_ARRAY_TASK_ID} complete (${task_ids[*]})"
