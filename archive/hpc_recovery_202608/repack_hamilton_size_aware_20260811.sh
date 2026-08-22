#!/bin/bash

# Replace only held Hamilton recovery work with short, size-aware allocations
# on the shared partition.  The original task IDs, seeds, output paths, and
# downstream analysis dependencies remain unchanged.

set -euo pipefail

repo="/home/fvkl37/durham-research-square-grid/research"
recovery_root="/home/fvkl37/hamilton-booster-20260809"
manifest_root="${recovery_root}/manifests"
wrapper="${recovery_root}/scripts/square_grid_packed.slurm"
map="${manifest_root}/size_aware_shared_job_map_20260811.txt"
lock="${recovery_root}/size-aware-shared-20260811.lock"

mkdir "${lock}" 2>/dev/null || {
  echo "Refusing to run: lock already exists at ${lock}" >&2
  exit 2
}
trap 'rmdir "${lock}" 2>/dev/null || true' EXIT

[[ ! -e "${map}" ]] || {
  echo "Refusing to overwrite existing job map: ${map}" >&2
  exit 2
}
test -r "${wrapper}" || { echo "Missing wrapper: ${wrapper}" >&2; exit 2; }
test -f "${repo}/.hamilton_env/Project.toml" || {
  echo "Missing Hamilton Julia environment under ${repo}" >&2
  exit 2
}

old_boosters=(
  18249839 18249841
  18250662 18250664 18250666 18250668 18250670
)
old_releases=(
  18249838 18249840 18249842
  18250663 18250665 18250667 18250669 18250671
)

# The old lognormal-ultra array has groups 1--13 complete, group 14 in flight,
# and groups 15--24 held.  Every other booster below must be wholly unstarted.
for job in "${old_boosters[@]}"; do
  if squeue -h -j "${job}" -t RUNNING,COMPLETING | grep -q .; then
    echo "Refusing to replace started booster ${job}" >&2
    exit 2
  fi
done

for job in "${old_boosters[@]}"; do
  pending=$(squeue -h -j "${job}" -t PENDING -o '%r' | sort -u)
  [[ "${pending}" == "JobHeldUser" ]] || {
    echo "Booster ${job} is not exclusively user-held (reasons: ${pending:-none})" >&2
    exit 2
  }
done

lognormal_state=$(sacct -X -n -P -j 18249837_14 --format=State | sed -n '1p')
case "${lognormal_state}" in
  RUNNING*|COMPLETED*) ;;
  *)
    echo "Unexpected state for protected in-flight job 18249837_14: ${lognormal_state}" >&2
    exit 2
    ;;
esac

completed_prefix=$(sacct -X -n -P -j 18249837 --format=JobID,State |
  awk -F'|' '$1 ~ /^18249837_([1-9]|1[0-3])$/ && $2 ~ /^COMPLETED/ {n++} END {print n+0}')
[[ "${completed_prefix}" -eq 13 ]] || {
  echo "Expected 13 completed lognormal-ultra groups, found ${completed_prefix}" >&2
  exit 2
}

create_piece() {
  local source="$1"
  local destination="$2"
  local minimum="$3"
  local maximum="$4"
  local expected="$5"

  test -r "${source}" || { echo "Missing source manifest: ${source}" >&2; exit 2; }
  awk -v lo="${minimum}" -v hi="${maximum}" \
    '$1 >= lo && $1 <= hi {print $1}' "${source}" > "${destination}"
  [[ "$(wc -l < "${destination}")" -eq "${expected}" ]] || {
    echo "Wrong task count in ${destination}" >&2
    exit 2
  }
  [[ "$(head -n 1 "${destination}")" -eq "${minimum}" ]] || {
    echo "Wrong first task in ${destination}" >&2
    exit 2
  }
  [[ "$(tail -n 1 "${destination}")" -eq "${maximum}" ]] || {
    echo "Wrong last task in ${destination}" >&2
    exit 2
  }
  sha256sum "${destination}" > "${destination}.sha256"
}

ultra_lognormal_source="${manifest_root}/original_18243079_tasks.txt"
ultra_uniform_source="${manifest_root}/original_18243081_tasks.txt"
ultra_undirected_source="${manifest_root}/original_18243083_tasks.txt"

create_piece "${ultra_lognormal_source}" \
  "${manifest_root}/sizeaware_18243079_L5120_tasks_071_120.txt" 71 120 50
create_piece "${ultra_uniform_source}" \
  "${manifest_root}/sizeaware_18243081_L6144_tasks_001_040.txt" 1 40 40
create_piece "${ultra_uniform_source}" \
  "${manifest_root}/sizeaware_18243081_L5120_tasks_041_120.txt" 41 120 80
create_piece "${ultra_undirected_source}" \
  "${manifest_root}/sizeaware_18243083_L6144_tasks_001_040.txt" 1 40 40
create_piece "${ultra_undirected_source}" \
  "${manifest_root}/sizeaware_18243083_L5120_tasks_041_120.txt" 41 120 80

high_originals=(18243074 18243076 18243078 18243080 18243082)
for original in "${high_originals[@]}"; do
  source="${manifest_root}/original_${original}_tasks.txt"
  if [[ "${original}" == 18243074 ]]; then
    source="${manifest_root}/original_18243074_tasks_repack11.txt"
    create_piece "${source}" \
      "${manifest_root}/sizeaware_${original}_L3072_tasks_158_450.txt" 158 450 293
  else
    create_piece "${source}" \
      "${manifest_root}/sizeaware_${original}_L4096_tasks_001_150.txt" 1 150 150
    create_piece "${source}" \
      "${manifest_root}/sizeaware_${original}_L3072_tasks_151_450.txt" 151 450 300
  fi
  create_piece "${source}" \
    "${manifest_root}/sizeaware_${original}_L2560_tasks_451_850.txt" 451 850 400
done

# Validate that every replacement manifest is an exact, non-overlapping
# partition of its corresponding pending source manifest.
validate_union() {
  local source="$1"
  shift
  local combined
  combined=$(mktemp "${manifest_root}/sizeaware_union.XXXXXX")
  sort -n "$@" > "${combined}"
  if ! cmp -s "${source}" "${combined}"; then
    echo "Replacement manifests do not exactly cover ${source}" >&2
    rm -f "${combined}"
    exit 2
  fi
  rm -f "${combined}"
}

expected_lognormal=$(mktemp "${manifest_root}/sizeaware_expected.XXXXXX")
sed -n '71,120p' "${ultra_lognormal_source}" > "${expected_lognormal}"
validate_union "${expected_lognormal}" \
  "${manifest_root}/sizeaware_18243079_L5120_tasks_071_120.txt"
rm -f "${expected_lognormal}"
validate_union "${ultra_uniform_source}" \
  "${manifest_root}/sizeaware_18243081_L6144_tasks_001_040.txt" \
  "${manifest_root}/sizeaware_18243081_L5120_tasks_041_120.txt"
validate_union "${ultra_undirected_source}" \
  "${manifest_root}/sizeaware_18243083_L6144_tasks_001_040.txt" \
  "${manifest_root}/sizeaware_18243083_L5120_tasks_041_120.txt"
validate_union "${manifest_root}/original_18243074_tasks_repack11.txt" \
  "${manifest_root}/sizeaware_18243074_L3072_tasks_158_450.txt" \
  "${manifest_root}/sizeaware_18243074_L2560_tasks_451_850.txt"
for original in 18243076 18243078 18243080 18243082; do
  validate_union "${manifest_root}/original_${original}_tasks.txt" \
    "${manifest_root}/sizeaware_${original}_L4096_tasks_001_150.txt" \
    "${manifest_root}/sizeaware_${original}_L3072_tasks_151_450.txt" \
    "${manifest_root}/sizeaware_${original}_L2560_tasks_451_850.txt"
done

: > "${map}"
printf 'created_at=%s\n' "$(date --iso-8601=seconds)" >> "${map}"
printf 'partition=shared|strategy=size_aware|wrapper=%s\n' "${wrapper}" >> "${map}"
printf 'protected_old_job=18249837_14|tasks=66-70|state=%s\n' \
  "${lognormal_state}" >> "${map}"
printf 'cancelled_old_boosters=18249837_[15-24] %s\n' \
  "${old_boosters[*]}" >> "${map}"
printf 'cancelled_old_releases=%s\n' "${old_releases[*]}" >> "${map}"

# Remove only held work.  The running lognormal group 14 is deliberately not
# part of this cancellation and remains a dependency of its replacement.
scancel '18249837_[15-24]'
for job in "${old_boosters[@]}" "${old_releases[@]}"; do
  if squeue -h -j "${job}" 2>/dev/null | grep -q .; then
    scancel "${job}"
  fi
done

submit_piece() {
  local label="$1"
  local task_file="$2"
  local group_size="$3"
  local allocation_memory="$4"
  local worker_memory="$5"
  local walltime="$6"
  local environment_model="$7"
  local distribution="$8"
  local parameter="$9"
  local output_dir="${10}"
  local config="${11}"
  local base_seed="${12}"
  local task_count groups job

  task_count=$(wc -l < "${task_file}")
  groups=$(( (task_count + group_size - 1) / group_size ))
  job=$(sbatch --parsable \
    --partition=shared --time="${walltime}" \
    --array="1-${groups}" --nodes=1 --ntasks="${group_size}" \
    --mem="${allocation_memory}" \
    --job-name="sq-${label}" \
    --export=ALL,REPOSITORY_ROOT="${repo}",TASK_FILE="${task_file}",GROUP_SIZE="${group_size}",WORKER_MEMORY="${worker_memory}",ENVIRONMENT_MODEL="${environment_model}",DISTRIBUTION="${distribution}",PARAMETER="${parameter}",OUTPUT_DIR="${output_dir}",CONFIG="${config}",BASE_SEED="${base_seed}" \
    "${wrapper}")
  printf 'piece=%s|tasks=%s|group_size=%s|groups=%s|mem=%s|worker_mem=%s|time=%s|job=%s|manifest=%s\n' \
    "${label}" "${task_count}" "${group_size}" "${groups}" \
    "${allocation_memory}" "${worker_memory}" "${walltime}" "${job}" \
    "${task_file}" | tee -a "${map}" >&2
  printf '%s' "${job}"
}

submit_release() {
  local original_job="$1"
  shift
  local dependency job
  dependency=$(IFS=:; echo "$*")
  job=$(sbatch --parsable \
    --dependency="afterok:${dependency}" \
    --partition=test --time=00:02:00 --nodes=1 --ntasks=1 --mem=256M \
    --job-name="release-${original_job}" \
    --output="${recovery_root}/logs/release-${original_job}-%j.out" \
    --error="${recovery_root}/logs/release-${original_job}-%j.err" \
    --wrap="scontrol release ${original_job}")
  printf 'release_original=%s|dependencies=%s|job=%s\n' \
    "${original_job}" "${dependency}" "${job}" | tee -a "${map}" >&2
}

# Ultra-size replacements: keep measured L=6144 packing at five workers and
# raise only L=5120 to seven workers with at least 20% measured RSS headroom.
lognormal_5120=$(submit_piece logn-5120 \
  "${manifest_root}/sizeaware_18243079_L5120_tasks_071_120.txt" \
  7 238G 34G 00:45:00 directed_site_iid lognormal 1.048147073968205 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_lognormal_var2/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261104)
submit_release 18243079 18249837_14 "${lognormal_5120}"

uniform_6144=$(submit_piece unif-6144 \
  "${manifest_root}/sizeaware_18243081_L6144_tasks_001_040.txt" \
  5 220G 44G 01:00:00 directed_site_iid uniform 2 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_uniform_0_2/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261105)
uniform_5120=$(submit_piece unif-5120 \
  "${manifest_root}/sizeaware_18243081_L5120_tasks_041_120.txt" \
  7 238G 34G 00:45:00 directed_site_iid uniform 2 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_uniform_0_2/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261105)
submit_release 18243081 "${uniform_6144}" "${uniform_5120}"

undirected_6144=$(submit_piece undir-6144 \
  "${manifest_root}/sizeaware_18243083_L6144_tasks_001_040.txt" \
  5 220G 44G 01:00:00 undirected_conductance gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/undirected_gamma_k05/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261106)
undirected_5120=$(submit_piece undir-5120 \
  "${manifest_root}/sizeaware_18243083_L5120_tasks_041_120.txt" \
  7 238G 34G 00:45:00 undirected_conductance gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/undirected_gamma_k05/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261106)
submit_release 18243083 "${undirected_6144}" "${undirected_5120}"

submit_high_model() {
  local original="$1"
  local short_label="$2"
  local environment_model="$3"
  local distribution="$4"
  local parameter="$5"
  local output_dir="$6"
  local base_seed="$7"
  local jobs=()

  if [[ "${original}" != 18243074 ]]; then
    jobs+=("$(submit_piece "${short_label}-4096" \
      "${manifest_root}/sizeaware_${original}_L4096_tasks_001_150.txt" \
      12 240G 20G 00:30:00 "${environment_model}" "${distribution}" "${parameter}" \
      "${output_dir}" aztec/configs/square_grid_high_l_extension.csv "${base_seed}")")
    l3072_manifest="${manifest_root}/sizeaware_${original}_L3072_tasks_151_450.txt"
  else
    l3072_manifest="${manifest_root}/sizeaware_${original}_L3072_tasks_158_450.txt"
  fi
  jobs+=("$(submit_piece "${short_label}-3072" "${l3072_manifest}" \
    20 240G 12G 00:20:00 "${environment_model}" "${distribution}" "${parameter}" \
    "${output_dir}" aztec/configs/square_grid_high_l_extension.csv "${base_seed}")")
  jobs+=("$(submit_piece "${short_label}-2560" \
    "${manifest_root}/sizeaware_${original}_L2560_tasks_451_850.txt" \
    30 240G 8G 00:20:00 "${environment_model}" "${distribution}" "${parameter}" \
    "${output_dir}" aztec/configs/square_grid_high_l_extension.csv "${base_seed}")")
  submit_release "${original}" "${jobs[@]}"
}

submit_high_model 18243074 dg05 directed_site_iid gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k05/high 20261102
submit_high_model 18243076 dg1 directed_site_iid gamma 1 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k1/high 20261103
submit_high_model 18243078 logn directed_site_iid lognormal 1.048147073968205 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_lognormal_var2/high 20261104
submit_high_model 18243080 unif directed_site_iid uniform 2 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_uniform_0_2/high 20261105
submit_high_model 18243082 undir undirected_conductance gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/undirected_gamma_k05/high 20261106

sha256sum "${map}" > "${map}.sha256"
echo "Size-aware shared-partition repack complete: ${map}"
