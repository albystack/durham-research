#!/bin/bash

# Replace only the five not-yet-started eight-worker high-size boosters from
# the 9 August recovery launch with eleven-worker versions.  Ultra-size and
# Aztec jobs are deliberately outside this script.

set -euo pipefail

repo="/home/fvkl37/durham-research-square-grid/research"
root="/home/fvkl37/hamilton-booster-20260809"
map="${root}/manifests/high_repack11_job_map.txt"

old_boosters=(18249843 18249845 18249847 18249849 18249851)
old_releases=(18249844 18249846 18249848 18249850 18249852)

# Do not create duplicate writers if scheduler capacity changed between the
# planning check and execution.
for job in "${old_boosters[@]}"; do
  if squeue -h -j "${job}" -t RUNNING,COMPLETING | grep -q .; then
    echo "Refusing to replace started booster ${job}" >&2
    exit 2
  fi
done

for job in "${old_boosters[@]}" "${old_releases[@]}"; do
  if squeue -h -j "${job}" 2>/dev/null | grep -q .; then
    scancel "${job}"
  fi
done

# Tasks 156 and 157 of the first high campaign were completed by the real-data
# packed smoke job.  Keep the original immutable manifest and derive a hashed
# repack manifest containing only tasks 158--850.
source_manifest="${root}/manifests/original_18243074_tasks.txt"
repack_manifest="${root}/manifests/original_18243074_tasks_repack11.txt"
sed -n '3,$p' "${source_manifest}" > "${repack_manifest}"
[[ "$(wc -l < "${repack_manifest}")" -eq 693 ]]
[[ "$(head -1 "${repack_manifest}")" == 158 ]]
[[ "$(tail -1 "${repack_manifest}")" == 850 ]]
sha256sum "${repack_manifest}" > "${repack_manifest}.sha256"

: > "${map}"
printf 'created_at=%s\n' "$(date --iso-8601=seconds)" >> "${map}"
printf 'superseded_boosters=%s\n' "${old_boosters[*]}" >> "${map}"
printf 'superseded_releases=%s\n' "${old_releases[*]}" >> "${map}"
printf 'group_size=11|allocation_memory=242G|worker_memory=22G\n' >> "${map}"

submit_release() {
  local original_job="$1"
  local booster_job="$2"
  sbatch --parsable \
    --dependency="afterok:${booster_job}" \
    --partition=test --time=00:02:00 --nodes=1 --ntasks=1 --mem=256M \
    --job-name="release-${original_job}" \
    --output="${root}/logs/release-${original_job}-%j.out" \
    --error="${root}/logs/release-${original_job}-%j.err" \
    --wrap="scontrol release ${original_job}"
}

submit_high() {
  local original_job="$1"
  local environment_model="$2"
  local distribution="$3"
  local parameter="$4"
  local output_dir="$5"
  local base_seed="$6"
  local task_file="${root}/manifests/original_${original_job}_tasks.txt"
  local count groups booster release

  if [[ "${original_job}" == 18243074 ]]; then
    task_file="${repack_manifest}"
  fi
  count=$(wc -l < "${task_file}")
  groups=$(( (count + 10) / 11 ))
  booster=$(sbatch --parsable \
    --array="1-${groups}" --nodes=1 --ntasks=11 --mem=242G \
    --export=ALL,REPOSITORY_ROOT="${repo}",TASK_FILE="${task_file}",GROUP_SIZE=11,WORKER_MEMORY=22G,ENVIRONMENT_MODEL="${environment_model}",DISTRIBUTION="${distribution}",PARAMETER="${parameter}",OUTPUT_DIR="${output_dir}",CONFIG=aztec/configs/square_grid_high_l_extension.csv,BASE_SEED="${base_seed}" \
    "${root}/scripts/square_grid_packed.slurm")
  release=$(submit_release "${original_job}" "${booster}")
  printf 'original=%s|tasks=%s|groups=%s|booster=%s|release=%s|manifest=%s\n' \
    "${original_job}" "${count}" "${groups}" "${booster}" "${release}" \
    "${task_file}" | tee -a "${map}"
}

submit_high 18243074 directed_site_iid gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k05/high 20261102
submit_high 18243076 directed_site_iid gamma 1 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k1/high 20261103
submit_high 18243078 directed_site_iid lognormal 1.048147073968205 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_lognormal_var2/high 20261104
submit_high 18243080 directed_site_iid uniform 2 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_uniform_0_2/high 20261105
submit_high 18243082 undirected_conductance gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/undirected_gamma_k05/high 20261106

sha256sum "${map}" > "${map}.sha256"
echo "High-size repack complete: ${map}"
