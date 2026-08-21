#!/bin/bash

# One-off, auditable launcher for the held 9 August 2026 Hamilton task
# manifests. Run from the remote repository root after both real-data packed
# smoke jobs have been submitted. Full boosters depend on both smoke jobs; each
# original array is released only after its corresponding booster succeeds.

set -euo pipefail

repo="/home/fvkl37/durham-research-square-grid/research"
root="/home/fvkl37/hamilton-booster-20260809"
square_smoke="18249791"
aztec_smoke="18249792"
smoke_dependency="afterok:${square_smoke}:${aztec_smoke}"
job_map="${root}/manifests/booster_job_map.txt"

mkdir -p "${root}/logs"
cd "${repo}"
: > "${job_map}"
printf 'created_at=%s\n' "$(date --iso-8601=seconds)" >> "${job_map}"
printf 'smoke_jobs=%s,%s\n' "${square_smoke}" "${aztec_smoke}" >> "${job_map}"

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

submit_square() {
  local original_job="$1"
  local group_size="$2"
  local allocation_memory="$3"
  local worker_memory="$4"
  local environment_model="$5"
  local distribution="$6"
  local parameter="$7"
  local output_dir="$8"
  local config="$9"
  local base_seed="${10}"
  local task_file="${root}/manifests/original_${original_job}_tasks.txt"
  local task_count groups booster release

  task_count=$(wc -l < "${task_file}")
  ((task_count > 0)) || { echo "Empty task manifest: ${task_file}" >&2; exit 2; }
  groups=$(( (task_count + group_size - 1) / group_size ))
  booster=$(sbatch --parsable \
    --dependency="${smoke_dependency}" \
    --array="1-${groups}" --nodes=1 --ntasks="${group_size}" \
    --mem="${allocation_memory}" \
    --export=ALL,REPOSITORY_ROOT="${repo}",TASK_FILE="${task_file}",GROUP_SIZE="${group_size}",WORKER_MEMORY="${worker_memory}",ENVIRONMENT_MODEL="${environment_model}",DISTRIBUTION="${distribution}",PARAMETER="${parameter}",OUTPUT_DIR="${output_dir}",CONFIG="${config}",BASE_SEED="${base_seed}" \
    "${root}/scripts/square_grid_packed.slurm")
  release=$(submit_release "${original_job}" "${booster}")
  printf 'original=%s|kind=square_grid|tasks=%s|group_size=%s|groups=%s|booster=%s|release=%s\n' \
    "${original_job}" "${task_count}" "${group_size}" "${groups}" \
    "${booster}" "${release}" | tee -a "${job_map}"
}

submit_aztec() {
  local original_job="$1"
  local output_dir="$2"
  local base_seed="$3"
  local alpha="$4"
  local beta="$5"
  local campaign_label="$6"
  local group_size=24
  local task_file="${root}/manifests/original_${original_job}_tasks.txt"
  local task_count groups booster release

  task_count=$(wc -l < "${task_file}")
  ((task_count > 0)) || { echo "Empty task manifest: ${task_file}" >&2; exit 2; }
  groups=$(( (task_count + group_size - 1) / group_size ))
  booster=$(sbatch --parsable \
    --dependency="${smoke_dependency}" \
    --array="1-${groups}" --nodes=1 --ntasks="${group_size}" --mem=48G \
    --export=ALL,REPOSITORY_ROOT="${repo}",TASK_FILE="${task_file}",GROUP_SIZE="${group_size}",WORKER_MEMORY=2G,CONFIG=aztec/configs/spatial_publication_gamma.csv,OUTPUT_DIR="${output_dir}",BASE_SEED="${base_seed}",ALPHA="${alpha}",BETA="${beta}",CAMPAIGN_LABEL="${campaign_label}" \
    "${root}/scripts/aztec_spatial_packed.slurm")
  release=$(submit_release "${original_job}" "${booster}")
  printf 'original=%s|kind=aztec|tasks=%s|group_size=%s|groups=%s|booster=%s|release=%s\n' \
    "${original_job}" "${task_count}" "${group_size}" "${groups}" \
    "${booster}" "${release}" | tee -a "${job_map}"
}

# Ultra-size square-grid arrays: five workers fit below 220 GB at the measured
# worst-case 40 GB RSS per size-6,144 process.
submit_square 18243073 5 220G 44G baseline gamma 1 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/baseline/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261101
submit_square 18243075 5 220G 44G directed_site_iid gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k05/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261102
submit_square 18243077 5 220G 44G directed_site_iid gamma 1 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k1/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261103
submit_square 18243079 5 220G 44G directed_site_iid lognormal 1.048147073968205 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_lognormal_var2/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261104
submit_square 18243081 5 220G 44G directed_site_iid uniform 2 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_uniform_0_2/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261105
submit_square 18243083 5 220G 44G undirected_conductance gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/undirected_gamma_k05/ultra \
  aztec/configs/square_grid_ultra_l_extension.csv 20261106

# High-size square-grid arrays: eleven 22 GB worker limits fit within the
# 256,000 MB Hamilton node allocation.  The measured production maximum is
# 17.41 GB, retaining more than 20% per-worker headroom.
submit_square 18243074 11 242G 22G directed_site_iid gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k05/high \
  aztec/configs/square_grid_high_l_extension.csv 20261102
submit_square 18243076 11 242G 22G directed_site_iid gamma 1 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_gamma_k1/high \
  aztec/configs/square_grid_high_l_extension.csv 20261103
submit_square 18243078 11 242G 22G directed_site_iid lognormal 1.048147073968205 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_lognormal_var2/high \
  aztec/configs/square_grid_high_l_extension.csv 20261104
submit_square 18243080 11 242G 22G directed_site_iid uniform 2 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/directed_uniform_0_2/high \
  aztec/configs/square_grid_high_l_extension.csv 20261105
submit_square 18243082 11 242G 22G undirected_conductance gamma 0.5 \
  /nobackup/fvkl37/square_grid_high_l_20260808/production/undirected_gamma_k05/high \
  aztec/configs/square_grid_high_l_extension.csv 20261106

submit_aztec 18243433 \
  /nobackup/fvkl37/aztec_gamma_parameters_20260808/weak_a0p8_b1p0 \
  20260912 0.8 1.0 aztec_gamma_a0p8_b1p0
submit_aztec 18243434 \
  /nobackup/fvkl37/aztec_gamma_parameters_20260808/symmetric_a0p25_b0p25 \
  20260913 0.25 0.25 aztec_gamma_a0p25_b0p25

sha256sum "${job_map}" > "${job_map}.sha256"
echo "Booster launch complete: ${job_map}"
