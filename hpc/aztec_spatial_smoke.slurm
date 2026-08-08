#!/bin/bash
#SBATCH --job-name=aztec-spatial-smoke
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || { echo "Run bash hpc/setup_hamilton_environment.sh first" >&2; exit 1; }
smoke_root=$(mktemp -d "${TMPDIR:-/tmp}/aztec-spatial-smoke.XXXXXX")
trap 'rm -rf "$smoke_root"' EXIT HUP INT TERM
gamma_dir="${smoke_root}/gamma"
uniform_dir="${smoke_root}/uniform"
analysis_dir="${smoke_root}/analysis"

"$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model gamma --config aztec/configs/spatial_smoke.csv \
  --output-dir "$gamma_dir" --base-seed 20260921 \
  --alpha 0.8 --beta 1.0 --campaign-label aztec_gamma_smoke \
  --task-id 1

# Resumption must validate and skip the existing atomic batch.
"$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model gamma --config aztec/configs/spatial_smoke.csv \
  --output-dir "$gamma_dir" --base-seed 20260921 \
  --alpha 0.8 --beta 1.0 --campaign-label aztec_gamma_smoke \
  --task-id 1

# A parameter change under the same label must be rejected.
if "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model gamma --config aztec/configs/spatial_smoke.csv \
  --output-dir "$gamma_dir" --base-seed 20260921 \
  --alpha 0.7 --beta 1.0 --campaign-label aztec_gamma_smoke \
  --task-id 1 >/dev/null 2>&1
then
  echo "parameter-mismatched Aztec batch was accepted" >&2
  exit 1
fi

# Finish the Gamma schedule, generate the control, and exercise analysis/plots.
"$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model gamma --config aztec/configs/spatial_smoke.csv \
  --output-dir "$gamma_dir" --base-seed 20260921 \
  --alpha 0.8 --beta 1.0 --campaign-label aztec_gamma_smoke
"$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model uniform --config aztec/configs/spatial_smoke.csv \
  --output-dir "$uniform_dir" --base-seed 20260922
"$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/analyze_spatial_campaign.jl \
  --gamma-results "$gamma_dir" --uniform-results "$uniform_dir" \
  --gamma-model aztec_gamma_smoke --uniform-model uniform \
  --output-dir "$analysis_dir" --bootstrap-reps 50 \
  --bootstrap-seed 20260923 --min-order 64 --holdout-orders 2
"$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/plot_spatial_campaign.jl \
  --analysis-dir "$analysis_dir" --output-dir "$analysis_dir" \
  --geometry-label "Aztec diamond smoke"

echo "Aztec spatial runner smoke passed"
