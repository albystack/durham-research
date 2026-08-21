#!/bin/bash
#SBATCH --job-name=aztec-summary
#SBATCH --partition=shared
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"
root="/nobackup/fvkl37/aztec_gamma_parameters_20260808"

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs "${root}/analysis/comparison"

arguments=()
for cutoff in 128 384 512; do
  arguments+=(--entry "original_min${cutoff}=${root}/analysis/original_a0p2_b0p25_min${cutoff}")
  arguments+=(--entry "strong_min${cutoff}=${root}/analysis/strong_min${cutoff}")
  arguments+=(--entry "weak_min${cutoff}=${root}/analysis/weak_min${cutoff}")
  arguments+=(--entry "symmetric_min${cutoff}=${root}/analysis/symmetric_min${cutoff}")
done

srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/summarize_spatial_gls.jl \
  "${arguments[@]}" \
  --output-dir "${root}/analysis/comparison"
