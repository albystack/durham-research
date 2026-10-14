#!/bin/bash
#SBATCH --job-name=aztec-robust-refine
#SBATCH --partition=shared
#SBATCH --array=0-5%2
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from verified isolated source snapshot}"
: "${REFINEMENT_OUTPUT:?Set isolated /nobackup output directory}"
[[ "$REFINEMENT_OUTPUT" == /nobackup/* ]] || exit 2
comparison_julia="${RESEARCH_JULIA:-/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia}"
orders=(16 32 48 64 96 128)
export JULIA_DEPOT_PATH="$REFINEMENT_OUTPUT/julia-depot-$SLURM_ARRAY_TASK_ID"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$JULIA_DEPOT_PATH"
srun "$comparison_julia" --startup-file=no --project=aztec aztec/scripts/run_robust_refinement.jl \
  --output "$REFINEMENT_OUTPUT/n${orders[$SLURM_ARRAY_TASK_ID]}" \
  --size "${orders[$SLURM_ARRAY_TASK_ID]}" --environments 2048 --seed 2026100102
