#!/bin/bash
#SBATCH --job-name=glauber-dimer-analysis
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
JULIA="/apps/developers/compilers/julia/1.10.4/1/default/bin.wrapped/julia"
export JULIA_NUM_THREADS=1
: "${PILOT_DIR:?Set PILOT_DIR to the completed /nobackup pilot directory}"
: "${ANALYSIS_DIR:?Set ANALYSIS_DIR under /nobackup}"

[[ "${PILOT_DIR}" == /nobackup/* && "${ANALYSIS_DIR}" == /nobackup/* ]] || {
  echo "PILOT_DIR and ANALYSIS_DIR must be under /nobackup on Hamilton8" >&2
  exit 2
}
cd "${SLURM_SUBMIT_DIR}"
test -f .hamilton_env/Project.toml || { echo "Run bash hpc/setup_hamilton_environment.sh first" >&2; exit 1; }
mkdir -p logs "${ANALYSIS_DIR}"
srun "$JULIA" --project=.hamilton_env --startup-file=no \
  aztec/scripts/analyze_glauber_square_grid_pilot.jl \
  --pilot-dir "${PILOT_DIR}" --analysis-dir "${ANALYSIS_DIR}"
