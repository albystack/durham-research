#!/bin/sh
set -eu

# Reproduce the August 2026 large-order extensions. Run from the repository
# root. Each Julia runner validates and skips complete batches, so this script
# is safe to resume after interruption.
JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-4}"
export JULIA_NUM_THREADS

# Complete orders 900 and 1000 from sample ID 513 through 1000.
julia --project=aztec --startup-file=no \
  aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_extension_stage2.csv \
  --output-dir aztec/output/gamma_height_extension_stage2 \
  --base-seed 20260801

# Add the three new single-height orders.
julia --project=aztec --startup-file=no \
  aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_new_sizes.csv \
  --output-dir aztec/output/gamma_height_new_sizes \
  --base-seed 20260801

# Generate two conditional tilings for every fresh disorder environment.
julia --project=aztec --startup-file=no \
  aztec/scripts/run_double_dimer_campaign.jl \
  --config aztec/configs/double_dimer_campaign.csv \
  --output-dir aztec/output/double_dimer_campaign \
  --base-seed 20260802
