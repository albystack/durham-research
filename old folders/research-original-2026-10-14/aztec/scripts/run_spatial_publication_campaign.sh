#!/bin/sh
set -eu

# Publication-focused experiment.  Both runners are resumable and validate
# every existing atomic batch before skipping it.
JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-8}"
export JULIA_NUM_THREADS

julia --project=aztec --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model gamma \
  --config aztec/configs/spatial_publication_gamma.csv \
  --output-dir aztec/output/spatial_publication/gamma \
  --fractions 1/32,1/16,1/8,1/4 \
  --base-seed 20260803 \
  --alpha 0.2 --beta 0.25

julia --project=aztec --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model uniform \
  --config aztec/configs/spatial_publication_uniform.csv \
  --output-dir aztec/output/spatial_publication/uniform \
  --fractions 1/32,1/16,1/8,1/4 \
  --base-seed 20260804

julia --project=aztec --startup-file=no \
  aztec/scripts/analyze_spatial_campaign.jl \
  --gamma-results aztec/output/spatial_publication/gamma \
  --uniform-results aztec/output/spatial_publication/uniform \
  --output-dir aztec/output/spatial_publication/analysis \
  --bootstrap-reps 10000 \
  --bootstrap-seed 20260805 \
  --min-order 128 \
  --holdout-orders 2

julia --project=aztec --startup-file=no \
  aztec/scripts/plot_spatial_campaign.jl \
  --analysis-dir aztec/output/spatial_publication/analysis \
  --output-dir aztec/output/spatial_publication/analysis
