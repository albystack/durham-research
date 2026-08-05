#!/bin/sh
set -eu

# Exercise every public command-line stage on tiny schedules. Unit tests cover
# the mathematical kernels; this script catches broken argument names, CSV
# schemas, atomic resumption, merge logic, analysis wiring, and SVG generation.
repository_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
smoke_root=$(mktemp -d "${TMPDIR:-/tmp}/aztec-smoke.XXXXXX")
trap 'rm -rf "$smoke_root"' EXIT HUP INT TERM

cd "$repository_root"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-2}"

height_batches="$smoke_root/height_batches"
height_csv="$smoke_root/height.csv"
height_analysis="$smoke_root/height_analysis"

julia --project=aztec --startup-file=no aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_smoke.csv \
  --output-dir "$height_batches" \
  --base-seed 101

# The same filenames with a different campaign seed must be rejected, not
# silently treated as complete. The failure is expected and therefore muted.
if julia --project=aztec --startup-file=no \
  aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_smoke.csv \
  --output-dir "$height_batches" \
  --base-seed 999 >/dev/null 2>&1
then
  echo "height runner accepted a batch generated with a different seed" >&2
  exit 1
fi

# A second identical invocation must validate and skip all complete batches.
julia --project=aztec --startup-file=no aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_smoke.csv \
  --output-dir "$height_batches" \
  --base-seed 101

julia --project=aztec --startup-file=no aztec/scripts/merge_height_batches.jl \
  --inputs "$height_batches" \
  --output "$height_csv"

julia --project=aztec --startup-file=no aztec/scripts/analyze_height_campaign.jl \
  --results-dir "$height_csv" \
  --output-dir "$height_analysis" \
  --bootstrap-reps 50 \
  --bootstrap-seed 102 \
  --min-order 4

julia --project=aztec --startup-file=no aztec/scripts/plot_height_campaign.jl \
  --analysis-dir "$height_analysis" \
  --output-dir "$height_analysis"

pair_batches="$smoke_root/pair_batches"
pair_csv="$smoke_root/pairs.csv"
pair_analysis="$smoke_root/pair_analysis"

julia --project=aztec --startup-file=no aztec/scripts/run_double_dimer_campaign.jl \
  --config aztec/configs/double_dimer_smoke.csv \
  --output-dir "$pair_batches" \
  --base-seed 103

if julia --project=aztec --startup-file=no \
  aztec/scripts/run_double_dimer_campaign.jl \
  --config aztec/configs/double_dimer_smoke.csv \
  --output-dir "$pair_batches" \
  --base-seed 999 >/dev/null 2>&1
then
  echo "double runner accepted a batch generated with a different seed" >&2
  exit 1
fi

julia --project=aztec --startup-file=no aztec/scripts/merge_double_dimer_batches.jl \
  --inputs "$pair_batches" \
  --output "$pair_csv"

julia --project=aztec --startup-file=no aztec/scripts/analyze_double_dimer_campaign.jl \
  --paired-results "$pair_csv" \
  --single-results "$height_csv" \
  --output-dir "$pair_analysis" \
  --bootstrap-reps 50 \
  --bootstrap-seed 104 \
  --min-order 4

julia --project=aztec --startup-file=no aztec/scripts/plot_double_dimer_campaign.jl \
  --analysis-dir "$pair_analysis" \
  --output-dir "$pair_analysis"

spatial_gamma="$smoke_root/spatial_gamma"
spatial_uniform="$smoke_root/spatial_uniform"
spatial_analysis="$smoke_root/spatial_analysis"

julia --project=aztec --startup-file=no aztec/scripts/run_spatial_campaign.jl \
  --model gamma \
  --config aztec/configs/spatial_smoke.csv \
  --output-dir "$spatial_gamma" \
  --base-seed 107

if julia --project=aztec --startup-file=no \
  aztec/scripts/run_spatial_campaign.jl \
  --model gamma \
  --config aztec/configs/spatial_smoke.csv \
  --output-dir "$spatial_gamma" \
  --base-seed 999 >/dev/null 2>&1
then
  echo "spatial runner accepted a batch generated with a different seed" >&2
  exit 1
fi

julia --project=aztec --startup-file=no aztec/scripts/run_spatial_campaign.jl \
  --model uniform \
  --config aztec/configs/spatial_smoke.csv \
  --output-dir "$spatial_uniform" \
  --base-seed 108

julia --project=aztec --startup-file=no aztec/scripts/analyze_spatial_campaign.jl \
  --gamma-results "$spatial_gamma" \
  --uniform-results "$spatial_uniform" \
  --output-dir "$spatial_analysis" \
  --bootstrap-reps 50 \
  --bootstrap-seed 109 \
  --min-order 64 \
  --holdout-orders 2

julia --project=aztec --startup-file=no aztec/scripts/plot_spatial_campaign.jl \
  --analysis-dir "$spatial_analysis" \
  --output-dir "$spatial_analysis"

# Also exercise the two illustrated one-off samplers at a tiny order.
julia --project=aztec --startup-file=no aztec/scripts/run_gamma_disordered.jl \
  --order 6 --seed 105 --output-dir "$smoke_root/gamma_example"
julia --project=aztec --startup-file=no aztec/scripts/run_random_weights.jl \
  --order 6 --seed 106 --output-dir "$smoke_root/uniform_example"

for artifact in \
  "$height_analysis/height_fits.txt" \
  "$height_analysis/center_height_variance.svg" \
  "$pair_analysis/double_dimer_fits.txt" \
  "$pair_analysis/variance_component_fits.txt" \
  "$pair_analysis/double_dimer_variance_fits.svg" \
  "$pair_analysis/disorder_covariance_fits.svg" \
  "$pair_analysis/variance_decomposition.svg" \
  "$spatial_analysis/spatial_summary.csv" \
  "$spatial_analysis/spatial_model_comparison.csv" \
  "$spatial_analysis/spatial_weighted_model_comparison.csv" \
  "$spatial_analysis/spatial_pooled_model_comparison.csv" \
  "$spatial_analysis/spatial_cutoff_sensitivity.csv" \
  "$spatial_analysis/spatial_analysis_report.md" \
  "$spatial_analysis/gamma_disorder_spatial_fits.svg" \
  "$spatial_analysis/uniform_control_spatial_fits.svg" \
  "$spatial_analysis/spatial_variance_decomposition.svg" \
  "$spatial_analysis/pooled_log2_coefficients.svg" \
  "$smoke_root/gamma_example/tiling_mathematica_style.svg" \
  "$smoke_root/uniform_example/tiling_mathematica_style.svg"
do
  test -s "$artifact"
done


square_baseline="$smoke_root/square_baseline"
square_directed="$smoke_root/square_directed"
square_baseline_csv="$smoke_root/square_baseline.csv"
square_directed_csv="$smoke_root/square_directed.csv"

julia --project=aztec --startup-file=no aztec/scripts/run_square_grid_campaign.jl \
  --environment-model baseline \
  --config aztec/configs/square_grid_smoke.csv \
  --output-dir "$square_baseline" \
  --base-seed 110

julia --project=aztec --startup-file=no aztec/scripts/run_square_grid_campaign.jl \
  --environment-model directed_site_iid \
  --distribution gamma \
  --parameter 0.5 \
  --config aztec/configs/square_grid_smoke.csv \
  --output-dir "$square_directed" \
  --base-seed 111

if julia --project=aztec --startup-file=no \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model directed_site_iid \
  --parameter 0.5 \
  --config aztec/configs/square_grid_smoke.csv \
  --output-dir "$square_directed" \
  --base-seed 999 >/dev/null 2>&1
then
  echo "square-grid runner accepted batches generated with a different seed" >&2
  exit 1
fi

julia --project=aztec --startup-file=no aztec/scripts/merge_square_grid_batches.jl \
  --inputs "$square_baseline" \
  --output "$square_baseline_csv"

julia --project=aztec --startup-file=no aztec/scripts/merge_square_grid_batches.jl \
  --inputs "$square_directed" \
  --output "$square_directed_csv"

test -s "$square_baseline_csv"
test -s "$square_directed_csv"
test -s "$square_baseline/campaign_metadata.txt"
test -s "$square_directed/campaign_metadata.txt"

echo "All Aztec and square-grid command-line smoke workflows passed."
