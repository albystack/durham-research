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
  "$smoke_root/gamma_example/tiling_mathematica_style.svg" \
  "$smoke_root/uniform_example/tiling_mathematica_style.svg"
do
  test -s "$artifact"
done

echo "All Aztec command-line smoke workflows passed."
