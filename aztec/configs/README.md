# Campaign schedules

Each CSV gives an Aztec order, the first deterministic sample ID, the number
of samples, and an atomic batch size. Files with the legacy three-column
header start sample IDs at one.

An order may appear only once in a config. To extend an existing campaign,
use a new output directory and choose `first_sample_id` immediately after the
last retained ID. When merging campaigns, sample IDs and generated seeds must
both be globally unique.

## Retained production data

- `gamma_height_pilot.csv`, `gamma_height_large_stage1.csv`, and
  `gamma_height_large_stage2.csv` produced the original single-height data.
- `gamma_height_extension_stage1.csv` and
  `gamma_height_extension_stage2.csv` brought orders 900 and 1,000 to 1,000
  observations each. `gamma_height_extension.csv` is the equivalent combined
  schedule for a fresh output directory.
- `gamma_height_new_sizes.csv` contains 512 observations at each of orders
  1,100, 1,200, and 1,300.
- `double_dimer_campaign.csv` produced the retained 28,304 paired observations.
- `double_dimer_large_continuation.csv` is the recommended next run: it adds
  another 256 pairs at each large order without reusing sample IDs.

## Fast checks

The `*_smoke.csv` schedules use matching orders 4, 8, and 12 so the complete
single/double variance-decomposition pipeline can be tested quickly. The
`*_benchmark.csv` schedules use deliberately distant sample IDs and are for
timing and memory measurements only; do not merge them into production data.
