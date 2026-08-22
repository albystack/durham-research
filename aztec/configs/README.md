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
- `spatial_publication_gamma.csv` and `spatial_publication_uniform.csv` are the
  publication-focused paired-replica schedules. They measure symmetric
  central-row increments at separations `L/32`, `L/16`, `L/8`, and `L/4`.
  The Gamma campaign shares one environment between replicas; the uniform
  campaign is the no-disorder control.
- `aztec_gamma_parameter_pilot.csv` is the three-task safety pilot used before
  launching an additional Gamma-parameter campaign on Hamilton.

## Fast checks

The original `*_smoke.csv` schedules use matching orders 4, 8, and 12.
`spatial_smoke.csv` uses five geometrically separated orders so the nested
quadratic-in-log fit and held-out prediction can also be tested. The
`*_benchmark.csv` schedules use deliberately distant sample IDs and are for
timing and memory measurements only; do not merge them into production data.

## Square-grid schedules

- `square_grid_smoke.csv`: two tiny sizes for CLI and Hamilton validation.
- `square_grid_pilot.csv`: 50 resumable tasks across `L=16,32,64,128,256` with 570 environments per model.
  Run it separately for the baseline and directed Gamma environments.
- `square_grid_robustness_pilot.csv`: six-size validation run for each new
  disorder law before production.
- `square_grid_robustness.csv`: nine sizes through `L=2048`, 38,600
  environments per disorder law, and 899 atomic tasks. Batch sizes use all
  four worker threads while keeping the largest tasks short. It is designed to
  reuse the completed 38,600-environment baseline while comparing directed
  Gamma `k=1`, directed mean-one lognormal `sigma=sqrt(log(3))`, directed
  `Uniform(0,2)`, and undirected Gamma `k=0.5` without pooling identities.
- `square_grid_high_l_pilot.csv`: one memory/runtime probe at each of
  `L=3072,4096,5120,6144` before releasing expensive extensions.
- `square_grid_high_l_extension.csv`: 850 independent high-size environments
  per law across `L=2560,3072,4096`.
- `square_grid_ultra_l_extension.csv`: 120 exploratory environments per law at
  `L=5120,6144`; these are held-out stress tests, not precision estimates.

## Direct Glauber schedules

- `glauber_square_grid_L2_calibration.csv`: exact-enumeration calibration at
  the smallest nontrivial size.
- `glauber_square_grid_mixing_ladder_*.csv`: successive ladder and run-length
  pilots; retained to document failed as well as successful calibration.
- `glauber_square_grid_tempered_pilot.csv` and
  `glauber_square_grid_tempered_stress.csv`: corrected parallel-tempering
  pilots at `L<=6`.
- `glauber_square_grid_tempered_L8_pilot.csv`: four-environment extension of
  the validated small-size schedule.
- `glauber_square_grid_production_smoke.csv`: tiny end-to-end check of the
  production schema and wrapper.
- `glauber_square_grid_production_core.csv`: Gamma production at `L=2--12`.
- `glauber_square_grid_production_frontier.csv`: isolated Gamma production at
  `L=16,20`, kept separate because of cost and mixing risk.
- `glauber_square_grid_production_control.csv`: all-one negative control.

Direct Glauber configs specify one environment per restart-safe task. Every
environment retains two independent extremal-start chains and their full
diagnostic records.
