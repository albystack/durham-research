# Direct square-grid Glauber environment blocks

These compact tables are the retained statistical input for the completed
frozen-edge weighted-dimer production campaign.

## Campaigns

| Arm | Sizes | Environments | Base seed |
|---|---|---:|---:|
| Gamma core | 2, 4, 6, 8, 10, 12 | 864 | 2026082201 |
| Gamma frontier | 16, 20 | 96 | 2026082201 |
| All-one control | 4, 6, 8, 10, 12, 16 | 352 | 2026082202 |

Each row of an `*_environment_blocks.csv` file is one independent frozen edge
environment. It contains the two conditionally independent chain means, the
within-environment conditional-variance estimate, ESS, extremal-start gap, and
exchange diagnostics. Each chain contributed 2,000 retained central-height
values to the raw analysis.

The `*_size_summary.csv` files were generated directly from the raw traces and
include the independent direct-annealed-variance identity check. The raw MCMC
traces are retained in external project storage because they are too large for
Git.

Regenerate the scaling tables with:

```bash
julia --project=aztec aztec/scripts/analyze_glauber_square_grid_scaling.jl \
  --gamma-blocks aztec/data/glauber_square_grid_20260822/gamma_core_environment_blocks.csv,aztec/data/glauber_square_grid_20260822/gamma_frontier_environment_blocks.csv \
  --control-blocks aztec/data/glauber_square_grid_20260822/control_environment_blocks.csv \
  --output-dir aztec/output/glauber_square_grid_20260822 \
  --bootstrap-reps 5000 \
  --bootstrap-seed 20260822 \
  --cutoffs 2,4,6,8 \
  --diagnostic-threshold 4
```

See [`../../results/glauber_square_grid_20260822/`](../../results/glauber_square_grid_20260822/)
for the retained MCMC output and interpretation. Determinantal moments for the
same 960 Gamma environments are in
[`../glauber_square_grid_kasteleyn_20260822/`](../glauber_square_grid_kasteleyn_20260822/).
