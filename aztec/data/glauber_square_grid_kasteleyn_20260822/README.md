# Determinantal square-grid environment moments

`kasteleyn_environment_moments.csv` contains the finite-volume conditional
central-height mean and variance for all 960 frozen Gamma environments from the
completed direct Glauber campaign. The environment seeds and size allocation
are identical to the original production:

| L | Environments |
|---:|---:|
| 2 | 64 |
| 4 | 192 |
| 6 | 192 |
| 8 | 192 |
| 10 | 128 |
| 12 | 96 |
| 16 | 64 |
| 20 | 32 |

The moments come from a finite bipartite Kasteleyn factorization, not MCMC.
Each row records solve and imaginary-residual diagnostics. One ill-conditioned
`L=12` environment automatically used a 256-bit factorization; all others used
`Float64`. `CHECKSUMS.sha256` covers the retained environment table.

Regenerate the raw restart-safe batches and merged analysis with:

```bash
julia --project=aztec aztec/scripts/run_glauber_kasteleyn_campaign.jl \
  --config aztec/configs/glauber_square_grid_kasteleyn_replay.csv \
  --output-dir aztec/output/glauber_square_grid_kasteleyn_replay_20260822 \
  --distribution gamma --parameter 0.5 --base-seed 2026082201

julia --project=aztec aztec/scripts/analyze_glauber_kasteleyn_campaign.jl \
  --input-root aztec/output/glauber_square_grid_kasteleyn_replay_20260822 \
  --output-dir aztec/output/glauber_square_grid_kasteleyn_analysis_20260822 \
  --mcmc-blocks aztec/data/glauber_square_grid_20260822/gamma_core_environment_blocks.csv,aztec/data/glauber_square_grid_20260822/gamma_frontier_environment_blocks.csv \
  --bootstrap-reps 5000 --bootstrap-seed 20260823 --cutoffs 2,4,6,8
```

The retained analysis is in
[`../../results/glauber_square_grid_kasteleyn_20260822/`](../../results/glauber_square_grid_kasteleyn_20260822/).
