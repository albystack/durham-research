# Command-line workflows

Scripts are grouped by pipeline stage through their filename prefixes:

- `run_*`: generate deterministic, restart-safe batches;
- `merge_*`: validate and combine compatible batches;
- `analyze_*`: compute environment-blocked summaries and model comparisons;
- `plot_*`: render deterministic figures from derived tables;
- `summarize_*`: consolidate several completed analyses.

## Main pipelines

| Model | Run | Analyze |
|---|---|---|
| Aztec central height | `run_height_campaign.jl` | `analyze_height_campaign.jl` |
| Aztec paired central height | `run_double_dimer_campaign.jl` | `analyze_double_dimer_campaign.jl` |
| Aztec paired spatial increments | `run_spatial_campaign.jl` | `analyze_spatial_campaign.jl` |
| Temperley square grid | `run_square_grid_campaign.jl` | `analyze_spatial_campaign.jl` |
| Direct square-grid Glauber | `run_glauber_square_grid_campaign.jl` | `analyze_glauber_square_grid_production.jl`, then `analyze_glauber_square_grid_scaling.jl` |
| Direct square-grid Kasteleyn | `run_glauber_kasteleyn_campaign.jl` | `analyze_glauber_kasteleyn_campaign.jl` |

Run a script with `--help` before using a new workflow. Campaign outputs should
go below `aztec/output/` locally or an explicit project-data directory on HPC.
Never merge batches from different base seeds or model labels.

The Glauber and Kasteleyn runners deliberately share the direct-model
environment seed namespace. This permits exact environment-by-environment
comparisons when their distribution, parameter, base seed, size, and
environment identifier agree.
