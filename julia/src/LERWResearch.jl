module LERWResearch

using CSV
using Dates
using Distributions
using JSON3
using Printf
using Random
using SHA
using StableRNGs
using Statistics

include("config.jl")
include("simulation.jl")
include("batch.jl")
include("analysis.jl")

export BatchConfig,
       SiteIIDEnvironment,
       TemporalIIDEnvironment,
       analyze_results,
       canonical_json,
       completed_result,
       distribution_spec,
       generate_config,
       generate_double_dimer_reproduction_config,
       generate_double_dimer_pilot_config,
       generate_temporal_iid_benchmark_config,
       generate_temporal_iid_gamma_length_config,
       generate_temporal_iid_pilot_config,
       generate_temporal_iid_smoke_config,
       generate_strict_annealed_reproduction_config,
       is_boundary,
       load_config_row,
       loop_erased_walk,
       loop_erased_walk_diagnostics,
       loop_erase,
       main_analyze,
       main_generate_config,
       main_run_batch,
       main_run_campaign,
       result_path,
       run_batch,
       run_campaign,
       stable_seed,
       temporal_weights!,
       winding

end
