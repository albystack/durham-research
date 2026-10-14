#!/usr/bin/env julia

"""Environment-blocked estimators for direct square-grid Glauber production traces."""

using Statistics

const DIAGNOSTIC_HEADER =
    "campaign_id,phase,L,environment_id,environment_seed,replica,start," *
    "replica_seed,mean,variance,integrated_autocorrelation_time," *
    "effective_sample_size,changed_rate,final_total_change_rate," *
    "burn_in_attempts,thin_attempts,trace_samples,swap_acceptance_rate," *
    "target_exchange_acceptance_rate,attempted_swaps_by_pair," *
    "accepted_swaps_by_pair"
const TRACE_HEADER =
    "campaign_id,phase,L,environment_id,environment_seed,replica,start," *
    "replica_seed,draw_index,attempt_time,center_height"

function parse_arguments(arguments)
    options = Dict("production-dir" => "", "analysis-dir" => "")
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        startswith(argument, "--") || error("unexpected argument: $argument")
        key = argument[3:end]
        haskey(options, key) || error("unknown option: --$key")
        index < length(arguments) || error("missing value after --$key")
        options[key] = arguments[index + 1]
        index += 2
    end
    isempty(options["production-dir"]) && error("--production-dir is required")
    isempty(options["analysis-dir"]) && error("--analysis-dir is required")
    return (production_dir=abspath(options["production-dir"]),
            analysis_dir=abspath(options["analysis-dir"]))
end

function metadata(path)
    values = Dict{String,String}()
    for line in eachline(path)
        key, value = split(chomp(line), '='; limit=2)
        values[key] = value
    end
    return values
end

function files_with_prefix(root, prefix)
    paths = String[]
    for (directory, _, files) in walkdir(root), file in files
        startswith(file, prefix) && endswith(file, ".csv") &&
            push!(paths, joinpath(directory, file))
    end
    return sort(paths)
end

function parse_integer_vector(value)
    value == "not_applicable" && return Int[]
    return parse.(Int, split(value, ';'))
end

function load_diagnostics(root)
    rows = Dict{Tuple{Int,Int,Int},NamedTuple}()
    for path in files_with_prefix(root, "diagnostic_")
        lines = readlines(path)
        strip(first(lines)) == DIAGNOSTIC_HEADER || error("unexpected diagnostic header: $path")
        for line in lines[2:end]
            fields = split(strip(line), ',')
            length(fields) == 21 || error("malformed diagnostic row: $path")
            L = parse(Int, fields[3])
            environment_id = parse(Int, fields[4])
            replica = parse(Int, fields[6])
            attempted = parse_integer_vector(fields[20])
            accepted = parse_integer_vector(fields[21])
            length(attempted) == length(accepted) || error("mismatched swap counts: $path")
            pair_rates = [accepted[index] / attempted[index] for index in eachindex(attempted)
                          if attempted[index] > 0]
            length(pair_rates) == length(attempted) || error("unattempted beta pair: $path")
            key = (L, environment_id, replica)
            haskey(rows, key) && error("duplicate diagnostic row for $key")
            rows[key] = (
                environment_seed=parse(UInt64, fields[5]), start=Symbol(fields[7]),
                mean=parse(Float64, fields[9]), variance=parse(Float64, fields[10]),
                ess=parse(Float64, fields[12]),
                minimum_pair_swap_acceptance=isempty(pair_rates) ? NaN : minimum(pair_rates),
                target_exchange_acceptance=parse(Float64, fields[19]),
            )
        end
    end
    isempty(rows) && error("no diagnostic files below $root")
    return rows
end

function load_traces(root)
    traces = Dict{Tuple{Int,Int,Int},Vector{Float64}}()
    seeds = Dict{Tuple{Int,Int},UInt64}()
    for path in files_with_prefix(root, "batch_")
        lines = readlines(path)
        strip(first(lines)) == TRACE_HEADER || error("unexpected trace header: $path")
        for line in lines[2:end]
            fields = split(strip(line), ',')
            length(fields) == 11 || error("malformed trace row: $path")
            fields[2] == "production" || error("non-production trace row: $path")
            L = parse(Int, fields[3])
            environment_id = parse(Int, fields[4])
            environment_seed = parse(UInt64, fields[5])
            replica = parse(Int, fields[6])
            seeds[(L, environment_id)] = environment_seed
            push!(get!(traces, (L, environment_id, replica), Float64[]),
                  parse(Float64, fields[11]))
        end
    end
    isempty(traces) && error("no trace files below $root")
    return traces, seeds
end

function environment_rows(diagnostics, traces, seeds)
    keys_by_environment = sort(collect(keys(seeds)))
    rows = NamedTuple[]
    for (L, environment_id) in keys_by_environment
        first_key = (L, environment_id, 1)
        second_key = (L, environment_id, 2)
        haskey(traces, first_key) && haskey(traces, second_key) ||
            error("missing paired trace for L=$L environment=$environment_id")
        haskey(diagnostics, first_key) && haskey(diagnostics, second_key) ||
            error("missing paired diagnostics for L=$L environment=$environment_id")
        first_trace = traces[first_key]
        second_trace = traces[second_key]
        length(first_trace) == length(second_trace) || error(
            "paired trace lengths differ for L=$L environment=$environment_id")
        first = diagnostics[first_key]
        second = diagnostics[second_key]
        first.environment_seed == second.environment_seed == seeds[(L, environment_id)] ||
            error("environment seed mismatch for L=$L environment=$environment_id")
        first.start == :max && second.start == :min || error(
            "expected max/min replicas for L=$L environment=$environment_id")
        standard_error = sqrt(first.variance / first.ess + second.variance / second.ess)
        gap = first.mean - second.mean
        difference = first_trace .- second_trace
        push!(rows, (
            L=L, environment_id=environment_id,
            environment_seed=seeds[(L, environment_id)], samples=length(first_trace),
            replica_1_mean=mean(first_trace), replica_2_mean=mean(second_trace),
            mean_height=(mean(first_trace) + mean(second_trace)) / 2,
            conditional_variance=var(difference; corrected=true) / 2,
            replica_1_variance=var(first_trace; corrected=true),
            replica_2_variance=var(second_trace; corrected=true),
            replica_1_ess=first.ess, replica_2_ess=second.ess,
            start_gap=gap,
            standardized_start_gap=standard_error == 0 ?
                (gap == 0 ? 0.0 : Inf) : abs(gap) / standard_error,
            minimum_pair_swap_acceptance=ifelse(
                isnan(first.minimum_pair_swap_acceptance),
                second.minimum_pair_swap_acceptance,
                isnan(second.minimum_pair_swap_acceptance) ?
                    first.minimum_pair_swap_acceptance :
                    min(first.minimum_pair_swap_acceptance,
                        second.minimum_pair_swap_acceptance)),
            target_exchange_acceptance=min(first.target_exchange_acceptance,
                                            second.target_exchange_acceptance),
            marginal_second_moment=(mean(abs2, first_trace) +
                                    mean(abs2, second_trace)) / 2,
        ))
    end
    return rows
end

function write_environment_rows(path, rows)
    open(path, "w") do io
        println(io, "L,environment_id,environment_seed,samples,replica_1_mean," *
                    "replica_2_mean,mean_height,conditional_variance," *
                    "replica_1_variance,replica_2_variance,replica_1_ess," *
                    "replica_2_ess,start_gap,standardized_start_gap," *
                    "minimum_pair_swap_acceptance,target_exchange_acceptance")
        for row in rows
            println(io, join((row.L, row.environment_id, row.environment_seed, row.samples,
                              row.replica_1_mean, row.replica_2_mean, row.mean_height,
                              row.conditional_variance, row.replica_1_variance,
                              row.replica_2_variance, row.replica_1_ess,
                              row.replica_2_ess, row.start_gap,
                              row.standardized_start_gap,
                              row.minimum_pair_swap_acceptance,
                              row.target_exchange_acceptance), ','))
        end
    end
end

function write_size_summary(path, rows)
    summaries = NamedTuple[]
    for L in sort(unique(row.L for row in rows))
        selected = filter(row -> row.L == L, rows)
        first_means = getproperty.(selected, :replica_1_mean)
        second_means = getproperty.(selected, :replica_2_mean)
        conditional = mean(getproperty.(selected, :conditional_variance))
        disorder = length(selected) > 1 ? cov(first_means, second_means; corrected=true) : NaN
        pooled_mean = mean(getproperty.(selected, :mean_height))
        direct_variance = mean(getproperty.(selected, :marginal_second_moment)) - pooled_mean^2
        push!(summaries, (
            L=L, environments=length(selected), pooled_mean=pooled_mean,
            conditional_component=conditional, disorder_component=disorder,
            total_component=conditional + disorder,
            direct_annealed_variance=direct_variance,
            maximum_absolute_start_gap=maximum(abs, getproperty.(selected, :start_gap)),
            maximum_standardized_start_gap=maximum(
                getproperty.(selected, :standardized_start_gap)),
            minimum_chain_ess=minimum(vcat(getproperty.(selected, :replica_1_ess),
                                           getproperty.(selected, :replica_2_ess))),
            minimum_pair_swap_acceptance=begin
                rates = filter(rate -> !isnan(rate),
                               getproperty.(selected, :minimum_pair_swap_acceptance))
                isempty(rates) ? NaN : minimum(rates)
            end,
            minimum_target_exchange_acceptance=minimum(
                getproperty.(selected, :target_exchange_acceptance)),
        ))
    end
    open(path, "w") do io
        println(io, "L,environments,pooled_mean,conditional_component,disorder_component," *
                    "total_component,direct_annealed_variance,maximum_absolute_start_gap," *
                    "maximum_standardized_start_gap,minimum_chain_ess," *
                    "minimum_pair_swap_acceptance,minimum_target_exchange_acceptance")
        for row in summaries
            println(io, join(values(row), ','))
        end
    end
    return summaries
end

function main(arguments)
    parsed = parse_arguments(arguments)
    metadata_path = joinpath(parsed.production_dir, "campaign_metadata.txt")
    isfile(metadata_path) || error("missing production metadata: $metadata_path")
    campaign_metadata = metadata(metadata_path)
    campaign_metadata["phase"] == "production" || error("input is not production output")
    diagnostics = load_diagnostics(parsed.production_dir)
    traces, seeds = load_traces(parsed.production_dir)
    rows = environment_rows(diagnostics, traces, seeds)
    mkpath(parsed.analysis_dir)
    write_environment_rows(joinpath(parsed.analysis_dir, "production_environment_blocks.csv"), rows)
    summaries = write_size_summary(joinpath(parsed.analysis_dir, "production_size_summary.csv"), rows)
    open(joinpath(parsed.analysis_dir, "ESTIMATORS.txt"), "w") do io
        println(io, "Independent resampling unit: frozen edge-weight environment.")
        println(io, "conditional_component = mean over environments of Var(H1-H2 | environment)/2.")
        println(io, "disorder_component = covariance across environments of the two independent chain means.")
        println(io, "total_component = conditional_component + disorder_component.")
        println(io, "Raw MCMC draws are retained for diagnostics but are not independent environment blocks.")
    end
    println("Wrote $(length(rows)) environment blocks and $(length(summaries)) size summaries")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
