#!/usr/bin/env julia

"""Summarize direct-dimer Glauber pilot traces without declaring convergence."""

using AztecDiamond
using AztecDiamond.GlauberSquareGrid
using Random
using Statistics

const DIAGNOSTIC_HEADER =
    "campaign_id,phase,L,environment_id,environment_seed,replica,start," *
    "replica_seed,mean,variance,integrated_autocorrelation_time," *
    "effective_sample_size,changed_rate,final_total_change_rate," *
    "burn_in_attempts,thin_attempts,trace_samples"

const TRACE_HEADER =
    "campaign_id,phase,L,environment_id,environment_seed,replica,start," *
    "replica_seed,draw_index,attempt_time,center_height"

function parse_arguments(arguments)
    options = Dict("pilot-dir" => "", "analysis-dir" => "")
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
    isempty(options["pilot-dir"]) && error("--pilot-dir is required")
    isempty(options["analysis-dir"]) && error("--analysis-dir is required")
    return (pilot_dir=abspath(options["pilot-dir"]), analysis_dir=abspath(options["analysis-dir"]))
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
        startswith(file, prefix) && endswith(file, ".csv") && push!(paths, joinpath(directory, file))
    end
    return sort(paths)
end

function load_diagnostics(root)
    rows = NamedTuple[]
    for path in files_with_prefix(root, "diagnostic_")
        lines = readlines(path)
        strip(first(lines)) == DIAGNOSTIC_HEADER || error("unexpected diagnostic header: $path")
        for line in lines[2:end]
            fields = split(strip(line), ',')
            length(fields) == 17 || error("malformed diagnostic row: $path")
            push!(rows, (
                L=parse(Int, fields[3]), environment_id=parse(Int, fields[4]),
                environment_seed=parse(UInt64, fields[5]), replica=parse(Int, fields[6]),
                start=Symbol(fields[7]), mean=parse(Float64, fields[9]), variance=parse(Float64, fields[10]),
                ess=parse(Float64, fields[12]), changed_rate=parse(Float64, fields[13]),
            ))
        end
    end
    isempty(rows) && error("no diagnostic files below $root")
    return rows
end

function load_traces(root)
    rows = NamedTuple[]
    for path in files_with_prefix(root, "batch_")
        lines = readlines(path)
        strip(first(lines)) == TRACE_HEADER || error("unexpected trace header: $path")
        for line in lines[2:end]
            fields = split(strip(line), ',')
            length(fields) == 11 || error("malformed trace row: $path")
            push!(rows, (
                L=parse(Int, fields[3]), environment_id=parse(Int, fields[4]),
                environment_seed=parse(UInt64, fields[5]), replica=parse(Int, fields[6]),
                height=parse(Int, fields[11]),
            ))
        end
    end
    return rows
end

function write_summary(path, diagnostics)
    by_environment = Dict{Tuple{Int,Int},Vector{NamedTuple}}()
    for row in diagnostics
        push!(get!(by_environment, (row.L, row.environment_id), NamedTuple[]), row)
    end
    summaries = NamedTuple[]
    for L in sort(unique(row.L for row in diagnostics))
        pairs = [rows for ((order, _), rows) in by_environment if order == L]
        all(length(rows) == 2 && Set(row.start for row in rows) == Set((:max, :min)) for rows in pairs) ||
            error("L=$L does not have one max/min chain pair per environment")
        gaps = Float64[]
        standardized_gaps = Float64[]
        ess_values = Float64[]
        changed_rates = Float64[]
        for pair in pairs
            maximum = only(filter(row -> row.start === :max, pair))
            minimum = only(filter(row -> row.start === :min, pair))
            gap = maximum.mean - minimum.mean
            push!(gaps, gap)
            standard_error = sqrt(maximum.variance / maximum.ess + minimum.variance / minimum.ess)
            push!(standardized_gaps, standard_error == 0 ?
                 (gap == 0 ? 0.0 : Inf) : abs(gap) / standard_error)
            append!(ess_values, (maximum.ess, minimum.ess))
            append!(changed_rates, (maximum.changed_rate, minimum.changed_rate))
        end
        push!(summaries, (
            L=L, environments=length(pairs), mean_start_gap=mean(gaps),
            maximum_absolute_start_gap=maximum(abs, gaps),
            maximum_standardized_start_gap=maximum(standardized_gaps),
            median_chain_ess=median(ess_values), minimum_chain_ess=minimum(ess_values),
            mean_changed_rate=mean(changed_rates),
        ))
    end
    open(path, "w") do io
        println(io, "L,environments,mean_start_gap,maximum_absolute_start_gap," *
                    "maximum_standardized_start_gap,median_chain_ess,minimum_chain_ess,mean_changed_rate")
        for row in summaries
            println(io, join((row.L, row.environments, row.mean_start_gap,
                              row.maximum_absolute_start_gap, row.maximum_standardized_start_gap,
                              row.median_chain_ess, row.minimum_chain_ess, row.mean_changed_rate), ','))
        end
    end
    return summaries
end

function write_exact_L2_check(path, pilot_metadata, traces)
    distribution = Symbol(pilot_metadata["distribution"])
    parameter = parse(Float64, pilot_metadata["parameter"])
    l2 = filter(row -> row.L == 2, traces)
    isempty(l2) && return nothing
    grouped = Dict{Tuple{Int,UInt64},Vector{Int}}()
    for row in l2
        push!(get!(grouped, (row.environment_id, row.environment_seed), Int[]), row.height)
    end
    open(path, "w") do io
        println(io, "L,environment_id,environment_seed,trace_samples,total_variation_distance")
        for ((environment_id, environment_seed), heights) in sort(collect(grouped); by=first)
            weights = distribution === :constant ? constant_edge_weights(2) :
                random_edge_weights(Xoshiro(environment_seed), 2;
                                    distribution=distribution, parameter=parameter)
            exact = exact_center_distribution(2, weights).masses
            observed = Dict(value => count(==(value), heights) / length(heights) for value in unique(heights))
            support = union(keys(exact), keys(observed))
            distance = 0.5 * sum(abs(get(exact, value, 0.0) - get(observed, value, 0.0))
                                 for value in support)
            println(io, join((2, environment_id, environment_seed, length(heights), distance), ','))
        end
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    metadata_path = joinpath(parsed.pilot_dir, "campaign_metadata.txt")
    isfile(metadata_path) || error("missing pilot metadata: $metadata_path")
    pilot_metadata = metadata(metadata_path)
    pilot_metadata["phase"] == "pilot" || error("input is not a Glauber pilot")
    diagnostics = load_diagnostics(parsed.pilot_dir)
    traces = load_traces(parsed.pilot_dir)
    mkpath(parsed.analysis_dir)
    summaries = write_summary(joinpath(parsed.analysis_dir, "pilot_mixing_summary.csv"), diagnostics)
    write_exact_L2_check(joinpath(parsed.analysis_dir, "pilot_L2_exact_check.csv"), pilot_metadata, traces)
    open(joinpath(parsed.analysis_dir, "REVIEW_REQUIRED.txt"), "w") do io
        println(io, "This pilot summary is diagnostic evidence, not an automatic equilibration certificate.")
        println(io, "Review the start gaps, ESS, L=2 exact check, Slurm resource use, and raw traces before any production submission.")
    end
    println("Wrote $(length(summaries)) L-level pilot summaries to $(parsed.analysis_dir)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
