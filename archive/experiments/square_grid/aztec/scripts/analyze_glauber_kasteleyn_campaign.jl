#!/usr/bin/env julia

module GlauberKasteleynAnalysis

using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "analyze_glauber_square_grid_scaling.jl"))
const Scaling = GlauberProductionScaling

const INPUT_HEADER =
    "campaign_id,L,environment_id,environment_seed,conditional_mean," *
    "conditional_variance,conditional_second_moment,log_partition," *
    "partition_phase_real,partition_phase_imaginary,matrix_order,precision_bits," *
    "crossed_edges," *
    "relative_solve_residual,maximum_probability_imaginary_residual," *
    "maximum_covariance_imaginary_residual,elapsed_seconds"

function print_help()
    println("""
    Analyze deterministic Kasteleyn environment moments and optionally compare
    them environment-by-environment with retained Glauber MCMC blocks.

    Options:
      --input-root PATH       completed Kasteleyn campaign root
      --output-dir PATH       analysis output directory
      --mcmc-blocks PATHS     optional comma-separated MCMC block CSVs
      --bootstrap-reps INT    default 5000
      --bootstrap-seed UINT   default 20260823
      --cutoffs INTS          default 2,4,6,8
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict(
        "input-root" => "",
        "output-dir" => "",
        "mcmc-blocks" => "",
        "bootstrap-reps" => "5000",
        "bootstrap-seed" => "20260823",
        "cutoffs" => "2,4,6,8",
    )
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
    isempty(options["input-root"]) && error("--input-root is required")
    isempty(options["output-dir"]) && error("--output-dir is required")
    repetitions = parse(Int, options["bootstrap-reps"])
    repetitions > 0 || error("--bootstrap-reps must be positive")
    cutoffs = sort(unique(parse.(Int, split(options["cutoffs"], ','))))
    !isempty(cutoffs) && all(>(0), cutoffs) || error("cutoffs must be positive")
    mcmc_paths = [abspath(strip(path)) for path in split(options["mcmc-blocks"], ',')
                  if !isempty(strip(path))]
    return (
        input_root=abspath(options["input-root"]),
        output_dir=abspath(options["output-dir"]),
        mcmc_paths=mcmc_paths,
        bootstrap_reps=repetitions,
        bootstrap_seed=parse(UInt64, options["bootstrap-seed"]),
        cutoffs=cutoffs,
    )
end

function batch_files(root)
    isdir(root) || error("Kasteleyn input root does not exist: $root")
    paths = String[]
    for (directory, _, files) in walkdir(root), file in files
        occursin(r"^batch_\d+\.csv$", file) || continue
        push!(paths, joinpath(directory, file))
    end
    sort!(paths)
    isempty(paths) && error("no Kasteleyn batch files found under $root")
    return paths
end

function load_rows(root)
    rows = NamedTuple[]
    seen = Set{Tuple{Int,Int,UInt64}}()
    for path in batch_files(root)
        lines = readlines(path)
        isempty(lines) && error("empty Kasteleyn batch: $path")
        strip(first(lines)) == INPUT_HEADER || error("unexpected header: $path")
        for (offset, line) in enumerate(lines[2:end])
            isempty(strip(line)) && continue
            fields = split(strip(line), ',')
            length(fields) == 17 || error("malformed row $(offset + 1): $path")
            row = (
                campaign_id=fields[1],
                L=parse(Int, fields[2]),
                environment_id=parse(Int, fields[3]),
                environment_seed=parse(UInt64, fields[4]),
                conditional_mean=parse(Float64, fields[5]),
                conditional_variance=parse(Float64, fields[6]),
                conditional_second_moment=parse(Float64, fields[7]),
                log_partition=parse(Float64, fields[8]),
                partition_phase_real=parse(Float64, fields[9]),
                partition_phase_imaginary=parse(Float64, fields[10]),
                matrix_order=parse(Int, fields[11]),
                precision_bits=parse(Int, fields[12]),
                crossed_edges=parse(Int, fields[13]),
                relative_solve_residual=parse(Float64, fields[14]),
                maximum_probability_imaginary_residual=parse(Float64, fields[15]),
                maximum_covariance_imaginary_residual=parse(Float64, fields[16]),
                elapsed_seconds=parse(Float64, fields[17]),
            )
            key = (row.L, row.environment_id, row.environment_seed)
            key in seen && error("duplicate Kasteleyn environment: $key")
            push!(seen, key)
            push!(rows, row)
        end
    end
    sort!(rows; by=row -> (row.L, row.environment_id))
    return rows
end

function analysis_row(row)
    return (
        replica_1_mean=row.conditional_mean,
        replica_2_mean=row.conditional_mean,
        conditional_variance=row.conditional_variance,
    )
end

function grouped_analysis_rows(rows)
    grouped = Dict{Int,Vector{NamedTuple}}()
    for row in rows
        push!(get!(grouped, row.L, NamedTuple[]), analysis_row(row))
    end
    for (L, selected) in grouped
        length(selected) >= 2 || error("Kasteleyn scaling needs at least two environments at L=$L")
    end
    return grouped
end

function write_merged(path, rows)
    open(path, "w") do io
        println(io, INPUT_HEADER)
        for row in rows
            println(io, join(values(row), ','))
        end
    end
end

function comparison_rows(exact_rows, mcmc_rows)
    exact = Dict((row.L, row.environment_seed) => row for row in exact_rows)
    comparisons = NamedTuple[]
    for row in mcmc_rows
        key = (row.L, row.environment_seed)
        haskey(exact, key) || continue
        reference = exact[key]
        row.environment_id == reference.environment_id || error(
            "environment ID mismatch for L=$(row.L), seed=$(row.environment_seed)")
        mcmc_mean = (row.replica_1_mean + row.replica_2_mean) / 2
        standard_error = sqrt(row.replica_1_variance / row.replica_1_ess +
                              row.replica_2_variance / row.replica_2_ess) / 2
        mean_error = mcmc_mean - reference.conditional_mean
        zero_standard_error_inconsistent = iszero(standard_error) && !iszero(mean_error)
        push!(comparisons, (
            L=row.L,
            environment_id=row.environment_id,
            environment_seed=row.environment_seed,
            exact_mean=reference.conditional_mean,
            mcmc_mean=mcmc_mean,
            mcmc_minus_exact_mean=mean_error,
            estimated_mcmc_mean_standard_error=standard_error,
            standardized_mean_error=standard_error > 0 ? mean_error / standard_error :
                                    (iszero(mean_error) ? 0.0 : NaN),
            zero_standard_error_inconsistent=zero_standard_error_inconsistent,
            exact_conditional_variance=reference.conditional_variance,
            mcmc_conditional_variance=row.conditional_variance,
            mcmc_minus_exact_conditional_variance=
                row.conditional_variance - reference.conditional_variance,
            extremal_start_gap=row.start_gap,
            standardized_extremal_start_gap=row.standardized_start_gap,
            minimum_chain_ess=min(row.replica_1_ess, row.replica_2_ess),
            minimum_pair_swap_acceptance=row.minimum_pair_swap_acceptance,
            target_exchange_acceptance=row.target_exchange_acceptance,
        ))
    end
    sort!(comparisons; by=row -> (row.L, row.environment_id))
    return comparisons
end

function write_comparisons(path, rows)
    isempty(rows) && return
    open(path, "w") do io
        println(io, join(keys(first(rows)), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
end

function comparison_summary(rows)
    grouped = Dict{Int,Vector{NamedTuple}}()
    for row in rows
        push!(get!(grouped, row.L, NamedTuple[]), row)
    end
    summaries = NamedTuple[]
    for (L, selected) in sort(collect(grouped); by=first)
        mean_errors = getproperty.(selected, :mcmc_minus_exact_mean)
        standardized = getproperty.(selected, :standardized_mean_error)
        finite_standardized = filter(isfinite, standardized)
        zero_standard_error_inconsistent = count(
            getproperty.(selected, :zero_standard_error_inconsistent))
        variance_errors = getproperty.(selected, :mcmc_minus_exact_conditional_variance)
        push!(summaries, (
            L=L,
            environments=length(selected),
            mean_mcmc_minus_exact=mean(mean_errors),
            mean_absolute_mean_error=mean(abs, mean_errors),
            root_mean_square_mean_error=sqrt(mean(abs2, mean_errors)),
            maximum_absolute_finite_standardized_mean_error=
                isempty(finite_standardized) ? NaN : maximum(abs, finite_standardized),
            count_absolute_standardized_mean_error_above_2=
                count(value -> abs(value) > 2, finite_standardized) +
                zero_standard_error_inconsistent,
            count_absolute_standardized_mean_error_above_4=
                count(value -> abs(value) > 4, finite_standardized) +
                zero_standard_error_inconsistent,
            count_zero_standard_error_inconsistent=zero_standard_error_inconsistent,
            mean_conditional_variance_error=mean(variance_errors),
            root_mean_square_conditional_variance_error=
                sqrt(mean(abs2, variance_errors)),
        ))
    end
    return summaries
end

function write_named_tuples(path, rows)
    isempty(rows) && return
    open(path, "w") do io
        println(io, join(keys(first(rows)), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
end

function write_numerical_summary(path, rows)
    open(path, "w") do io
        println(io, "environments=$(length(rows))")
        println(io, "maximum_relative_solve_residual=$(maximum(
            getproperty.(rows, :relative_solve_residual)))")
        println(io, "maximum_probability_imaginary_residual=$(maximum(
            getproperty.(rows, :maximum_probability_imaginary_residual)))")
        println(io, "maximum_covariance_imaginary_residual=$(maximum(
            getproperty.(rows, :maximum_covariance_imaginary_residual)))")
        println(io, "high_precision_fallback_environments=$(count(
            row -> row.precision_bits > 53, rows))")
        println(io, "maximum_elapsed_seconds=$(maximum(getproperty.(rows, :elapsed_seconds)))")
    end
end

function main(arguments=ARGS)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return print_help()
    rows = load_rows(parsed.input_root)
    grouped = grouped_analysis_rows(rows)
    rng = Xoshiro(parsed.bootstrap_seed)
    analysis = (
        model="kasteleyn_gamma",
        filter="all",
        grouped=grouped,
        draws=Scaling.bootstrap_size_draws(rng, grouped, parsed.bootstrap_reps),
    )

    mkpath(parsed.output_dir)
    write_merged(joinpath(parsed.output_dir, "kasteleyn_environment_moments.csv"), rows)
    Scaling.write_size_components(
        joinpath(parsed.output_dir, "kasteleyn_component_summary.csv"), [analysis])
    scaling_rows = Scaling.model_comparisons(
        [analysis], parsed.cutoffs, parsed.bootstrap_reps)
    Scaling.write_model_comparisons(
        joinpath(parsed.output_dir, "kasteleyn_scaling_comparison.csv"), scaling_rows)
    write_numerical_summary(joinpath(parsed.output_dir, "NUMERICAL_DIAGNOSTICS.txt"), rows)

    if !isempty(parsed.mcmc_paths)
        mcmc = Scaling.load_environment_blocks(parsed.mcmc_paths)
        comparisons = comparison_rows(rows, mcmc)
        isempty(comparisons) && error("no shared Kasteleyn/MCMC environments")
        write_comparisons(joinpath(parsed.output_dir, "kasteleyn_mcmc_comparison.csv"),
                          comparisons)
        write_named_tuples(joinpath(parsed.output_dir, "kasteleyn_mcmc_size_summary.csv"),
                           comparison_summary(comparisons))
    end

    open(joinpath(parsed.output_dir, "ANALYSIS_METHOD.txt"), "w") do io
        println(io, "Independent block: one frozen edge-weight environment within one L.")
        println(io, "Conditional moments: finite dense Kasteleyn factorization; no MCMC draws.")
        println(io, "Disorder component: sample variance of exact conditional means.")
        println(io, "Total component: mean exact conditional variance plus disorder component.")
        println(io, "Bootstrap: whole environments within size; $(parsed.bootstrap_reps) repetitions.")
        println(io, "Bootstrap seed: $(parsed.bootstrap_seed).")
        println(io, "Nested models: a+b*log(L) versus a+b*log(L)+c*log(L)^2.")
        println(io, "Positive delta BIC favors the quadratic extension.")
        println(io, "MCMC comparison, when supplied, joins identical L and environment seeds.")
    end
    println("Kasteleyn campaign analysis complete")
    println("  environments: $(length(rows))")
    println("  sizes:        $(join(sort(collect(keys(grouped))), ','))")
    println("  output:       $(parsed.output_dir)")
end

end # module

abspath(PROGRAM_FILE) == abspath(@__FILE__) && GlauberKasteleynAnalysis.main()
