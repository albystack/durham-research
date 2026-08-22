#!/usr/bin/env julia

module GlauberProductionScaling

using LinearAlgebra
using Printf
using Random
using Statistics

const ENVIRONMENT_HEADER =
    "L,environment_id,environment_seed,samples,replica_1_mean," *
    "replica_2_mean,mean_height,conditional_variance," *
    "replica_1_variance,replica_2_variance,replica_1_ess," *
    "replica_2_ess,start_gap,standardized_start_gap," *
    "minimum_pair_swap_acceptance,target_exchange_acceptance"

function print_help()
    println("""
    Fit size scaling for environment-blocked square-grid Glauber production data.

    Primary comparison:
      H0: a + b log(L)
      H1: a + b log(L) + c (log(L))^2

    Each bootstrap resamples whole frozen environments within a size.  The two
    chain means and conditional-variance estimate always remain in one block.

    Options:
      --gamma-blocks PATHS       comma-separated environment-block CSVs
      --control-blocks PATHS     comma-separated environment-block CSVs
      --output-dir PATH          required
      --bootstrap-reps INT       default 5000
      --bootstrap-seed UINT      default 20260822
      --cutoffs INTS             comma-separated lower-size cutoffs (default 2,4,6,8)
      --diagnostic-threshold X   |standardized start gap| sensitivity cutoff (default 4)
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict(
        "gamma-blocks" => "",
        "control-blocks" => "",
        "output-dir" => "",
        "bootstrap-reps" => "5000",
        "bootstrap-seed" => "20260822",
        "cutoffs" => "2,4,6,8",
        "diagnostic-threshold" => "4.0",
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
    isempty(options["gamma-blocks"]) && error("--gamma-blocks is required")
    isempty(options["control-blocks"]) && error("--control-blocks is required")
    isempty(options["output-dir"]) && error("--output-dir is required")
    paths(value) = [abspath(strip(path)) for path in split(value, ',') if !isempty(strip(path))]
    cutoffs = sort(unique(parse.(Int, split(options["cutoffs"], ','))))
    repetitions = parse(Int, options["bootstrap-reps"])
    threshold = parse(Float64, options["diagnostic-threshold"])
    repetitions > 0 || error("--bootstrap-reps must be positive")
    !isempty(cutoffs) && all(>(0), cutoffs) || error("cutoffs must be positive")
    threshold > 0 || error("--diagnostic-threshold must be positive")
    return (
        gamma_paths=paths(options["gamma-blocks"]),
        control_paths=paths(options["control-blocks"]),
        output_dir=abspath(options["output-dir"]),
        bootstrap_reps=repetitions,
        bootstrap_seed=parse(UInt64, options["bootstrap-seed"]),
        cutoffs=cutoffs,
        diagnostic_threshold=threshold,
    )
end

function load_environment_blocks(paths)
    rows = NamedTuple[]
    seen = Set{Tuple{Int,UInt64}}()
    for path in paths
        isfile(path) || error("environment-block CSV does not exist: $path")
        lines = readlines(path)
        isempty(lines) && error("empty environment-block CSV: $path")
        strip(first(lines)) == ENVIRONMENT_HEADER || error("unexpected header: $path")
        for (offset, line) in enumerate(lines[2:end])
            isempty(strip(line)) && continue
            fields = split(strip(line), ',')
            length(fields) == 16 || error("malformed row $(offset + 1): $path")
            row = (
                L=parse(Int, fields[1]),
                environment_id=parse(Int, fields[2]),
                environment_seed=parse(UInt64, fields[3]),
                samples=parse(Int, fields[4]),
                replica_1_mean=parse(Float64, fields[5]),
                replica_2_mean=parse(Float64, fields[6]),
                mean_height=parse(Float64, fields[7]),
                conditional_variance=parse(Float64, fields[8]),
                replica_1_ess=parse(Float64, fields[11]),
                replica_2_ess=parse(Float64, fields[12]),
                start_gap=parse(Float64, fields[13]),
                standardized_start_gap=parse(Float64, fields[14]),
                minimum_pair_swap_acceptance=parse(Float64, fields[15]),
                target_exchange_acceptance=parse(Float64, fields[16]),
            )
            key = (row.L, row.environment_seed)
            key in seen && error("duplicate environment block L=$(row.L), seed=$(row.environment_seed)")
            push!(seen, key)
            push!(rows, row)
        end
    end
    isempty(rows) && error("no environment blocks loaded")
    return rows
end

function grouped_rows(rows; diagnostic_threshold=nothing)
    grouped = Dict{Int,Vector{NamedTuple}}()
    for row in rows
        if !isnothing(diagnostic_threshold) &&
                row.standardized_start_gap > diagnostic_threshold
            continue
        end
        push!(get!(grouped, row.L, NamedTuple[]), row)
    end
    for (L, selected) in grouped
        length(selected) >= 2 || error("fewer than two environment blocks remain at L=$L")
    end
    return grouped
end

function component_statistics(rows)
    length(rows) >= 2 || error("component estimates need at least two environments")
    first_means = getproperty.(rows, :replica_1_mean)
    second_means = getproperty.(rows, :replica_2_mean)
    conditional = mean(getproperty.(rows, :conditional_variance))
    disorder = cov(first_means, second_means; corrected=true)
    return (conditional=conditional, disorder=disorder, total=conditional + disorder)
end

function bootstrap_size_draws(rng, grouped, repetitions)
    draws = Dict{Int,NamedTuple}()
    for L in sort(collect(keys(grouped)))
        rows = grouped[L]
        conditional = Vector{Float64}(undef, repetitions)
        disorder = similar(conditional)
        total = similar(conditional)
        for repetition in 1:repetitions
            resampled = [rows[rand(rng, eachindex(rows))] for _ in eachindex(rows)]
            estimate = component_statistics(resampled)
            conditional[repetition] = estimate.conditional
            disorder[repetition] = estimate.disorder
            total[repetition] = estimate.total
        end
        draws[L] = (conditional=conditional, disorder=disorder, total=total)
    end
    return draws
end

function percentile(values, probability)
    sorted = sort(values)
    position = 1 + (length(sorted) - 1) * probability
    lower, upper = floor(Int, position), ceil(Int, position)
    lower == upper && return sorted[lower]
    fraction = position - lower
    return (1 - fraction) * sorted[lower] + fraction * sorted[upper]
end

function regression(x, y; quadratic=false, weights=nothing)
    design = quadratic ? hcat(ones(length(x)), x, x .^ 2) : hcat(ones(length(x)), x)
    if isnothing(weights)
        coefficients = design \ y
        rss = sum(abs2, y - design * coefficients)
    else
        length(weights) == length(y) || throw(DimensionMismatch("weights and response differ"))
        all(>(0), weights) || error("regression weights must be positive")
        root_weights = sqrt.(weights)
        coefficients = (design .* root_weights) \ (y .* root_weights)
        rss = sum(weights .* abs2.(y - design * coefficients))
    end
    count = length(y)
    parameters = size(design, 2)
    bic = count * log(max(rss / count, eps(Float64))) + parameters * log(count)
    return (coefficients=coefficients, rss=rss, bic=bic)
end

function predict(fit, x; quadratic=false)
    quadratic && return fit.coefficients[1] + fit.coefficients[2] * x +
                        fit.coefficients[3] * x^2
    return fit.coefficients[1] + fit.coefficients[2] * x
end

function leave_one_out_rmse(x, y; quadratic=false, weights=nothing)
    errors = Float64[]
    for held_out in eachindex(x)
        training = [index for index in eachindex(x) if index != held_out]
        training_weights = isnothing(weights) ? nothing : weights[training]
        fit = regression(x[training], y[training]; quadratic=quadratic,
                         weights=training_weights)
        push!(errors, y[held_out] - predict(fit, x[held_out]; quadratic=quadratic))
    end
    return sqrt(mean(abs2, errors))
end

function fit_comparison(x, y; weights=nothing)
    length(x) >= 4 || error("model comparison requires at least four sizes")
    log_fit = regression(x, y; weights=weights)
    quadratic_fit = regression(x, y; quadratic=true, weights=weights)
    return (
        log=log_fit,
        quadratic=quadratic_fit,
        delta_bic=log_fit.bic - quadratic_fit.bic,
        log_loocv_rmse=leave_one_out_rmse(x, y; weights=weights),
        quadratic_loocv_rmse=leave_one_out_rmse(
            x, y; quadratic=true, weights=weights),
    )
end

function diagnostic_summary(rows)
    gaps = getproperty.(rows, :start_gap)
    standardized = getproperty.(rows, :standardized_start_gap)
    gap_t = length(rows) > 1 && std(gaps; corrected=true) > 0 ?
        mean(gaps) / (std(gaps; corrected=true) / sqrt(length(rows))) : 0.0
    rates = filter(rate -> !isnan(rate),
                   getproperty.(rows, :minimum_pair_swap_acceptance))
    targets = filter(rate -> !isnan(rate),
                     getproperty.(rows, :target_exchange_acceptance))
    return (
        environments=length(rows),
        mean_signed_start_gap=mean(gaps),
        signed_start_gap_t=gap_t,
        maximum_absolute_start_gap=maximum(abs, gaps),
        count_standardized_above_2=count(>(2), standardized),
        count_standardized_above_3=count(>(3), standardized),
        count_standardized_above_4=count(>(4), standardized),
        maximum_standardized_start_gap=maximum(standardized),
        minimum_chain_ess=minimum(vcat(getproperty.(rows, :replica_1_ess),
                                       getproperty.(rows, :replica_2_ess))),
        minimum_pair_swap_acceptance=isempty(rates) ? NaN : minimum(rates),
        minimum_target_exchange_acceptance=isempty(targets) ? NaN : minimum(targets),
    )
end

function write_diagnostics(path, datasets)
    open(path, "w") do io
        println(io, "model,L,environments,mean_signed_start_gap,signed_start_gap_t," *
                    "maximum_absolute_start_gap,count_standardized_above_2," *
                    "count_standardized_above_3,count_standardized_above_4," *
                    "maximum_standardized_start_gap,minimum_chain_ess," *
                    "minimum_pair_swap_acceptance,minimum_target_exchange_acceptance")
        for model in sort(collect(keys(datasets))),
                (L, rows) in sort(collect(grouped_rows(datasets[model])); by=first)
            summary = diagnostic_summary(rows)
            println(io, join((model, L, values(summary)...), ','))
        end
    end
end

function write_size_components(path, analyses)
    open(path, "w") do io
        println(io, "model,filter,L,environments,component,estimate,bootstrap_low,bootstrap_high")
        for analysis in analyses
            for L in sort(collect(keys(analysis.grouped)))
                estimate = component_statistics(analysis.grouped[L])
                for component in (:conditional, :disorder, :total)
                    values_drawn = getproperty(analysis.draws[L], component)
                    println(io, join((analysis.model, analysis.filter, L,
                        length(analysis.grouped[L]), component, getproperty(estimate, component),
                        percentile(values_drawn, 0.025), percentile(values_drawn, 0.975)), ','))
                end
            end
        end
    end
end

function model_comparisons(analyses, cutoffs, repetitions)
    rows = NamedTuple[]
    for analysis in analyses
        sizes = sort(collect(keys(analysis.grouped)))
        points = Dict(L => component_statistics(analysis.grouped[L]) for L in sizes)
        seen_size_sets = Set{Tuple}()
        for cutoff in cutoffs
            selected = filter(>=(cutoff), sizes)
            length(selected) >= 4 || continue
            selected_key = Tuple(selected)
            selected_key in seen_size_sets && continue
            push!(seen_size_sets, selected_key)
            x = log.(Float64.(selected))
            for component in (:conditional, :disorder, :total)
                y = [getproperty(points[L], component) for L in selected]
                raw_weights = [
                    inv(max(var(getproperty(analysis.draws[L], component); corrected=true),
                            eps(Float64))) for L in selected
                ]
                normalized_weights = raw_weights ./ mean(raw_weights)
                for (method, weights) in (("unweighted", nothing),
                                          ("inverse_bootstrap_variance", normalized_weights))
                    comparison = fit_comparison(x, y; weights=weights)
                    coefficients = Vector{Float64}(undef, repetitions)
                    deltas = similar(coefficients)
                    for repetition in 1:repetitions
                        drawn = [getproperty(analysis.draws[L], component)[repetition]
                                 for L in selected]
                        fitted = fit_comparison(x, drawn; weights=weights)
                        coefficients[repetition] = fitted.quadratic.coefficients[3]
                        deltas[repetition] = fitted.delta_bic
                    end
                    push!(rows, (
                        model=analysis.model,
                        filter=analysis.filter,
                        cutoff=minimum(selected),
                        sizes=join(selected, ';'),
                        component=String(component),
                        method=method,
                        log_intercept=comparison.log.coefficients[1],
                        log_slope=comparison.log.coefficients[2],
                        quadratic_intercept=comparison.quadratic.coefficients[1],
                        quadratic_log_slope=comparison.quadratic.coefficients[2],
                        log2_coefficient=comparison.quadratic.coefficients[3],
                        coefficient_low=percentile(coefficients, 0.025),
                        coefficient_high=percentile(coefficients, 0.975),
                        coefficient_positive_fraction=mean(coefficients .> 0),
                        delta_bic=comparison.delta_bic,
                        delta_bic_low=percentile(deltas, 0.025),
                        delta_bic_high=percentile(deltas, 0.975),
                        quadratic_bic_win_fraction=mean(deltas .> 0),
                        log_loocv_rmse=comparison.log_loocv_rmse,
                        quadratic_loocv_rmse=comparison.quadratic_loocv_rmse,
                    ))
                end
            end
        end
    end
    return rows
end

function write_model_comparisons(path, rows)
    open(path, "w") do io
        println(io, "model,filter,cutoff,sizes,component,method,log_intercept," *
                    "log_slope,quadratic_intercept,quadratic_log_slope," *
                    "log2_coefficient,coefficient_low,coefficient_high," *
                    "coefficient_positive_fraction,delta_bic,delta_bic_low," *
                    "delta_bic_high,quadratic_bic_win_fraction,log_loocv_rmse," *
                    "quadratic_loocv_rmse")
        for row in rows
            println(io, join(values(row), ','))
        end
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return print_help()
    datasets = Dict(
        "gamma" => load_environment_blocks(parsed.gamma_paths),
        "control" => load_environment_blocks(parsed.control_paths),
    )
    analyses = NamedTuple[]
    rng = Xoshiro(parsed.bootstrap_seed)
    for model in ("gamma", "control")
        all_rows = grouped_rows(datasets[model])
        filtered_rows = grouped_rows(
            datasets[model]; diagnostic_threshold=parsed.diagnostic_threshold)
        for (filter_name, groups) in (("all", all_rows),
                                      ("start_gap_le_$(parsed.diagnostic_threshold)",
                                       filtered_rows))
            push!(analyses, (
                model=model,
                filter=filter_name,
                grouped=groups,
                draws=bootstrap_size_draws(rng, groups, parsed.bootstrap_reps),
            ))
        end
    end
    mkpath(parsed.output_dir)
    write_diagnostics(joinpath(parsed.output_dir, "production_diagnostic_summary.csv"), datasets)
    write_size_components(joinpath(parsed.output_dir, "production_component_summary.csv"), analyses)
    comparisons = model_comparisons(analyses, parsed.cutoffs, parsed.bootstrap_reps)
    write_model_comparisons(
        joinpath(parsed.output_dir, "production_scaling_comparison.csv"), comparisons)
    open(joinpath(parsed.output_dir, "ANALYSIS_METHOD.txt"), "w") do io
        println(io, "Independent bootstrap block: one frozen edge-weight environment within one L.")
        println(io, "Paired replica means and conditional variance remain in the same block.")
        println(io, "Primary data filter: all environments.")
        println(io, "Diagnostic sensitivity only: standardized extremal-start gap <= $(parsed.diagnostic_threshold).")
        println(io, "Nested models: a+b*log(L) versus a+b*log(L)+c*log(L)^2.")
        println(io, "Positive delta BIC favors the quadratic extension.")
        println(io, "Prediction metric: leave-one-size-out RMSE; lower is better.")
        println(io, "Weighted sensitivity uses fixed inverse environment-bootstrap variance weights.")
        println(io, "Bootstrap repetitions: $(parsed.bootstrap_reps). Seed: $(parsed.bootstrap_seed).")
    end
    gamma_count = length(datasets["gamma"])
    control_count = length(datasets["control"])
    println("Glauber production scaling analysis complete")
    println("  gamma environments:   $gamma_count")
    println("  control environments: $control_count")
    println("  comparisons:          $(length(comparisons))")
    println("  output:               $(parsed.output_dir)")
end

end # module

abspath(PROGRAM_FILE) == abspath(@__FILE__) && GlauberProductionScaling.main(ARGS)
