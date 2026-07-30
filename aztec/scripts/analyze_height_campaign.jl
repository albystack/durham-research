#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using Statistics

function parse_arguments(arguments)
    options = Dict{String,String}(
        "results-dir" =>
            joinpath(@__DIR__, "..", "data", "height", "center_height_samples.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "gamma_height_analysis"),
        "bootstrap-reps" => "2000",
        "bootstrap-seed" => "20260729",
        "min-order" => "24",
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
    return (
        results_dir=abspath(options["results-dir"]),
        output_dir=abspath(options["output-dir"]),
        bootstrap_reps=parse(Int, options["bootstrap-reps"]),
        bootstrap_seed=parse(UInt64, options["bootstrap-seed"]),
        min_order=parse(Int, options["min-order"]),
    )
end

function load_results(results_path)
    grouped = Dict{Int,Vector{Int}}()
    seen_ids = Dict{Int,Set{Int}}()
    files = String[]
    if isfile(results_path)
        push!(files, results_path)
    elseif isdir(results_path)
        for (root, _, names) in walkdir(results_path)
            for name in names
                startswith(name, "batch_") && endswith(name, ".csv") || continue
                push!(files, joinpath(root, name))
            end
        end
    else
        error("results path does not exist: $results_path")
    end
    isempty(files) && error("no batch CSV files found under $results_path")

    for path in sort(files)
        lines = readlines(path)
        first(lines) ==
        "order,sample_id,seed,center_row,center_column,center_height" ||
            error("unexpected header in $path")
        for line in lines[2:end]
            fields = split(line, ',')
            length(fields) == 6 || error("malformed row in $path")
            order = parse(Int, fields[1])
            sample_id = parse(Int, fields[2])
            height = parse(Int, fields[6])
            ids = get!(seen_ids, order, Set{Int}())
            sample_id in ids && error("duplicate sample id $sample_id at order $order")
            push!(ids, sample_id)
            push!(get!(grouped, order, Int[]), height)
        end
    end
    return grouped
end

function percentile(sorted_values, probability)
    isempty(sorted_values) && error("cannot take percentile of empty data")
    position = 1 + (length(sorted_values) - 1) * probability
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return sorted_values[lower]
    fraction = position - lower
    return (1 - fraction) * sorted_values[lower] + fraction * sorted_values[upper]
end

function linear_fit(x, y)
    design = hcat(ones(length(x)), x)
    coefficients = design \ y
    fitted = design * coefficients
    residuals = y - fitted
    rss = sum(abs2, residuals)
    count = length(y)
    bic = count * log(max(rss / count, eps(Float64))) + 2 * log(count)
    return (
        intercept=coefficients[1],
        slope=coefficients[2],
        fitted=fitted,
        rss=rss,
        bic=bic,
    )
end

function exponent_fit(orders, variances)
    x = log.(log.(Float64.(orders)))
    y = log.(variances)
    fit = linear_fit(x, y)
    return (
        prefactor=exp(fit.intercept),
        exponent=fit.slope,
        rss_log_scale=fit.rss,
        bic_log_scale=fit.bic,
    )
end

function bootstrap_statistics(rng, grouped, orders, repetitions)
    variance_draws = Dict(order => Vector{Float64}(undef, repetitions) for order in orders)
    exponent_draws = Vector{Float64}(undef, repetitions)
    delta_bic_draws = Vector{Float64}(undef, repetitions)
    log_orders = log.(Float64.(orders))
    for repetition in 1:repetitions
        variances = Vector{Float64}(undef, length(orders))
        for (index, order) in enumerate(orders)
            values = grouped[order]
            count_values = length(values)
            resampled = Vector{Float64}(undef, count_values)
            for draw in eachindex(resampled)
                resampled[draw] = values[rand(rng, 1:count_values)]
            end
            estimate = var(resampled; corrected=true)
            variance_draws[order][repetition] = estimate
            variances[index] = estimate
        end
        exponent_draws[repetition] = exponent_fit(orders, variances).exponent
        log_fit = linear_fit(log_orders, variances)
        log2_fit = linear_fit(log_orders .^ 2, variances)
        delta_bic_draws[repetition] = log_fit.bic - log2_fit.bic
    end
    return variance_draws, exponent_draws, delta_bic_draws
end

function write_summary(path, grouped, orders, variance_draws)
    open(path, "w") do io
        println(
            io,
            "order,n_samples,mean_height,variance_height,std_height,",
            "variance_bootstrap_low,variance_bootstrap_high",
        )
        for order in orders
            values = grouped[order]
            draws = sort(variance_draws[order])
            @printf(
                io,
                "%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                order,
                length(values),
                mean(values),
                var(values; corrected=true),
                std(values; corrected=true),
                percentile(draws, 0.025),
                percentile(draws, 0.975),
            )
        end
    end
end

function write_fits(
    path,
    fit_orders,
    fit_variances,
    exponent,
    exponent_draws,
    delta_bic_draws,
)
    log_fit = linear_fit(log.(Float64.(fit_orders)), fit_variances)
    log2_fit = linear_fit(log.(Float64.(fit_orders)) .^ 2, fit_variances)
    sorted_exponents = sort(exponent_draws)
    sorted_delta_bic = sort(delta_bic_draws)
    open(path, "w") do io
        println(io, "fit_min_order=$(minimum(fit_orders))")
        println(io, "fit_max_order=$(maximum(fit_orders))")
        println(io, "fit_points=$(length(fit_orders))")
        println(io, "log_intercept=$(log_fit.intercept)")
        println(io, "log_slope=$(log_fit.slope)")
        println(io, "log_rss=$(log_fit.rss)")
        println(io, "log_bic=$(log_fit.bic)")
        println(io, "log2_intercept=$(log2_fit.intercept)")
        println(io, "log2_slope=$(log2_fit.slope)")
        println(io, "log2_rss=$(log2_fit.rss)")
        println(io, "log2_bic=$(log2_fit.bic)")
        println(io, "delta_bic_log_minus_log2=$(log_fit.bic - log2_fit.bic)")
        println(
            io,
            "delta_bic_bootstrap_low=$(percentile(sorted_delta_bic, 0.025))",
        )
        println(
            io,
            "delta_bic_bootstrap_high=$(percentile(sorted_delta_bic, 0.975))",
        )
        println(
            io,
            "bootstrap_fraction_favoring_log2=$(mean(delta_bic_draws .> 0))",
        )
        println(io, "power_prefactor=$(exponent.prefactor)")
        println(io, "power_exponent=$(exponent.exponent)")
        println(io, "power_exponent_bootstrap_low=$(percentile(sorted_exponents, 0.025))")
        println(io, "power_exponent_bootstrap_high=$(percentile(sorted_exponents, 0.975))")
    end
    return log_fit, log2_fit
end

function write_curve(path, fit_orders, log_fit, log2_fit, exponent)
    open(path, "w") do io
        println(io, "order,log_fit,log2_fit,power_fit")
        for order in fit_orders
            log_order = log(Float64(order))
            @printf(
                io,
                "%d,%.10g,%.10g,%.10g\n",
                order,
                log_fit.intercept + log_fit.slope * log_order,
                log2_fit.intercept + log2_fit.slope * log_order^2,
                exponent.prefactor * log_order^exponent.exponent,
            )
        end
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    parsed.bootstrap_reps > 0 || error("--bootstrap-reps must be positive")
    grouped = load_results(parsed.results_dir)
    orders = sort(collect(keys(grouped)))
    all(length(grouped[order]) >= 2 for order in orders) ||
        error("each order requires at least two samples")
    fit_orders = filter(>=(parsed.min_order), orders)
    length(fit_orders) >= 3 || error("at least three fitted orders are required")
    fit_variances = [var(grouped[order]; corrected=true) for order in fit_orders]
    all(>(0), fit_variances) || error("all fitted variances must be positive")

    rng = Xoshiro(parsed.bootstrap_seed)
    variance_draws, exponent_draws, delta_bic_draws = bootstrap_statistics(
        rng,
        grouped,
        fit_orders,
        parsed.bootstrap_reps,
    )
    for order in setdiff(orders, fit_orders)
        values = grouped[order]
        draws = Vector{Float64}(undef, parsed.bootstrap_reps)
        for repetition in eachindex(draws)
            resampled = [values[rand(rng, eachindex(values))] for _ in eachindex(values)]
            draws[repetition] = var(resampled; corrected=true)
        end
        variance_draws[order] = draws
    end

    mkpath(parsed.output_dir)
    summary_path = joinpath(parsed.output_dir, "height_summary.csv")
    fits_path = joinpath(parsed.output_dir, "height_fits.txt")
    curve_path = joinpath(parsed.output_dir, "height_fit_curves.csv")
    write_summary(summary_path, grouped, orders, variance_draws)
    exponent = exponent_fit(fit_orders, fit_variances)
    log_fit, log2_fit = write_fits(
        fits_path,
        fit_orders,
        fit_variances,
        exponent,
        exponent_draws,
        delta_bic_draws,
    )
    write_curve(curve_path, fit_orders, log_fit, log2_fit, exponent)

    sorted_exponents = sort(exponent_draws)
    println("Height analysis complete")
    println("  samples:        $(sum(length(sample_values) for sample_values in values(grouped)))")
    println("  orders:         $(join(orders, ", "))")
    @printf(
        "  exponent p:     %.4f (95%% bootstrap %.4f to %.4f)\n",
        exponent.exponent,
        percentile(sorted_exponents, 0.025),
        percentile(sorted_exponents, 0.975),
    )
    @printf(
        "  delta BIC:      %.3f (positive favors log^2)\n",
        log_fit.bic - log2_fit.bic,
    )
    @printf(
        "  log^2 wins:     %.1f%% of bootstrap resamples\n",
        100 * mean(delta_bic_draws .> 0),
    )
    println("  output:         $(parsed.output_dir)")
end

main(ARGS)
