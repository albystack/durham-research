#!/usr/bin/env julia

"""Analyse matched temporal-IID and ordinary-UST height-variance campaigns.

Both models use the marginal estimator V(r).  Replica covariance is retained
only as an independent-replica control; it is never named disorder covariance.
"""
module TemporalSpatialAnalysisCore
include(joinpath(@__DIR__, "analyze_spatial_campaign.jl"))
end

using Printf
using Random

const Core = TemporalSpatialAnalysisCore
const TEMPORAL_ESTIMATOR = "temporal_marginal_variance"
const ORDINARY_ESTIMATOR = "ordinary_ust_marginal_variance"

function parse_arguments(arguments)
    options = Dict{String,String}(
        "temporal-results" => joinpath(@__DIR__, "..", "output", "temporal_square_grid_smoke"),
        "temporal-model" => "square_grid__temporal_iid__gamma__p_0p5",
        "ordinary-results" => "",
        "ordinary-model" => "square_grid__ordinary_ust",
        "results" => "",
        "model" => "",
        "output-dir" => joinpath(@__DIR__, "..", "output", "temporal_square_grid_analysis"),
        "bootstrap-reps" => "500",
        "bootstrap-seed" => "20260816",
        "min-order" => "32",
        "holdout-orders" => "2",
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument in ("-h", "--help")
            println("Usage: analyze_temporal_square_grid_campaign.jl --temporal-results PATH --temporal-model LABEL [--ordinary-results PATH --ordinary-model LABEL] --output-dir PATH [--bootstrap-reps N --bootstrap-seed N --min-order L --holdout-orders N]")
            return nothing
        end
        startswith(argument, "--") || error("unexpected argument: $argument")
        key = argument[3:end]
        haskey(options, key) || error("unknown option: --$key")
        index < length(arguments) || error("missing value after --$key")
        options[key] = arguments[index + 1]
        index += 2
    end
    !isempty(options["results"]) && (options["temporal-results"] = options["results"])
    !isempty(options["model"]) && (options["temporal-model"] = options["model"])
    return (
        temporal_paths=Core.parse_paths(options["temporal-results"]),
        temporal_model=options["temporal-model"],
        ordinary_paths=isempty(options["ordinary-results"]) ? nothing :
                       Core.parse_paths(options["ordinary-results"]),
        ordinary_model=options["ordinary-model"],
        output_dir=abspath(options["output-dir"]),
        bootstrap_reps=parse(Int, options["bootstrap-reps"]),
        bootstrap_seed=parse(UInt64, options["bootstrap-seed"]),
        min_order=parse(Int, options["min-order"]),
        holdout_orders=parse(Int, options["holdout-orders"]),
    )
end

function summarise(data, repetitions, rng)
    points = Dict{Tuple{Int,Int,Int},NamedTuple}()
    draws = Dict{Tuple{Int,Int,Int},NamedTuple}()
    for (key, values) in data.groups
        length(values) >= 2 || error("each spatial group needs at least two samples")
        points[key] = Core.paired_statistics(values)
        draws[key] = Core.bootstrap_group(rng, values, repetitions)
    end
    return (points=points, draws=draws)
end

function selected_keys(data, fraction, min_order)
    sort([key for key in keys(data.groups) if key[1:2] == fraction && key[3] >= min_order]; by=last)
end

function fit_row(data, points, draws, fraction, min_order, repetitions, holdout_orders)
    keys_for_fraction = selected_keys(data, fraction, min_order)
    length(keys_for_fraction) >= holdout_orders + 3 || error("too few fitted orders at fraction $fraction")
    x = log.([Float64(data.separations[key]) for key in keys_for_fraction])
    y = [points[key].marginal for key in keys_for_fraction]
    comparison = Core.fit_comparison(x, y, holdout_orders)
    bootstrap = Core.bootstrap_fit(Dict(key => draws[key] for key in keys_for_fraction),
                                   keys_for_fraction, data.separations, repetitions,
                                   holdout_orders, :marginal)
    coefficient_ci, bic_ci = Core.interval(bootstrap.coefficient), Core.interval(bootstrap.delta_bic)
    return (
        fraction_num=fraction[1], fraction_den=fraction[2],
        min_order=first(keys_for_fraction)[3], max_order=last(keys_for_fraction)[3],
        fit_points=length(keys_for_fraction),
        log_intercept=comparison.log.coefficients[1], log_slope=comparison.log.coefficients[2],
        log2_coefficient=comparison.quadratic.coefficients[3],
        log2_coefficient_low=coefficient_ci.low, log2_coefficient_high=coefficient_ci.high,
        delta_bic=comparison.delta_bic, delta_bic_low=bic_ci.low, delta_bic_high=bic_ci.high,
        log_holdout_rmse=comparison.log_holdout_rmse,
        quadratic_holdout_rmse=comparison.quadratic_holdout_rmse,
    )
end

function write_summary(path, datasets)
    open(path, "w") do io
        println(io, "model,estimator,fraction_num,fraction_den,order,separation,n,marginal_variance,variance_low,variance_high,independent_replica_covariance,covariance_low,covariance_high,half_difference_variance,half_difference_low,half_difference_high")
        for (label, estimator, data, summary) in datasets
            for key in sort(collect(keys(data.groups)); by=key -> (key[1] / key[2], key[3]))
                point, draw = summary.points[key], summary.draws[key]
                variance_ci, covariance_ci, difference_ci = Core.interval(draw.marginal), Core.interval(draw.disorder), Core.interval(draw.conditional)
                @printf(io, "%s,%s,%d,%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                    label, estimator, key..., data.separations[key], length(data.groups[key]),
                    point.marginal, variance_ci.low, variance_ci.high, point.disorder,
                    covariance_ci.low, covariance_ci.high, point.conditional,
                    difference_ci.low, difference_ci.high)
            end
        end
    end
end

function write_fraction_fits(path, label, estimator, data, summary, repetitions, min_order, holdout_orders)
    rows = [fit_row(data, summary.points, summary.draws, fraction, min_order, repetitions, holdout_orders)
            for fraction in data.fractions]
    open(path, "w") do io
        println(io, "model,estimator,fraction_num,fraction_den,min_order,max_order,fit_points,log_intercept,log_slope,log2_coefficient,log2_coefficient_low,log2_coefficient_high,delta_bic_log_minus_quadratic,delta_bic_low,delta_bic_high,log_holdout_rmse,quadratic_holdout_rmse")
        for row in rows
            @printf(io, "%s,%s,%d,%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n", label, estimator, values(row)...)
        end
    end
end

function pooled_fit(data, summary, repetitions, min_order, holdout_orders, rng)
    orders = sort(unique(key[3] for key in keys(data.groups) if key[3] >= min_order))
    length(orders) >= holdout_orders + 3 || error("too few orders for pooled fit")
    values = Dict(key => summary.points[key].marginal for key in keys(data.groups) if key[3] >= min_order)
    comparison = Core.pooled_comparison(data, orders, values, holdout_orders)
    coefficient_draws, bic_draws = Vector{Float64}(undef, repetitions), Vector{Float64}(undef, repetitions)
    for repetition in 1:repetitions
        boot_values = Dict{Tuple{Int,Int,Int},Float64}()
        for order in orders
            reference = (data.fractions[1]..., order)
            count = length(data.groups[reference])
            indices = rand(rng, 1:count, count)
            for fraction in data.fractions
                key = (fraction..., order)
                boot_values[key] = Core.paired_statistics(data.groups[key], indices).marginal
            end
        end
        boot = Core.pooled_comparison(data, orders, boot_values, holdout_orders)
        coefficient_draws[repetition], bic_draws[repetition] = boot.common_coefficient, boot.delta_bic
    end
    coefficient_ci, bic_ci = Core.interval(coefficient_draws), Core.interval(bic_draws)
    return (min_order=first(orders), max_order=last(orders), fit_points=length(orders) * length(data.fractions),
            common_log2_coefficient=comparison.common_coefficient,
            coefficient_low=coefficient_ci.low, coefficient_high=coefficient_ci.high,
            delta_bic=comparison.delta_bic, delta_bic_low=bic_ci.low, delta_bic_high=bic_ci.high,
            log_holdout_rmse=comparison.null_holdout_rmse,
            quadratic_holdout_rmse=comparison.quadratic_holdout_rmse)
end

function write_pooled_fit(path, label, estimator, row)
    open(path, "w") do io
        println(io, "model,estimator,min_order,max_order,fit_points,common_log2_coefficient,coefficient_low,coefficient_high,delta_bic_log_minus_quadratic,delta_bic_low,delta_bic_high,log_holdout_rmse,quadratic_holdout_rmse")
        @printf(io, "%s,%s,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n", label, estimator, values(row)...)
    end
end

function write_cutoff_sensitivity(path, label, estimator, data, summary, repetitions, min_order, holdout_orders)
    cutoffs = sort(unique(key[3] for key in keys(data.groups) if key[3] >= min_order))
    open(path, "w") do io
        println(io, "model,estimator,fraction_num,fraction_den,min_order,max_order,fit_points,log2_coefficient,log2_coefficient_low,log2_coefficient_high,delta_bic_log_minus_quadratic,log_holdout_rmse,quadratic_holdout_rmse")
        for cutoff in cutoffs, fraction in data.fractions
            length(selected_keys(data, fraction, cutoff)) >= holdout_orders + 3 || continue
            row = fit_row(data, summary.points, summary.draws, fraction, cutoff, repetitions, holdout_orders)
            @printf(io, "%s,%s,%d,%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n", label, estimator, row.fraction_num, row.fraction_den, row.min_order, row.max_order, row.fit_points, row.log2_coefficient, row.log2_coefficient_low, row.log2_coefficient_high, row.delta_bic, row.log_holdout_rmse, row.quadratic_holdout_rmse)
        end
    end
end

function write_differences(path, temporal_data, temporal_summary, ordinary_data, ordinary_summary)
    temporal_data.fractions == ordinary_data.fractions || error("temporal and ordinary campaigns use different rho values")
    Set(keys(temporal_data.groups)) == Set(keys(ordinary_data.groups)) ||
        error("temporal and ordinary campaigns use different L/rho schedule")
    open(path, "w") do io
        println(io, "fraction_num,fraction_den,order,separation,temporal_minus_ust,low,high")
        for key in sort(collect(keys(temporal_data.groups)); by=key -> (key[1] / key[2], key[3]))
            point = temporal_summary.points[key].marginal - ordinary_summary.points[key].marginal
            draws = temporal_summary.draws[key].marginal .- ordinary_summary.draws[key].marginal
            ci = Core.interval(draws)
            @printf(io, "%d,%d,%d,%d,%.10g,%.10g,%.10g\n", key..., temporal_data.separations[key], point, ci.low, ci.high)
        end
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return
    parsed.bootstrap_reps > 0 || error("bootstrap repetitions must be positive")
    parsed.holdout_orders > 0 || error("holdout orders must be positive")
    temporal_data = Core.load_results(parsed.temporal_paths, parsed.temporal_model)
    temporal_summary = summarise(temporal_data, parsed.bootstrap_reps, Xoshiro(parsed.bootstrap_seed))
    datasets = [(parsed.temporal_model, TEMPORAL_ESTIMATOR, temporal_data, temporal_summary)]
    ordinary_data = ordinary_summary = nothing
    if !isnothing(parsed.ordinary_paths)
        ordinary_data = Core.load_results(parsed.ordinary_paths, parsed.ordinary_model)
        ordinary_summary = summarise(ordinary_data, parsed.bootstrap_reps, Xoshiro(parsed.bootstrap_seed ⊻ 0x6f7264696e617279))
        push!(datasets, (parsed.ordinary_model, ORDINARY_ESTIMATOR, ordinary_data, ordinary_summary))
    end
    mkpath(parsed.output_dir)
    write_summary(joinpath(parsed.output_dir, "temporal_ust_marginal_variance_summary.csv"), datasets)
    for (index, (label, estimator, data, summary)) in enumerate(datasets)
        tag = index == 1 ? "temporal" : "ordinary_ust"
        write_fraction_fits(joinpath(parsed.output_dir, "$(tag)_marginal_variance_model_comparison.csv"), label, estimator, data, summary, parsed.bootstrap_reps, parsed.min_order, parsed.holdout_orders)
        pooled = pooled_fit(data, summary, parsed.bootstrap_reps, parsed.min_order, parsed.holdout_orders, Xoshiro(parsed.bootstrap_seed ⊻ UInt64(index)))
        write_pooled_fit(joinpath(parsed.output_dir, "$(tag)_marginal_variance_pooled_model_comparison.csv"), label, estimator, pooled)
        write_cutoff_sensitivity(joinpath(parsed.output_dir, "$(tag)_marginal_variance_cutoff_sensitivity.csv"), label, estimator, data, summary, parsed.bootstrap_reps, parsed.min_order, parsed.holdout_orders)
    end
    if !isnothing(ordinary_data)
        write_differences(joinpath(parsed.output_dir, "temporal_minus_ust_marginal_variance.csv"), temporal_data, temporal_summary, ordinary_data, ordinary_summary)
    end
    open(joinpath(parsed.output_dir, "analysis_metadata.txt"), "w") do io
        println(io, "temporal_estimator=$(TEMPORAL_ESTIMATOR)")
        println(io, "ordinary_estimator=$(isnothing(ordinary_data) ? "not_supplied" : ORDINARY_ESTIMATOR)")
        println(io, "replica_covariance_role=independence_control_not_disorder_covariance")
        println(io, "bootstrap_reps=$(parsed.bootstrap_reps)")
        println(io, "bootstrap_seed=$(parsed.bootstrap_seed)")
        println(io, "min_order=$(parsed.min_order)")
        println(io, "holdout_orders=$(parsed.holdout_orders)")
    end
    println("Temporal/ordinary-UST marginal spatial analysis complete: $(parsed.output_dir)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
