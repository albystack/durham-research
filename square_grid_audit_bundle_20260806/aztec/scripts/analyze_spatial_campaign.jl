#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using Statistics

const HEADER =
    "model,order,sample_id,seed,fraction_num,fraction_den,separation," *
    "left_column,right_column,increment_1,increment_2,difference"

function print_help()
    println("""
    Analyse Gamma spatial increments and the uniform control.

    The primary nested comparison is

      H0: a + b log(r)
      H1: a + b log(r) + c (log(r))^2.

    Bootstrap resampling keeps the two replicas paired within each independent
    environment.  The reported coefficient interval is for c.

    Options:
      --gamma-results PATHS   files/directories, comma separated
      --uniform-results PATHS files/directories, comma separated
      --output-dir PATH
      --bootstrap-reps INT    default 2000
      --bootstrap-seed UINT
      --min-order INT         default 128
      --holdout-orders INT    largest orders reserved for prediction (default 2)
      -h, --help
    """)
end

parse_paths(value) = [abspath(strip(path)) for path in split(value, ',') if !isempty(strip(path))]

function parse_arguments(arguments)
    if any(argument -> argument in ("-h", "--help"), arguments)
        print_help()
        return nothing
    end
    options = Dict{String,String}(
        "gamma-results" => joinpath(@__DIR__, "..", "output", "spatial_gamma_smoke"),
        "uniform-results" => joinpath(@__DIR__, "..", "output", "spatial_uniform_smoke"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "spatial_analysis"),
        "bootstrap-reps" => "2000",
        "bootstrap-seed" => "20260805",
        "min-order" => "128",
        "holdout-orders" => "2",
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
        gamma_paths=parse_paths(options["gamma-results"]),
        uniform_paths=parse_paths(options["uniform-results"]),
        output_dir=abspath(options["output-dir"]),
        bootstrap_reps=parse(Int, options["bootstrap-reps"]),
        bootstrap_seed=parse(UInt64, options["bootstrap-seed"]),
        min_order=parse(Int, options["min-order"]),
        holdout_orders=parse(Int, options["holdout-orders"]),
    )
end

function collect_files(paths)
    files = String[]
    for path in paths
        if isfile(path)
            push!(files, path)
        elseif isdir(path)
            for (root, _, names) in walkdir(path), name in names
                startswith(name, "batch_") && endswith(name, ".csv") || continue
                push!(files, joinpath(root, name))
            end
        else
            error("results path does not exist: $path")
        end
    end
    unique!(files)
    sort!(files)
    isempty(files) && error("no batch CSV files found")
    return files
end

function load_results(paths, expected_model)
    # group key: (fraction numerator, fraction denominator, order)
    groups = Dict{Tuple{Int,Int,Int},Vector{NTuple{3,Int}}}()
    seen = Set{Tuple{Int,Int,Int,Int}}()
    seed_by_sample = Dict{Tuple{Int,Int},UInt64}()
    fractions_by_sample = Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}()
    separation_by_group = Dict{Tuple{Int,Int,Int},Int}()
    sample_keys_by_group = Dict{Tuple{Int,Int,Int},Vector{Tuple{Int,UInt64}}}()
    for path in collect_files(paths)
        lines = readlines(path)
        isempty(lines) && error("empty CSV: $path")
        strip(first(lines)) == HEADER || error("unexpected header in $path")
        for (offset, line) in enumerate(lines[2:end])
            isempty(strip(line)) && continue
            fields = split(strip(line), ',')
            length(fields) == 12 || error("malformed row $(offset + 1) in $path")
            fields[1] == expected_model || error("unexpected model in $path")
            order = parse(Int, fields[2])
            sample_id = parse(Int, fields[3])
            seed = parse(UInt64, fields[4])
            numerator = parse(Int, fields[5])
            denominator = parse(Int, fields[6])
            separation = parse(Int, fields[7])
            left_column = parse(Int, fields[8])
            right_column = parse(Int, fields[9])
            increment_1 = parse(Int, fields[10])
            increment_2 = parse(Int, fields[11])
            difference = parse(Int, fields[12])
            difference == increment_1 - increment_2 || error("wrong difference in $path")
            right_column - left_column == separation || error("wrong separation in $path")
            row_key = (order, sample_id, numerator, denominator)
            row_key in seen && error("duplicate spatial row $row_key")
            push!(seen, row_key)
            sample_key = (order, sample_id)
            if haskey(seed_by_sample, sample_key)
                seed_by_sample[sample_key] == seed || error("inconsistent seed for $sample_key")
            else
                seed_by_sample[sample_key] = seed
            end
            push!(get!(fractions_by_sample, sample_key, Set{Tuple{Int,Int}}()), (numerator, denominator))
            group_key = (numerator, denominator, order)
            if haskey(separation_by_group, group_key)
                separation_by_group[group_key] == separation ||
                    error("inconsistent separation for $group_key")
            else
                separation_by_group[group_key] = separation
            end
            push!(get!(groups, group_key, NTuple{3,Int}[]), (increment_1, increment_2, difference))
            push!(
                get!(sample_keys_by_group, group_key, Tuple{Int,UInt64}[]),
                (sample_id, seed),
            )
        end
    end
    expected_fractions = first(values(fractions_by_sample))
    all(fractions == expected_fractions for fractions in values(fractions_by_sample)) ||
        error("not every sample contains the same fractions")
    fractions = sort(collect(expected_fractions); by=fraction -> fraction[1] / fraction[2])
    orders = sort(unique(key[3] for key in keys(groups)))
    for order in orders
        reference = sample_keys_by_group[(fractions[1]..., order)]
        for fraction in fractions[2:end]
            sample_keys_by_group[(fraction..., order)] == reference ||
                error("fraction rows are not aligned by environment at order $order")
        end
    end
    return (
        groups=groups,
        separations=separation_by_group,
        fractions=fractions,
        sample_keys=sample_keys_by_group,
    )
end

function paired_statistics(values, indices=nothing)
    selected = isnothing(indices) ? eachindex(values) : indices
    count_values = length(selected)
    count_values >= 2 || error("paired statistics require at least two samples")
    total_1 = total_2 = total_d = 0.0
    square_1 = square_2 = square_d = product = 0.0
    @inbounds for index in selected
        value_1, value_2, difference = values[index]
        total_1 += value_1
        total_2 += value_2
        total_d += difference
        square_1 += value_1^2
        square_2 += value_2^2
        square_d += difference^2
        product += value_1 * value_2
    end
    denominator = count_values - 1
    variance_1 = (square_1 - total_1^2 / count_values) / denominator
    variance_2 = (square_2 - total_2^2 / count_values) / denominator
    variance_difference = (square_d - total_d^2 / count_values) / denominator
    covariance = (product - total_1 * total_2 / count_values) / denominator
    marginal = (variance_1 + variance_2) / 2
    conditional = variance_difference / 2
    isapprox(marginal - covariance, conditional; atol=1e-9) ||
        error("finite-sample replica identity failed")
    return (
        mean_1=total_1 / count_values,
        mean_2=total_2 / count_values,
        marginal=marginal,
        conditional=conditional,
        disorder=covariance,
    )
end

function bootstrap_group(rng, values, repetitions)
    count_values = length(values)
    marginal = Vector{Float64}(undef, repetitions)
    conditional = similar(marginal)
    disorder = similar(marginal)
    indices = Vector{Int}(undef, count_values)
    for repetition in 1:repetitions
        rand!(rng, indices, 1:count_values)
        statistics = paired_statistics(values, indices)
        marginal[repetition] = statistics.marginal
        conditional[repetition] = statistics.conditional
        disorder[repetition] = statistics.disorder
    end
    return (marginal=marginal, conditional=conditional, disorder=disorder)
end

function percentile(sorted_values, probability)
    position = 1 + (length(sorted_values) - 1) * probability
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return sorted_values[lower]
    return sorted_values[lower] + (position - lower) * (sorted_values[upper] - sorted_values[lower])
end

interval(values) = begin
    sorted_values = sort(values)
    (low=percentile(sorted_values, 0.025), high=percentile(sorted_values, 0.975))
end

function matrix_regression(design, y; weights=nothing)
    length(y) == size(design, 1) || throw(DimensionMismatch("design and response differ"))
    if isnothing(weights)
        coefficients = design \ y
    else
        length(weights) == length(y) || throw(DimensionMismatch("weights and response differ"))
        all(>(0), weights) || error("regression weights must be positive")
        square_root_weights = sqrt.(weights)
        coefficients = (design .* square_root_weights) \ (y .* square_root_weights)
    end
    fitted = design * coefficients
    residual_sum = if isnothing(weights)
        sum(abs2, y - fitted)
    else
        sum(weights .* abs2.(y - fitted))
    end
    count_values = length(y)
    parameter_count = size(design, 2)
    bic = count_values * log(max(residual_sum / count_values, eps(Float64))) +
          parameter_count * log(count_values)
    return (coefficients=coefficients, fitted=fitted, rss=residual_sum, bic=bic)
end

function regression(x, y; quadratic=false, weights=nothing)
    design = quadratic ? hcat(ones(length(x)), x, x .^ 2) : hcat(ones(length(x)), x)
    return matrix_regression(design, y; weights=weights)
end

predict(fit, x; quadratic=false) = quadratic ?
    fit.coefficients[1] .+ fit.coefficients[2] .* x .+ fit.coefficients[3] .* x .^ 2 :
    fit.coefficients[1] .+ fit.coefficients[2] .* x

function fit_comparison(x, y, holdout_orders; weights=nothing)
    log_fit = regression(x, y; weights=weights)
    quadratic_fit = regression(x, y; quadratic=true, weights=weights)
    training_count = length(x) - holdout_orders
    training_count >= 3 || error("not enough training sizes for quadratic held-out fit")
    training = 1:training_count
    testing = (training_count + 1):length(x)
    training_weights = isnothing(weights) ? nothing : weights[training]
    train_log = regression(x[training], y[training]; weights=training_weights)
    train_quadratic = regression(
        x[training],
        y[training];
        quadratic=true,
        weights=training_weights,
    )
    log_rmse = sqrt(mean(abs2, y[testing] - predict(train_log, x[testing])))
    quadratic_rmse = sqrt(mean(abs2, y[testing] - predict(train_quadratic, x[testing]; quadratic=true)))
    return (
        log=log_fit,
        quadratic=quadratic_fit,
        delta_bic=log_fit.bic - quadratic_fit.bic,
        log_holdout_rmse=log_rmse,
        quadratic_holdout_rmse=quadratic_rmse,
    )
end

function bootstrap_fit(
    draws,
    keys,
    separations,
    repetitions,
    holdout_orders,
    component;
    weights=nothing,
)
    x = log.([Float64(separations[key]) for key in keys])
    coefficient = Vector{Float64}(undef, repetitions)
    delta_bic = similar(coefficient)
    for repetition in 1:repetitions
        y = [getproperty(draws[key], component)[repetition] for key in keys]
        comparison = fit_comparison(x, y, holdout_orders; weights=weights)
        coefficient[repetition] = comparison.quadratic.coefficients[3]
        delta_bic[repetition] = comparison.delta_bic
    end
    return (coefficient=coefficient, delta_bic=delta_bic)
end

function write_summary(path, datasets, point_statistics, bootstrap_draws)
    open(path, "w") do io
        println(
            io,
            "model,fraction_num,fraction_den,order,separation,n,mean_1,mean_2," *
            "marginal,marginal_low,marginal_high,conditional,conditional_low," *
            "conditional_high,disorder,disorder_low,disorder_high",
        )
        for model in ("gamma", "uniform")
            data = datasets[model]
            keys_sorted = sort(collect(keys(data.groups)); by=key -> (key[1] / key[2], key[3]))
            for key in keys_sorted
                numerator, denominator, order = key
                point = point_statistics[(model, key...)]
                draws = bootstrap_draws[(model, key...)]
                marginal_interval = interval(draws.marginal)
                conditional_interval = interval(draws.conditional)
                disorder_interval = interval(draws.disorder)
                @printf(
                    io,
                    "%s,%d,%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                    model, numerator, denominator, order, data.separations[key],
                    length(data.groups[key]), point.mean_1, point.mean_2,
                    point.marginal, marginal_interval.low, marginal_interval.high,
                    point.conditional, conditional_interval.low, conditional_interval.high,
                    point.disorder, disorder_interval.low, disorder_interval.high,
                )
            end
        end
    end
end

function write_fits(
    path,
    prediction_path,
    datasets,
    points,
    draws,
    repetitions,
    min_order,
    holdout_orders,
)
    fit_rows = NamedTuple[]
    open(prediction_path, "w") do prediction_io
        println(
            prediction_io,
            "model,fraction_num,fraction_den,component,order,separation,observed,log_fit,log_plus_log2_fit",
        )
        for model in ("gamma", "uniform")
            data = datasets[model]
            for fraction in data.fractions
                numerator, denominator = fraction
                keys_for_fraction = sort(
                    [key for key in keys(data.groups) if key[1:2] == fraction && key[3] >= min_order];
                    by=last,
                )
                length(keys_for_fraction) >= holdout_orders + 3 ||
                    error("too few fitted orders for $model fraction $fraction")
                x = log.([Float64(data.separations[key]) for key in keys_for_fraction])
                for component in (:marginal, :conditional, :disorder)
                    y = [getproperty(points[(model, key...)], component) for key in keys_for_fraction]
                    comparison = fit_comparison(x, y, holdout_orders)
                    bootstrap = bootstrap_fit(
                        Dict(key => draws[(model, key...)] for key in keys_for_fraction),
                        keys_for_fraction,
                        data.separations,
                        repetitions,
                        holdout_orders,
                        component,
                    )
                    coefficient_interval = interval(bootstrap.coefficient)
                    delta_interval = interval(bootstrap.delta_bic)
                    push!(
                        fit_rows,
                        (
                            model=model,
                            fraction_num=numerator,
                            fraction_den=denominator,
                            component=String(component),
                            min_order=minimum(key[3] for key in keys_for_fraction),
                            max_order=maximum(key[3] for key in keys_for_fraction),
                            fit_points=length(keys_for_fraction),
                            log_intercept=comparison.log.coefficients[1],
                            log_slope=comparison.log.coefficients[2],
                            quadratic_intercept=comparison.quadratic.coefficients[1],
                            quadratic_log_slope=comparison.quadratic.coefficients[2],
                            log2_coefficient=comparison.quadratic.coefficients[3],
                            log2_coefficient_low=coefficient_interval.low,
                            log2_coefficient_high=coefficient_interval.high,
                            coefficient_positive_fraction=mean(bootstrap.coefficient .> 0),
                            delta_bic=comparison.delta_bic,
                            delta_bic_low=delta_interval.low,
                            delta_bic_high=delta_interval.high,
                            quadratic_bic_win_fraction=mean(bootstrap.delta_bic .> 0),
                            log_holdout_rmse=comparison.log_holdout_rmse,
                            quadratic_holdout_rmse=comparison.quadratic_holdout_rmse,
                        ),
                    )
                    log_predictions = predict(comparison.log, x)
                    quadratic_predictions = predict(comparison.quadratic, x; quadratic=true)
                    for (key, observed, log_value, quadratic_value) in zip(
                        keys_for_fraction,
                        y,
                        log_predictions,
                        quadratic_predictions,
                    )
                        @printf(
                            prediction_io,
                            "%s,%d,%d,%s,%d,%d,%.10g,%.10g,%.10g\n",
                            model, numerator, denominator, component, key[3],
                            data.separations[key], observed, log_value, quadratic_value,
                        )
                    end
                end
            end
        end
    end
    open(path, "w") do io
        println(
            io,
            "model,fraction_num,fraction_den,component,min_order,max_order,fit_points," *
            "log_intercept,log_slope,quadratic_intercept,quadratic_log_slope," *
            "log2_coefficient,log2_coefficient_low,log2_coefficient_high," *
            "coefficient_positive_fraction,delta_bic_log_minus_quadratic," *
            "delta_bic_low,delta_bic_high,quadratic_bic_win_fraction," *
            "log_holdout_rmse,quadratic_holdout_rmse",
        )
        for row in fit_rows
            @printf(
                io,
                "%s,%d,%d,%s,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                values(row)...,
            )
        end
    end
    return fit_rows
end

function write_weighted_fits(
    path,
    datasets,
    points,
    draws,
    repetitions,
    min_order,
    holdout_orders,
)
    rows = NamedTuple[]
    for model in ("gamma", "uniform")
        data = datasets[model]
        for fraction in data.fractions
            selected = sort(
                [key for key in keys(data.groups) if key[1:2] == fraction && key[3] >= min_order];
                by=last,
            )
            x = log.([Float64(data.separations[key]) for key in selected])
            for component in (:marginal, :conditional, :disorder)
                y = [getproperty(points[(model, key...)], component) for key in selected]
                # The reciprocal bootstrap variance is a plug-in precision
                # weight for each independently estimated size-level statistic.
                weights = [
                    inv(var(getproperty(draws[(model, key...)], component); corrected=true))
                    for key in selected
                ]
                comparison = fit_comparison(x, y, holdout_orders; weights=weights)
                bootstrap = bootstrap_fit(
                    Dict(key => draws[(model, key...)] for key in selected),
                    selected,
                    data.separations,
                    repetitions,
                    holdout_orders,
                    component;
                    weights=weights,
                )
                coefficient_interval = interval(bootstrap.coefficient)
                delta_interval = interval(bootstrap.delta_bic)
                push!(rows, (
                    model=model,
                    fraction_num=fraction[1],
                    fraction_den=fraction[2],
                    component=String(component),
                    log2_coefficient=comparison.quadratic.coefficients[3],
                    log2_coefficient_low=coefficient_interval.low,
                    log2_coefficient_high=coefficient_interval.high,
                    coefficient_positive_fraction=mean(bootstrap.coefficient .> 0),
                    delta_bic=comparison.delta_bic,
                    delta_bic_low=delta_interval.low,
                    delta_bic_high=delta_interval.high,
                    quadratic_bic_win_fraction=mean(bootstrap.delta_bic .> 0),
                    log_holdout_rmse=comparison.log_holdout_rmse,
                    quadratic_holdout_rmse=comparison.quadratic_holdout_rmse,
                ))
            end
        end
    end
    open(path, "w") do io
        println(
            io,
            "model,fraction_num,fraction_den,component,log2_coefficient," *
            "log2_coefficient_low,log2_coefficient_high,coefficient_positive_fraction," *
            "delta_bic_log_minus_quadratic,delta_bic_low,delta_bic_high," *
            "quadratic_bic_win_fraction,log_holdout_rmse,quadratic_holdout_rmse",
        )
        for row in rows
            @printf(
                io,
                "%s,%d,%d,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                values(row)...,
            )
        end
    end
    return rows
end

function pooled_matrices(data, selected_orders, component_values; precision_weights=nothing)
    row_count = length(data.fractions) * length(selected_orders)
    null_design = zeros(row_count, 2 * length(data.fractions))
    quadratic_design = zeros(row_count, 2 * length(data.fractions) + 1)
    response = zeros(row_count)
    weights = isnothing(precision_weights) ? nothing : zeros(row_count)
    row = 1
    for (fraction_index, fraction) in enumerate(data.fractions), order in selected_orders
        key = (fraction..., order)
        log_separation = log(Float64(data.separations[key]))
        null_design[row, 2 * fraction_index - 1] = 1
        null_design[row, 2 * fraction_index] = log_separation
        quadratic_design[row, 1:end-1] .= null_design[row, :]
        quadratic_design[row, end] = log_separation^2
        response[row] = component_values[key]
        !isnothing(weights) && (weights[row] = precision_weights[key])
        row += 1
    end
    return null_design, quadratic_design, response, weights
end

function pooled_comparison(
    data,
    selected_orders,
    component_values,
    holdout_orders;
    precision_weights=nothing,
)
    null_design, quadratic_design, response, weights = pooled_matrices(
        data,
        selected_orders,
        component_values;
        precision_weights=precision_weights,
    )
    null_fit = matrix_regression(null_design, response; weights=weights)
    quadratic_fit = matrix_regression(quadratic_design, response; weights=weights)

    training_orders = selected_orders[1:(end - holdout_orders)]
    testing_orders = selected_orders[(end - holdout_orders + 1):end]
    train_null, train_quadratic, train_response, train_weights = pooled_matrices(
        data,
        training_orders,
        component_values;
        precision_weights=precision_weights,
    )
    test_null, test_quadratic, test_response, _ = pooled_matrices(
        data,
        testing_orders,
        component_values;
        precision_weights=precision_weights,
    )
    trained_null = matrix_regression(train_null, train_response; weights=train_weights)
    trained_quadratic = matrix_regression(
        train_quadratic,
        train_response;
        weights=train_weights,
    )
    null_rmse = sqrt(mean(abs2, test_response - test_null * trained_null.coefficients))
    quadratic_rmse = sqrt(
        mean(abs2, test_response - test_quadratic * trained_quadratic.coefficients),
    )
    return (
        common_coefficient=quadratic_fit.coefficients[end],
        delta_bic=null_fit.bic - quadratic_fit.bic,
        null_holdout_rmse=null_rmse,
        quadratic_holdout_rmse=quadratic_rmse,
    )
end

function write_pooled_fits(
    path,
    datasets,
    points,
    independent_draws,
    repetitions,
    min_order,
    holdout_orders,
    rng,
)
    output_rows = NamedTuple[]
    components = (:marginal, :conditional, :disorder)
    for model in ("gamma", "uniform")
        data = datasets[model]
        selected_orders = sort(
            unique(key[3] for key in keys(data.groups) if key[3] >= min_order),
        )
        length(selected_orders) >= holdout_orders + 3 ||
            error("too few orders for pooled $model fit")

        point_values = Dict(
            component => Dict(
                key => getproperty(points[(model, key...)], component)
                for key in keys(data.groups) if key[3] >= min_order
            )
            for component in components
        )
        precision_weights = Dict(
            component => Dict(
                key => inv(
                    var(
                        getproperty(independent_draws[(model, key...)], component);
                        corrected=true,
                    ),
                )
                for key in keys(data.groups) if key[3] >= min_order
            )
            for component in components
        )
        coefficient_draws = Dict(
            (component, method) => Vector{Float64}(undef, repetitions)
            for component in components for method in (:unweighted, :weighted)
        )
        delta_draws = Dict(
            (component, method) => Vector{Float64}(undef, repetitions)
            for component in components for method in (:unweighted, :weighted)
        )

        # One index sample per order is shared across all fractions.  This is
        # the environment-clustered step that preserves their strong dependence.
        for repetition in 1:repetitions
            bootstrap_values = Dict(
                component => Dict{Tuple{Int,Int,Int},Float64}()
                for component in components
            )
            for order in selected_orders
                reference_key = (data.fractions[1]..., order)
                sample_count = length(data.groups[reference_key])
                indices = rand(rng, 1:sample_count, sample_count)
                for fraction in data.fractions
                    key = (fraction..., order)
                    statistics = paired_statistics(data.groups[key], indices)
                    for component in components
                        bootstrap_values[component][key] = getproperty(statistics, component)
                    end
                end
            end
            for component in components, method in (:unweighted, :weighted)
                weights = method == :weighted ? precision_weights[component] : nothing
                comparison = pooled_comparison(
                    data,
                    selected_orders,
                    bootstrap_values[component],
                    holdout_orders;
                    precision_weights=weights,
                )
                coefficient_draws[(component, method)][repetition] =
                    comparison.common_coefficient
                delta_draws[(component, method)][repetition] = comparison.delta_bic
            end
        end

        for component in components, method in (:unweighted, :weighted)
            weights = method == :weighted ? precision_weights[component] : nothing
            comparison = pooled_comparison(
                data,
                selected_orders,
                point_values[component],
                holdout_orders;
                precision_weights=weights,
            )
            coefficient_interval = interval(coefficient_draws[(component, method)])
            delta_interval = interval(delta_draws[(component, method)])
            push!(output_rows, (
                model=model,
                component=String(component),
                method=String(method),
                common_log2_coefficient=comparison.common_coefficient,
                coefficient_low=coefficient_interval.low,
                coefficient_high=coefficient_interval.high,
                coefficient_positive_fraction=mean(
                    coefficient_draws[(component, method)] .> 0,
                ),
                delta_bic=comparison.delta_bic,
                delta_bic_low=delta_interval.low,
                delta_bic_high=delta_interval.high,
                quadratic_bic_win_fraction=mean(delta_draws[(component, method)] .> 0),
                log_holdout_rmse=comparison.null_holdout_rmse,
                quadratic_holdout_rmse=comparison.quadratic_holdout_rmse,
            ))
        end
    end

    open(path, "w") do io
        println(
            io,
            "model,component,method,common_log2_coefficient,coefficient_low," *
            "coefficient_high,coefficient_positive_fraction," *
            "delta_bic_log_minus_quadratic,delta_bic_low,delta_bic_high," *
            "quadratic_bic_win_fraction,log_holdout_rmse,quadratic_holdout_rmse",
        )
        for row in output_rows
            @printf(
                io,
                "%s,%s,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                values(row)...,
            )
        end
    end
    return output_rows
end

function write_report(
    path,
    fit_rows,
    weighted_rows,
    pooled_rows,
    datasets,
    repetitions,
    min_order,
    holdout_orders,
)
    gamma_rows = [
        row for row in fit_rows
        if row.model == "gamma" && row.component == "disorder"
    ]
    uniform_rows = [
        row for row in fit_rows
        if row.model == "uniform" && row.component == "marginal"
    ]
    pooled_primary = [
        row for row in pooled_rows
        if (row.model == "gamma" && row.component in ("disorder", "conditional")) ||
           (row.model == "uniform" && row.component == "marginal")
    ]
    pooled_order(row) = if row.model == "gamma" && row.component == "disorder"
        1
    elseif row.model == "gamma"
        2
    else
        3
    end
    sort!(pooled_primary; by=row -> (pooled_order(row), row.method == "unweighted" ? 1 : 2))
    gamma_environment_count =
        sum(length(group) for group in values(datasets["gamma"].groups)) ÷
        length(datasets["gamma"].fractions)
    uniform_pair_count =
        sum(length(group) for group in values(datasets["uniform"].groups)) ÷
        length(datasets["uniform"].fractions)
    open(path, "w") do io
        println(io, "# Spatial height-increment analysis")
        println(io)
        println(io, "Independent Gamma environments: **$gamma_environment_count**.")
        println(io, "Independent uniform replica pairs: **$uniform_pair_count**.")
        println(io, "Bootstrap repetitions: **$repetitions**; fitted orders start at **$min_order**; the largest **$holdout_orders** orders test prediction.")
        println(io)
        println(io, "## Joint environment-clustered result")
        println(io)
        println(io, "The four separations share environments. The pooled model therefore gives each separation its own intercept and log slope, imposes one common log-squared coefficient, and resamples whole environments jointly. The weighted sensitivity fit uses inverse bootstrap-variance precision weights.")
        println(io)
        println(io, "| observable | method | common c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |")
        println(io, "|---|---|---:|---:|---:|---:|---:|")
        for row in pooled_primary
            observable = if row.model == "gamma" && row.component == "disorder"
                "Gamma disorder covariance"
            elseif row.model == "gamma"
                "Gamma conditional"
            else
                "Uniform marginal control"
            end
            @printf(
                io,
                "| %s | %s | %.4g | [%.4g, %.4g] | %.3f | %.3f | %.3f / %.3f |\n",
                observable, row.method, row.common_log2_coefficient,
                row.coefficient_low, row.coefficient_high,
                row.coefficient_positive_fraction, row.delta_bic,
                row.log_holdout_rmse, row.quadratic_holdout_rmse,
            )
        end
        println(io)
        println(io, "## Primary Gamma disorder result")
        println(io)
        println(io, "The nested comparison is `a + b log(r)` versus `a + b log(r) + c(log(r))^2`. Positive delta BIC favors the quadratic extension.")
        println(io)
        println(io, "| separation | c | 95% bootstrap interval | P(c>0) | delta BIC | held-out RMSE log / quadratic |")
        println(io, "|---:|---:|---:|---:|---:|---:|")
        for row in sort(gamma_rows; by=row -> row.fraction_num / row.fraction_den)
            @printf(
                io,
                "| %d/%d | %.4g | [%.4g, %.4g] | %.3f | %.3f | %.3f / %.3f |\n",
                row.fraction_num, row.fraction_den, row.log2_coefficient,
                row.log2_coefficient_low, row.log2_coefficient_high,
                row.coefficient_positive_fraction, row.delta_bic,
                row.log_holdout_rmse, row.quadratic_holdout_rmse,
            )
        end
        println(io)
        println(io, "## Uniform control")
        println(io)
        println(io, "| separation | c | 95% bootstrap interval | delta BIC | held-out RMSE log / quadratic |")
        println(io, "|---:|---:|---:|---:|---:|")
        for row in sort(uniform_rows; by=row -> row.fraction_num / row.fraction_den)
            @printf(
                io,
                "| %d/%d | %.4g | [%.4g, %.4g] | %.3f | %.3f / %.3f |\n",
                row.fraction_num, row.fraction_den, row.log2_coefficient,
                row.log2_coefficient_low, row.log2_coefficient_high,
                row.delta_bic, row.log_holdout_rmse, row.quadratic_holdout_rmse,
            )
        end
        println(io)
        println(io, "These finite-size comparisons are numerical evidence, not an asymptotic proof. Dependence among separation fractions is respected by interpreting each fraction as a robustness check rather than four independent experiments.")
        println(io)
        println(io, "The full uncertainty-weighted per-separation sensitivity results are in `spatial_weighted_model_comparison.csv`; they are not used to select a preferred separation after seeing the data.")
    end
end

function write_cutoff_sensitivity(
    path,
    datasets,
    points,
    draws,
    repetitions,
    min_order,
    holdout_orders,
)
    open(path, "w") do io
        println(
            io,
            "model,fraction_num,fraction_den,component,min_order,max_order,fit_points," *
            "log2_coefficient,log2_coefficient_low,log2_coefficient_high," *
            "delta_bic_log_minus_quadratic,delta_bic_low,delta_bic_high",
        )
        for (model, component) in (("gamma", :disorder), ("uniform", :marginal))
            data = datasets[model]
            for fraction in data.fractions
                all_keys = sort(
                    [key for key in keys(data.groups) if key[1:2] == fraction && key[3] >= min_order];
                    by=last,
                )
                # Five points are the minimum for a three-parameter curve plus
                # a meaningful residual comparison.  Successively discarding
                # the smallest sizes exposes finite-size sensitivity.
                for first_index in 1:max(1, length(all_keys) - 4)
                    selected = all_keys[first_index:end]
                    length(selected) >= 5 || continue
                    x = log.([Float64(data.separations[key]) for key in selected])
                    y = [getproperty(points[(model, key...)], component) for key in selected]
                    comparison = fit_comparison(x, y, min(holdout_orders, length(selected) - 3))
                    bootstrap = bootstrap_fit(
                        Dict(key => draws[(model, key...)] for key in selected),
                        selected,
                        data.separations,
                        repetitions,
                        min(holdout_orders, length(selected) - 3),
                        component,
                    )
                    coefficient_interval = interval(bootstrap.coefficient)
                    delta_interval = interval(bootstrap.delta_bic)
                    @printf(
                        io,
                        "%s,%d,%d,%s,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                        model, fraction[1], fraction[2], component, selected[1][3],
                        selected[end][3], length(selected),
                        comparison.quadratic.coefficients[3], coefficient_interval.low,
                        coefficient_interval.high, comparison.delta_bic,
                        delta_interval.low, delta_interval.high,
                    )
                end
            end
        end
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return
    parsed.bootstrap_reps > 0 || error("bootstrap repetitions must be positive")
    parsed.holdout_orders > 0 || error("holdout orders must be positive")
    datasets = Dict(
        "gamma" => load_results(parsed.gamma_paths, "gamma"),
        "uniform" => load_results(parsed.uniform_paths, "uniform"),
    )
    datasets["gamma"].fractions == datasets["uniform"].fractions ||
        error("Gamma and uniform campaigns must use the same fractions")

    rng = Xoshiro(parsed.bootstrap_seed)
    points = Dict{Tuple,NamedTuple}()
    draws = Dict{Tuple,NamedTuple}()
    for model in ("gamma", "uniform")
        for (key, values) in datasets[model].groups
            length(values) >= 2 || error("each spatial group needs at least two samples")
            full_key = (model, key...)
            points[full_key] = paired_statistics(values)
            draws[full_key] = bootstrap_group(rng, values, parsed.bootstrap_reps)
        end
    end

    mkpath(parsed.output_dir)
    summary_path = joinpath(parsed.output_dir, "spatial_summary.csv")
    fits_path = joinpath(parsed.output_dir, "spatial_model_comparison.csv")
    weighted_fits_path = joinpath(
        parsed.output_dir,
        "spatial_weighted_model_comparison.csv",
    )
    pooled_fits_path = joinpath(parsed.output_dir, "spatial_pooled_model_comparison.csv")
    predictions_path = joinpath(parsed.output_dir, "spatial_fit_curves.csv")
    cutoff_path = joinpath(parsed.output_dir, "spatial_cutoff_sensitivity.csv")
    report_path = joinpath(parsed.output_dir, "spatial_analysis_report.md")
    write_summary(summary_path, datasets, points, draws)
    fit_rows = write_fits(
        fits_path,
        predictions_path,
        datasets,
        points,
        draws,
        parsed.bootstrap_reps,
        parsed.min_order,
        parsed.holdout_orders,
    )
    weighted_rows = write_weighted_fits(
        weighted_fits_path,
        datasets,
        points,
        draws,
        parsed.bootstrap_reps,
        parsed.min_order,
        parsed.holdout_orders,
    )
    pooled_rows = write_pooled_fits(
        pooled_fits_path,
        datasets,
        points,
        draws,
        parsed.bootstrap_reps,
        parsed.min_order,
        parsed.holdout_orders,
        Xoshiro(parsed.bootstrap_seed ⊻ 0x706f6f6c65645f62),
    )
    write_cutoff_sensitivity(
        cutoff_path,
        datasets,
        points,
        draws,
        parsed.bootstrap_reps,
        parsed.min_order,
        parsed.holdout_orders,
    )
    write_report(
        report_path,
        fit_rows,
        weighted_rows,
        pooled_rows,
        datasets,
        parsed.bootstrap_reps,
        parsed.min_order,
        parsed.holdout_orders,
    )
    println("Spatial analysis complete")
    println("  output: $(parsed.output_dir)")
    println("  report: $report_path")
end

main(ARGS)
