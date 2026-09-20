#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using Statistics

const PAIR_HEADER =
    "order,sample_id,seed,center_row,center_column,height_1,height_2,height_difference"
const SINGLE_HEADER =
    "order,sample_id,seed,center_row,center_column,center_height"

# For a shared environment E and conditionally i.i.d. heights H1,H2:
#
#   Var(H1-H2)/2 = E[Var(H | E)]
#   Cov(H1,H2)   = Var(E[H | E]).
#
# The first identity isolates tiling noise; the second directly measures the
# disorder contribution.  Bootstrap resampling always keeps (H1,H2) paired.
# As in the single-height analysis, BIC is a heuristic unweighted comparison
# of two affine curves fitted to size-level variance/covariance estimates.

function print_help()
    println("""
    Analyse double-dimer center-height differences and variance decomposition.

    Usage:
      julia --project=aztec aztec/scripts/analyze_double_dimer_campaign.jl [options]

    Options:
      --paired-results PATHS  paired CSV files/directories, comma separated
      --single-results PATHS  single-height CSV files/directories, comma separated
      --output-dir PATH       output directory
      --bootstrap-reps INT    within-size bootstrap repetitions (default: 2000)
      --bootstrap-seed UINT   deterministic bootstrap seed
      --min-order INT         smallest size included in fits (default: 24)
      -h, --help              show this message
    """)
end

function parse_paths(value)
    return [abspath(strip(path)) for path in split(value, ',') if !isempty(strip(path))]
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict{String,String}(
        "paired-results" => joinpath(@__DIR__, "..", "data", "double_dimer", "pairs.csv"),
        "single-results" =>
            joinpath(@__DIR__, "..", "data", "height", "center_height_samples.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "double_dimer_analysis"),
        "bootstrap-reps" => "2000",
        "bootstrap-seed" => "20260802",
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
    paired_paths = parse_paths(options["paired-results"])
    single_paths = parse_paths(options["single-results"])
    isempty(paired_paths) && error("--paired-results cannot be empty")
    isempty(single_paths) && error("--single-results cannot be empty")
    return (
        paired_paths=paired_paths,
        single_paths=single_paths,
        output_dir=abspath(options["output-dir"]),
        bootstrap_reps=parse(Int, options["bootstrap-reps"]),
        bootstrap_seed=parse(UInt64, options["bootstrap-seed"]),
        min_order=parse(Int, options["min-order"]),
    )
end

function collect_files(paths; batch_prefix="batch_")
    files = String[]
    for path in paths
        if isfile(path)
            push!(files, path)
        elseif isdir(path)
            for (root, _, names) in walkdir(path)
                for name in names
                    startswith(name, batch_prefix) && endswith(name, ".csv") || continue
                    push!(files, joinpath(root, name))
                end
            end
        else
            error("results path does not exist: $path")
        end
    end
    unique!(files)
    sort!(files)
    isempty(files) && error("no result CSV files found")
    return files
end

function load_pairs(paths)
    grouped = Dict{Int,Vector{NTuple{3,Int}}}()
    keys_seen = Set{Tuple{Int,Int}}()
    seeds_seen = Set{UInt64}()
    for path in collect_files(paths)
        lines = readlines(path)
        isempty(lines) && error("empty CSV: $path")
        strip(first(lines)) == PAIR_HEADER || error("unexpected paired header in $path")
        for (offset, line) in enumerate(lines[2:end])
            line_number = offset + 1
            isempty(strip(line)) && continue
            fields = split(strip(line), ',')
            length(fields) == 8 || error("malformed row $line_number in $path")
            order = parse(Int, fields[1])
            sample_id = parse(Int, fields[2])
            seed = parse(UInt64, fields[3])
            parse(Int, fields[4]) == order + 1 || error("wrong center row in $path")
            parse(Int, fields[5]) == fld(order, 2) + 1 || error("wrong center column in $path")
            height_1 = parse(Int, fields[6])
            height_2 = parse(Int, fields[7])
            difference = parse(Int, fields[8])
            difference == height_1 - height_2 || error("wrong difference in $path")
            key = (order, sample_id)
            key in keys_seen && error("duplicate paired key $key")
            seed in seeds_seen && error("duplicate paired seed $seed")
            push!(keys_seen, key)
            push!(seeds_seen, seed)
            push!(get!(grouped, order, NTuple{3,Int}[]), (height_1, height_2, difference))
        end
    end
    return grouped
end

function load_single(paths)
    grouped = Dict{Int,Vector{Int}}()
    keys_seen = Set{Tuple{Int,Int}}()
    seeds_seen = Set{UInt64}()
    for path in collect_files(paths; batch_prefix="batch_")
        lines = readlines(path)
        isempty(lines) && error("empty CSV: $path")
        strip(first(lines)) == SINGLE_HEADER || error("unexpected single header in $path")
        for (offset, line) in enumerate(lines[2:end])
            line_number = offset + 1
            isempty(strip(line)) && continue
            fields = split(strip(line), ',')
            length(fields) == 6 || error("malformed row $line_number in $path")
            order = parse(Int, fields[1])
            sample_id = parse(Int, fields[2])
            seed = parse(UInt64, fields[3])
            parse(Int, fields[4]) == order + 1 || error("wrong center row in $path")
            parse(Int, fields[5]) == fld(order, 2) + 1 || error("wrong center column in $path")
            key = (order, sample_id)
            key in keys_seen && error("duplicate single key $key")
            seed in seeds_seen && error("duplicate single seed $seed")
            push!(keys_seen, key)
            push!(seeds_seen, seed)
            push!(get!(grouped, order, Int[]), parse(Int, fields[6]))
        end
    end
    return grouped
end

function percentile(sorted_values, probability)
    isempty(sorted_values) && error("cannot take percentile of empty data")
    0 <= probability <= 1 || error("percentile probability must lie in [0, 1]")
    position = 1 + (length(sorted_values) - 1) * probability
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return sorted_values[lower]
    fraction = position - lower
    return (1 - fraction) * sorted_values[lower] + fraction * sorted_values[upper]
end

function linear_fit(x, y)
    length(x) == length(y) || throw(DimensionMismatch("x and y lengths differ"))
    length(x) >= 2 || error("linear fit needs at least two observations")
    design = hcat(ones(length(x)), x)
    coefficients = design \ y
    fitted = design * coefficients
    rss = sum(abs2, y - fitted)
    count = length(y)
    bic = count * log(max(rss / count, eps(Float64))) + 2 * log(count)
    return (intercept=coefficients[1], slope=coefficients[2], fitted=fitted, rss=rss, bic=bic)
end

function exponent_fit(orders, variances)
    all(>(1), orders) || error("power fit requires orders greater than one")
    all(>(0), variances) || error("power fit requires positive values")
    fit = linear_fit(log.(log.(Float64.(orders))), log.(variances))
    return (prefactor=exp(fit.intercept), exponent=fit.slope)
end

function variance_from_indices(values, indices)
    count_values = length(indices)
    count_values >= 2 || error("variance requires at least two draws")
    total = 0.0
    total_squared = 0.0
    @inbounds for index in indices
        value = values[index]
        total += value
        total_squared += value * value
    end
    return (total_squared - total^2 / count_values) / (count_values - 1)
end

function covariance_from_indices(values, indices)
    count_values = length(indices)
    count_values >= 2 || error("covariance requires at least two pairs")
    total_1 = 0.0
    total_2 = 0.0
    total_product = 0.0
    @inbounds for index in indices
        height_1, height_2, _ = values[index]
        total_1 += height_1
        total_2 += height_2
        total_product += height_1 * height_2
    end
    return (total_product - total_1 * total_2 / count_values) / (count_values - 1)
end

function bootstrap_statistics(rng, pairs, singles, orders, fit_orders, repetitions)
    pair_variance_draws = Dict(order => Vector{Float64}(undef, repetitions) for order in orders)
    conditional_draws = Dict(order => Vector{Float64}(undef, repetitions) for order in orders)
    covariance_draws = Dict(order => Vector{Float64}(undef, repetitions) for order in orders)
    total_draws = Dict(order => Vector{Float64}(undef, repetitions) for order in orders)
    disorder_draws = Dict(order => Vector{Float64}(undef, repetitions) for order in orders)
    # Difference variance may be zero in a tiny bootstrap resample.  Affine
    # fits still work, but the log-scale power fit does not.
    exponent_draws = Float64[]
    delta_bic_draws = Vector{Float64}(undef, repetitions)
    log_orders = log.(Float64.(fit_orders))

    # Differences do not change between bootstrap repetitions, so build them
    # once rather than allocating 19 new vectors on every repetition.
    differences_by_order = Dict(
        order => [value[3] for value in pairs[order]]
        for order in orders
    )

    for repetition in 1:repetitions
        fitted_variances = Float64[]
        for order in orders
            pair_values = pairs[order]
            differences = differences_by_order[order]
            pair_indices = rand(rng, 1:length(differences), length(differences))
            pair_variance = variance_from_indices(differences, pair_indices)
            conditional = pair_variance / 2
            paired_covariance = covariance_from_indices(pair_values, pair_indices)

            single_values = singles[order]
            single_indices = rand(rng, 1:length(single_values), length(single_values))
            total = variance_from_indices(single_values, single_indices)
            disorder = total - conditional

            pair_variance_draws[order][repetition] = pair_variance
            conditional_draws[order][repetition] = conditional
            covariance_draws[order][repetition] = paired_covariance
            total_draws[order][repetition] = total
            disorder_draws[order][repetition] = disorder
            order in fit_orders && push!(fitted_variances, pair_variance)
        end
        all(>(0), fitted_variances) &&
            push!(exponent_draws, exponent_fit(fit_orders, fitted_variances).exponent)
        log_fit = linear_fit(log_orders, fitted_variances)
        log2_fit = linear_fit(log_orders .^ 2, fitted_variances)
        delta_bic_draws[repetition] = log_fit.bic - log2_fit.bic
    end
    return (
        pair_variance=pair_variance_draws,
        conditional=conditional_draws,
        covariance=covariance_draws,
        total=total_draws,
        disorder=disorder_draws,
        exponent=exponent_draws,
        delta_bic=delta_bic_draws,
    )
end

interval(draws) = begin
    sorted = sort(draws)
    (low=percentile(sorted, 0.025), high=percentile(sorted, 0.975))
end

function write_summary(path, pairs, singles, orders, draws)
    open(path, "w") do io
        println(
            io,
            "order,n_pairs,mean_h1,mean_h2,mean_difference,variance_h1,variance_h2,",
            "covariance_h1_h2,covariance_low,covariance_high,correlation_h1_h2,",
            "paired_marginal_variance,variance_difference,",
            "variance_difference_low,variance_difference_high,conditional_tiling_variance,",
            "conditional_tiling_low,conditional_tiling_high,n_single,total_single_variance,",
            "total_single_low,total_single_high,disorder_variance,disorder_low,disorder_high",
        )
        for order in orders
            values = pairs[order]
            h1 = [value[1] for value in values]
            h2 = [value[2] for value in values]
            differences = [value[3] for value in values]
            single_values = singles[order]
            variance_difference = var(differences; corrected=true)
            conditional = variance_difference / 2
            total = var(single_values; corrected=true)
            disorder = total - conditional
            pair_interval = interval(draws.pair_variance[order])
            conditional_interval = interval(draws.conditional[order])
            total_interval = interval(draws.total[order])
            disorder_interval = interval(draws.disorder[order])
            covariance_value = cov(h1, h2; corrected=true)
            covariance_interval = interval(draws.covariance[order])
            paired_marginal = (var(h1; corrected=true) + var(h2; corrected=true)) / 2

            # This equality holds algebraically even for finite samples when
            # all variances/covariances use the same corrected denominator.
            # Retaining the check catches column-order and pairing mistakes.
            isapprox(paired_marginal - covariance_value, conditional; atol=1e-10) ||
                error("finite-sample variance identity failed at order $order")
            @printf(
                io,
                "%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,",
                order, length(values), mean(h1), mean(h2), mean(differences),
                var(h1; corrected=true), var(h2; corrected=true), covariance_value,
                covariance_interval.low, covariance_interval.high, cor(h1, h2),
                paired_marginal, variance_difference,
            )
            @printf(
                io,
                "%.10g,%.10g,%.10g,%.10g,%.10g,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                pair_interval.low, pair_interval.high, conditional,
                conditional_interval.low, conditional_interval.high, length(single_values),
                total, total_interval.low, total_interval.high, disorder,
                disorder_interval.low, disorder_interval.high,
            )
        end
    end
end

function write_fits(path, fit_orders, fit_variances, draws)
    log_orders = log.(Float64.(fit_orders))
    log_fit = linear_fit(log_orders, fit_variances)
    log2_fit = linear_fit(log_orders .^ 2, fit_variances)
    exponent = exponent_fit(fit_orders, fit_variances)
    exponent_interval = interval(draws.exponent)
    delta_interval = interval(draws.delta_bic)
    open(path, "w") do io
        println(io, "observable=variance(height_1-height_2)")
        println(io, "fit_min_order=$(minimum(fit_orders))")
        println(io, "fit_max_order=$(maximum(fit_orders))")
        println(io, "fit_points=$(length(fit_orders))")
        println(io, "bootstrap_reps=$(length(draws.delta_bic))")
        println(io, "log_intercept=$(log_fit.intercept)")
        println(io, "log_slope=$(log_fit.slope)")
        println(io, "log_rss=$(log_fit.rss)")
        println(io, "log_bic=$(log_fit.bic)")
        println(io, "log2_intercept=$(log2_fit.intercept)")
        println(io, "log2_slope=$(log2_fit.slope)")
        println(io, "log2_rss=$(log2_fit.rss)")
        println(io, "log2_bic=$(log2_fit.bic)")
        println(io, "delta_bic_log_minus_log2=$(log_fit.bic - log2_fit.bic)")
        println(io, "delta_bic_bootstrap_low=$(delta_interval.low)")
        println(io, "delta_bic_bootstrap_high=$(delta_interval.high)")
        println(io, "bootstrap_fraction_favoring_log2=$(mean(draws.delta_bic .> 0))")
        println(io, "power_prefactor=$(exponent.prefactor)")
        println(io, "power_exponent=$(exponent.exponent)")
        println(io, "power_exponent_bootstrap_low=$(exponent_interval.low)")
        println(io, "power_exponent_bootstrap_high=$(exponent_interval.high)")
        println(io, "power_positive_bootstrap_reps=$(length(draws.exponent))")
    end
    return log_fit, log2_fit, exponent
end

function write_curves(path, fit_orders, log_fit, log2_fit, exponent)
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

function component_bootstrap(draws_by_order, fit_orders)
    repetitions = length(draws_by_order[first(fit_orders)])
    exponent_draws = Float64[]
    delta_bic_draws = Vector{Float64}(undef, repetitions)
    log_orders = log.(Float64.(fit_orders))
    for repetition in 1:repetitions
        values = [draws_by_order[order][repetition] for order in fit_orders]
        log_fit = linear_fit(log_orders, values)
        log2_fit = linear_fit(log_orders .^ 2, values)
        delta_bic_draws[repetition] = log_fit.bic - log2_fit.bic
        # The affine fits remain defined if a noisy covariance/remainder draw
        # is negative.  The logarithmic power fit does not, so record its
        # interval only from positive replicates and report the retained count.
        all(>(0), values) &&
            push!(exponent_draws, exponent_fit(fit_orders, values).exponent)
    end
    return (exponent=exponent_draws, delta_bic=delta_bic_draws)
end

function write_component_fits(path, pairs, singles, fit_orders, draws)
    pair_variances = Dict(
        order => var([value[3] for value in pairs[order]]; corrected=true)
        for order in fit_orders
    )
    total_values = Dict(order => var(singles[order]; corrected=true) for order in fit_orders)
    conditional_values = Dict(order => pair_variances[order] / 2 for order in fit_orders)
    disorder_values = Dict(
        order => total_values[order] - conditional_values[order]
        for order in fit_orders
    )
    components = (
        (name="total_single", point=total_values, bootstrap=draws.total),
        (name="conditional_tiling", point=conditional_values, bootstrap=draws.conditional),
        (
            name="disorder_paired_covariance",
            point=Dict(
                order => cov(
                    [value[1] for value in pairs[order]],
                    [value[2] for value in pairs[order]];
                    corrected=true,
                )
                for order in fit_orders
            ),
            bootstrap=draws.covariance,
        ),
        (name="disorder", point=disorder_values, bootstrap=draws.disorder),
    )
    log_orders = log.(Float64.(fit_orders))

    open(path, "w") do io
        println(io, "fit_min_order=$(minimum(fit_orders))")
        println(io, "fit_max_order=$(maximum(fit_orders))")
        println(io, "fit_points=$(length(fit_orders))")
        println(io, "bootstrap_reps=$(length(draws.delta_bic))")
        for component in components
            values = [component.point[order] for order in fit_orders]
            all(>(0), values) || error("$(component.name) has a non-positive fitted value")
            log_fit = linear_fit(log_orders, values)
            log2_fit = linear_fit(log_orders .^ 2, values)
            exponent = exponent_fit(fit_orders, values)
            bootstrap = component_bootstrap(component.bootstrap, fit_orders)
            exponent_interval = interval(bootstrap.exponent)
            delta_interval = interval(bootstrap.delta_bic)
            prefix = component.name
            println(io, "$(prefix)_log_intercept=$(log_fit.intercept)")
            println(io, "$(prefix)_log_slope=$(log_fit.slope)")
            println(io, "$(prefix)_log_bic=$(log_fit.bic)")
            println(io, "$(prefix)_log2_intercept=$(log2_fit.intercept)")
            println(io, "$(prefix)_log2_slope=$(log2_fit.slope)")
            println(io, "$(prefix)_log2_bic=$(log2_fit.bic)")
            println(io, "$(prefix)_delta_bic_log_minus_log2=$(log_fit.bic - log2_fit.bic)")
            println(io, "$(prefix)_delta_bic_bootstrap_low=$(delta_interval.low)")
            println(io, "$(prefix)_delta_bic_bootstrap_high=$(delta_interval.high)")
            println(
                io,
                "$(prefix)_bootstrap_fraction_favoring_log2=" *
                "$(mean(bootstrap.delta_bic .> 0))",
            )
            println(io, "$(prefix)_power_prefactor=$(exponent.prefactor)")
            println(io, "$(prefix)_power_exponent=$(exponent.exponent)")
            println(io, "$(prefix)_power_exponent_bootstrap_low=$(exponent_interval.low)")
            println(io, "$(prefix)_power_exponent_bootstrap_high=$(exponent_interval.high)")
            println(io, "$(prefix)_positive_exponent_bootstrap_reps=$(length(bootstrap.exponent))")
        end
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    if isnothing(parsed)
        print_help()
        return
    end
    parsed.bootstrap_reps > 0 || error("bootstrap repetitions must be positive")
    pairs = load_pairs(parsed.paired_paths)
    singles = load_single(parsed.single_paths)
    orders = sort(collect(keys(pairs)))
    missing_single = setdiff(orders, collect(keys(singles)))
    isempty(missing_single) || error("single-height data missing orders: $missing_single")
    all(length(pairs[order]) >= 2 && length(singles[order]) >= 2 for order in orders) ||
        error("each order needs at least two paired and single observations")
    fit_orders = filter(>=(parsed.min_order), orders)
    length(fit_orders) >= 3 || error("at least three fitted orders are required")

    rng = Xoshiro(parsed.bootstrap_seed)
    draws = bootstrap_statistics(
        rng,
        pairs,
        singles,
        orders,
        fit_orders,
        parsed.bootstrap_reps,
    )
    isempty(draws.exponent) &&
        error("no bootstrap replicate had positive difference variances")
    fit_variances = [
        var([value[3] for value in pairs[order]]; corrected=true)
        for order in fit_orders
    ]
    all(>(0), fit_variances) || error("all fitted difference variances must be positive")

    mkpath(parsed.output_dir)
    summary_path = joinpath(parsed.output_dir, "double_dimer_summary.csv")
    fits_path = joinpath(parsed.output_dir, "double_dimer_fits.txt")
    curves_path = joinpath(parsed.output_dir, "double_dimer_fit_curves.csv")
    component_fits_path = joinpath(parsed.output_dir, "variance_component_fits.txt")
    write_summary(summary_path, pairs, singles, orders, draws)
    log_fit, log2_fit, exponent = write_fits(fits_path, fit_orders, fit_variances, draws)
    write_curves(curves_path, fit_orders, log_fit, log2_fit, exponent)
    write_component_fits(component_fits_path, pairs, singles, fit_orders, draws)

    exponent_interval = interval(draws.exponent)
    println("Double-dimer analysis complete")
    println("  pairs:          $(sum(length(group) for group in Base.values(pairs)))")
    println("  orders:         $(join(orders, ", "))")
    @printf(
        "  exponent p:     %.4f (95%% bootstrap %.4f to %.4f)\n",
        exponent.exponent,
        exponent_interval.low,
        exponent_interval.high,
    )
    @printf(
        "  delta BIC:      %.3f (positive favors log^2)\n",
        log_fit.bic - log2_fit.bic,
    )
    @printf("  log^2 wins:     %.1f%% of bootstraps\n", 100 * mean(draws.delta_bic .> 0))
    println("  output:         $(parsed.output_dir)")
end

main(ARGS)
