function sample_variance(values)::Float64
    length(values) >= 2 || return NaN
    return var(Float64.(values); corrected=true)
end

function csv_rows(path::AbstractString)
    return NamedTuple.(CSV.File(path; types=String))
end

function read_config(path::AbstractString)::Vector{BatchConfig}
    ids = [parse(Int, row.task_id) for row in CSV.File(path; types=String)]
    return [load_config_row(path, id) for id in ids]
end

function validate_and_collect(tasks, results_dir::AbstractString; allow_incomplete=false)
    raw = NamedTuple[]
    validation = NamedTuple[]
    for task in tasks
        path = result_path(results_dir, task)
        expected = expected_observations(task)
        rows = isfile(path) ? csv_rows(path) : NamedTuple[]
        ok_rows = [row for row in rows if row.status == "ok"]
        append!(raw, ok_rows)
        complete = isfile(path) && completed_result(path, task)
        push!(validation, (
            task_id=task.task_id, environment_model=task.environment_model,
            distribution=task.distribution,
            distribution_params=params_json(task), L=task.L, batch_id=task.batch_id,
            path=path, exists=isfile(path), expected_rows=expected,
            row_count=length(rows), ok_count=length(ok_rows),
            failed_count=length(rows) - length(ok_rows), complete=complete,
        ))
    end
    incomplete = [row.task_id for row in validation if !row.complete]
    !allow_incomplete && !isempty(incomplete) &&
        throw(ErrorException("incomplete or missing tasks: $(join(incomplete[1:min(end, 10)], ", "))"))
    return raw, validation
end

function assert_no_duplicates(raw)
    seen = Set{Tuple{String,String,String,String,String,String}}()
    for row in raw
        key = (row.distribution, row.distribution_params, row.L, row.batch_id,
               row.environment_id, row.walk_id)
        key in seen && throw(ErrorException("duplicate observation key: $key"))
        push!(seen, key)
    end
end

function assert_direction_counts(raw)
    for row in raw
        hasproperty(row, :north_steps) || continue # batch_v6 compatibility
        counts = parse.(Int, (row.north_steps, row.east_steps,
                              row.south_steps, row.west_steps))
        sum(counts) == parse(Int, row.raw_walk_length) ||
            throw(ErrorException(
                "raw direction counts do not sum to raw_walk_length for task " *
                "$(row.task_id), environment $(row.environment_id)"))
    end
    return nothing
end

function assert_strict_annealed_rows(raw)
    all(row -> row.walk_id == "0", raw) ||
        throw(ErrorException("strict annealed results must have walk_id=0 for every row"))
    seeds = getproperty.(raw, :environment_seed)
    length(seeds) == length(unique(seeds)) ||
        throw(ErrorException("strict annealed results reuse at least one environment seed"))
    return nothing
end

function assert_temporal_rows(raw)
    all(row -> row.environment_model == "temporal_iid", raw) ||
        throw(ErrorException("temporal analysis received non-temporal result rows"))
    all(row -> row.walk_id == "0", raw) ||
        throw(ErrorException("temporal results must have walk_id=0 for every row"))
    environment_seeds = getproperty.(raw, :environment_seed)
    length(environment_seeds) == length(unique(environment_seeds)) ||
        throw(ErrorException("temporal results reuse at least one weight-stream seed"))
    direction_seeds = getproperty.(raw, :direction_seed)
    length(direction_seeds) == length(unique(direction_seeds)) ||
        throw(ErrorException("temporal results reuse at least one direction-stream seed"))
    assert_direction_counts(raw)
    return nothing
end

function double_dimer_key(row)
    return (row.distribution, row.distribution_params, row.L, row.task_id,
            row.batch_id, row.environment_id)
end

"Validate two conditionally independent walks and collapse each pair to W1-W2."
function make_double_dimer_pairs(raw)
    groups = Dict{NTuple{6,String},Vector{Any}}()
    for row in raw
        push!(get!(groups, double_dimer_key(row), Any[]), row)
    end

    environment_seeds = String[]
    walk_seeds = String[]
    pairs = NamedTuple[]
    for (_, rows) in sort!(collect(groups); by=first)
        length(rows) == 2 || throw(ErrorException(
            "double-dimer environment must contain exactly two successful walks: " *
            string(double_dimer_key(first(rows)))))
        sort!(rows; by=row -> parse(Int, row.walk_id))
        getproperty.(rows, :walk_id) == ["0", "1"] || throw(ErrorException(
            "double-dimer walk IDs must be 0 and 1: " * string(double_dimer_key(first(rows)))))
        rows[1].environment_seed == rows[2].environment_seed || throw(ErrorException(
            "double-dimer pair does not share one environment seed"))
        rows[1].walk_seed != rows[2].walk_seed || throw(ErrorException(
            "double-dimer pair reuses a walk seed"))

        first_row, second_row = rows
        winding_1 = parse(Int, first_row.winding)
        winding_2 = parse(Int, second_row.winding)
        push!(environment_seeds, first_row.environment_seed)
        append!(walk_seeds, (first_row.walk_seed, second_row.walk_seed))
        push!(pairs, (
            distribution=first_row.distribution,
            distribution_params=first_row.distribution_params,
            L=first_row.L,
            task_id=first_row.task_id,
            batch_id=first_row.batch_id,
            environment_id=first_row.environment_id,
            environment_seed=first_row.environment_seed,
            walk_seed_1=first_row.walk_seed,
            walk_seed_2=second_row.walk_seed,
            winding_1=winding_1,
            winding_2=winding_2,
            winding_difference=winding_1 - winding_2,
            # The generic fitting/bootstrap code treats this as its observable.
            winding=string(winding_1 - winding_2),
            loop_erased_path_length_1=parse(Int, first_row.loop_erased_path_length),
            loop_erased_path_length_2=parse(Int, second_row.loop_erased_path_length),
            raw_walk_length_1=parse(Int, first_row.raw_walk_length),
            raw_walk_length_2=parse(Int, second_row.raw_walk_length),
            exit_x_1=parse(Int, first_row.exit_x), exit_y_1=parse(Int, first_row.exit_y),
            exit_x_2=parse(Int, second_row.exit_x), exit_y_2=parse(Int, second_row.exit_y),
            runtime_1=parse(Float64, first_row.runtime),
            runtime_2=parse(Float64, second_row.runtime),
        ))
    end
    length(environment_seeds) == length(unique(environment_seeds)) ||
        throw(ErrorException("double-dimer results reuse an environment across pairs"))
    length(walk_seeds) == length(unique(walk_seeds)) ||
        throw(ErrorException("double-dimer results reuse at least one walk seed"))
    return pairs
end

function group_observations(raw)
    groups = Dict{Tuple{String,String,Int},Vector{Any}}()
    for row in raw
        key = (row.distribution, row.distribution_params, parse(Int, row.L))
        push!(get!(groups, key, Any[]), row)
    end
    return groups
end

function group_environments(rows)
    groups = Dict{String,Vector{Float64}}()
    for row in rows
        push!(get!(groups, row.environment_id, Float64[]), parse(Float64, row.winding))
    end
    return groups
end

"Cluster-jackknife SE for the annealed variance, leaving out one environment."
function cluster_variance_se(by_environment)::Float64
    items = collect(values(by_environment))
    count = length(items)
    count >= 2 || return NaN
    if all(values -> length(values) == 1, items)
        count >= 3 || return NaN
        observations = first.(items)
        total = sum(observations)
        total_squares = sum(abs2, observations)
        remaining = count - 1
        estimates = [
            ((total_squares - value^2) - (total - value)^2 / remaining) /
            (remaining - 1)
            for value in observations
        ]
        centre = mean(estimates)
        return sqrt((count - 1) / count *
                    sum((estimate - centre)^2 for estimate in estimates))
    end
    estimates = Float64[]
    for omitted in eachindex(items)
        values = Float64[]
        for index in eachindex(items)
            index == omitted || append!(values, items[index])
        end
        push!(estimates, sample_variance(values))
    end
    centre = mean(estimates)
    return sqrt((count - 1) / count * sum((value - centre)^2 for value in estimates))
end

function aggregate_direction_counts(rows)
    isempty(rows) && return nothing
    hasproperty(first(rows), :north_steps) || return nothing # batch_v6 compatibility
    totals = ntuple(index -> sum(parse(Int, getproperty(row,
        (:north_steps, :east_steps, :south_steps, :west_steps)[index]))
        for row in rows), 4)
    return totals
end

function make_summary(raw; annealed_only::Bool=false)
    summary = NamedTuple[]
    for ((distribution, params, L), rows) in sort!(collect(group_observations(raw)); by=first)
        windings = parse.(Float64, getproperty.(rows, :winding))
        by_environment = group_environments(rows)
        env_variances = [sample_variance(values) for values in values(by_environment)
                         if length(values) >= 2]
        env_means = [mean(values) for values in values(by_environment)]
        raw_lengths = parse.(Float64, getproperty.(rows, :raw_walk_length))
        path_lengths = parse.(Float64, getproperty.(rows, :loop_erased_path_length))
        runtimes = parse.(Float64, getproperty.(rows, :runtime))
        exit_x = parse.(Float64, getproperty.(rows, :exit_x)) ./ L
        exit_y = parse.(Float64, getproperty.(rows, :exit_y)) ./ L
        counts = length.(values(by_environment))
        direction_totals = aggregate_direction_counts(rows)
        total_raw_steps = direction_totals === nothing ? missing : sum(direction_totals)
        direction_frequencies = direction_totals === nothing ?
            (missing, missing, missing, missing) :
            Tuple(count / total_raw_steps for count in direction_totals)
        quenched = isempty(env_variances) ? NaN : mean(env_variances)
        quenched_se = length(env_variances) >= 2 ? std(env_variances) / sqrt(length(env_variances)) : NaN
        common = (
            environment_model=first(rows).environment_model,
            distribution=distribution, distribution_params=params, L=L,
            log_L=log(L), log_log_L=log(log(L)), observations=length(rows),
            environments=length(by_environment), walks_per_environment_min=minimum(counts),
            walks_per_environment_max=maximum(counts), mean_winding=mean(windings),
            annealed_variance=sample_variance(windings),
            annealed_variance_se=cluster_variance_se(by_environment),
            mean_exit_x_over_L=mean(exit_x), mean_exit_y_over_L=mean(exit_y),
            mean_raw_walk_length=mean(raw_lengths),
            mean_loop_erased_path_length=mean(path_lengths), mean_runtime=mean(runtimes),
            total_raw_steps=total_raw_steps,
            total_north_steps=direction_totals === nothing ? missing : direction_totals[1],
            total_east_steps=direction_totals === nothing ? missing : direction_totals[2],
            total_south_steps=direction_totals === nothing ? missing : direction_totals[3],
            total_west_steps=direction_totals === nothing ? missing : direction_totals[4],
            north_step_frequency=direction_frequencies[1],
            east_step_frequency=direction_frequencies[2],
            south_step_frequency=direction_frequencies[3],
            west_step_frequency=direction_frequencies[4],
        )
        if annealed_only
            push!(summary, common)
        else
            push!(summary, merge(common, (
                quenched_variance=quenched, quenched_variance_se=quenched_se,
                environment_mean_variance=sample_variance(env_means),
            )))
        end
    end
    return summary
end

function resample_rows(rng::AbstractRNG, rows)
    return [rows[rand(rng, eachindex(rows))] for _ in eachindex(rows)]
end

function temporal_comparison_statistics(rows, baseline_rows)
    windings = parse.(Float64, getproperty.(rows, :winding))
    baseline_windings = parse.(Float64, getproperty.(baseline_rows, :winding))
    counts = aggregate_direction_counts(rows)
    baseline_counts = aggregate_direction_counts(baseline_rows)
    frequencies = Tuple(count / sum(counts) for count in counts)
    baseline_frequencies = Tuple(count / sum(baseline_counts) for count in baseline_counts)
    direction_differences = frequencies .- baseline_frequencies
    return (
        variance_ratio=sample_variance(windings) / sample_variance(baseline_windings),
        mean_winding_difference=mean(windings) - mean(baseline_windings),
        direction_differences=direction_differences,
    )
end

function temporal_baseline_comparisons(raw; bootstrap_reps::Int=0,
                                       bootstrap_seed::UInt64=0x000000000134d5ef)
    groups = group_observations(raw)
    rng = StableRNG(bootstrap_seed)
    output = NamedTuple[]
    for ((distribution, params, L), rows) in sort!(collect(groups); by=first)
        distribution == "baseline" && continue
        baseline_key = ("baseline", "{}", L)
        haskey(groups, baseline_key) || throw(ErrorException(
            "temporal baseline is missing at L=$L"))
        baseline_rows = groups[baseline_key]
        estimate = temporal_comparison_statistics(rows, baseline_rows)
        variance_ratios = Float64[]
        winding_differences = Float64[]
        direction_differences = [Float64[] for _ in 1:4]
        for _ in 1:bootstrap_reps
            candidate = temporal_comparison_statistics(
                resample_rows(rng, rows), resample_rows(rng, baseline_rows))
            isfinite(candidate.variance_ratio) &&
                push!(variance_ratios, candidate.variance_ratio)
            push!(winding_differences, candidate.mean_winding_difference)
            for direction in 1:4
                push!(direction_differences[direction],
                      candidate.direction_differences[direction])
            end
        end
        interval(values) = isempty(values) ? (NaN, NaN) :
            (quantile(values, 0.025), quantile(values, 0.975))
        ratio_ci = interval(variance_ratios)
        winding_ci = interval(winding_differences)
        direction_cis = interval.(direction_differences)
        separated = ratio_ci[1] > 1 || ratio_ci[2] < 1 ||
            winding_ci[1] > 0 || winding_ci[2] < 0 ||
            any(ci -> ci[1] > 0 || ci[2] < 0, direction_cis)
        differences = estimate.direction_differences
        push!(output, (
            environment_model="temporal_iid",
            distribution=distribution, distribution_params=params, L=L,
            observations=length(rows), baseline_observations=length(baseline_rows),
            variance_ratio=estimate.variance_ratio,
            variance_ratio_ci_low=ratio_ci[1], variance_ratio_ci_high=ratio_ci[2],
            mean_winding_difference=estimate.mean_winding_difference,
            mean_winding_difference_ci_low=winding_ci[1],
            mean_winding_difference_ci_high=winding_ci[2],
            north_frequency_difference=differences[1],
            north_frequency_difference_ci_low=direction_cis[1][1],
            north_frequency_difference_ci_high=direction_cis[1][2],
            east_frequency_difference=differences[2],
            east_frequency_difference_ci_low=direction_cis[2][1],
            east_frequency_difference_ci_high=direction_cis[2][2],
            south_frequency_difference=differences[3],
            south_frequency_difference_ci_low=direction_cis[3][1],
            south_frequency_difference_ci_high=direction_cis[3][2],
            west_frequency_difference=differences[4],
            west_frequency_difference_ci_low=direction_cis[4][1],
            west_frequency_difference_ci_high=direction_cis[4][2],
            any_unadjusted_95pct_ci_excludes_null=separated,
        ))
    end
    return output
end

function temporal_direction_diagnostics(raw; bootstrap_reps::Int=0,
                                        bootstrap_seed::UInt64=0x00000000013504f7)
    rng = StableRNG(bootstrap_seed)
    output = NamedTuple[]
    for ((distribution, params, L), rows) in
            sort!(collect(group_observations(raw)); by=first)
        totals = aggregate_direction_counts(rows)
        frequencies = Tuple(count / sum(totals) for count in totals)
        bootstrap_frequencies = [Float64[] for _ in 1:4]
        for _ in 1:bootstrap_reps
            candidate_totals = aggregate_direction_counts(resample_rows(rng, rows))
            candidate_total = sum(candidate_totals)
            for direction in 1:4
                push!(bootstrap_frequencies[direction],
                      candidate_totals[direction] / candidate_total)
            end
        end
        interval(values) = isempty(values) ? (NaN, NaN) :
            (quantile(values, 0.025), quantile(values, 0.975))
        intervals = interval.(bootstrap_frequencies)
        consistent = all(ci -> ci[1] <= 0.25 <= ci[2], intervals)
        push!(output, (
            environment_model="temporal_iid",
            distribution=distribution, distribution_params=params, L=L,
            observations=length(rows), total_raw_steps=sum(totals),
            north_step_frequency=frequencies[1],
            north_step_frequency_ci_low=intervals[1][1],
            north_step_frequency_ci_high=intervals[1][2],
            east_step_frequency=frequencies[2],
            east_step_frequency_ci_low=intervals[2][1],
            east_step_frequency_ci_high=intervals[2][2],
            south_step_frequency=frequencies[3],
            south_step_frequency_ci_low=intervals[3][1],
            south_step_frequency_ci_high=intervals[3][2],
            west_step_frequency=frequencies[4],
            west_step_frequency_ci_low=intervals[4][1],
            west_step_frequency_ci_high=intervals[4][2],
            max_abs_frequency_deviation=maximum(abs.(frequencies .- 0.25)),
            all_unadjusted_95pct_cis_include_one_quarter=consistent,
        ))
    end
    return output
end

function make_double_dimer_summary(pairs)
    summary = NamedTuple[]
    for ((distribution, params, L), rows) in sort!(collect(group_observations(pairs)); by=first)
        winding_1 = Float64.(getproperty.(rows, :winding_1))
        winding_2 = Float64.(getproperty.(rows, :winding_2))
        differences = Float64.(getproperty.(rows, :winding_difference))
        variance_1 = sample_variance(winding_1)
        variance_2 = sample_variance(winding_2)
        covariance = length(rows) >= 2 ? cov(winding_1, winding_2; corrected=true) : NaN
        difference_variance = sample_variance(differences)
        by_environment = group_environments(rows)
        pooled = vcat(winding_1, winding_2)
        identity_rhs = variance_1 + variance_2 - 2covariance
        exit_x = vcat(Float64.(getproperty.(rows, :exit_x_1)),
                      Float64.(getproperty.(rows, :exit_x_2))) ./ L
        exit_y = vcat(Float64.(getproperty.(rows, :exit_y_1)),
                      Float64.(getproperty.(rows, :exit_y_2))) ./ L
        raw_lengths = vcat(Float64.(getproperty.(rows, :raw_walk_length_1)),
                           Float64.(getproperty.(rows, :raw_walk_length_2)))
        path_lengths = vcat(Float64.(getproperty.(rows, :loop_erased_path_length_1)),
                            Float64.(getproperty.(rows, :loop_erased_path_length_2)))
        runtimes = vcat(Float64.(getproperty.(rows, :runtime_1)),
                        Float64.(getproperty.(rows, :runtime_2)))
        push!(summary, (
            distribution=distribution, distribution_params=params, L=L,
            log_L=log(L), log_log_L=log(log(L)), pairs=length(rows),
            observations=length(rows), walks=2length(rows), environments=length(rows),
            walks_per_environment_min=2, walks_per_environment_max=2,
            mean_winding_1=mean(winding_1), mean_winding_2=mean(winding_2),
            mean_winding_difference=mean(differences),
            variance_winding_1=variance_1, variance_winding_2=variance_2,
            covariance_winding_1_winding_2=covariance,
            correlation_winding_1_winding_2=cor(winding_1, winding_2),
            pooled_single_winding_variance=sample_variance(pooled),
            double_dimer_variance=difference_variance,
            double_dimer_variance_se=cluster_variance_se(by_environment),
            variance_ratio_to_twice_pooled=difference_variance / (2sample_variance(pooled)),
            variance_identity_rhs=identity_rhs,
            variance_identity_residual=difference_variance - identity_rhs,
            # Compatibility aliases used by the common scaling analysis.
            mean_winding=mean(differences), annealed_variance=difference_variance,
            annealed_variance_se=cluster_variance_se(by_environment),
            mean_exit_x_over_L=mean(exit_x), mean_exit_y_over_L=mean(exit_y),
            mean_raw_walk_length=mean(raw_lengths),
            mean_loop_erased_path_length=mean(path_lengths), mean_runtime=mean(runtimes),
        ))
    end
    return summary
end

function linear_fit(xs, ys)
    n = length(xs)
    n >= 2 || return (intercept=NaN, slope=NaN, slope_se=NaN, sse=NaN, r2=NaN, aic=NaN, bic=NaN)
    xbar, ybar = mean(xs), mean(ys)
    sxx = sum((x - xbar)^2 for x in xs)
    sxx > 0 || return (intercept=NaN, slope=NaN, slope_se=NaN, sse=NaN, r2=NaN, aic=NaN, bic=NaN)
    slope = sum((x - xbar) * (y - ybar) for (x, y) in zip(xs, ys)) / sxx
    intercept = ybar - slope * xbar
    residuals = ys .- (intercept .+ slope .* xs)
    sse = sum(abs2, residuals)
    total = sum((y - ybar)^2 for y in ys)
    r2 = total == 0 ? NaN : 1 - sse / total
    slope_se = n > 2 ? sqrt((sse / (n - 2)) / sxx) : NaN
    safe_sse = max(sse, 1e-300)
    aic = n * log(safe_sse / n) + 4
    bic = n * log(safe_sse / n) + 2 * log(n)
    return (; intercept, slope, slope_se, sse, r2, aic, bic)
end

function grouped_summary(summary; min_L=nothing)
    groups = Dict{Tuple{String,String},Vector{Any}}()
    for row in summary
        min_L !== nothing && row.L < min_L && continue
        push!(get!(groups, (row.distribution, row.distribution_params), Any[]), row)
    end
    return groups
end

function fit_exponents(summary, variance_kind::Symbol; min_L=nothing)
    column = variance_kind === :annealed ? :annealed_variance : :quenched_variance
    fits = NamedTuple[]
    for ((distribution, params), rows) in sort!(collect(grouped_summary(summary; min_L)); by=first)
        ordered = sort(rows; by=row -> row.L)
        ordered = [row for row in ordered if isfinite(getproperty(row, column)) && getproperty(row, column) > 0]
        xs = getproperty.(ordered, :log_log_L)
        ys = log.(getproperty.(ordered, column))
        fit = linear_fit(xs, ys)
        push!(fits, (
            distribution=distribution, distribution_params=params,
            variance_kind=String(variance_kind), fit_min_L=min_L,
            n_sizes=length(ordered), L_values=join(getproperty.(ordered, :L), " "),
            intercept_log_C=fit.intercept, p=fit.slope, p_se=fit.slope_se,
            p_ci_low=fit.slope - 1.96 * fit.slope_se,
            p_ci_high=fit.slope + 1.96 * fit.slope_se,
            p_bootstrap_ci_low=NaN, p_bootstrap_ci_high=NaN,
            r2=fit.r2, sse=fit.sse, aic=fit.aic, bic=fit.bic,
        ))
    end
    return fits
end

"Fit E[LERW length] = C L^d on log-log axes."
function fit_path_length_exponents(summary; min_L=nothing)
    fits = NamedTuple[]
    for ((distribution, params), rows) in
            sort!(collect(grouped_summary(summary; min_L)); by=first)
        ordered = sort(rows; by=row -> row.L)
        ordered = [row for row in ordered
                   if isfinite(row.mean_loop_erased_path_length) &&
                      row.mean_loop_erased_path_length > 0]
        xs = Float64.(getproperty.(ordered, :log_L))
        ys = log.(Float64.(getproperty.(ordered, :mean_loop_erased_path_length)))
        fit = linear_fit(xs, ys)
        push!(fits, (
            distribution=distribution, distribution_params=params,
            fit_min_L=min_L, n_sizes=length(ordered),
            L_values=join(getproperty.(ordered, :L), " "),
            intercept_log_C=fit.intercept, length_exponent=fit.slope,
            length_exponent_se=fit.slope_se,
            length_exponent_ci_low=fit.slope - 1.96 * fit.slope_se,
            length_exponent_ci_high=fit.slope + 1.96 * fit.slope_se,
            length_exponent_bootstrap_ci_low=NaN,
            length_exponent_bootstrap_ci_high=NaN,
            r2=fit.r2, sse=fit.sse, aic=fit.aic, bic=fit.bic,
        ))
    end
    return fits
end

function group_path_lengths_by_environment(rows)
    groups = Dict{String,Vector{Float64}}()
    for row in rows
        push!(get!(groups, row.environment_id, Float64[]),
              parse(Float64, row.loop_erased_path_length))
    end
    return groups
end

function bootstrap_path_length_intervals!(
    fits, raw, summary, reps::Int, seed::UInt64,
)
    reps > 1 || return fits
    rng = StableRNG(seed)
    raw_groups = group_observations(raw)
    for index in eachindex(fits)
        fit_row = fits[index]
        model_summary = [row for row in summary
                         if row.distribution == fit_row.distribution &&
                            row.distribution_params == fit_row.distribution_params &&
                            (fit_row.fit_min_L === nothing ||
                             row.L >= fit_row.fit_min_L)]
        length(model_summary) >= 3 || continue
        exponents = Float64[]
        for _ in 1:reps
            pseudo = NamedTuple[]
            for row in model_summary
                observations = raw_groups[
                    (row.distribution, row.distribution_params, row.L)]
                environments = collect(values(
                    group_path_lengths_by_environment(observations)))
                sampled = Float64[]
                for _ in eachindex(environments)
                    append!(sampled, environments[rand(rng, eachindex(environments))])
                end
                push!(pseudo, merge(row, (
                    mean_loop_erased_path_length=mean(sampled),)))
            end
            candidate = fit_path_length_exponents(pseudo)
            if !isempty(candidate) && isfinite(candidate[1].length_exponent)
                push!(exponents, candidate[1].length_exponent)
            end
        end
        if !isempty(exponents)
            fits[index] = merge(fit_row, (
                length_exponent_bootstrap_ci_low=quantile(exponents, 0.025),
                length_exponent_bootstrap_ci_high=quantile(exponents, 0.975),
            ))
        end
    end
    return fits
end

function path_length_pointwise_rows(summary)
    output = NamedTuple[]
    for row in summary
        mean_length = row.mean_loop_erased_path_length
        push!(output, (
            distribution=row.distribution,
            distribution_params=row.distribution_params,
            L=row.L, observations=row.observations,
            mean_loop_erased_path_length=mean_length,
            log_mean_length_over_log_L=log(mean_length) / row.log_L,
            mean_length_over_L_to_5_over_4=mean_length / row.L^1.25,
        ))
    end
    return output
end

function path_length_local_exponent_rows(summary)
    output = NamedTuple[]
    for ((distribution, params), rows) in
            sort!(collect(grouped_summary(summary)); by=first)
        ordered = sort(rows; by=row -> row.L)
        for (left, right) in zip(ordered[1:end-1], ordered[2:end])
            exponent = (
                log(right.mean_loop_erased_path_length) -
                log(left.mean_loop_erased_path_length)
            ) / (right.log_L - left.log_L)
            push!(output, (
                distribution=distribution, distribution_params=params,
                L_left=left.L, L_right=right.L,
                L_mid_geometric=sqrt(left.L * right.L),
                length_exponent_local=exponent,
            ))
        end
    end
    return output
end

"Compare additive `a+b log L` and `a+b(log L)^2` finite-size models."
function fit_scaling_models(summary)
    return fit_scaling_models(summary, (:annealed, :quenched))
end

function fit_scaling_models(summary, variance_kinds)
    rows_out = NamedTuple[]
    for variance_kind in variance_kinds
        column = variance_kind === :annealed ? :annealed_variance : :quenched_variance
        for ((distribution, params), rows) in sort!(collect(grouped_summary(summary)); by=first)
            ordered = sort(rows; by=row -> row.L)
            ys = Float64[getproperty(row, column) for row in ordered]
            for (model, power) in (("a_plus_b_log_L", 1), ("a_plus_b_log_L_squared", 2))
                xs = [row.log_L^power for row in ordered]
                fit = linear_fit(xs, ys)
                push!(rows_out, (
                    distribution=distribution, distribution_params=params,
                    variance_kind=String(variance_kind), model=model,
                    intercept=fit.intercept, coefficient=fit.slope,
                    n_sizes=length(ordered), sse=fit.sse, r2=fit.r2,
                    aic=fit.aic, bic=fit.bic,
                ))
            end
        end
    end
    return rows_out
end

function bootstrap_intervals!(fits, raw, summary, variance_kind::Symbol, reps::Int, seed::UInt64)
    reps > 1 || return fits
    rng = StableRNG(seed)
    raw_groups = group_observations(raw)
    for index in eachindex(fits)
        fit_row = fits[index]
        slopes = Float64[]
        model_summary = [row for row in summary if row.distribution == fit_row.distribution &&
                         row.distribution_params == fit_row.distribution_params &&
                         (fit_row.fit_min_L === nothing || row.L >= fit_row.fit_min_L)]
        length(model_summary) >= 3 || continue
        for _ in 1:reps
            pseudo = NamedTuple[]
            for row in model_summary
                observations = raw_groups[(row.distribution, row.distribution_params, row.L)]
                environments = collect(values(group_environments(observations)))
                sampled = Float64[]
                sampled_variances = Float64[]
                for _ in eachindex(environments)
                    values = environments[rand(rng, eachindex(environments))]
                    append!(sampled, values) # cluster bootstrap: preserve walks within environment
                    length(values) >= 2 && push!(sampled_variances, sample_variance(values))
                end
                value = variance_kind === :annealed ? sample_variance(sampled) : mean(sampled_variances)
                replacement = variance_kind === :annealed ?
                    (; annealed_variance=value) : (; quenched_variance=value)
                push!(pseudo, merge(row, replacement))
            end
            candidate = fit_exponents(pseudo, variance_kind)
            !isempty(candidate) && isfinite(candidate[1].p) && push!(slopes, candidate[1].p)
        end
        if !isempty(slopes)
            sort!(slopes)
            low = quantile(slopes, 0.025)
            high = quantile(slopes, 0.975)
            fits[index] = merge(fit_row, (; p_bootstrap_ci_low=low, p_bootstrap_ci_high=high))
        end
    end
    return fits
end

function pointwise_rows(summary)
    return pointwise_rows(summary, (:annealed, :quenched))
end

function pointwise_rows(summary, variance_kinds)
    output = NamedTuple[]
    for row in summary, kind in variance_kinds
        value = kind === :annealed ? row.annealed_variance : row.quenched_variance
        push!(output, (
            distribution=row.distribution, distribution_params=row.distribution_params,
            variance_kind=String(kind), L=row.L, variance=value,
            log_variance=value > 0 ? log(value) : NaN, log_log_L=row.log_log_L,
            log_variance_over_log_log_L=value > 0 ? log(value) / row.log_log_L : NaN,
        ))
    end
    return output
end

function local_exponent_rows(summary)
    return local_exponent_rows(summary, (:annealed, :quenched))
end

function local_exponent_rows(summary, variance_kinds)
    output = NamedTuple[]
    for kind in variance_kinds
        column = kind === :annealed ? :annealed_variance : :quenched_variance
        for ((distribution, params), rows) in sort!(collect(grouped_summary(summary)); by=first)
            ordered = sort(rows; by=row -> row.L)
            for (left, right) in zip(ordered[1:end-1], ordered[2:end])
                a, b = getproperty(left, column), getproperty(right, column)
                denominator = right.log_log_L - left.log_log_L
                p = a > 0 && b > 0 ? (log(b) - log(a)) / denominator : NaN
                push!(output, (
                    distribution=distribution, distribution_params=params,
                    variance_kind=String(kind), L_left=left.L, L_right=right.L,
                    L_mid_geometric=sqrt(left.L * right.L), p_local=p,
                ))
            end
        end
    end
    return output
end

function write_table(path::AbstractString, rows)
    mkpath(dirname(path))
    if isempty(rows)
        open(path, "w") do _ end
    else
        columns = Symbol[]
        seen = Set{Symbol}()
        for row in rows, column in keys(row)
            if column ∉ seen
                push!(columns, column)
                push!(seen, column)
            end
        end
        names = Tuple(columns)
        normalized = [NamedTuple{names}(ntuple(index -> begin
                          column = names[index]
                          hasproperty(row, column) ? getproperty(row, column) : missing
                      end, length(names))) for row in rows]
        CSV.write(path, normalized;
                  transform=(_column, value) -> something(value, missing))
    end
end

function analyze_results(config::AbstractString, results_dir::AbstractString,
                         output_dir::AbstractString; fit_min_L=nothing,
                         bootstrap_reps::Int=0, bootstrap_seed::Integer=20260623,
                         allow_incomplete::Bool=false, strict_annealed::Bool=false,
                         double_dimer::Bool=false)
    tasks = read_config(config)
    strict_annealed && double_dimer && throw(ArgumentError(
        "--strict-annealed and --double-dimer are mutually exclusive"))
    strict_annealed && require_strict_annealed(tasks)
    double_dimer && require_double_dimer(tasks)
    environment_models = unique(getproperty.(tasks, :environment_model))
    length(environment_models) == 1 || throw(ArgumentError(
        "analyze one environment model at a time; found: " *
        join(sort(environment_models), ", ")))
    temporal = only(environment_models) == "temporal_iid"
    temporal && require_temporal_iid(tasks)
    temporal && double_dimer && throw(ArgumentError(
        "temporal_iid observations are independent annealed walks, not double-dimer pairs"))
    raw, validation = validate_and_collect(tasks, results_dir; allow_incomplete)
    assert_no_duplicates(raw)
    assert_direction_counts(raw)
    strict_annealed && assert_strict_annealed_rows(raw)
    temporal && assert_temporal_rows(raw)
    pairs = double_dimer ? make_double_dimer_pairs(raw) : NamedTuple[]
    analysis_rows = double_dimer ? pairs : raw
    summary = double_dimer ? make_double_dimer_summary(pairs) :
              make_summary(raw; annealed_only=temporal)
    path_length_fits = double_dimer ? NamedTuple[] :
        fit_path_length_exponents(summary; min_L=fit_min_L)
    !double_dimer && bootstrap_path_length_intervals!(
        path_length_fits, raw, summary, bootstrap_reps,
        UInt64(bootstrap_seed + 15485863))
    annealed = fit_exponents(summary, :annealed; min_L=fit_min_L)
    bootstrap_intervals!(annealed, analysis_rows, summary, :annealed,
                         bootstrap_reps, UInt64(bootstrap_seed))
    variance_kinds = (strict_annealed || double_dimer || temporal) ?
        (:annealed,) : (:annealed, :quenched)
    fits = if strict_annealed || double_dimer || temporal
        annealed
    else
        quenched = fit_exponents(summary, :quenched; min_L=fit_min_L)
        bootstrap_intervals!(quenched, raw, summary, :quenched, bootstrap_reps,
                             UInt64(bootstrap_seed + 7919))
        vcat(annealed, quenched)
    end

    write_table(joinpath(output_dir, "validation.csv"), validation)
    write_table(joinpath(output_dir, "combined_raw.csv"), raw)
    double_dimer && write_table(joinpath(output_dir, "double_dimer_pairs.csv"), pairs)
    write_table(joinpath(output_dir, "summary.csv"), summary)
    write_table(joinpath(output_dir, "loglog_fits.csv"), fits)
    if !double_dimer
        write_table(joinpath(output_dir, "path_length_fits.csv"),
                    path_length_fits)
        write_table(joinpath(output_dir, "path_length_pointwise.csv"),
                    path_length_pointwise_rows(summary))
        write_table(joinpath(output_dir,
                             "path_length_local_effective_exponents.csv"),
                    path_length_local_exponent_rows(summary))
    end
    write_table(joinpath(output_dir, "scaling_model_comparison.csv"),
                fit_scaling_models(summary, variance_kinds))
    write_table(joinpath(output_dir, "pointwise_ratios.csv"),
                pointwise_rows(summary, variance_kinds))
    write_table(joinpath(output_dir, "local_effective_exponents.csv"),
                local_exponent_rows(summary, variance_kinds))
    has_temporal_baseline = temporal &&
        any(row -> row.distribution == "baseline", summary)
    has_temporal_baseline &&
        write_table(joinpath(output_dir, "temporal_baseline_comparisons.csv"),
        temporal_baseline_comparisons(raw; bootstrap_reps,
            bootstrap_seed=UInt64(bootstrap_seed + 104729)))
    temporal && write_table(joinpath(output_dir, "temporal_direction_diagnostics.csv"),
        temporal_direction_diagnostics(raw; bootstrap_reps,
            bootstrap_seed=UInt64(bootstrap_seed + 130363)))
    return (; observations=length(analysis_rows), raw_walks=length(raw), tasks=length(tasks),
            summary_rows=length(summary), fits=length(fits))
end

function main_analyze(args=ARGS)::Int
    options = parse_cli(args)
    config = require_option(options, "config")
    results_dir = get(options, "results-dir", "results")
    output_dir = get(options, "output-dir", "analysis")
    fit_min_L = haskey(options, "fit-min-L") ? parse(Int, options["fit-min-L"]) : nothing
    report = analyze_results(config, results_dir, output_dir;
        fit_min_L=fit_min_L,
        bootstrap_reps=parse(Int, get(options, "bootstrap-reps", "0")),
        bootstrap_seed=parse(Int, get(options, "bootstrap-seed", "20260623")),
        allow_incomplete=haskey(options, "allow-incomplete"),
        strict_annealed=haskey(options, "strict-annealed"),
        double_dimer=haskey(options, "double-dimer"))
    unit = haskey(options, "double-dimer") ? "pairs" : "observations"
    println("Analysed $(report.observations) $unit ($(report.raw_walks) raw walks) " *
            "from $(report.tasks) tasks")
    println("Summary rows: $(report.summary_rows); fit rows: $(report.fits)")
    println("Output: $output_dir")
    return 0
end
