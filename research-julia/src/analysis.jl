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
        complete = isfile(path) && completed_result(path, expected)
        push!(validation, (
            task_id=task.task_id, distribution=task.distribution,
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

function make_summary(raw)
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
        quenched = isempty(env_variances) ? NaN : mean(env_variances)
        quenched_se = length(env_variances) >= 2 ? std(env_variances) / sqrt(length(env_variances)) : NaN
        push!(summary, (
            distribution=distribution, distribution_params=params, L=L,
            log_L=log(L), log_log_L=log(log(L)), observations=length(rows),
            environments=length(by_environment), walks_per_environment_min=minimum(counts),
            walks_per_environment_max=maximum(counts), mean_winding=mean(windings),
            annealed_variance=sample_variance(windings),
            annealed_variance_se=cluster_variance_se(by_environment),
            quenched_variance=quenched, quenched_variance_se=quenched_se,
            environment_mean_variance=sample_variance(env_means),
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

"Compare additive `a+b log L` and `a+b(log L)^2` finite-size models."
function fit_scaling_models(summary)
    rows_out = NamedTuple[]
    for variance_kind in (:annealed, :quenched)
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
    output = NamedTuple[]
    for row in summary, kind in (:annealed, :quenched)
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
    output = NamedTuple[]
    for kind in (:annealed, :quenched)
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
                         allow_incomplete::Bool=false)
    tasks = read_config(config)
    raw, validation = validate_and_collect(tasks, results_dir; allow_incomplete)
    assert_no_duplicates(raw)
    summary = make_summary(raw)
    annealed = fit_exponents(summary, :annealed; min_L=fit_min_L)
    quenched = fit_exponents(summary, :quenched; min_L=fit_min_L)
    bootstrap_intervals!(annealed, raw, summary, :annealed, bootstrap_reps, UInt64(bootstrap_seed))
    bootstrap_intervals!(quenched, raw, summary, :quenched, bootstrap_reps,
                         UInt64(bootstrap_seed + 7919))
    fits = vcat(annealed, quenched)

    write_table(joinpath(output_dir, "validation.csv"), validation)
    write_table(joinpath(output_dir, "combined_raw.csv"), raw)
    write_table(joinpath(output_dir, "summary.csv"), summary)
    write_table(joinpath(output_dir, "loglog_fits.csv"), fits)
    write_table(joinpath(output_dir, "scaling_model_comparison.csv"), fit_scaling_models(summary))
    write_table(joinpath(output_dir, "pointwise_ratios.csv"), pointwise_rows(summary))
    write_table(joinpath(output_dir, "local_effective_exponents.csv"), local_exponent_rows(summary))
    return (; observations=length(raw), tasks=length(tasks), summary_rows=length(summary), fits=length(fits))
end

function main_analyze(args=ARGS)::Int
    options = parse_cli(args)
    config = require_option(options, "config")
    results_dir = get(options, "results-dir", "results_hpc")
    output_dir = get(options, "output-dir", "analysis_hpc")
    fit_min_L = haskey(options, "fit-min-L") ? parse(Int, options["fit-min-L"]) : nothing
    report = analyze_results(config, results_dir, output_dir;
        fit_min_L=fit_min_L,
        bootstrap_reps=parse(Int, get(options, "bootstrap-reps", "0")),
        bootstrap_seed=parse(Int, get(options, "bootstrap-seed", "20260623")),
        allow_incomplete=haskey(options, "allow-incomplete"))
    println("Analysed $(report.observations) observations from $(report.tasks) tasks")
    println("Summary rows: $(report.summary_rows); fit rows: $(report.fits)")
    println("Output: $output_dir")
    return 0
end
