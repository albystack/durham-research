const SCHEMA_VERSION = "hpc_batch_v5_julia"
const ENVIRONMENT_MODEL = "site_iid"

struct CTimespec
    seconds::Clong
    nanoseconds::Clong
end

"Monotonic elapsed-time source independent of wall-clock adjustments."
function monotonic_seconds()::Float64
    clock_id = Sys.isapple() ? Cint(6) : Cint(1) # CLOCK_MONOTONIC on macOS/Linux
    result = Ref{CTimespec}()
    status = ccall(:clock_gettime, Cint, (Cint, Ref{CTimespec}), clock_id, result)
    status == 0 || throw(SystemError("clock_gettime", Libc.errno()))
    value = result[]
    return Float64(value.seconds) + Float64(value.nanoseconds) * 1.0e-9
end

struct BatchConfig
    task_id::Int
    distribution::String
    distribution_params::Dict{String,Float64}
    L::Int
    batch_id::Int
    num_environments::Int
    walks_per_environment::Int
    base_seed::UInt64
    environment_model::String
end

params_json(task::BatchConfig) = canonical_json(task.distribution_params)

function load_config_row(config_path::AbstractString, task_id::Integer)::BatchConfig
    for row in CSV.File(config_path; types=String)
        parse(Int, row.task_id) == task_id || continue
        raw_environment_model = hasproperty(row, :environment_model) ? row.environment_model : ""
        environment_model = isempty(raw_environment_model) ? ENVIRONMENT_MODEL : raw_environment_model
        environment_model == ENVIRONMENT_MODEL ||
            throw(ArgumentError("Julia runner supports only site_iid"))
        task = BatchConfig(
            parse(Int, row.task_id), row.distribution, parse_params(row.distribution_params),
            parse(Int, row.L), parse(Int, row.batch_id), parse(Int, row.num_environments),
            parse(Int, row.walks_per_environment), parse(UInt64, row.base_seed),
            environment_model,
        )
        distribution_spec(task.distribution, task.distribution_params)
        return task
    end
    throw(ArgumentError("task_id $task_id not found in $config_path"))
end

function safe_token(value)::String
    text = value isa AbstractFloat ? @sprintf("%g", value) : string(value)
    text = replace(text, "." => "p", "-" => "m")
    return strip(replace(text, r"[^A-Za-z0-9_]+" => "_"), '_')
end

function safe_distribution_name(distribution::String, params::AbstractDict)::String
    isempty(params) && return safe_token(distribution)
    suffix = join((string(safe_token(key), "_", safe_token(params[key]))
                   for key in sort!(collect(keys(params)))), "_")
    return string(safe_token(distribution), "_", suffix)
end

function result_path(output_dir::AbstractString, task::BatchConfig)::String
    return joinpath(output_dir, safe_token(task.environment_model),
                    safe_distribution_name(task.distribution, task.distribution_params),
                    @sprintf("L_%04d", task.L), @sprintf("batch_%04d.csv", task.batch_id))
end

expected_observations(task::BatchConfig) = task.num_environments * task.walks_per_environment

function completed_result(path::AbstractString, expected_rows::Integer)::Bool
    isfile(path) || return false
    rows = collect(CSV.File(path; types=String))
    return length(rows) == expected_rows && all(row -> row.status == "ok", rows)
end

function task_seed(task::BatchConfig, parts...)::UInt64
    return stable_seed(task.base_seed, task.environment_model, task.distribution,
                       params_json(task), task.L, task.batch_id, parts...)
end

function git_version()::String
    project_root = normpath(joinpath(@__DIR__, ".."))
    try
        git_commit = strip(read(`git -C $project_root rev-parse --short HEAD`, String))
        source_paths = [joinpath(@__DIR__, name) for name in
                        ("LERWResearch.jl", "config.jl", "simulation.jl", "batch.jl",
                         "analysis.jl")]
        append!(source_paths, [joinpath(project_root, "Project.toml"),
                               joinpath(project_root, "Manifest.toml")])
        source_text = join((string(path, '\n', read(path, String))
                            for path in source_paths), "\n")
        source_hash = bytes2hex(sha256(source_text))[1:12]
        return string(git_commit, "+source.", source_hash)
    catch
        return "unknown"
    end
end

function result_row(task, code_version, model, parameter, environment_id, walk_id,
                    environment_seed, walk_seed, winding_value, path_length,
                    raw_length, exit_x, exit_y, runtime, sampled_sites, status,
                    error_type, error_message, started_at, finished_at)
    return (
        schema_version=SCHEMA_VERSION,
        code_version=code_version,
        julia_version=string(VERSION),
        task_id=task.task_id,
        batch_id=task.batch_id,
        distribution=task.distribution,
        distribution_params=params_json(task),
        environment_model=task.environment_model,
        model=String(model),
        model_parameter=parameter === nothing ? missing : parameter,
        L=task.L,
        environment_id=environment_id,
        walk_id=walk_id,
        environment_seed=environment_seed,
        walk_seed=walk_seed,
        winding=winding_value,
        loop_erased_path_length=path_length,
        raw_walk_length=raw_length,
        exit_location=exit_x === missing ? missing : string(exit_x, ",", exit_y),
        exit_x=exit_x,
        exit_y=exit_y,
        runtime=runtime,
        sampled_site_count=sampled_sites,
        status=status,
        error_type=error_type,
        error_message=error_message,
        started_at_utc=started_at,
        finished_at_utc=finished_at,
    )
end

function run_batch(task::BatchConfig; max_steps::Union{Nothing,Integer}=nothing)
    model, parameter = distribution_spec(task.distribution, task.distribution_params)
    version = git_version() # warm the command cache before timing observations
    rows_by_environment = Vector{Vector{NamedTuple}}(undef, task.num_environments)

    Threads.@threads for environment_offset in 0:(task.num_environments - 1)
        environment_id = task.batch_id * task.num_environments + environment_offset
        environment_seed = task_seed(task, "environment", environment_id)
        cache_capacity = model === :baseline ? 1 : site_cache_capacity(task.L)
        environment = SiteIIDEnvironment(environment_seed, model, parameter; cache_capacity)
        environment_rows = NamedTuple[]

        for walk_id in 0:(task.walks_per_environment - 1)
            walk_seed = task_seed(task, "walk", environment_id, walk_id)
            walk_rng = StableRNG(walk_seed)
            started_at = string(now(UTC))
            start = monotonic_seconds()
            try
                path, raw_steps = loop_erased_walk(task.L, environment, walk_rng;
                                                   max_steps=max_steps)
                elapsed = monotonic_seconds() - start
                exit_x, exit_y = unpack_point(path[end])
                row = result_row(task, version, model, parameter, environment_id, walk_id,
                                 environment_seed, walk_seed, winding(path), length(path) - 1,
                                 raw_steps, exit_x, exit_y, elapsed, length(environment.weights),
                                 "ok", missing, missing, started_at, string(now(UTC)))
                push!(environment_rows, row)
            catch error
                elapsed = monotonic_seconds() - start
                row = result_row(task, version, model, parameter, environment_id, walk_id,
                                 environment_seed, walk_seed, missing, missing, missing,
                                 missing, missing, elapsed, length(environment.weights),
                                 "failed", string(typeof(error)), sprint(showerror, error),
                                 started_at, string(now(UTC)))
                push!(environment_rows, row)
            end
        end
        rows_by_environment[environment_offset + 1] = environment_rows
    end
    return reduce(vcat, rows_by_environment; init=NamedTuple[])
end

function write_csv_atomic(path::AbstractString, rows)
    mkpath(dirname(path))
    temporary = joinpath(dirname(path), string(".", basename(path), ".tmp.", getpid()))
    CSV.write(temporary, rows; transform=(_column, value) -> something(value, missing))
    mv(temporary, path; force=true)
end

function parse_cli(args)::Dict{String,String}
    options = Dict{String,String}()
    index = 1
    while index <= length(args)
        startswith(args[index], "--") || throw(ArgumentError("unexpected argument: $(args[index])"))
        key = args[index][3:end]
        if key in ("force", "rerun-failed", "allow-incomplete", "no-plots")
            options[key] = "true"
            index += 1
        else
            index == length(args) && throw(ArgumentError("missing value for --$key"))
            options[key] = args[index + 1]
            index += 2
        end
    end
    return options
end

function require_option(options, key)
    haskey(options, key) || throw(ArgumentError("missing required option --$key"))
    return options[key]
end

function main_run_batch(args=ARGS)::Int
    options = parse_cli(args)
    task_id = parse(Int, require_option(options, "task-id"))
    config = require_option(options, "config")
    output_dir = get(options, "output-dir", "results_hpc")
    task = load_config_row(config, task_id)
    path = result_path(output_dir, task)
    expected = expected_observations(task)
    force = haskey(options, "force")
    rerun_failed = haskey(options, "rerun-failed")

    if isfile(path) && !force
        if completed_result(path, expected)
            println("Task $task_id already complete: $path")
            return 0
        elseif !rerun_failed
            println(stderr, "Existing output is incomplete; use --rerun-failed or --force: $path")
            return 3
        end
    end

    max_steps = haskey(options, "max-steps") ? parse(Int, options["max-steps"]) : nothing
    start = monotonic_seconds()
    rows = run_batch(task; max_steps=max_steps)
    write_csv_atomic(path, rows)
    seconds = monotonic_seconds() - start
    failures = count(row -> row.status != "ok", rows)
    @printf("Task %d wrote %d rows to %s in %.2fs; failures=%d\n",
            task_id, length(rows), path, seconds, failures)
    return failures == 0 ? 0 : 2
end

function generate_config(output::AbstractString, distributions, sizes;
                         batches::Int=1, num_environments::Int=2,
                         walks_per_environment::Int=3,
                         base_seed::Integer=20260623)
    rows = NamedTuple[]
    task_id = 0
    for spec in distributions
        distribution, params = parse_distribution_argument(spec)
        distribution_spec(distribution, params)
        for L in sizes, batch_id in 0:(batches - 1)
            push!(rows, (
                task_id=task_id,
                environment_model=ENVIRONMENT_MODEL,
                distribution=distribution,
                distribution_params=canonical_json(params),
                L=L,
                batch_id=batch_id,
                num_environments=num_environments,
                walks_per_environment=walks_per_environment,
                base_seed=stable_seed(base_seed, task_id),
            ))
            task_id += 1
        end
    end
    mkpath(dirname(output))
    CSV.write(output, rows)
    return rows
end

function main_generate_config(args=ARGS)::Int
    options = parse_cli(args)
    output = get(options, "output", "configs/hpc_test.csv")
    preset = get(options, "preset", "")
    distributions = if preset == "hpc_full"
        ["baseline", "gamma:0.5", "gamma:1.0", "gamma:2.0", "exponential",
         "lognormal:0.5", "lognormal:1.0", "pareto:2.0", "pareto:3.0",
         "uniform:0.7", "beta:0.5", "weibull:0.7", "inverse_gamma:2.2",
         "bernoulli:0.8", "triangular:0.8"]
    elseif haskey(options, "distributions")
        split(options["distributions"], ',')
    else
        ["baseline", "gamma:1.0"]
    end
    sizes = preset == "hpc_full" ? [32, 64, 128, 256, 512, 1024] :
            parse.(Int, split(get(options, "sizes", "64,128"), ','))
    rows = generate_config(output, distributions, sizes;
        batches=parse(Int, get(options, "batches", "1")),
        num_environments=parse(Int, get(options, "num-environments", "2")),
        walks_per_environment=parse(Int, get(options, "walks-per-environment", "3")),
        base_seed=parse(Int, get(options, "base-seed", "20260623")))
    println("Wrote $(length(rows)) tasks to $output")
    println("Slurm array range: 0-$(length(rows) - 1)")
    return 0
end
