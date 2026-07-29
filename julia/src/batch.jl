const SCHEMA_VERSION = "batch_v7_julia"
const ENVIRONMENT_MODEL = "site_iid"
const ENVIRONMENT_MODELS = Set(["site_iid", "temporal_iid"])
const TEMPORAL_PILOT_DISTRIBUTIONS = [
    "baseline", "gamma:0.5", "lognormal:1.0", "pareto:2.0",
]
const FULL_DISTRIBUTIONS = [
    "baseline", "gamma:0.5", "gamma:1.0", "gamma:2.0", "exponential",
    "lognormal:0.5", "lognormal:1.0", "pareto:2.0", "pareto:3.0",
    "uniform:0.7", "beta:0.5", "weibull:0.7", "inverse_gamma:2.2",
    "bernoulli:0.8", "triangular:0.8",
]

# These six cases had 5,000 observations at L=32,...,1024 in the reported data;
# the remaining cases had 2,000. All cases had 5,000 at L=16 and 500 at the
# two largest common sizes. Only baseline and gamma(shape=1) extended to L=8192.
const HIGH_REPLICATION_DISTRIBUTIONS = Set([
    "baseline", "gamma:0.5", "lognormal:1.0", "pareto:2.0", "uniform:0.7",
    "weibull:0.7",
])

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

function validate_environment_model(environment_model::AbstractString)::String
    model = String(environment_model)
    model in ENVIRONMENT_MODELS ||
        throw(ArgumentError("unknown environment model: $model"))
    return model
end

function load_config_row(config_path::AbstractString, task_id::Integer)::BatchConfig
    for row in CSV.File(config_path; types=String)
        parse(Int, row.task_id) == task_id || continue
        raw_environment_model = hasproperty(row, :environment_model) ? row.environment_model : ""
        environment_model = isempty(raw_environment_model) ? ENVIRONMENT_MODEL : raw_environment_model
        validate_environment_model(environment_model)
        task = BatchConfig(
            parse(Int, row.task_id), row.distribution, parse_params(row.distribution_params),
            parse(Int, row.L), parse(Int, row.batch_id), parse(Int, row.num_environments),
            parse(Int, row.walks_per_environment), parse(UInt64, row.base_seed),
            environment_model,
        )
        task.num_environments > 0 || throw(ArgumentError("num_environments must be positive"))
        task.walks_per_environment > 0 ||
            throw(ArgumentError("walks_per_environment must be positive"))
        task.environment_model == "temporal_iid" && task.walks_per_environment != 1 &&
            throw(ArgumentError(
                "temporal_iid requires walks_per_environment=1: each walk is one " *
                "independent temporal realisation"))
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

function require_strict_annealed(tasks)
    invalid = [task.task_id for task in tasks if task.walks_per_environment != 1]
    isempty(invalid) || throw(ArgumentError(
        "strict annealed sampling requires walks_per_environment=1; invalid task IDs: " *
        join(invalid[1:min(end, 10)], ", ")))
    return nothing
end

"Require independent, one-walk temporal realisations and no quenched grouping."
function require_temporal_iid(tasks)
    wrong_model = [task.task_id for task in tasks
                   if task.environment_model != "temporal_iid"]
    isempty(wrong_model) || throw(ArgumentError(
        "temporal campaign contains non-temporal task IDs: " *
        join(wrong_model[1:min(end, 10)], ", ")))
    invalid = [task.task_id for task in tasks if task.walks_per_environment != 1]
    isempty(invalid) || throw(ArgumentError(
        "temporal_iid requires walks_per_environment=1; invalid task IDs: " *
        join(invalid[1:min(end, 10)], ", ")))
    return nothing
end

"Require exactly one independent pair of walks in every sampled environment."
function require_double_dimer(tasks)
    invalid = [task.task_id for task in tasks
               if task.environment_model != "site_iid" ||
                  task.walks_per_environment != 2]
    isempty(invalid) || throw(ArgumentError(
        "double-dimer sampling requires site_iid and walks_per_environment=2; " *
        "invalid task IDs: " *
        join(invalid[1:min(end, 10)], ", ")))
    return nothing
end

function completed_result(path::AbstractString, expected_rows::Integer)::Bool
    isfile(path) || return false
    rows = collect(CSV.File(path; types=String))
    return length(rows) == expected_rows && all(row -> row.status == "ok", rows)
end

"Restart-safe completion check that also rejects output from a different task."
function completed_result(path::AbstractString, task::BatchConfig)::Bool
    isfile(path) || return false
    rows = collect(CSV.File(path; types=String))
    length(rows) == expected_observations(task) || return false
    expected_params = params_json(task)
    return all(rows) do row
        row.status == "ok" &&
        row.task_id == string(task.task_id) &&
        row.batch_id == string(task.batch_id) &&
        row.distribution == task.distribution &&
        row.distribution_params == expected_params &&
        row.environment_model == task.environment_model &&
        row.L == string(task.L)
    end
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
                    raw_length, exit_x, exit_y, runtime, weight_seed, direction_seed,
                    sampled_sites, sampled_weight_vectors, direction_counts, status,
                    error_type, error_message, started_at, finished_at)
    north_steps, east_steps, south_steps, west_steps = direction_counts
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
        weight_seed=weight_seed,
        direction_seed=direction_seed,
        winding=winding_value,
        loop_erased_path_length=path_length,
        raw_walk_length=raw_length,
        exit_location=exit_x === missing ? missing : string(exit_x, ",", exit_y),
        exit_x=exit_x,
        exit_y=exit_y,
        runtime=runtime,
        sampled_site_count=sampled_sites,
        sampled_weight_vector_count=sampled_weight_vectors,
        north_steps=north_steps,
        east_steps=east_steps,
        south_steps=south_steps,
        west_steps=west_steps,
        status=status,
        error_type=error_type,
        error_message=error_message,
        started_at_utc=started_at,
        finished_at_utc=finished_at,
    )
end

function run_batch(task::BatchConfig; max_steps::Union{Nothing,Integer}=nothing)
    validate_environment_model(task.environment_model)
    task.environment_model == "temporal_iid" && task.walks_per_environment != 1 &&
        throw(ArgumentError(
            "temporal_iid requires walks_per_environment=1: each walk is one " *
            "independent temporal realisation"))
    model, parameter = distribution_spec(task.distribution, task.distribution_params)
    version = git_version()
    rows_by_environment = Vector{Vector{NamedTuple}}(undef, task.num_environments)

    Threads.@threads for environment_offset in 0:(task.num_environments - 1)
        environment_id = task.batch_id * task.num_environments + environment_offset
        environment_seed = task_seed(task, "environment", environment_id)
        environment = if task.environment_model == "site_iid"
            cache_capacity = model === :baseline ? 1 : site_cache_capacity(task.L)
            SiteIIDEnvironment(environment_seed, model, parameter; cache_capacity)
        else
            TemporalIIDEnvironment(environment_seed, model, parameter)
        end
        environment_rows = NamedTuple[]

        for walk_id in 0:(task.walks_per_environment - 1)
            walk_seed = task_seed(task, "walk", environment_id, walk_id)
            direction_rng = StableRNG(walk_seed)
            started_at = string(now(UTC))
            start = monotonic_seconds()
            try
                path, raw_steps, direction_counts = loop_erased_walk_diagnostics(
                    task.L, environment, direction_rng; max_steps)
                elapsed = monotonic_seconds() - start
                exit_x, exit_y = unpack_point(path[end])
                sampled_sites = environment isa SiteIIDEnvironment ?
                    length(environment.weights) : missing
                sampled_weight_vectors = environment isa TemporalIIDEnvironment ?
                    environment.sampled_weight_vectors : missing
                row = result_row(task, version, model, parameter, environment_id, walk_id,
                                 environment_seed, walk_seed, winding(path), length(path) - 1,
                                 raw_steps, exit_x, exit_y, elapsed, environment.weight_seed,
                                 walk_seed, sampled_sites, sampled_weight_vectors,
                                 direction_counts,
                                 "ok", missing, missing, started_at, string(now(UTC)))
                push!(environment_rows, row)
            catch error
                elapsed = monotonic_seconds() - start
                sampled_sites = environment isa SiteIIDEnvironment ?
                    length(environment.weights) : missing
                sampled_weight_vectors = environment isa TemporalIIDEnvironment ?
                    environment.sampled_weight_vectors : missing
                row = result_row(task, version, model, parameter, environment_id, walk_id,
                                 environment_seed, walk_seed, missing, missing, missing,
                                 missing, missing, elapsed, environment.weight_seed, walk_seed,
                                 sampled_sites, sampled_weight_vectors,
                                 (missing, missing, missing, missing),
                                 "failed", string(typeof(error)), sprint(showerror, error),
                                 started_at, string(now(UTC)))
                push!(environment_rows, row)
            end
        end
        rows_by_environment[environment_offset + 1] = environment_rows
    end
    return reduce(vcat, rows_by_environment; init=NamedTuple[])
end

function run_campaign(config_path::AbstractString, output_dir::AbstractString;
                      strict_annealed::Bool=false, double_dimer::Bool=false,
                      rerun_failed::Bool=false,
                      start_task::Union{Nothing,Int}=nothing,
                      end_task::Union{Nothing,Int}=nothing,
                      max_steps::Union{Nothing,Int}=nothing)
    tasks = read_config(config_path)
    strict_annealed && double_dimer && throw(ArgumentError(
        "--strict-annealed and --double-dimer are mutually exclusive"))
    strict_annealed && require_strict_annealed(tasks)
    double_dimer && require_double_dimer(tasks)
    any(task -> task.environment_model == "temporal_iid", tasks) &&
        require_temporal_iid(tasks)
    selected = [task for task in tasks
                if (start_task === nothing || task.task_id >= start_task) &&
                   (end_task === nothing || task.task_id <= end_task)]
    isempty(selected) && throw(ArgumentError("no tasks selected"))

    pending = BatchConfig[]
    complete = 0
    for task in selected
        path = result_path(output_dir, task)
        if completed_result(path, task)
            complete += 1
        elseif isfile(path) && !rerun_failed
            throw(ErrorException(
                "incomplete output for task $(task.task_id); rerun with --rerun-failed: $path"))
        else
            push!(pending, task)
        end
    end

    println("Campaign tasks: $(length(selected)); already complete: $complete; pending: $(length(pending))")
    flush(stdout)
    campaign_start = monotonic_seconds()
    observations = 0
    for (index, task) in enumerate(pending)
        start = monotonic_seconds()
        rows = run_batch(task; max_steps)
        path = result_path(output_dir, task)
        write_csv_atomic(path, rows)
        failures = count(row -> row.status != "ok", rows)
        failures == 0 || throw(ErrorException(
            "task $(task.task_id) wrote $failures failed observations to $path"))
        observations += length(rows)
        elapsed = monotonic_seconds() - start
        total_done = complete + index
        @printf("[%d/%d] task=%d distribution=%s L=%d rows=%d seconds=%.2f\n",
                total_done, length(selected), task.task_id,
                safe_distribution_name(task.distribution, task.distribution_params),
                task.L, length(rows), elapsed)
        flush(stdout)
        task.L >= 2048 && GC.gc(false)
    end
    elapsed = monotonic_seconds() - campaign_start
    @printf("Campaign complete: tasks=%d new_observations=%d seconds=%.2f\n",
            length(selected), observations, elapsed)
    return (; tasks=length(selected), already_complete=complete,
            new_tasks=length(pending), new_observations=observations, elapsed)
end

function main_run_campaign(args=ARGS)::Int
    options = parse_cli(args)
    config = require_option(options, "config")
    output_dir = get(options, "output-dir", "results_strict_annealed")
    report = run_campaign(config, output_dir;
        strict_annealed=haskey(options, "strict-annealed"),
        double_dimer=haskey(options, "double-dimer"),
        rerun_failed=haskey(options, "rerun-failed"),
        start_task=haskey(options, "start-task") ? parse(Int, options["start-task"]) : nothing,
        end_task=haskey(options, "end-task") ? parse(Int, options["end-task"]) : nothing,
        max_steps=haskey(options, "max-steps") ? parse(Int, options["max-steps"]) : nothing)
    return 0
end

function write_csv_atomic(path::AbstractString, rows)
    mkpath(dirname(path))
    temporary = joinpath(dirname(path), string(".", basename(path), ".tmp.", getpid()))
    try
        CSV.write(temporary, rows; transform=(_column, value) -> something(value, missing))
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary)
    end
end

function parse_cli(args)::Dict{String,String}
    options = Dict{String,String}()
    index = 1
    while index <= length(args)
        startswith(args[index], "--") || throw(ArgumentError("unexpected argument: $(args[index])"))
        key = args[index][3:end]
        if key in ("force", "rerun-failed", "allow-incomplete", "strict-annealed",
                   "double-dimer")
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
    output_dir = get(options, "output-dir", "results")
    task = load_config_row(config, task_id)
    haskey(options, "strict-annealed") && haskey(options, "double-dimer") &&
        throw(ArgumentError("--strict-annealed and --double-dimer are mutually exclusive"))
    haskey(options, "strict-annealed") && require_strict_annealed([task])
    haskey(options, "double-dimer") && require_double_dimer([task])
    path = result_path(output_dir, task)
    expected = expected_observations(task)
    force = haskey(options, "force")
    rerun_failed = haskey(options, "rerun-failed")

    if isfile(path) && !force
        if completed_result(path, task)
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
                         base_seed::Integer=20260623,
                         environment_model::AbstractString=ENVIRONMENT_MODEL)
    environment_model = validate_environment_model(environment_model)
    environment_model == "temporal_iid" && walks_per_environment != 1 &&
        throw(ArgumentError("temporal_iid requires walks_per_environment=1"))
    rows = NamedTuple[]
    task_id = 0
    for spec in distributions
        distribution, params = parse_distribution_argument(spec)
        distribution_spec(distribution, params)
        for L in sizes, batch_id in 0:(batches - 1)
            push!(rows, (
                task_id=task_id,
                environment_model=environment_model,
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

function generate_temporal_iid_pilot_config(output::AbstractString;
                                            batch_walks::Int=100,
                                            base_seed::Integer=20260726)
    batch_walks > 0 || throw(ArgumentError("batch_walks must be positive"))
    size_counts = Pair{Int,Int}[
        16 => 2_000, 32 => 2_000, 64 => 2_000, 128 => 2_000, 256 => 2_000,
        512 => 1_000, 1024 => 500,
    ]
    rows = NamedTuple[]
    task_id = 0
    for spec in TEMPORAL_PILOT_DISTRIBUTIONS
        distribution, params = parse_distribution_argument(spec)
        for (L, walks) in size_counts
            walks % batch_walks == 0 || throw(ArgumentError(
                "$walks walks at L=$L are not divisible by batch_walks=$batch_walks"))
            for batch_id in 0:(walks ÷ batch_walks - 1)
                push!(rows, (
                    task_id=task_id,
                    environment_model="temporal_iid",
                    distribution=distribution,
                    distribution_params=canonical_json(params),
                    L=L,
                    batch_id=batch_id,
                    num_environments=batch_walks,
                    walks_per_environment=1,
                    base_seed=stable_seed(base_seed, "temporal_iid_pilot", task_id),
                ))
                task_id += 1
            end
        end
    end
    sum(row.num_environments for row in rows) == 46_000 ||
        error("internal error: temporal pilot must contain exactly 46,000 walks")
    mkpath(dirname(output))
    CSV.write(output, rows)
    return rows
end

"Reuse the Gamma arm of the temporal pilot and extend it beyond L=5000."
function generate_temporal_iid_gamma_length_config(
    output::AbstractString; extension_walks::Int=100,
    base_seed::Integer=20260726,
)
    extension_walks > 1 || throw(ArgumentError(
        "extension_walks must exceed one so uncertainty can be estimated"))
    rows = NamedTuple[]
    task_id = 0
    size_counts = Pair{Int,Int}[
        16 => 2_000, 32 => 2_000, 64 => 2_000, 128 => 2_000, 256 => 2_000,
        512 => 1_000, 1024 => 500,
    ]

    # Preserve the task IDs and seeds from temporal_iid_pilot.csv, so the
    # already-computed Gamma batches are restart-safe and directly reusable.
    for spec in TEMPORAL_PILOT_DISTRIBUTIONS
        distribution, params = parse_distribution_argument(spec)
        for (L, walks) in size_counts
            for batch_id in 0:(walks ÷ 100 - 1)
                if spec == "gamma:0.5"
                    push!(rows, (
                        task_id=task_id,
                        environment_model="temporal_iid",
                        distribution=distribution,
                        distribution_params=canonical_json(params),
                        L=L,
                        batch_id=batch_id,
                        num_environments=100,
                        walks_per_environment=1,
                        base_seed=stable_seed(
                            base_seed, "temporal_iid_pilot", task_id),
                    ))
                end
                task_id += 1
            end
        end
    end

    distribution, params = parse_distribution_argument("gamma:0.5")
    for L in (2048, 4096, 5000, 8192)
        push!(rows, (
            task_id=task_id,
            environment_model="temporal_iid",
            distribution=distribution,
            distribution_params=canonical_json(params),
            L=L,
            batch_id=0,
            num_environments=extension_walks,
            walks_per_environment=1,
            base_seed=stable_seed(
                base_seed, "temporal_iid_gamma_length", task_id),
        ))
        task_id += 1
    end
    mkpath(dirname(output))
    CSV.write(output, rows)
    return rows
end

function generate_temporal_iid_smoke_config(output::AbstractString;
                                            walks_per_task::Int=100,
                                            base_seed::Integer=20260726)
    walks_per_task > 1 || throw(ArgumentError("walks_per_task must exceed one"))
    return generate_config(output, TEMPORAL_PILOT_DISTRIBUTIONS, (16, 32, 64);
        batches=1, num_environments=walks_per_task, walks_per_environment=1,
        base_seed=stable_seed(base_seed, "temporal_iid_smoke"),
        environment_model="temporal_iid")
end

function generate_temporal_iid_benchmark_config(output::AbstractString;
                                                walks_per_task::Int=100,
                                                base_seed::Integer=20260726)
    walks_per_task > 0 || throw(ArgumentError("walks_per_task must be positive"))
    return generate_config(output, ("pareto:2.0",), (256, 512, 1024);
        batches=1, num_environments=walks_per_task, walks_per_environment=1,
        base_seed=stable_seed(base_seed, "temporal_iid_benchmark"),
        environment_model="temporal_iid")
end

function generate_strict_annealed_reproduction_config(output::AbstractString;
                                                       batch_observations::Int=100,
                                                       base_seed::Integer=20260713)
    batch_observations > 0 || throw(ArgumentError("batch_observations must be positive"))
    rows = NamedTuple[]
    task_id = 0
    for spec in FULL_DISTRIBUTIONS
        distribution, params = parse_distribution_argument(spec)
        size_counts = Pair{Int,Int}[16 => 5_000]
        middle_count = spec in HIGH_REPLICATION_DISTRIBUTIONS ? 5_000 : 2_000
        append!(size_counts, [L => middle_count for L in (32, 64, 128, 256, 512, 1024)])
        append!(size_counts, [2048 => 500, 4096 => 500])
        spec in ("baseline", "gamma:1.0") && push!(size_counts, 8192 => 500)

        for (L, observations) in size_counts
            observations % batch_observations == 0 || throw(ArgumentError(
                "$observations observations at L=$L are not divisible by " *
                "batch_observations=$batch_observations"))
            for batch_id in 0:(observations ÷ batch_observations - 1)
                push!(rows, (
                    task_id=task_id,
                    environment_model=ENVIRONMENT_MODEL,
                    distribution=distribution,
                    distribution_params=canonical_json(params),
                    L=L,
                    batch_id=batch_id,
                    num_environments=batch_observations,
                    walks_per_environment=1,
                    base_seed=stable_seed(base_seed, "strict_annealed", task_id),
                ))
                task_id += 1
            end
        end
    end
    mkpath(dirname(output))
    CSV.write(output, rows)
    return rows
end

"Generate the matched double-dimer campaign: one two-walk pair per environment."
function generate_double_dimer_reproduction_config(output::AbstractString;
                                                     batch_environments::Int=100,
                                                     base_seed::Integer=20260718)
    batch_environments > 0 || throw(ArgumentError("batch_environments must be positive"))
    rows = NamedTuple[]
    task_id = 0
    for spec in FULL_DISTRIBUTIONS
        distribution, params = parse_distribution_argument(spec)
        size_counts = Pair{Int,Int}[16 => 5_000]
        middle_count = spec in HIGH_REPLICATION_DISTRIBUTIONS ? 5_000 : 2_000
        append!(size_counts, [L => middle_count for L in (32, 64, 128, 256, 512, 1024)])
        append!(size_counts, [2048 => 500, 4096 => 500])
        spec in ("baseline", "gamma:1.0") && push!(size_counts, 8192 => 500)

        for (L, environments) in size_counts
            environments % batch_environments == 0 || throw(ArgumentError(
                "$environments environments at L=$L are not divisible by " *
                "batch_environments=$batch_environments"))
            for batch_id in 0:(environments ÷ batch_environments - 1)
                push!(rows, (
                    task_id=task_id,
                    environment_model=ENVIRONMENT_MODEL,
                    distribution=distribution,
                    distribution_params=canonical_json(params),
                    L=L,
                    batch_id=batch_id,
                    num_environments=batch_environments,
                    walks_per_environment=2,
                    base_seed=stable_seed(base_seed, "double_dimer", task_id),
                ))
                task_id += 1
            end
        end
    end
    mkpath(dirname(output))
    CSV.write(output, rows)
    return rows
end

"Generate one lightweight all-cases task per size for an end-to-end pilot."
function generate_double_dimer_pilot_config(output::AbstractString;
                                            environments_per_size::Int=100,
                                            base_seed::Integer=20260718)
    environments_per_size > 1 || throw(ArgumentError(
        "environments_per_size must exceed one so a variance is defined"))
    rows = NamedTuple[]
    task_id = 0
    for spec in FULL_DISTRIBUTIONS
        distribution, params = parse_distribution_argument(spec)
        sizes = collect(16 .* (2 .^ (0:8)))
        spec in ("baseline", "gamma:1.0") && push!(sizes, 8192)
        for L in sizes
            push!(rows, (
                task_id=task_id,
                environment_model=ENVIRONMENT_MODEL,
                distribution=distribution,
                distribution_params=canonical_json(params),
                L=L,
                batch_id=0,
                num_environments=environments_per_size,
                walks_per_environment=2,
                base_seed=stable_seed(base_seed, "double_dimer_pilot", task_id),
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
    output = get(options, "output", "configs/test.csv")
    preset = get(options, "preset", "")
    if preset == "strict_annealed_reproduction"
        rows = generate_strict_annealed_reproduction_config(output;
            batch_observations=parse(Int, get(options, "batch-observations", "100")),
            base_seed=parse(Int, get(options, "base-seed", "20260713")))
        println("Wrote $(length(rows)) strict-annealed tasks to $output")
        println("Slurm array range: 0-$(length(rows) - 1)")
        return 0
    elseif preset == "double_dimer_reproduction"
        rows = generate_double_dimer_reproduction_config(output;
            batch_environments=parse(Int, get(options, "batch-environments", "100")),
            base_seed=parse(Int, get(options, "base-seed", "20260718")))
        println("Wrote $(length(rows)) double-dimer tasks to $output")
        println("Each environment contains exactly two independent walks")
        println("Slurm array range: 0-$(length(rows) - 1)")
        return 0
    elseif preset == "double_dimer_pilot"
        rows = generate_double_dimer_pilot_config(output;
            environments_per_size=parse(Int, get(options, "pilot-environments", "100")),
            base_seed=parse(Int, get(options, "base-seed", "20260718")))
        println("Wrote $(length(rows)) double-dimer pilot tasks to $output")
        println("Slurm array range: 0-$(length(rows) - 1)")
        return 0
    elseif preset == "temporal_iid_pilot"
        rows = generate_temporal_iid_pilot_config(output;
            batch_walks=parse(Int, get(options, "batch-walks", "100")),
            base_seed=parse(Int, get(options, "base-seed", "20260726")))
        println("Wrote $(length(rows)) temporal-i.i.d. pilot tasks and 46,000 walks to $output")
        println("Each temporal realisation contains exactly one walk")
        return 0
    elseif preset == "temporal_iid_smoke"
        rows = generate_temporal_iid_smoke_config(output;
            walks_per_task=parse(Int, get(options, "walks-per-task", "100")),
            base_seed=parse(Int, get(options, "base-seed", "20260726")))
        println("Wrote $(length(rows)) temporal-i.i.d. smoke tasks to $output")
        return 0
    elseif preset == "temporal_iid_gamma_length"
        rows = generate_temporal_iid_gamma_length_config(output;
            extension_walks=parse(Int, get(options, "extension-walks", "100")),
            base_seed=parse(Int, get(options, "base-seed", "20260726")))
        println("Wrote $(length(rows)) temporal Gamma length-scaling tasks to $output")
        println("The campaign extends through L=8192, including L=5000")
        return 0
    elseif preset == "temporal_iid_benchmark"
        rows = generate_temporal_iid_benchmark_config(output;
            walks_per_task=parse(Int, get(options, "walks-per-task", "100")),
            base_seed=parse(Int, get(options, "base-seed", "20260726")))
        println("Wrote $(length(rows)) temporal-i.i.d. benchmark tasks to $output")
        return 0
    end

    distributions = if haskey(options, "distributions")
        split(options["distributions"], ',')
    else
        ["baseline", "gamma:1.0"]
    end
    sizes = parse.(Int, split(get(options, "sizes", "64,128"), ','))
    rows = generate_config(output, distributions, sizes;
        batches=parse(Int, get(options, "batches", "1")),
        num_environments=parse(Int, get(options, "num-environments", "2")),
        walks_per_environment=parse(Int, get(options, "walks-per-environment", "3")),
        base_seed=parse(Int, get(options, "base-seed", "20260623")),
        environment_model=get(options, "environment-model", ENVIRONMENT_MODEL))
    println("Wrote $(length(rows)) tasks to $output")
    println("Slurm array range: 0-$(length(rows) - 1)")
    return 0
end
