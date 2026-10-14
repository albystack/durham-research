#!/usr/bin/env julia

"""
Restart-safe finite-volume Kasteleyn evaluator for the direct weighted square-
grid dimer model. Each output row is one independent frozen environment and
contains deterministic conditional central-height moments. No Markov-chain
draws or replica pairing are involved.
"""

using AztecDiamond.GlauberSquareGrid
using Printf
using Random
using SHA

const OUTPUT_HEADER =
    "campaign_id,L,environment_id,environment_seed,conditional_mean," *
    "conditional_variance,conditional_second_moment,log_partition," *
    "partition_phase_real,partition_phase_imaginary,matrix_order,precision_bits," *
    "crossed_edges," *
    "relative_solve_residual,maximum_probability_imaginary_residual," *
    "maximum_covariance_imaginary_residual,elapsed_seconds"

function print_help()
    println("""
    Evaluate central-height moments with a finite Kasteleyn matrix.

    Usage:
      julia --project=aztec aztec/scripts/run_glauber_kasteleyn_campaign.jl [options]

    Options:
      --config PATH          CSV schedule (default: Kasteleyn pilot config)
      --output-dir PATH      restart-safe batch output directory
      --distribution NAME    constant, gamma, lognormal, or uniform
      --parameter FLOAT      Gamma shape or lognormal sigma (default: 0.5)
      --base-seed UINT       deterministic environment seed namespace
      --task-id INT          one expanded task; 0 runs all tasks
      --list-tasks           print the expanded task table and exit
      -h, --help             show this help

    The dense factorization is a validation implementation. Benchmark memory
    and runtime before extending the schedule beyond the supplied pilot sizes.
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return (help=true,)
    options = Dict(
        "config" => joinpath(@__DIR__, "..", "configs",
                             "glauber_square_grid_kasteleyn_pilot.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output",
                                 "glauber_square_grid_kasteleyn_pilot"),
        "distribution" => "gamma",
        "parameter" => "0.5",
        "base-seed" => "2026082201",
        "task-id" => "0",
    )
    list_tasks = false
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--list-tasks"
            list_tasks = true
            index += 1
            continue
        end
        startswith(argument, "--") || error("unexpected argument: $argument")
        key = argument[3:end]
        haskey(options, key) || error("unknown option: --$key")
        index < length(arguments) || error("missing value after --$key")
        options[key] = arguments[index + 1]
        index += 2
    end

    distribution = Symbol(lowercase(options["distribution"]))
    distribution in (:constant, :gamma, :lognormal, :uniform) || error(
        "--distribution must be constant, gamma, lognormal, or uniform")
    parameter = parse(Float64, options["parameter"])
    isfinite(parameter) && parameter > 0 || error("--parameter must be positive")
    task_id = parse(Int, options["task-id"])
    task_id >= 0 || error("--task-id must be nonnegative")
    return (
        help=false,
        config=abspath(options["config"]),
        output_dir=abspath(options["output-dir"]),
        distribution=distribution,
        parameter=parameter,
        base_seed=parse(UInt64, options["base-seed"]),
        task_id=task_id,
        list_tasks=list_tasks,
    )
end

function read_config(path)
    lines = readlines(path)
    isempty(lines) && error("empty config: $path")
    expected = "L,first_environment_id,environments,batch_size"
    strip(first(lines)) == expected || error("unexpected Kasteleyn config header in $path")
    rows = NamedTuple[]
    for (offset, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        length(fields) == 4 || error("invalid config row $(offset + 1)")
        row = (
            L=parse(Int, fields[1]),
            first_environment_id=parse(Int, fields[2]),
            environments=parse(Int, fields[3]),
            batch_size=parse(Int, fields[4]),
        )
        row.L > 0 || error("L must be positive")
        row.first_environment_id > 0 || error("first_environment_id must be positive")
        row.environments > 0 || error("environments must be positive")
        row.batch_size > 0 || error("batch_size must be positive")
        push!(rows, row)
    end
    isempty(rows) && error("config has no data rows")
    length(unique(row.L for row in rows)) == length(rows) || error("L values must be unique")
    return rows
end

function expand_tasks(rows)
    tasks = NamedTuple[]
    for row in rows
        first_id = row.first_environment_id
        last_id = first_id + row.environments - 1
        batch_id = 0
        while first_id <= last_id
            batch_id += 1
            final_id = min(first_id + row.batch_size - 1, last_id)
            push!(tasks, merge(row, (
                task_id=length(tasks) + 1,
                batch_id=batch_id,
                first_id=first_id,
                last_id=final_id,
            )))
            first_id = final_id + 1
        end
    end
    return tasks
end

@inline function splitmix64(value::UInt64)
    value += 0x9e3779b97f4a7c15
    value = (value ⊻ (value >> 30)) * 0xbf58476d1ce4e5b9
    value = (value ⊻ (value >> 27)) * 0x94d049bb133111eb
    return value ⊻ (value >> 31)
end

"Match the environment namespace used by the direct Glauber campaign runner."
function environment_seed(base_seed::UInt64, L::Int, environment_id::Int)
    key = splitmix64(base_seed ⊻ (UInt64(L) * 0xd6e8feb86659fd93) ⊻
                     (UInt64(environment_id) * 0xa5a3564e27f8862f))
    return splitmix64(key ⊻ 0x656e7669726f6e6d)
end

parameter_token(value::Float64) = replace(@sprintf("%.12g", value), "-" => "m", "." => "p")
campaign_id(parsed) = "glauber_square_grid_kasteleyn__$(parsed.distribution)__p_$(parameter_token(parsed.parameter))"

function weights_for_environment(parsed, L, seed)
    parsed.distribution === :constant && return constant_edge_weights(L)
    return random_edge_weights(Random.Xoshiro(seed), L;
                               distribution=parsed.distribution, parameter=parsed.parameter)
end

function batch_paths(parsed, task)
    directory = joinpath(parsed.output_dir, @sprintf("L_%04d", task.L))
    stem = @sprintf("batch_%04d", task.batch_id)
    return (
        data=joinpath(directory, stem * ".csv"),
        execution=joinpath(directory, "execution_$(lpad(task.batch_id, 4, '0')).txt"),
    )
end

function write_atomic(writer::Function, path::AbstractString)
    mkpath(dirname(path))
    temporary = path * ".$(getpid()).tmp"
    try
        open(temporary, "w") do io
            writer(io)
        end
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
end

function valid_existing_batch(path, parsed, task)
    isfile(path) || return false
    lines = readlines(path)
    length(lines) == 2 + task.last_id - task.first_id || return false
    strip(first(lines)) == OUTPUT_HEADER || return false
    for (line, environment_id) in zip(lines[2:end], task.first_id:task.last_id)
        fields = split(strip(line), ',')
        length(fields) == 17 || return false
        fields[1] == campaign_id(parsed) || return false
        parse(Int, fields[2]) == task.L || return false
        parse(Int, fields[3]) == environment_id || return false
        parse(UInt64, fields[4]) == environment_seed(
            parsed.base_seed, task.L, environment_id) || return false
        floating = parse.(Float64, fields[[5, 6, 7, 8, 9, 10, 14, 15, 16, 17]])
        all(isfinite, floating) || return false
        mean_value, variance, second_moment = floating[1:3]
        variance >= 0 || return false
        isapprox(second_moment, variance + mean_value^2;
                 rtol=1e-12, atol=1e-12) || return false
        parse(Int, fields[11]) == 2 * task.L^2 || return false
        parse(Int, fields[12]) in (53, 256) || return false
        parse(Int, fields[13]) == task.L || return false
        floating[7] <= 1e-10 || return false
        floating[10] >= 0 || return false
    end
    return true
end

function result_row(parsed, task, environment_id)
    seed = environment_seed(parsed.base_seed, task.L, environment_id)
    weights = weights_for_environment(parsed, task.L, seed)
    result = nothing
    elapsed = @elapsed result = center_height_moments_kasteleyn(weights)
    precision_bits = 53
    if result.relative_solve_residual > 1e-10
        high_precision_elapsed = @elapsed result = setprecision(256) do
            center_height_moments_kasteleyn(weights; number_type=BigFloat)
        end
        elapsed += high_precision_elapsed
        precision_bits = 256
    end
    return (
        campaign_id=campaign_id(parsed),
        L=task.L,
        environment_id=environment_id,
        environment_seed=seed,
        conditional_mean=result.mean,
        conditional_variance=result.variance,
        conditional_second_moment=result.variance + result.mean^2,
        log_partition=result.log_partition,
        partition_phase_real=real(result.partition_phase),
        partition_phase_imaginary=imag(result.partition_phase),
        matrix_order=result.matrix_order,
        precision_bits=precision_bits,
        crossed_edges=result.crossed_edges,
        relative_solve_residual=result.relative_solve_residual,
        maximum_probability_imaginary_residual=
            result.maximum_probability_imaginary_residual,
        maximum_covariance_imaginary_residual=
            result.maximum_covariance_imaginary_residual,
        elapsed_seconds=elapsed,
    )
end

function write_rows(path, rows)
    write_atomic(path) do io
        println(io, OUTPUT_HEADER)
        for row in rows
            println(io, join((
                row.campaign_id,
                row.L,
                row.environment_id,
                row.environment_seed,
                row.conditional_mean,
                row.conditional_variance,
                row.conditional_second_moment,
                row.log_partition,
                row.partition_phase_real,
                row.partition_phase_imaginary,
                row.matrix_order,
                row.precision_bits,
                row.crossed_edges,
                row.relative_solve_residual,
                row.maximum_probability_imaginary_residual,
                row.maximum_covariance_imaginary_residual,
                row.elapsed_seconds,
            ), ','))
        end
    end
end

file_hash(path) = isfile(path) ? bytes2hex(sha256(read(path))) : "unavailable"

function metadata_text(parsed, rows, tasks)
    repository_root = normpath(joinpath(@__DIR__, "..", ".."))
    schedule = join((join((row.L, row.first_environment_id, row.environments,
                           row.batch_size), ':') for row in rows), ',')
    return join([
        "schema_version=1",
        "campaign_id=$(campaign_id(parsed))",
        "geometry=square_grid_direct_weighted_dimer",
        "boundary=tileable_extremal_height_boundary_v1",
        "observable=conditional_central_face_height_mean_and_variance",
        "method=dense_finite_bipartite_kasteleyn_factorization",
        "precision_policy=Float64 with automatic 256-bit fallback when relative solve residual exceeds 1e-10",
        "horizontal_kasteleyn_phase=1",
        "vertical_kasteleyn_phase=im",
        "distribution=$(parsed.distribution)",
        "parameter=$(parsed.parameter)",
        "parameterization=$(parsed.distribution === :gamma ? "Gamma(shape=k,scale=1/k), mean one" : parsed.distribution === :lognormal ? "LogNormal(-sigma^2/2,sigma), mean one" : parsed.distribution === :uniform ? "Uniform(0,2), mean one" : "all weights one")",
        "base_seed=$(parsed.base_seed)",
        "config=$(parsed.config)",
        "config_sha256=$(file_hash(parsed.config))",
        "project_sha256=$(file_hash(joinpath(repository_root, "aztec", "Project.toml")))",
        "manifest_sha256=$(file_hash(joinpath(repository_root, "aztec", "Manifest.toml")))",
        "kernel_sha256=$(file_hash(joinpath(repository_root, "aztec", "src", "GlauberSquareGrid.jl")))",
        "runner_sha256=$(file_hash(abspath(@__FILE__)))",
        "schedule=$schedule",
        "expanded_tasks=$(length(tasks))",
        "julia_version=$(VERSION)",
        "seed_identity=identical environment namespace to run_glauber_square_grid_campaign.jl",
        "numerical_scope=dense validation implementation; benchmark before larger sizes",
    ], '\n') * "\n"
end

function ensure_metadata(parsed, rows, tasks)
    mkpath(parsed.output_dir)
    expected = metadata_text(parsed, rows, tasks)
    path = joinpath(parsed.output_dir, "campaign_metadata.txt")
    if isfile(path)
        read(path, String) == expected || error("existing campaign metadata differs: $path")
    else
        write_atomic(path) do io
            print(io, expected)
        end
    end
end

function run_task(parsed, task)
    paths = batch_paths(parsed, task)
    artifacts_exist = isfile(paths.data) || isfile(paths.execution)
    if artifacts_exist
        isfile(paths.data) && isfile(paths.execution) || error(
            "incomplete existing task artifacts: $(paths.data)")
        valid_existing_batch(paths.data, parsed, task) || error(
            "existing batch failed restart validation: $(paths.data)")
        println("skip existing $(paths.data)")
        return
    end
    started = time()
    rows = [result_row(parsed, task, environment_id)
            for environment_id in task.first_id:task.last_id]
    write_rows(paths.data, rows)
    write_atomic(paths.execution) do io
        println(io, "task_id=$(task.task_id)")
        println(io, "L=$(task.L)")
        println(io, "first_environment_id=$(task.first_id)")
        println(io, "last_environment_id=$(task.last_id)")
        println(io, "elapsed_seconds=$(time() - started)")
    end
    @printf("completed Kasteleyn L=%d task=%d environments=%d:%d in %.2fs\n",
            task.L, task.task_id, task.first_id, task.last_id, time() - started)
end

function main(arguments=ARGS)
    parsed = parse_arguments(arguments)
    if parsed.help
        print_help()
        return
    end
    rows = read_config(parsed.config)
    tasks = expand_tasks(rows)
    if parsed.list_tasks
        println("task_id,L,batch_id,first_environment_id,last_environment_id")
        for task in tasks
            println(join((task.task_id, task.L, task.batch_id,
                          task.first_id, task.last_id), ','))
        end
        return
    end
    parsed.task_id <= length(tasks) || error(
        "task $(parsed.task_id) exceeds expanded task count $(length(tasks))")
    ensure_metadata(parsed, rows, tasks)
    selected = parsed.task_id == 0 ? tasks : [tasks[parsed.task_id]]
    for task in selected
        run_task(parsed, task)
    end
end

main()
