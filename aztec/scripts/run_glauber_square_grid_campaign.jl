#!/usr/bin/env julia

"""
Restart-safe pilot runner for the direct square-grid weighted-dimer Glauber
sampler.  Every environment owns two independent chains and one atomic output
batch.  In `pilot` mode the chains start from the maximal and minimal height
configurations; production sampling is intentionally not enabled here until
the pilot's mixing evidence has been reviewed.
"""

using AztecDiamond
using AztecDiamond.GlauberSquareGrid
using Base.Threads
using Printf
using Random
using SHA

const TRACE_HEADER =
    "campaign_id,phase,L,environment_id,environment_seed,replica,start," *
    "replica_seed,draw_index,attempt_time,center_height"
const DIAGNOSTIC_HEADER =
    "campaign_id,phase,L,environment_id,environment_seed,replica,start," *
    "replica_seed,mean,variance,integrated_autocorrelation_time," *
    "effective_sample_size,changed_rate,final_total_change_rate," *
    "burn_in_attempts,thin_attempts,trace_samples"

function print_help()
    println("""
    Run a resumable direct square-grid dimer Glauber pilot.

    Usage:
      julia --project=aztec aztec/scripts/run_glauber_square_grid_campaign.jl [options]

    Options:
      --config PATH              pilot CSV; see aztec/configs/glauber_square_grid_pilot.csv
      --output-dir PATH          directory for atomic raw batches and metadata
      --distribution NAME        constant, gamma, lognormal, or uniform (default: gamma)
      --parameter FLOAT          Gamma shape or lognormal sigma (default: 0.5)
      --base-seed UINT           public deterministic campaign seed
      --task-id INT              one expanded batch task; 0 runs all (default: 0)
      --algorithm NAME           accelerated, literal, or tempered (default: accelerated)
      --tempering-betas CSV      ordered auxiliary inverse temperatures (default: 0:0.1:1)
      --swap-interval-attempts N attempts per replica between exchange rounds (default: 10000)
      --list-tasks               print task table and exit
      -h, --help                 show this help

    The pilot always compares independent chains from the maximum and minimum
    height configurations in one frozen edge environment.  It is not a
    production size-scaling runner.
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return (help=true,)
    options = Dict(
        "config" => joinpath(@__DIR__, "..", "configs", "glauber_square_grid_pilot.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "glauber_square_grid_pilot"),
        "distribution" => "gamma",
        "parameter" => "0.5",
        "base-seed" => "20260821",
        "task-id" => "0",
        "algorithm" => "accelerated",
        "tempering-betas" => "0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
        "swap-interval-attempts" => "10000",
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
    parameter > 0 || error("--parameter must be positive")
    algorithm = Symbol(lowercase(options["algorithm"]))
    algorithm in (:accelerated, :literal, :tempered) || error(
        "--algorithm must be accelerated, literal, or tempered")
    tempering_betas = parse.(Float64, split(options["tempering-betas"], ','))
    length(tempering_betas) >= 2 || error("--tempering-betas needs at least two values")
    issorted(tempering_betas) && all(beta -> 0 <= beta <= 1, tempering_betas) &&
        last(tempering_betas) == 1.0 || error(
            "--tempering-betas must be sorted in [0,1] and end at 1")
    swap_interval_attempts = parse(Int, options["swap-interval-attempts"])
    swap_interval_attempts > 0 || error("--swap-interval-attempts must be positive")
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
        algorithm=algorithm,
        tempering_betas=tempering_betas,
        swap_interval_attempts=swap_interval_attempts,
        list_tasks=list_tasks,
    )
end

function read_config(path)
    lines = readlines(path)
    isempty(lines) && error("empty config: $path")
    expected = "L,first_environment_id,environments,batch_size,burn_in_attempts,thin_attempts,trace_samples"
    strip(first(lines)) == expected || error("unexpected Glauber config header in $path")
    rows = NamedTuple[]
    for (offset, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        length(fields) == 7 || error("invalid config row $(offset + 1)")
        row = (
            L=parse(Int, fields[1]),
            first_environment_id=parse(Int, fields[2]),
            environments=parse(Int, fields[3]),
            batch_size=parse(Int, fields[4]),
            burn_in_attempts=parse(Int, fields[5]),
            thin_attempts=parse(Int, fields[6]),
            trace_samples=parse(Int, fields[7]),
        )
        row.L >= 1 || error("L must be positive")
        row.first_environment_id > 0 || error("first_environment_id must be positive")
        row.environments > 0 || error("environments must be positive")
        row.batch_size > 0 || error("batch_size must be positive")
        row.burn_in_attempts > 0 || error("burn_in_attempts must be positive")
        row.thin_attempts > 0 || error("thin_attempts must be positive")
        row.trace_samples >= 20 || error("pilot trace_samples must be at least 20")
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

function derive_seeds(base_seed::UInt64, L::Int, environment_id::Int)
    key = splitmix64(base_seed ⊻ (UInt64(L) * 0xd6e8feb86659fd93) ⊻
                     (UInt64(environment_id) * 0xa5a3564e27f8862f))
    return (
        environment=splitmix64(key ⊻ 0x656e7669726f6e6d),
        maximum=splitmix64(key ⊻ 0x6d61785f63686169),
        minimum=splitmix64(key ⊻ 0x6d696e5f63686169),
    )
end

parameter_token(value::Float64) = replace(@sprintf("%.12g", value), "-" => "m", "." => "p")
campaign_id(parsed) = "glauber_square_grid__$(parsed.distribution)__p_$(parameter_token(parsed.parameter))"

function weights_for_environment(parsed, L, seed)
    parsed.distribution === :constant && return constant_edge_weights(L)
    return random_edge_weights(Random.Xoshiro(seed), L;
                               distribution=parsed.distribution, parameter=parsed.parameter)
end

function run_chain(rng, parsed, task, weights, start)
    keyword = (
        start=start,
        burn_in_attempts=task.burn_in_attempts,
        thin_attempts=task.thin_attempts,
        samples=task.trace_samples,
    )
    if parsed.algorithm === :accelerated
        return sample_center_height_chain_accelerated(rng, task.L, weights; keyword...)
    elseif parsed.algorithm === :tempered
        return sample_center_height_chain_parallel_tempering(
            rng, task.L, weights; keyword..., betas=parsed.tempering_betas,
            swap_interval_attempts=parsed.swap_interval_attempts)
    end
    return sample_center_height_chain(rng, task.L, weights; keyword...)
end

function trace_rows(parsed, task, environment_id, seeds, replica, start, result)
    return [(
        campaign_id=campaign_id(parsed), phase="pilot", L=task.L,
        environment_id=environment_id, environment_seed=seeds.environment,
        replica=replica, start=String(start),
        replica_seed=replica === 1 ? seeds.maximum : seeds.minimum,
        draw_index=index,
        attempt_time=task.burn_in_attempts + index * task.thin_attempts,
        center_height=value,
    ) for (index, value) in enumerate(result.heights)]
end

function diagnostic_row(parsed, task, environment_id, seeds, replica, start, result)
    changed_rate = parsed.algorithm === :tempered ?
                   result.diagnostics.target_changed_rate : result.diagnostics.changed_rate
    auxiliary_rate = parsed.algorithm === :tempered ?
                     result.diagnostics.swap_acceptance_rate :
                     get(result.diagnostics, :final_total_change_rate, NaN)
    return (
        campaign_id=campaign_id(parsed), phase="pilot", L=task.L,
        environment_id=environment_id, environment_seed=seeds.environment,
        replica=replica, start=String(start),
        replica_seed=replica === 1 ? seeds.maximum : seeds.minimum,
        mean=result.mean, variance=result.variance,
        integrated_autocorrelation_time=result.diagnostics.integrated_autocorrelation_time,
        effective_sample_size=result.diagnostics.effective_sample_size,
        changed_rate=changed_rate,
        final_total_change_rate=auxiliary_rate,
        burn_in_attempts=task.burn_in_attempts, thin_attempts=task.thin_attempts,
        trace_samples=task.trace_samples,
    )
end

function batch_paths(parsed, task)
    directory = joinpath(parsed.output_dir, @sprintf("L_%04d", task.L))
    stem = @sprintf("batch_%04d", task.batch_id)
    return (trace=joinpath(directory, stem * ".csv"),
            diagnostics=joinpath(directory, "diagnostic_$(lpad(task.batch_id, 4, '0')).csv"),
            execution=joinpath(directory, "execution_$(lpad(task.batch_id, 4, '0')).txt"))
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

function valid_existing_batch(paths, task)
    isfile(paths.trace) && isfile(paths.diagnostics) || return false
    trace_lines = readlines(paths.trace)
    diagnostic_lines = readlines(paths.diagnostics)
    environments = task.last_id - task.first_id + 1
    length(trace_lines) == 1 + 2 * environments * task.trace_samples || return false
    length(diagnostic_lines) == 1 + 2 * environments || return false
    return strip(first(trace_lines)) == TRACE_HEADER && strip(first(diagnostic_lines)) == DIAGNOSTIC_HEADER
end

function write_batch(paths, traces, diagnostics)
    write_atomic(paths.trace) do io
        println(io, TRACE_HEADER)
        for rows in traces, row in rows
            println(io, join((row.campaign_id, row.phase, row.L, row.environment_id,
                              row.environment_seed, row.replica, row.start, row.replica_seed,
                              row.draw_index, row.attempt_time, row.center_height), ','))
        end
    end
    write_atomic(paths.diagnostics) do io
        println(io, DIAGNOSTIC_HEADER)
        for rows in diagnostics, row in rows
            println(io, join((row.campaign_id, row.phase, row.L, row.environment_id,
                              row.environment_seed, row.replica, row.start, row.replica_seed,
                              row.mean, row.variance, row.integrated_autocorrelation_time,
                              row.effective_sample_size, row.changed_rate,
                              row.final_total_change_rate, row.burn_in_attempts,
                              row.thin_attempts, row.trace_samples), ','))
        end
    end
end

function file_hash(path)
    return isfile(path) ? bytes2hex(sha256(read(path))) : "unavailable"
end

function metadata_text(parsed, rows, tasks)
    repository_root = normpath(joinpath(@__DIR__, "..", ".."))
    schedule = join((join((row.L, row.first_environment_id, row.environments, row.batch_size,
                           row.burn_in_attempts, row.thin_attempts, row.trace_samples), ':')
                     for row in rows), ',')
    return join([
        "schema_version=1",
        "campaign_id=$(campaign_id(parsed))",
        "phase=pilot",
        "geometry=square_grid_direct_weighted_dimer",
        "observable=central_face_height_trace",
        "boundary=Sunil_supplied_tileable_extremal_height_boundary",
        "kernel=literal_random_face_heat_bath",
        "accelerator=$(parsed.algorithm)",
        "accelerator_note=$(parsed.algorithm === :accelerated ? "exactly skips self-loops at fixed attempted-update time" : parsed.algorithm === :tempered ? "parallel tempering with exact beta=1 target; final_total_change_rate column stores swap acceptance" : "literal random-face updates")",
        "tempering_betas=$(join(parsed.tempering_betas, ','))",
        "tempering_swap_interval_attempts=$(parsed.swap_interval_attempts)",
        "pairing=two independent chains in one frozen edge environment; starts=max,min",
        "distribution=$(parsed.distribution)",
        "parameter=$(parsed.parameter)",
        "parameterization=$(parsed.distribution === :gamma ? "Gamma(shape=k,scale=1/k), mean one" : parsed.distribution === :lognormal ? "LogNormal(-sigma^2/2,sigma), mean one" : parsed.distribution === :uniform ? "Uniform(0,2), mean one" : "all weights one")",
        "base_seed=$(parsed.base_seed)",
        "config=$(parsed.config)",
        "config_sha256=$(file_hash(parsed.config))",
        "project_sha256=$(file_hash(joinpath(repository_root, "aztec", "Project.toml")))",
        "manifest_sha256=$(file_hash(joinpath(repository_root, "aztec", "Manifest.toml")))",
        "sampler_sha256=$(file_hash(joinpath(repository_root, "aztec", "src", "GlauberSquareGrid.jl")))",
        "runner_sha256=$(file_hash(abspath(@__FILE__)))",
        "schedule=$schedule",
        "expanded_tasks=$(length(tasks))",
        "julia_version=$(VERSION)",
        "threads=$(nthreads())",
        "seed_identity=base_seed+L+environment_id with separate environment/max/min namespaces",
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
    if valid_existing_batch(paths, task)
        println("skip existing $(paths.trace)")
        return
    end
    (isfile(paths.trace) || isfile(paths.diagnostics)) && error(
        "incomplete or inconsistent existing output for task $(task.task_id)")
    count = task.last_id - task.first_id + 1
    traces = Vector{Vector{NamedTuple}}(undef, count)
    diagnostics = Vector{Vector{NamedTuple}}(undef, count)
    started = time()
    @threads :dynamic for offset in 1:count
        environment_id = task.first_id + offset - 1
        seeds = derive_seeds(parsed.base_seed, task.L, environment_id)
        weights = weights_for_environment(parsed, task.L, seeds.environment)
        maximum = run_chain(Random.Xoshiro(seeds.maximum), parsed, task, weights, :max)
        minimum = run_chain(Random.Xoshiro(seeds.minimum), parsed, task, weights, :min)
        traces[offset] = vcat(
            trace_rows(parsed, task, environment_id, seeds, 1, :max, maximum),
            trace_rows(parsed, task, environment_id, seeds, 2, :min, minimum),
        )
        diagnostics[offset] = [
            diagnostic_row(parsed, task, environment_id, seeds, 1, :max, maximum),
            diagnostic_row(parsed, task, environment_id, seeds, 2, :min, minimum),
        ]
    end
    write_batch(paths, traces, diagnostics)
    write_atomic(paths.execution) do io
        println(io, "task_id=$(task.task_id)")
        println(io, "elapsed_seconds=$(time() - started)")
        println(io, "hostname=$(get(ENV, "HOSTNAME", "unknown"))")
        println(io, "slurm_job_id=$(get(ENV, "SLURM_JOB_ID", "not_applicable"))")
        println(io, "slurm_array_task_id=$(get(ENV, "SLURM_ARRAY_TASK_ID", "not_applicable"))")
        println(io, "threads=$(nthreads())")
        println(io, "julia_version=$(VERSION)")
    end
    @printf("completed Glauber pilot L=%d task=%d environments=%d:%d in %.2fs\n",
            task.L, task.task_id, task.first_id, task.last_id, time() - started)
end

function main(arguments)
    parsed = parse_arguments(arguments)
    parsed.help && (print_help(); return)
    rows = read_config(parsed.config)
    tasks = expand_tasks(rows)
    if parsed.list_tasks
        println("task_id,L,batch_id,first_environment_id,last_environment_id")
        for task in tasks
            println(join((task.task_id, task.L, task.batch_id, task.first_id, task.last_id), ','))
        end
        return
    end
    parsed.task_id <= length(tasks) || error("--task-id exceeds $(length(tasks))")
    ensure_metadata(parsed, rows, tasks)
    selected = parsed.task_id == 0 ? tasks : [tasks[parsed.task_id]]
    foreach(task -> run_task(parsed, task), selected)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
