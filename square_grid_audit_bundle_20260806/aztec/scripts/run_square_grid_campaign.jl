#!/usr/bin/env julia

using AztecDiamond
using AztecDiamond.SquareGrid
using Base.Threads
using Printf

# This schema intentionally matches analyze_spatial_campaign.jl. The first
# field is an analysis role: `gamma` for a random Gamma environment and
# `uniform` for the no-disorder baseline. Exact square-grid model semantics are
# frozen in campaign_metadata.txt and in the diagnostic files.
const HEADER =
    "model,order,sample_id,seed,fraction_num,fraction_den,separation," *
    "left_column,right_column,increment_1,increment_2,difference"
const DIAGNOSTIC_HEADER =
    "order,sample_id,sample_seed,environment_seed,replica_1_seed,replica_2_seed," *
    "raw_steps_1,raw_steps_2,branches_1,branches_2"

function print_help()
    println("""
    Run a resumable paired square-grid Temperley-dimer height campaign.

    Usage:
      JULIA_NUM_THREADS=8 julia --project=aztec \\
        aztec/scripts/run_square_grid_campaign.jl [options]

    Options:
      --environment-model NAME   baseline, directed_site_iid, or
                                 undirected_conductance
                                 (default: directed_site_iid)
      --distribution NAME        gamma (default: gamma)
      --parameter FLOAT          mean-one Gamma shape (default: 0.5)
      --config PATH              campaign CSV
      --output-dir PATH          atomic batch output directory
      --fractions LIST           rational fractions of the full side 2L
                                 (default: 1/32,1/16,1/8,1/4)
      --base-seed UINT           public campaign seed
      --task-id INT              run one expanded batch task; 0 runs all
      --max-steps-per-walk INT   optional Wilson safety cap
      --list-tasks               print the expanded task table and exit
      -h, --help                 show this message

    The output CSV is compatible with analyze_spatial_campaign.jl. Run one
    baseline campaign and one disordered campaign in separate directories.
    """)
end

function parse_fraction(value)
    fields = split(strip(value), '/')
    length(fields) == 2 || error("fraction must have NUM/DEN form: $value")
    numerator = parse(Int, fields[1])
    denominator = parse(Int, fields[2])
    0 < numerator < denominator || error("fraction must lie in (0,1): $value")
    divisor = gcd(numerator, denominator)
    return (num=numerator ÷ divisor, den=denominator ÷ divisor)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) &&
        return (help=true,)
    options = Dict{String,String}(
        "environment-model" => "directed_site_iid",
        "distribution" => "gamma",
        "parameter" => "0.5",
        "config" => joinpath(@__DIR__, "..", "configs", "square_grid_smoke.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "square_grid_directed_smoke"),
        "fractions" => "1/32,1/16,1/8,1/4",
        "base-seed" => "20260805",
        "task-id" => "0",
        "max-steps-per-walk" => "",
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
    environment_model = Symbol(lowercase(options["environment-model"]))
    environment_model in (:baseline, :directed_site_iid, :undirected_conductance) ||
        error("invalid --environment-model")
    distribution = Symbol(lowercase(options["distribution"]))
    distribution == :gamma || error("the current campaign runner supports gamma only")
    parameter = parse(Float64, options["parameter"])
    parameter > 0 || error("--parameter must be positive")
    fractions = parse_fraction.(split(options["fractions"], ','))
    length(unique(fractions)) == length(fractions) || error("fractions must be unique")
    sort!(fractions; by=fraction -> fraction.num / fraction.den)
    task_id = parse(Int, options["task-id"])
    task_id >= 0 || error("--task-id must be nonnegative")
    max_steps = isempty(options["max-steps-per-walk"]) ? nothing :
                parse(Int, options["max-steps-per-walk"])
    !isnothing(max_steps) && max_steps <= 0 && error("Wilson cap must be positive")
    return (
        help=false,
        environment_model=environment_model,
        distribution=distribution,
        parameter=parameter,
        config=abspath(options["config"]),
        output_dir=abspath(options["output-dir"]),
        fractions=fractions,
        base_seed=parse(UInt64, options["base-seed"]),
        task_id=task_id,
        max_steps_per_walk=max_steps,
        list_tasks=list_tasks,
    )
end

function read_config(path)
    lines = readlines(path)
    isempty(lines) && error("empty config: $path")
    strip(first(lines)) == "order,first_sample_id,samples,batch_size" ||
        error("unexpected config header in $path")
    rows = NamedTuple[]
    for (offset, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        length(fields) == 4 || error("invalid config row $(offset + 1)")
        row = (
            order=parse(Int, fields[1]),
            first_sample_id=parse(Int, fields[2]),
            samples=parse(Int, fields[3]),
            batch_size=parse(Int, fields[4]),
        )
        row.order >= 2 || error("order L must be at least two")
        row.first_sample_id > 0 || error("first sample ID must be positive")
        row.samples > 0 || error("samples must be positive")
        row.batch_size > 0 || error("batch size must be positive")
        push!(rows, row)
    end
    isempty(rows) && error("config has no data rows")
    orders = [row.order for row in rows]
    length(unique(orders)) == length(orders) || error("orders must be unique")
    return rows
end

function expand_tasks(config_rows)
    tasks = NamedTuple[]
    for row in config_rows
        first_sample = row.first_sample_id
        last_config_sample = row.first_sample_id + row.samples - 1
        batch_id = 0
        while first_sample <= last_config_sample
            batch_id += 1
            last_sample = min(first_sample + row.batch_size - 1, last_config_sample)
            push!(tasks, (
                task_id=length(tasks) + 1,
                order=row.order,
                batch_id=batch_id,
                first_sample=first_sample,
                last_sample=last_sample,
            ))
            first_sample = last_sample + 1
        end
    end
    return tasks
end

function splitmix64(value::UInt64)
    value += 0x9e3779b97f4a7c15
    value = (value ⊻ (value >> 30)) * 0xbf58476d1ce4e5b9
    value = (value ⊻ (value >> 27)) * 0x94d049bb133111eb
    return value ⊻ (value >> 31)
end

function model_salt(model::Symbol)
    model === :baseline && return UInt64(0x626173656c696e65)
    model === :directed_site_iid && return UInt64(0x6469726563746564)
    return UInt64(0x756e646972656374)
end

function sample_seed(parsed, order, sample_id)
    parameter_bits = reinterpret(UInt64, parsed.parameter)
    key = parsed.base_seed ⊻ model_salt(parsed.environment_model) ⊻
          parameter_bits ⊻
          (UInt64(order) * 0xd6e8feb86659fd93) ⊻
          (UInt64(sample_id) * 0xa5a3564e27f8862f)
    return splitmix64(key)
end

# Fractions refer to the full square side 2L. This gives distinct probe
# separations 1,2,4,8 already at L=16.
function separations_for_order(order, fractions)
    separations = [
        clamp(round(Int, 2 * order * fraction.num / fraction.den), 1, order - 1)
        for fraction in fractions
    ]
    length(unique(separations)) == length(separations) ||
        error("fraction separations collide at L=$order")
    return separations
end

analysis_model(parsed) = parsed.environment_model === :baseline ? "uniform" : "gamma"

function batch_path(output_dir, order, batch_id)
    order_dir = joinpath(output_dir, @sprintf("L_%04d", order))
    return joinpath(order_dir, @sprintf("batch_%04d.csv", batch_id))
end

diagnostic_path(output_dir, order, batch_id) =
    joinpath(dirname(batch_path(output_dir, order, batch_id)),
             @sprintf("diagnostic_%04d.csv", batch_id))

function valid_existing_batch(data_path, diagnostics_path, parsed, task)
    isfile(data_path) && isfile(diagnostics_path) || return false
    lines = readlines(data_path)
    diagnostic_lines = readlines(diagnostics_path)
    sample_count = task.last_sample - task.first_sample + 1
    length(lines) == sample_count * length(parsed.fractions) + 1 || return false
    length(diagnostic_lines) == sample_count + 1 || return false
    strip(first(lines)) == HEADER || return false
    strip(first(diagnostic_lines)) == DIAGNOSTIC_HEADER || return false
    line_index = 2
    separations = separations_for_order(task.order, parsed.fractions)
    expected_model = analysis_model(parsed)
    for sample_id in task.first_sample:task.last_sample
        expected_seed = sample_seed(parsed, task.order, sample_id)
        for (fraction, separation) in zip(parsed.fractions, separations)
            fields = split(strip(lines[line_index]), ',')
            length(fields) == 12 || return false
            try
                fields[1] == expected_model || return false
                parse(Int, fields[2]) == task.order || return false
                parse(Int, fields[3]) == sample_id || return false
                parse(UInt64, fields[4]) == expected_seed || return false
                parse(Int, fields[5]) == fraction.num || return false
                parse(Int, fields[6]) == fraction.den || return false
                parse(Int, fields[7]) == separation || return false
                parse(Int, fields[9]) - parse(Int, fields[8]) == separation || return false
                first = parse(Int, fields[10])
                second = parse(Int, fields[11])
                parse(Int, fields[12]) == first - second || return false
            catch
                return false
            end
            line_index += 1
        end
    end
    return true
end

function write_atomic(writer, path)
    mkpath(dirname(path))
    temporary_path = path * ".$(getpid()).tmp"
    try
        open(temporary_path, "w") do io
            writer(io)
        end
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function write_batch(data_path, diagnostics_path, rows_by_sample, diagnostics)
    write_atomic(data_path) do io
        println(io, HEADER)
        for rows in rows_by_sample, row in rows
            println(
                io,
                row.model, ',', row.order, ',', row.sample_id, ',', row.seed, ',',
                row.fraction_num, ',', row.fraction_den, ',', row.separation, ',',
                row.left_x, ',', row.right_x, ',', row.increment_1, ',',
                row.increment_2, ',', row.difference,
            )
        end
    end
    write_atomic(diagnostics_path) do io
        println(io, DIAGNOSTIC_HEADER)
        for row in diagnostics
            println(
                io,
                row.order, ',', row.sample_id, ',', row.sample_seed, ',',
                row.environment_seed, ',', row.replica_1_seed, ',', row.replica_2_seed, ',',
                row.raw_steps_1, ',', row.raw_steps_2, ',',
                row.branches_1, ',', row.branches_2,
            )
        end
    end
end

function run_task(parsed, task)
    data_path = batch_path(parsed.output_dir, task.order, task.batch_id)
    diagnostics_path = diagnostic_path(parsed.output_dir, task.order, task.batch_id)
    if valid_existing_batch(data_path, diagnostics_path, parsed, task)
        println("skip existing $data_path")
        return
    end
    if isfile(data_path) || isfile(diagnostics_path)
        error("existing square-grid batch is incomplete or inconsistent: $data_path")
    end

    separations = separations_for_order(task.order, parsed.fractions)
    sample_count = task.last_sample - task.first_sample + 1
    rows_by_sample = Vector{Vector{NamedTuple}}(undef, sample_count)
    diagnostics = Vector{NamedTuple}(undef, sample_count)
    started = time()
    @threads :dynamic for offset in 1:sample_count
        sample_id = task.first_sample + offset - 1
        seed = sample_seed(parsed, task.order, sample_id)
        result = sample_spatial_increment_pair(
            seed,
            task.order,
            separations;
            environment_model=parsed.environment_model,
            distribution=parsed.distribution,
            parameter=parsed.parameter,
            max_steps_per_walk=parsed.max_steps_per_walk,
        )
        rows_by_sample[offset] = [
            (
                model=analysis_model(parsed),
                order=task.order,
                sample_id=sample_id,
                seed=seed,
                fraction_num=fraction.num,
                fraction_den=fraction.den,
                separation=row.separation,
                left_x=row.left_x,
                right_x=row.right_x,
                increment_1=row.increment_1,
                increment_2=row.increment_2,
                difference=row.difference,
            )
            for (fraction, row) in zip(parsed.fractions, result.rows)
        ]
        diagnostics[offset] = merge(
            (order=task.order, sample_id=sample_id, sample_seed=seed),
            result.diagnostics,
        )
    end
    write_batch(data_path, diagnostics_path, rows_by_sample, diagnostics)
    @printf(
        "completed square-grid %s L=%d task=%d batch=%d samples=%d:%d in %.2fs\n",
        parsed.environment_model,
        task.order,
        task.task_id,
        task.batch_id,
        task.first_sample,
        task.last_sample,
        time() - started,
    )
    GC.gc()
end

function metadata_text(parsed, config_rows, tasks)
    fractions = join((string(f.num, '/', f.den) for f in parsed.fractions), ',')
    return join([
        "geometry=square_grid_temperley",
        "observable=central_horizontal_dimer_height_increment",
        "environment_model=$(parsed.environment_model)",
        "analysis_model=$(analysis_model(parsed))",
        "distribution=$(parsed.distribution)",
        "parameter=$(parsed.parameter)",
        "parameterization=Gamma(shape=k,scale=1/k), mean one",
        "boundary=wired_square",
        "outer_face=lower_left_cell_-L_-L",
        "fractions_of_full_side=$fractions",
        "base_seed=$(parsed.base_seed)",
        "config=$(parsed.config)",
        "config_rows=$(length(config_rows))",
        "expanded_tasks=$(length(tasks))",
        "julia_version=$(VERSION)",
        "threads=$(nthreads())",
        "height_note=exact Temperley matching flow across central cut; deterministic reference flow cancels",
    ], '\n') * "\n"
end

function ensure_metadata(parsed, config_rows, tasks)
    mkpath(parsed.output_dir)
    path = joinpath(parsed.output_dir, "campaign_metadata.txt")
    expected = metadata_text(parsed, config_rows, tasks)
    if isfile(path)
        read(path, String) == expected || error("existing campaign metadata differs: $path")
        return
    end
    write_atomic(path) do io
        print(io, expected)
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    parsed.help && (print_help(); return)
    config_rows = read_config(parsed.config)
    tasks = expand_tasks(config_rows)
    if parsed.list_tasks
        println("task_id,order,batch_id,first_sample,last_sample")
        for task in tasks
            println(
                task.task_id, ',', task.order, ',', task.batch_id, ',',
                task.first_sample, ',', task.last_sample)
        end
        return
    end
    parsed.task_id <= length(tasks) ||
        error("--task-id $(parsed.task_id) exceeds $(length(tasks)) tasks")
    ensure_metadata(parsed, config_rows, tasks)
    selected = parsed.task_id == 0 ? tasks : [tasks[parsed.task_id]]
    for task in selected
        run_task(parsed, task)
    end
end

main(ARGS)
