#!/usr/bin/env julia

using AztecDiamond
using Base.Threads
using Printf

const HEADER =
    "model,order,sample_id,seed,fraction_num,fraction_den,separation," *
    "left_column,right_column,increment_1,increment_2,difference"

function print_help()
    println("""
    Run a resumable spatial height-increment experiment.

    For the Gamma model, each row group uses one fresh environment and two
    conditionally independent tilings.  The uniform model uses two independent
    uniform tilings and is the no-disorder control.

    Usage:
      JULIA_NUM_THREADS=8 julia --project=aztec \\
        aztec/scripts/run_spatial_campaign.jl [options]

    Options:
      --model NAME        gamma or uniform (default: gamma)
      --config PATH       CSV schedule
      --output-dir PATH   directory for atomic batch CSVs
      --fractions LIST    comma-separated rational separations
                          (default: 1/32,1/16,1/8,1/4)
      --base-seed UINT    campaign seed
      --alpha FLOAT       Gamma a-shape (default: 0.2)
      --beta FLOAT        Gamma b-shape (default: 0.25)
      -h, --help          show this message
    """)
end

function parse_fraction(value)
    fields = split(strip(value), '/')
    length(fields) == 2 || error("fraction must have NUM/DEN form: $value")
    numerator = parse(Int, fields[1])
    denominator = parse(Int, fields[2])
    0 < numerator < denominator || error("fraction must lie strictly in (0,1): $value")
    divisor = gcd(numerator, denominator)
    return (num=numerator ÷ divisor, den=denominator ÷ divisor)
end

function parse_arguments(arguments)
    if any(argument -> argument in ("-h", "--help"), arguments)
        print_help()
        return nothing
    end
    options = Dict{String,String}(
        "model" => "gamma",
        "config" => joinpath(@__DIR__, "..", "configs", "spatial_smoke.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "spatial_gamma_smoke"),
        "fractions" => "1/32,1/16,1/8,1/4",
        "base-seed" => "20260803",
        "alpha" => "0.2",
        "beta" => "0.25",
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
    model = lowercase(options["model"])
    model in ("gamma", "uniform") || error("--model must be gamma or uniform")
    fractions = parse_fraction.(split(options["fractions"], ','))
    length(unique(fractions)) == length(fractions) || error("fractions must be unique")
    sort!(fractions; by=fraction -> fraction.num / fraction.den)
    return (
        model=model,
        config=abspath(options["config"]),
        output_dir=abspath(options["output-dir"]),
        fractions=fractions,
        base_seed=parse(UInt64, options["base-seed"]),
        alpha=parse(Float64, options["alpha"]),
        beta=parse(Float64, options["beta"]),
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
        row.order > 1 || error("order must exceed one")
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

function splitmix64(value::UInt64)
    z = value + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

function sample_seed(base_seed, model, order, sample_id)
    model_salt = model == "gamma" ? 0x67616d6d615f7370 : 0x756e69666f726d73
    key =
        base_seed ⊻ model_salt ⊻
        (UInt64(order) * 0xd6e8feb86659fd93) ⊻
        (UInt64(sample_id) * 0xa5a3564e27f8862f)
    return splitmix64(key)
end

function separations_for_order(order, fractions)
    separations = [
        clamp(round(Int, order * fraction.num / fraction.den), 1, order - 1)
        for fraction in fractions
    ]
    length(unique(separations)) == length(separations) ||
        error("fraction separations collide at order $order")
    return separations
end

function batch_path(output_dir, order, batch_id)
    order_dir = joinpath(output_dir, @sprintf("L_%04d", order))
    return joinpath(order_dir, @sprintf("batch_%04d.csv", batch_id))
end

function valid_existing_batch(path, parsed, order, first_sample, last_sample)
    isfile(path) || return false
    lines = readlines(path)
    expected_rows = (last_sample - first_sample + 1) * length(parsed.fractions)
    length(lines) == expected_rows + 1 || return false
    strip(first(lines)) == HEADER || return false
    line_index = 2
    separations = separations_for_order(order, parsed.fractions)
    for sample_id in first_sample:last_sample
        expected_seed = sample_seed(parsed.base_seed, parsed.model, order, sample_id)
        for (fraction, separation) in zip(parsed.fractions, separations)
            fields = split(strip(lines[line_index]), ',')
            length(fields) == 12 || return false
            try
                fields[1] == parsed.model || return false
                parse(Int, fields[2]) == order || return false
                parse(Int, fields[3]) == sample_id || return false
                parse(UInt64, fields[4]) == expected_seed || return false
                parse(Int, fields[5]) == fraction.num || return false
                parse(Int, fields[6]) == fraction.den || return false
                parse(Int, fields[7]) == separation || return false
                increment_1 = parse(Int, fields[10])
                increment_2 = parse(Int, fields[11])
                parse(Int, fields[12]) == increment_1 - increment_2 || return false
            catch
                return false
            end
            line_index += 1
        end
    end
    return true
end

function write_batch(path, rows_by_sample)
    mkpath(dirname(path))
    temporary_path = path * ".tmp"
    try
        open(temporary_path, "w") do io
            println(io, HEADER)
            for rows in rows_by_sample, row in rows
                println(
                    io,
                    row.model, ',', row.order, ',', row.sample_id, ',', row.seed, ',',
                    row.fraction_num, ',', row.fraction_den, ',', row.separation, ',',
                    row.left_column, ',', row.right_column, ',', row.increment_1, ',',
                    row.increment_2, ',', row.difference,
                )
            end
        end
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function run_batch(parsed, order, batch_id, first_sample, last_sample)
    output_path = batch_path(parsed.output_dir, order, batch_id)
    if valid_existing_batch(output_path, parsed, order, first_sample, last_sample)
        println("skip existing $output_path")
        return
    end
    isfile(output_path) && error(
        "existing batch is incomplete or inconsistent: $output_path. " *
        "Move that file before resuming.",
    )

    separations = separations_for_order(order, parsed.fractions)
    sample_count = last_sample - first_sample + 1
    rows_by_sample = Vector{Vector{NamedTuple}}(undef, sample_count)
    started = time()
    @threads :dynamic for offset in 1:sample_count
        sample_id = first_sample + offset - 1
        seed = sample_seed(parsed.base_seed, parsed.model, order, sample_id)
        increments = if parsed.model == "gamma"
            sample_gamma_spatial_increment_pair(
                seed,
                order,
                separations;
                alpha=parsed.alpha,
                beta=parsed.beta,
            )
        else
            sample_uniform_spatial_increment_pair(seed, order, separations)
        end
        rows_by_sample[offset] = [
            (
                model=parsed.model,
                order=order,
                sample_id=sample_id,
                seed=seed,
                fraction_num=fraction.num,
                fraction_den=fraction.den,
                separation=result.separation,
                left_column=result.left_column,
                right_column=result.right_column,
                increment_1=result.increment_1,
                increment_2=result.increment_2,
                difference=result.difference,
            )
            for (fraction, result) in zip(parsed.fractions, increments)
        ]
    end
    write_batch(output_path, rows_by_sample)
    @printf(
        "completed %s spatial L=%d batch=%d samples=%d:%d in %.2fs\n",
        parsed.model,
        order,
        batch_id,
        first_sample,
        last_sample,
        time() - started,
    )
    GC.gc()
end

function write_metadata(path, parsed, config_rows)
    temporary_path = path * ".tmp"
    shared_environment = parsed.model == "gamma"
    fraction_text = join((string(f.num, '/', f.den) for f in parsed.fractions), ',')
    try
        open(temporary_path, "w") do io
            println(io, "model=$(parsed.model)")
            println(io, "observable=symmetric central-row spatial height increments")
            println(io, "replicas=2")
            println(io, "gamma_shared_environment=$shared_environment")
            println(io, "fractions=$fraction_text")
            println(io, "alpha=$(parsed.alpha)")
            println(io, "beta=$(parsed.beta)")
            println(io, "base_seed=$(parsed.base_seed)")
            println(io, "rng=Xoshiro with SplitMix64-derived per-sample seeds")
            println(io, "julia_version=$(VERSION)")
            println(io, "julia_threads=$(nthreads())")
            println(io, "config=$(parsed.config)")
            println(io, "orders=$(join((row.order for row in config_rows), ','))")
            println(io, "planned_environments=$(sum(row.samples for row in config_rows))")
        end
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    isnothing(parsed) && return
    parsed.alpha > 0 || error("--alpha must be positive")
    parsed.beta > 0 || error("--beta must be positive")
    config_rows = read_config(parsed.config)
    for row in config_rows
        separations_for_order(row.order, parsed.fractions)
    end
    mkpath(parsed.output_dir)
    println(
        "Starting $(parsed.model) spatial campaign with $(nthreads()) threads and ",
        "$(sum(row.samples for row in config_rows)) planned samples",
    )
    for row in config_rows
        first_sample = row.first_sample_id
        final_sample = row.first_sample_id + row.samples - 1
        batch_id = 1
        while first_sample <= final_sample
            last_sample = min(first_sample + row.batch_size - 1, final_sample)
            run_batch(parsed, row.order, batch_id, first_sample, last_sample)
            first_sample = last_sample + 1
            batch_id += 1
        end
    end
    write_metadata(joinpath(parsed.output_dir, "campaign_metadata.txt"), parsed, config_rows)
    println("Spatial campaign complete: $(parsed.output_dir)")
end

main(ARGS)
