#!/usr/bin/env julia

using AztecDiamond
using Base.Threads
using Printf

const HEADER =
    "order,sample_id,seed,center_row,center_column,height_1,height_2,height_difference"

# One CSV row is one disorder environment.  The two heights in that row share
# only that environment; all creation coins are drawn independently.

function print_help()
    println("""
    Run a resumable Gamma-disordered double-dimer center-height campaign.

    Each observation uses one fresh Gamma environment and two independent
    tilings conditional on that shared environment.

    Usage:
      JULIA_NUM_THREADS=4 julia --project=aztec \\
        aztec/scripts/run_double_dimer_campaign.jl [options]

    Options:
      --config PATH       CSV sample schedule
      --output-dir PATH   directory for atomic batch CSVs
      --base-seed UINT    campaign seed (default: 20260802)
      --alpha FLOAT       Gamma shape for a weights (default: 0.2)
      --beta FLOAT        Gamma shape for b weights (default: 0.25)
      -h, --help          show this message
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict{String,String}(
        "config" => joinpath(@__DIR__, "..", "configs", "double_dimer_smoke.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "double_dimer_smoke"),
        "base-seed" => "20260802",
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
    return (
        config=abspath(options["config"]),
        output_dir=abspath(options["output-dir"]),
        base_seed=parse(UInt64, options["base-seed"]),
        alpha=parse(Float64, options["alpha"]),
        beta=parse(Float64, options["beta"]),
    )
end

function read_config(path)
    lines = readlines(path)
    isempty(lines) && error("empty config: $path")
    header = strip(first(lines))
    legacy = header == "order,samples,batch_size"
    extended = header == "order,first_sample_id,samples,batch_size"
    legacy || extended || error("unexpected config header in $path: $header")

    rows = NamedTuple{
        (:order, :first_sample_id, :samples, :batch_size),
        Tuple{Int,Int,Int,Int},
    }[]
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        row = if legacy
            length(fields) == 3 || error("invalid config row $line_number")
            (
                order=parse(Int, fields[1]),
                first_sample_id=1,
                samples=parse(Int, fields[2]),
                batch_size=parse(Int, fields[3]),
            )
        else
            length(fields) == 4 || error("invalid config row $line_number")
            (
                order=parse(Int, fields[1]),
                first_sample_id=parse(Int, fields[2]),
                samples=parse(Int, fields[3]),
                batch_size=parse(Int, fields[4]),
            )
        end
        row.order > 0 || error("order must be positive on row $line_number")
        row.first_sample_id > 0 || error("first sample ID must be positive")
        row.samples > 0 || error("samples must be positive on row $line_number")
        row.batch_size > 0 || error("batch size must be positive on row $line_number")
        push!(rows, row)
    end
    isempty(rows) && error("config has no data rows")
    orders = [row.order for row in rows]
    length(unique(orders)) == length(orders) ||
        error("config must contain at most one row per order")
    return rows
end

function splitmix64(value::UInt64)
    # Deterministically hash public sample identifiers into well-separated
    # Xoshiro seeds.  This makes batching, threading, and resumption irrelevant
    # to the random stream assigned to a sample.
    z = value + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

function sample_seed(base_seed::UInt64, order::Int, sample_id::Int)
    key =
        base_seed ⊻
        (UInt64(order) * 0xd6e8feb86659fd93) ⊻
        (UInt64(sample_id) * 0xa5a3564e27f8862f)
    return splitmix64(key)
end

function batch_path(output_dir, order, batch_id)
    order_dir = joinpath(output_dir, @sprintf("L_%04d", order))
    return joinpath(order_dir, @sprintf("batch_%04d.csv", batch_id))
end

function valid_existing_batch(
    path,
    base_seed,
    expected_order,
    first_sample,
    last_sample,
)
    isfile(path) || return false
    lines = readlines(path)
    length(lines) == last_sample - first_sample + 2 || return false
    strip(first(lines)) == HEADER || return false
    for (line, expected_id) in zip(lines[2:end], first_sample:last_sample)
        fields = split(strip(line), ',')
        length(fields) == 8 || return false
        try
            parse(Int, fields[1]) == expected_order || return false
            parse(Int, fields[2]) == expected_id || return false
            parse(UInt64, fields[3]) ==
                sample_seed(base_seed, expected_order, expected_id) || return false
            parse(Int, fields[4]) == expected_order + 1 || return false
            parse(Int, fields[5]) == fld(expected_order, 2) + 1 || return false
            height_1 = parse(Int, fields[6])
            height_2 = parse(Int, fields[7])
            parse(Int, fields[8]) == height_1 - height_2 || return false
        catch
            return false
        end
    end
    return true
end

function write_batch(path, results)
    mkpath(dirname(path))
    temporary_path = path * ".tmp"
    try
        open(temporary_path, "w") do io
            println(io, HEADER)
            for result in results
                println(
                    io,
                    result.order, ',', result.sample_id, ',', result.seed, ',',
                    result.center_row, ',', result.center_column, ',',
                    result.height_1, ',', result.height_2, ',', result.height_difference,
                )
            end
        end
        # Atomic rename prevents an interrupted write from looking complete.
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function write_metadata(path, parsed, config_rows)
    temporary_path = path * ".tmp"
    try
        open(temporary_path, "w") do io
            println(io, "model=double dimer in biased Gamma-disordered Aztec diamond")
            println(io, "observable=center-height difference H1-H2")
            println(
                io,
                "one_sample=fresh environment plus two conditionally independent tilings",
            )
            println(io, "variance_identity=Var(H1-H2)/2=E[Var(H|environment)]")
            println(io, "covariance_identity=Cov(H1,H2)=Var(E[H|environment])")
            println(io, "center_index=(L+1, floor(L/2)+1)")
            println(io, "alpha=$(parsed.alpha)")
            println(io, "beta=$(parsed.beta)")
            println(io, "gamma_scale=1")
            println(io, "base_seed=$(parsed.base_seed)")
            println(io, "rng=Xoshiro with SplitMix64-derived per-sample seeds")
            println(io, "julia_version=$(VERSION)")
            println(io, "julia_threads=$(nthreads())")
            println(io, "config=$(parsed.config)")
            println(io, "orders=$(join((row.order for row in config_rows), ','))")
            println(io, "planned_pairs=$(sum(row.samples for row in config_rows))")
        end
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
end

function run_batch(parsed, order, batch_id, first_sample, last_sample)
    output_path = batch_path(parsed.output_dir, order, batch_id)
    if valid_existing_batch(
        output_path,
        parsed.base_seed,
        order,
        first_sample,
        last_sample,
    )
        println("skip existing $output_path")
        return
    end
    isfile(output_path) && error(
        "existing batch is incomplete or invalid: $output_path. " *
        "Move only that file before rerunning.",
    )

    count_samples = last_sample - first_sample + 1
    center = center_face_index(order)
    results = Vector{
        NamedTuple{
            (:order, :sample_id, :seed, :center_row, :center_column,
             :height_1, :height_2, :height_difference),
            Tuple{Int,Int,UInt64,Int,Int,Int,Int,Int},
        },
    }(undef, count_samples)

    started = time()
    # Results are written by array index, so thread scheduling cannot reorder
    # rows.  The deterministic per-sample seed makes thread count irrelevant.
    @threads :dynamic for offset in 1:count_samples
        sample_id = first_sample + offset - 1
        seed = sample_seed(parsed.base_seed, order, sample_id)
        pair = sample_gamma_center_height_pair(
            seed,
            order;
            alpha=parsed.alpha,
            beta=parsed.beta,
        )
        results[offset] = (
            order=order,
            sample_id=sample_id,
            seed=seed,
            center_row=center.row,
            center_column=center.column,
            height_1=pair.height_1,
            height_2=pair.height_2,
            height_difference=pair.difference,
        )
    end
    write_batch(output_path, results)
    @printf(
        "completed double L=%d batch=%d samples=%d:%d in %.2fs\n",
        order,
        batch_id,
        first_sample,
        last_sample,
        time() - started,
    )
    GC.gc()
end

function main(arguments)
    parsed = parse_arguments(arguments)
    if isnothing(parsed)
        print_help()
        return
    end
    parsed.alpha > 0 || error("--alpha must be positive")
    parsed.beta > 0 || error("--beta must be positive")
    config_rows = read_config(parsed.config)
    mkpath(parsed.output_dir)
    println(
        "Starting double-dimer campaign with $(nthreads()) threads and ",
        "$(sum(row.samples for row in config_rows)) planned pairs",
    )

    for row in config_rows
        batch_id = 1
        first_sample = row.first_sample_id
        final_sample = row.first_sample_id + row.samples - 1
        while first_sample <= final_sample
            last_sample = min(first_sample + row.batch_size - 1, final_sample)
            run_batch(parsed, row.order, batch_id, first_sample, last_sample)
            first_sample = last_sample + 1
            batch_id += 1
        end
    end
    write_metadata(
        joinpath(parsed.output_dir, "campaign_metadata.txt"),
        parsed,
        config_rows,
    )
    println("Double-dimer campaign complete: $(parsed.output_dir)")
end

main(ARGS)
