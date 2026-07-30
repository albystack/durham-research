#!/usr/bin/env julia

using Random
using Printf
using Base.Threads

include(joinpath(@__DIR__, "..", "src", "AztecDiamond.jl"))
using .AztecDiamond

function parse_arguments(arguments)
    options = Dict{String,String}(
        "config" => joinpath(@__DIR__, "..", "configs", "gamma_height_pilot.csv"),
        "output-dir" => joinpath(@__DIR__, "..", "output", "gamma_height_pilot"),
        "base-seed" => "20260729",
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
    strip(lines[1]) == "order,samples,batch_size" ||
        error("unexpected config header in $path")
    rows = NamedTuple{(:order, :samples, :batch_size),Tuple{Int,Int,Int}}[]
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        length(fields) == 3 || error("invalid config row $line_number")
        row = (
            order=parse(Int, fields[1]),
            samples=parse(Int, fields[2]),
            batch_size=parse(Int, fields[3]),
        )
        row.order > 0 && row.samples > 0 && row.batch_size > 0 ||
            error("config values must be positive on row $line_number")
        push!(rows, row)
    end
    isempty(rows) && error("config has no data rows")
    return rows
end

function splitmix64(value::UInt64)
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

function simulate_center_height(
    order::Int,
    seed::UInt64,
    alpha::Float64,
    beta::Float64,
)
    rng = Xoshiro(seed)
    weights = gamma_disordered_weights(rng, order; alpha=alpha, beta=beta)
    choices = gamma_disordered_creation_choices(rng, weights.a, weights.b)
    weights = nothing
    tiling = sample_tiling_from_choices(choices)
    count(tiling) == order * (order + 1) ||
        error("invalid domino count for order=$order seed=$seed")
    return center_height(tiling)
end

function batch_path(output_dir, order, batch_id)
    order_dir = joinpath(output_dir, @sprintf("L_%04d", order))
    return joinpath(order_dir, @sprintf("batch_%04d.csv", batch_id))
end

function valid_existing_batch(path, expected_rows)
    isfile(path) || return false
    lines = readlines(path)
    length(lines) == expected_rows + 1 || return false
    return strip(first(lines)) ==
           "order,sample_id,seed,center_row,center_column,center_height"
end

function write_batch(path, results)
    mkpath(dirname(path))
    temporary_path = path * ".tmp"
    open(temporary_path, "w") do io
        println(io, "order,sample_id,seed,center_row,center_column,center_height")
        for result in results
            println(
                io,
                result.order,
                ',',
                result.sample_id,
                ',',
                result.seed,
                ',',
                result.center_row,
                ',',
                result.center_column,
                ',',
                result.center_height,
            )
        end
    end
    mv(temporary_path, path; force=true)
end

function write_campaign_metadata(path, parsed, config_rows)
    open(path, "w") do io
        println(io, "model=biased Gamma-disordered Aztec diamond")
        println(io, "observable=center face height")
        println(io, "one_sample=fresh random environment plus one conditional dimer cover")
        println(io, "height_convention=Sunil Chhita email 2026-07-29")
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
        println(io, "planned_samples=$(sum(row.samples for row in config_rows))")
    end
end

function run_batch(parsed, order, batch_id, first_sample, last_sample)
    count_samples = last_sample - first_sample + 1
    output_path = batch_path(parsed.output_dir, order, batch_id)
    if valid_existing_batch(output_path, count_samples)
        println("skip existing $output_path")
        return
    end

    center = center_face_index(order)
    results = Vector{
        NamedTuple{
            (:order, :sample_id, :seed, :center_row, :center_column, :center_height),
            Tuple{Int,Int,UInt64,Int,Int,Int},
        },
    }(undef, count_samples)

    started = time()
    @threads :dynamic for offset in 1:count_samples
        sample_id = first_sample + offset - 1
        seed = sample_seed(parsed.base_seed, order, sample_id)
        height = simulate_center_height(order, seed, parsed.alpha, parsed.beta)
        results[offset] = (
            order=order,
            sample_id=sample_id,
            seed=seed,
            center_row=center.row,
            center_column=center.column,
            center_height=height,
        )
    end
    write_batch(output_path, results)
    elapsed = time() - started
    @printf(
        "completed L=%d batch=%d samples=%d:%d in %.2fs\n",
        order,
        batch_id,
        first_sample,
        last_sample,
        elapsed,
    )
    GC.gc()
end

function main(arguments)
    parsed = parse_arguments(arguments)
    parsed.alpha > 0 || error("--alpha must be positive")
    parsed.beta > 0 || error("--beta must be positive")
    config_rows = read_config(parsed.config)
    mkpath(parsed.output_dir)
    write_campaign_metadata(
        joinpath(parsed.output_dir, "campaign_metadata.txt"),
        parsed,
        config_rows,
    )

    println(
        "Starting Gamma height campaign with $(nthreads()) Julia threads and ",
        "$(sum(row.samples for row in config_rows)) planned samples",
    )
    for row in config_rows
        batch_id = 1
        first_sample = 1
        while first_sample <= row.samples
            last_sample = min(first_sample + row.batch_size - 1, row.samples)
            run_batch(
                parsed,
                row.order,
                batch_id,
                first_sample,
                last_sample,
            )
            first_sample = last_sample + 1
            batch_id += 1
        end
    end
    println("Campaign complete: $(parsed.output_dir)")
end

main(ARGS)
