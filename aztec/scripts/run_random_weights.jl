#!/usr/bin/env julia

using Random
using Printf
using SHA
using Statistics
using AztecDiamond

function print_help()
    println("""
    Generate and validate one illustrated arbitrary-weight Aztec diamond.

    Usage:
      julia --project=aztec aztec/scripts/run_random_weights.jl [options]

    Options:
      --order INT        Aztec-diamond order (default: 200)
      --seed UINT        Xoshiro seed
      --output-dir PATH  output directory
      -h, --help         show this message
    """)
end

function parse_arguments(arguments)
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    options = Dict{String,String}(
        "order" => "200",
        "seed" => "20260728",
        "output-dir" => joinpath(@__DIR__, "..", "output", "random_uniform_n200_seed20260728"),
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
        order=parse(Int, options["order"]),
        seed=parse(UInt64, options["seed"]),
        output_dir=abspath(options["output-dir"]),
    )
end

function file_sha256(path)
    return bytes2hex(open(sha256, path))
end

function write_metadata(
    path,
    arguments,
    weights,
    validation,
    counts,
    weights_path,
    tiling_path,
    figure_path,
    probability_min,
    probability_max,
    elapsed_probability,
    elapsed_sampling,
)
    open(path, "w") do io
        println(io, "model=random-weight Aztec diamond")
        println(io, "algorithm=domino shuffling from supplied simulatorfinal.jl (2022)")
        println(io, "order=$(arguments.order)")
        println(io, "matrix_side=$(2 * arguments.order)")
        println(io, "weight_distribution=iid Uniform(0,1)")
        println(io, "rng=Xoshiro")
        println(io, "seed=$(arguments.seed)")
        println(io, "rng_stream=weights first, then domino-creation draws")
        println(io, "julia_version=$(VERSION)")
        println(io, "weight_min=$(@sprintf("%.17g", minimum(weights)))")
        println(io, "weight_max=$(@sprintf("%.17g", maximum(weights)))")
        println(io, "weight_mean=$(@sprintf("%.17g", mean(weights)))")
        println(
            io,
            "weight_variance_population=" *
            "$(@sprintf("%.17g", var(weights; corrected=false)))",
        )
        println(io, "valid=$(validation.valid)")
        println(io, "dominoes=$(validation.dominoes)")
        println(io, "expected_dominoes=$(validation.expected_dominoes)")
        println(io, "covered_cells=$(validation.covered_cells)")
        println(io, "expected_cells=$(validation.expected_cells)")
        println(io, "overlaps=$(validation.overlaps)")
        println(io, "missing_cells=$(validation.missing_cells)")
        println(io, "outside_cells=$(validation.outside_cells)")
        println(io, "green=$(counts.green)")
        println(io, "blue=$(counts.blue)")
        println(io, "yellow=$(counts.yellow)")
        println(io, "red=$(counts.red)")
        println(io, "creation_probability_min=$(@sprintf("%.17g", probability_min))")
        println(io, "creation_probability_max=$(@sprintf("%.17g", probability_max))")
        println(io, "probability_seconds=$(@sprintf("%.6f", elapsed_probability))")
        println(io, "sampling_seconds=$(@sprintf("%.6f", elapsed_sampling))")
        println(io, "weights_file=$(basename(weights_path))")
        println(io, "weights_sha256=$(file_sha256(weights_path))")
        println(io, "tiling_file=$(basename(tiling_path))")
        println(io, "tiling_sha256=$(file_sha256(tiling_path))")
        println(io, "figure_file=$(basename(figure_path))")
        println(io, "figure_sha256=$(file_sha256(figure_path))")
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    if isnothing(parsed)
        print_help()
        return
    end
    parsed.order > 0 || error("--order must be positive")
    mkpath(parsed.output_dir)

    rng = Xoshiro(parsed.seed)
    weights = random_uniform_weights(rng, parsed.order)

    start_probability = time_ns()
    # The generic recurrence stores every Float64 probability table and is
    # intended for modest illustrated examples.  The Gamma production sampler
    # has a separate rolling-buffer implementation for large orders.
    probabilities = creation_probabilities(weights)
    elapsed_probability = (time_ns() - start_probability) / 1e9

    start_sampling = time_ns()
    tiling = sample_tiling(rng, probabilities)
    elapsed_sampling = (time_ns() - start_sampling) / 1e9

    validation = validate_tiling(tiling)
    validation.valid || error("generated tiling did not pass validation: $validation")
    counts = orientation_counts(tiling)
    probability_min = minimum(minimum, probabilities)
    probability_max = maximum(maximum, probabilities)

    weights_path = joinpath(parsed.output_dir, "weights_uniform_0_1.txt")
    tiling_path = joinpath(parsed.output_dir, "tiling_matrix.txt")
    figure_path = joinpath(parsed.output_dir, "tiling_mathematica_style.svg")
    metadata_path = joinpath(parsed.output_dir, "run_metadata.txt")

    write_table(weights_path, weights)
    write_table(tiling_path, tiling)
    write_svg(figure_path, tiling)
    write_metadata(
        metadata_path,
        parsed,
        weights,
        validation,
        counts,
        weights_path,
        tiling_path,
        figure_path,
        probability_min,
        probability_max,
        elapsed_probability,
        elapsed_sampling,
    )

    println("Random-weight Aztec diamond complete")
    println("  order:                 $(parsed.order)")
    println("  seed:                  $(parsed.seed)")
    println("  weight distribution:   iid Uniform(0,1)")
    println("  valid tiling:          $(validation.valid)")
    println("  dominoes:              $(validation.dominoes)")
    println("  orientation counts:    $counts")
    println("  creation p range:      [$probability_min, $probability_max]")
    println("  output directory:      $(parsed.output_dir)")
end

main(ARGS)
