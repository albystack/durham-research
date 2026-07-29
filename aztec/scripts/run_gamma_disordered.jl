#!/usr/bin/env julia

using Random
using Printf
using SHA
using Statistics

include(joinpath(@__DIR__, "..", "src", "AztecDiamond.jl"))
using .AztecDiamond

function parse_arguments(arguments)
    options = Dict{String,String}(
        "order" => "200",
        "seed" => "20260728",
        "alpha" => "0.2",
        "beta" => "0.25",
        "output-dir" => joinpath(
            @__DIR__,
            "..",
            "output",
            "gamma_disordered_alpha_0p2_beta_0p25_n200_seed20260728",
        ),
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
        alpha=parse(Float64, options["alpha"]),
        beta=parse(Float64, options["beta"]),
        output_dir=abspath(options["output-dir"]),
    )
end

file_sha256(path) = bytes2hex(open(sha256, path))

function print_statistics(io, name, values)
    println(io, "$(name)_min=$(@sprintf("%.17g", minimum(values)))")
    println(io, "$(name)_max=$(@sprintf("%.17g", maximum(values)))")
    println(io, "$(name)_mean=$(@sprintf("%.17g", mean(values)))")
    println(
        io,
        "$(name)_variance_population=$(@sprintf("%.17g", var(values; corrected=false)))",
    )
end

function write_metadata(
    path,
    arguments,
    weights,
    validation,
    counts,
    a_path,
    b_path,
    tiling_path,
    figure_path,
    probability_min,
    probability_max,
    elapsed_probability,
    elapsed_sampling,
)
    open(path, "w") do io
        println(io, "model=biased Gamma-disordered Aztec diamond")
        println(io, "reference=Duits and Van Peski, arXiv:2512.03033, Definition 1.1")
        println(io, "algorithm=equations (1.22)-(1.23) followed by domino shuffling")
        println(io, "order=$(arguments.order)")
        println(io, "matrix_side=$(2 * arguments.order)")
        println(io, "a_distribution=Gamma(shape=$(arguments.alpha), scale=1)")
        println(io, "b_distribution=Gamma(shape=$(arguments.beta), scale=1)")
        println(io, "other_edge_weights=gauge-fixed to 1")
        println(io, "rng=Xoshiro")
        println(io, "seed=$(arguments.seed)")
        println(io, "rng_stream=a weights, then b weights, then domino-creation draws")
        println(io, "julia_version=$(VERSION)")
        print_statistics(io, "a_weight", weights.a)
        print_statistics(io, "b_weight", weights.b)
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
        println(io, "a_weights_file=$(basename(a_path))")
        println(io, "a_weights_sha256=$(file_sha256(a_path))")
        println(io, "b_weights_file=$(basename(b_path))")
        println(io, "b_weights_sha256=$(file_sha256(b_path))")
        println(io, "tiling_file=$(basename(tiling_path))")
        println(io, "tiling_sha256=$(file_sha256(tiling_path))")
        println(io, "figure_file=$(basename(figure_path))")
        println(io, "figure_sha256=$(file_sha256(figure_path))")
    end
end

function main(arguments)
    parsed = parse_arguments(arguments)
    parsed.order > 0 || error("--order must be positive")
    parsed.alpha > 0 || error("--alpha must be positive")
    parsed.beta > 0 || error("--beta must be positive")
    mkpath(parsed.output_dir)

    rng = Xoshiro(parsed.seed)
    weights = gamma_disordered_weights(
        rng,
        parsed.order;
        alpha=parsed.alpha,
        beta=parsed.beta,
    )

    start_probability = time_ns()
    probabilities = gamma_disordered_probabilities(weights.a, weights.b)
    elapsed_probability = (time_ns() - start_probability) / 1e9

    start_sampling = time_ns()
    tiling = sample_tiling(rng, probabilities)
    elapsed_sampling = (time_ns() - start_sampling) / 1e9

    validation = validate_tiling(tiling)
    validation.valid || error("generated tiling did not pass validation: $validation")
    counts = orientation_counts(tiling)
    probability_min = minimum(minimum, probabilities)
    probability_max = maximum(maximum, probabilities)

    a_path = joinpath(parsed.output_dir, "a_weights_gamma.txt")
    b_path = joinpath(parsed.output_dir, "b_weights_gamma.txt")
    tiling_path = joinpath(parsed.output_dir, "tiling_matrix.txt")
    figure_path = joinpath(parsed.output_dir, "tiling_mathematica_style.svg")
    metadata_path = joinpath(parsed.output_dir, "run_metadata.txt")

    write_table(a_path, weights.a)
    write_table(b_path, weights.b)
    write_table(tiling_path, tiling)
    write_svg(figure_path, tiling)
    write_metadata(
        metadata_path,
        parsed,
        weights,
        validation,
        counts,
        a_path,
        b_path,
        tiling_path,
        figure_path,
        probability_min,
        probability_max,
        elapsed_probability,
        elapsed_sampling,
    )

    println("Gamma-disordered Aztec diamond complete")
    println("  order:                 $(parsed.order)")
    println("  seed:                  $(parsed.seed)")
    println("  a weights:             Gamma($(parsed.alpha), 1)")
    println("  b weights:             Gamma($(parsed.beta), 1)")
    println("  valid tiling:          $(validation.valid)")
    println("  dominoes:              $(validation.dominoes)")
    println("  orientation counts:    $counts")
    println("  creation p range:      [$probability_min, $probability_max]")
    println("  output directory:      $(parsed.output_dir)")
end

main(ARGS)
