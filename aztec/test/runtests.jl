using Test
using Random
using Statistics

include(joinpath(@__DIR__, "..", "src", "AztecDiamond.jl"))
using .AztecDiamond

@testset "weighted Aztec diamond sampler" begin
    for order in 1:12
        rng = Xoshiro(10_000 + order)
        weights = random_uniform_weights(rng, order)
        probabilities = creation_probabilities(weights)
        tiling = sample_tiling(rng, probabilities)
        validation = validate_tiling(tiling)

        @test length(probabilities) == order
        @test size(first(probabilities)) == (1, 1)
        @test size(last(probabilities)) == (order, order)
        @test all(matrix -> all(p -> 0 <= p <= 1, matrix), probabilities)
        @test validation.valid
        @test validation.dominoes == order * (order + 1)
        @test sum(values(orientation_counts(tiling))) == order * (order + 1)
        center = center_face_index(order)
        @test center_height(tiling) == height_function(tiling)[center.row, center.column]

    end
end

@testset "uniform order-one probabilities" begin
    probabilities = creation_probabilities(ones(2, 2))
    @test probabilities[1][1, 1] == 0.5
end

@testset "Gamma-disordered model" begin
    for order in 1:12
        rng = Xoshiro(20_000 + order)
        weights = gamma_disordered_weights(rng, order; alpha=0.2, beta=0.25)
        probabilities = gamma_disordered_probabilities(weights.a, weights.b)
        tiling = sample_tiling(rng, probabilities)
        validation = validate_tiling(tiling)

        @test all(>(0), weights.a)
        @test all(>(0), weights.b)
        @test length(probabilities) == order
        @test size(first(probabilities)) == (1, 1)
        @test size(last(probabilities)) == (order, order)
        @test all(matrix -> all(p -> 0 <= p <= 1, matrix), probabilities)
        @test validation.valid
        center = center_face_index(order)
        @test center_height(tiling) == height_function(tiling)[center.row, center.column]

        choice_rng = Xoshiro(30_000 + order)
        choices = gamma_disordered_creation_choices(
            choice_rng,
            weights.a,
            weights.b,
        )
        choice_tiling = sample_tiling_from_choices(choices)
        choice_validation = validate_tiling(choice_tiling)
        @test length(choices) == order
        @test all(level -> size(choices[level]) == (level, level), 1:order)
        @test choice_validation.valid
        choice_center = center_face_index(order)
        @test center_height(choice_tiling) ==
              height_function(choice_tiling)[choice_center.row, choice_center.column]
    end
end

@testset "Gamma shape/scale convention" begin
    rng = Xoshiro(42)
    draws = random_gamma_weights(rng, 500, 400; shape=0.5, scale=2.0)
    @test all(>(0), draws)
    @test isapprox(mean(draws), 1.0; atol=0.015)
    @test isapprox(var(draws; corrected=false), 2.0; atol=0.06)
end
