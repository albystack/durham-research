using Test
using Random
using Statistics
using AztecDiamond

# The tests are ordered from generic shuffling, through the paper-specific
# Gamma recurrence, to production wrappers.  Small orders permit full geometric
# validation; production runs use the same code paths at larger orders.

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

@testset "order-two shuffle enumerates uniform tilings" begin
    # There are 2^(2*3/2) = 8 order-two Aztec tilings.  Enumerating all five
    # potential creation bits produces 32 bit patterns; unused bits should make
    # every valid tiling occur exactly four times, never biasing the result.
    tiling_counts = Dict{Tuple,Int}()
    for mask in 0:31
        choices = [falses(1, 1), falses(2, 2)]
        bit = 0
        for level in 1:2, index in eachindex(choices[level])
            choices[level][index] = !iszero(mask & (1 << bit))
            bit += 1
        end
        tiling = sample_tiling_from_choices(choices)
        @test validate_tiling(tiling).valid
        key = Tuple(vec(tiling))
        tiling_counts[key] = get(tiling_counts, key, 0) + 1
    end
    @test length(tiling_counts) == 8
    @test all(==(4), values(tiling_counts))
end

@testset "Duits--Van Peski recurrence reference values" begin
    # A hand-computable order-two example guards equation (1.22) independently
    # of the random sampler.  The order-one values use a'[1,1] and b'[1,1].
    a = [2.0 3.0; 5.0 7.0]
    b = [11.0 13.0; 17.0 19.0]
    probabilities = gamma_disordered_probabilities(a, b)

    reduced_a = (2 / (2 + 11)) * (5 + 17)
    reduced_b = (13 / (3 + 13)) * (7 + 19)
    @test probabilities[2] == b ./ (a .+ b)
    @test probabilities[1][1, 1] == reduced_b / (reduced_a + reduced_b)
end

@testset "pre-drawn Gamma coins match reduced probabilities" begin
    # Both public implementations must encode exactly the same recurrence and
    # b/(a+b) choice convention.  Replaying an identical scalar RNG stream in
    # descending level order makes this an exact, bit-for-bit comparison.
    environment_rng = Xoshiro(12_345)
    weights = gamma_disordered_weights(environment_rng, 9; alpha=0.2, beta=0.25)
    probabilities = gamma_disordered_probabilities(weights.a, weights.b)

    coin_seed = 54_321
    actual = gamma_disordered_creation_choices(
        Xoshiro(coin_seed),
        weights.a,
        weights.b,
    )
    expected = Vector{BitMatrix}(undef, 9)
    expected_rng = Xoshiro(coin_seed)
    for level in 9:-1:1
        expected[level] = falses(level, level)
        for index in eachindex(expected[level], probabilities[level])
            expected[level][index] =
                rand(expected_rng) < probabilities[level][index]
        end
    end
    @test actual == expected
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


@testset "public argument validation" begin
    rng = Xoshiro(1)
    @test_throws ArgumentError random_uniform_weights(rng, 0)
    @test_throws ArgumentError random_gamma_weights(rng, 2, 2; shape=0.0)
    @test_throws ArgumentError center_face_index(0)
    @test_throws ArgumentError sample_tiling(rng, [fill(0.5, 2, 2)])
    @test_throws ArgumentError sample_tiling_from_choices([falses(2, 2)])
end


@testset "Gamma center-height reproducibility" begin
    for order in (4, 8, 12)
        seed = UInt64(50_000 + order)
        first_height = sample_gamma_center_height(
            seed,
            order;
            alpha=0.2,
            beta=0.25,
        )
        second_height = sample_gamma_center_height(
            seed,
            order;
            alpha=0.2,
            beta=0.25,
        )
        @test first_height == second_height
    end
end

@testset "Gamma center-height wrapper agrees with manual pipeline" begin
    for order in (3, 6, 9)
        seed = UInt64(60_000 + order)
        expected = sample_gamma_center_height(
            seed,
            order;
            alpha=0.2,
            beta=0.25,
        )

        rng = Xoshiro(seed)
        weights = gamma_disordered_weights(
            rng,
            order;
            alpha=0.2,
            beta=0.25,
        )
        choices = gamma_disordered_creation_choices(
            rng,
            weights.a,
            weights.b,
        )
        tiling = sample_tiling_from_choices(choices)
        @test expected == center_height(tiling)
    end
end

@testset "independent Gamma center-height seeds" begin
    heights = [
        sample_gamma_center_height(
            UInt64(70_000 + sample_id),
            16;
            alpha=0.2,
            beta=0.25,
        )
        for sample_id in 1:20
    ]
    @test length(heights) == 20
    @test length(unique(heights)) > 1
end

@testset "double-dimer center-height pairs" begin
    differences = Int[]
    for order in (4, 8, 12)
        seed = UInt64(80_000 + order)
        first_result = sample_gamma_center_height_pair(
            seed,
            order;
            alpha=0.2,
            beta=0.25,
        )
        second_result = sample_gamma_center_height_pair(
            seed,
            order;
            alpha=0.2,
            beta=0.25,
        )
        @test first_result == second_result
        @test first_result.difference ==
              first_result.height_1 - first_result.height_2
        push!(differences, first_result.difference)

        rng = Xoshiro(seed)
        weights = gamma_disordered_weights(rng, order; alpha=0.2, beta=0.25)
        choices = gamma_disordered_creation_choice_pair(rng, weights.a, weights.b)
        first_tiling = sample_tiling_from_choices(choices.first)
        second_tiling = sample_tiling_from_choices(choices.second)
        @test first_result.height_1 == center_height(first_tiling)
        @test first_result.height_2 == center_height(second_tiling)
        @test validate_tiling(first_tiling).valid
        @test validate_tiling(second_tiling).valid
    end
    @test any(!=(0), differences)
end

@testset "multiple conditional copies share only the environment" begin
    rng = Xoshiro(90_000)
    weights = gamma_disordered_weights(rng, 12; alpha=0.2, beta=0.25)
    copies = gamma_disordered_creation_choices(rng, weights.a, weights.b, 3)
    @test length(copies) == 3
    @test all(length(copy) == 12 for copy in copies)
    @test any(copies[1][level] != copies[2][level] for level in 1:12)
    @test_throws ArgumentError gamma_disordered_creation_choices(
        rng,
        weights.a,
        weights.b,
        0,
    )
end

@testset "finite-sample double-dimer variance identity" begin
    # This is the algebra used by the analysis script.  Corrected sample
    # variance and covariance must all use the same n-1 denominator.
    height_1 = [10, 13, 11, 15, 12]
    height_2 = [9, 12, 14, 13, 10]
    difference = height_1 .- height_2
    paired_marginal = (var(height_1) + var(height_2)) / 2
    @test isapprox(
        var(difference) / 2,
        paired_marginal - cov(height_1, height_2);
        atol=1e-12,
    )
end


@testset "direct spatial height increments" begin
    rng = Xoshiro(100_001)
    order = 14
    tiling = sample_tiling_from_choices(uniform_creation_choices(rng, order))
    full_heights = height_function(tiling)
    for separation in (1, 2, 4, 7, 13)
        columns = symmetric_face_columns(order, separation)
        expected =
            full_heights[order + 1, columns.right] -
            full_heights[order + 1, columns.left]
        @test symmetric_height_increment(tiling, separation) == expected
        @test face_height(tiling, order + 1, columns.left) ==
              full_heights[order + 1, columns.left]
        @test face_height(tiling, order + 1, columns.right) ==
              full_heights[order + 1, columns.right]
    end
    @test_throws ArgumentError symmetric_face_columns(order, 0)
    @test_throws ArgumentError symmetric_face_columns(order, order)
end

@testset "Gamma and uniform spatial replica wrappers" begin
    separations = [2, 4, 8]
    gamma_first = sample_gamma_spatial_increment_pair(
        UInt64(110_001),
        16,
        separations;
        alpha=0.2,
        beta=0.25,
    )
    gamma_second = sample_gamma_spatial_increment_pair(
        UInt64(110_001),
        16,
        separations;
        alpha=0.2,
        beta=0.25,
    )
    @test gamma_first == gamma_second
    @test [row.separation for row in gamma_first] == separations
    @test all(
        row.difference == row.increment_1 - row.increment_2
        for row in gamma_first
    )

    uniform_first = sample_uniform_spatial_increment_pair(
        UInt64(120_001),
        16,
        separations,
    )
    uniform_second = sample_uniform_spatial_increment_pair(
        UInt64(120_001),
        16,
        separations,
    )
    @test uniform_first == uniform_second
    @test [row.separation for row in uniform_first] == separations
    @test all(
        row.difference == row.increment_1 - row.increment_2
        for row in uniform_first
    )
end
