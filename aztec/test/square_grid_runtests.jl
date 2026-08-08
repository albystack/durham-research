const SG = AztecDiamond.SquareGrid

module SquareGridRunnerHarness
include(joinpath(@__DIR__, "..", "scripts", "run_square_grid_campaign.jl"))
end

module SpatialAnalyzerHarness
include(joinpath(@__DIR__, "..", "scripts", "analyze_spatial_campaign.jl"))
end

const SGR = SquareGridRunnerHarness
const SGA = SpatialAnalyzerHarness

@testset "square-grid campaign identities and provenance" begin
    common = [
        "--config", joinpath(@__DIR__, "..", "configs", "square_grid_smoke.csv"),
        "--base-seed", "12345",
    ]
    baseline_gamma = SGR.parse_arguments(vcat(
        ["--environment-model", "baseline", "--distribution", "gamma", "--parameter", "0.5"],
        common,
    ))
    baseline_lognormal = SGR.parse_arguments(vcat(
        ["--environment-model", "baseline", "--distribution", "lognormal", "--parameter", "1.2"],
        common,
    ))
    directed_gamma = SGR.parse_arguments(vcat(
        ["--environment-model", "directed_site_iid", "--distribution", "gamma", "--parameter", "0.5"],
        common,
    ))
    directed_lognormal = SGR.parse_arguments(vcat(
        ["--environment-model", "directed_site_iid", "--distribution", "lognormal", "--parameter", "0.5"],
        common,
    ))
    directed_uniform = SGR.parse_arguments(vcat(
        ["--environment-model", "directed_site_iid", "--distribution", "uniform", "--parameter", "2"],
        common,
    ))
    undirected_gamma = SGR.parse_arguments(vcat(
        ["--environment-model", "undirected_conductance", "--distribution", "gamma", "--parameter", "0.5"],
        common,
    ))
    @test SGR.campaign_id(baseline_gamma) == "square_grid__baseline"
    @test SGR.campaign_id(baseline_gamma) == SGR.campaign_id(baseline_lognormal)
    @test SGR.sample_seed(baseline_gamma, 64, 7) ==
          SGR.sample_seed(baseline_lognormal, 64, 7)
    identities = SGR.campaign_id.((
        directed_gamma, directed_lognormal, directed_uniform, undirected_gamma))
    @test length(unique(identities)) == 4
    seeds = SGR.sample_seed.((
        directed_gamma, directed_lognormal, directed_uniform, undirected_gamma), 64, 7)
    @test length(unique(seeds)) == 4
    rows = SGR.read_config(directed_gamma.config)
    metadata = SGR.metadata_text(directed_gamma, rows, SGR.expand_tasks(rows))
    @test occursin("distribution=gamma", metadata)
    @test occursin("config_sha256=", metadata)
    @test occursin("square_grid_source_sha256=", metadata)
    @test occursin("schedule_order_first_samples_batch=", metadata)
    @test occursin("seed_identity=", metadata)
    manifest = SGR.manifest_text(metadata)
    @test occursin("campaign_id = \"square_grid__directed_site_iid__gamma__p_0p5\"", manifest)
end

@testset "square-grid analysis validation and block GLS" begin
    mktempdir() do directory
        valid_path = joinpath(directory, "batch_valid.csv")
        open(valid_path, "w") do io
            println(io, SGA.HEADER)
            for order in (16, 32), sample_id in 1:3, (numerator, denominator) in ((1, 4), (1, 2))
                separation = round(Int, 2 * order * numerator / denominator)
                println(io, join((
                    "test_campaign", order, sample_id, UInt64(order * 100 + sample_id),
                    numerator, denominator, separation, 0, separation,
                    sample_id, -sample_id, 2sample_id,
                ), ','))
            end
        end
        loaded = SGA.load_results([valid_path], "test_campaign")
        @test loaded.fractions == [(1, 4), (1, 2)]
        @test sort(unique(key[3] for key in keys(loaded.groups))) == [16, 32]
        @test SGA.load_results([valid_path], "legacy:test_campaign").fractions ==
              loaded.fractions
        @test_throws ErrorException SGA.load_results([valid_path], "wrong_campaign")

        incomplete_path = joinpath(directory, "batch_incomplete.csv")
        lines = readlines(valid_path)
        open(incomplete_path, "w") do io
            for (index, line) in enumerate(lines)
                index == 3 && continue
                println(io, line)
            end
        end
        @test_throws ErrorException SGA.load_results([incomplete_path], "test_campaign")
    end

    fractions = [(1, 32), (1, 16), (1, 8), (1, 4)]
    orders = [32, 64, 128, 256, 512, 1024]
    separations = Dict(
        (fraction..., order) => round(Int, 2 * order * fraction[1] / fraction[2])
        for fraction in fractions for order in orders
    )
    synthetic_data = (fractions=fractions, separations=separations)
    coefficient = 0.25
    values = Dict{Tuple{Int,Int,Int},Float64}()
    for (fraction_index, fraction) in enumerate(fractions), order in orders
        x = log(Float64(separations[(fraction..., order)]))
        values[(fraction..., order)] = 0.3fraction_index + 0.1fraction_index * x + coefficient * x^2
    end
    covariance_by_order = Dict(order => Matrix{Float64}(SGA.I, 4, 4) for order in orders)
    comparison = SGA.pooled_gls_comparison(
        synthetic_data, orders, values, 2, covariance_by_order)
    @test isapprox(comparison.common_coefficient, coefficient; atol=1e-10)
    @test comparison.quadratic_holdout_rmse < 1e-9
    @test comparison.null_holdout_rmse > comparison.quadratic_holdout_rmse
end

@testset "square-grid deterministic environments" begin
    L = 8
    directed = SG.DirectedSiteIIDEnvironment(0x1234; distribution=:gamma, parameter=0.5)
    first = SG.edge_weight(directed, 0, 0, SG.EAST, L)
    @test first == SG.edge_weight(directed, 0, 0, SG.EAST, L)
    @test first != SG.edge_weight(directed, 1, 0, SG.WEST, L)

    undirected = SG.UndirectedConductanceEnvironment(
        0x5678; distribution=:gamma, parameter=0.5)
    forward = SG.edge_weight(undirected, 0, 0, SG.EAST, L)
    reverse_weight = SG.edge_weight(undirected, 1, 0, SG.WEST, L)
    @test forward == reverse_weight

    # Query order must not alter any environment value.
    keys = [
        (-2, 1, SG.NORTH),
        (0, 0, SG.EAST),
        (3, -1, SG.SOUTH),
        (-1, -2, SG.WEST),
    ]
    values_forward = [SG.edge_weight(directed, x, y, d, L) for (x, y, d) in keys]
    values_reverse = Dict(
        (x, y, d) => SG.edge_weight(directed, x, y, d, L)
        for (x, y, d) in reverse(keys)
    )
    @test all(
        values_forward[index] == values_reverse[keys[index]]
        for index in eachindex(keys)
    )
end

@testset "square-grid disorder laws and materialization" begin
    gamma = SG.DirectedSiteIIDEnvironment(0x7711; distribution=:gamma, parameter=0.5)
    lognormal = SG.DirectedSiteIIDEnvironment(
        0x7712; distribution=:lognormal, parameter=0.7)
    uniform = SG.DirectedSiteIIDEnvironment(
        0x7713; distribution=:uniform, parameter=2.0)
    gamma_draws = [SG.edge_weight(gamma, index, -index, SG.NORTH, 8) for index in 1:20_000]
    lognormal_draws = [
        SG.edge_weight(lognormal, index, -index, SG.EAST, 8) for index in 1:20_000
    ]
    uniform_draws = [
        SG.edge_weight(uniform, index, -index, SG.SOUTH, 8) for index in 1:20_000
    ]
    @test all(>(0), gamma_draws)
    @test all(>(0), lognormal_draws)
    @test all(>(0), uniform_draws)
    @test isapprox(mean(gamma_draws), 1.0; atol=0.06)
    @test isapprox(var(gamma_draws), 2.0; atol=0.22)
    @test isapprox(mean(lognormal_draws), 1.0; atol=0.04)
    @test isapprox(var(lognormal_draws), exp(0.7^2) - 1; atol=0.08)
    @test isapprox(mean(uniform_draws), 1.0; atol=0.03)
    @test isapprox(var(uniform_draws), 1 / 3; atol=0.03)

    materialized = SG.materialize_environment(gamma, 6)
    first = SG.sample_full_tree(gamma, 6, Xoshiro(77_013))
    second = SG.sample_full_tree(materialized, 6, Xoshiro(77_013))
    @test first.parent_direction == second.parent_direction
    @test first.raw_steps == second.raw_steps
end

@testset "square-grid categorical transition law" begin
    environment = SG.FixedEnvironment((1.0, 2.0, 3.0, 4.0))
    rng = Xoshiro(91_001)
    counts = zeros(Int, 4)
    samples = 50_000
    for _ in 1:samples
        direction = SG.choose_direction(environment, 0, 0, 1, rng)
        counts[Int(direction)] += 1
    end
    observed = counts ./ samples
    expected = [0.1, 0.2, 0.3, 0.4]
    @test maximum(abs.(observed .- expected)) < 0.012
end

@testset "square-grid Wilson tree" begin
    for L in 1:5
        tree = SG.sample_full_tree(SG.BaselineEnvironment(), L, Xoshiro(92_000 + L))
        validation = SG.validate_tree(tree)
        @test validation.valid
        @test length(tree.parent_direction) == (2L - 1)^2
        @test all(direction -> direction in (SG.NORTH, SG.EAST, SG.SOUTH, SG.WEST),
                  tree.parent_direction)
        @test tree.raw_steps >= length(tree.parent_direction)
    end

    first = SG.sample_full_tree(SG.BaselineEnvironment(), 5, Xoshiro(92_100))
    second = SG.sample_full_tree(SG.BaselineEnvironment(), 5, Xoshiro(92_100))
    @test first.parent_direction == second.parent_direction
    @test first.raw_steps == second.raw_steps
end

@testset "generalized Temperley matching" begin
    for L in 1:5
        tree = SG.sample_full_tree(SG.BaselineEnvironment(), L, Xoshiro(93_000 + L))
        matching = SG.build_temperley_matching(tree)
        validation = SG.validate_matching(matching)
        @test validation.valid
        @test validation.edges == SG.primal_edge_count(L)
        @test validation.vertex_matches == (2L - 1)^2
        @test validation.face_matches == (2L)^2 - 1
        @test count(==(UInt8(1)), matching.edge_owner) == (2L - 1)^2
        @test count(==(UInt8(2)), matching.edge_owner) == (2L)^2 - 1
    end
end

@testset "exhaustive order-two Temperley incidences" begin
    valid_tree_count = 0
    for code in 0:(4^9 - 1)
        value = code
        parents = Vector{UInt8}(undef, 9)
        for index in eachindex(parents)
            parents[index] = UInt8(mod(value, 4) + 1)
            value = div(value, 4)
        end
        tree = SG.TreeSample(2, parents, 0, 0)
        SG.validate_tree(tree).valid || continue
        valid_tree_count += 1
        matching = SG.build_temperley_matching(tree)
        @test SG.validate_matching(matching).valid
        probes = SG.symmetric_probe_vertices(2, 1)
        direct_path = [(x, 0) for x in (2 * probes.left.x):(2 * probes.right.x)]
        @test SG.dimer_height_increment_along_path(matching, direct_path) ==
              SG.spatial_height_increment(matching, 1)
    end
    @test valid_tree_count == 100_352
end

@testset "square-grid central dimer height" begin
    for seed in 94_001:94_020
        tree = SG.sample_full_tree(SG.BaselineEnvironment(), 8, Xoshiro(seed))
        matching = SG.build_temperley_matching(tree)
        for separation in (1, 2, 4, 8)
            probes = SG.symmetric_probe_vertices(8, separation)
            @test probes.right.x - probes.left.x == separation
            increment = SG.spatial_height_increment(matching, separation)
            @test -separation <= increment <= separation

            left = 2 * probes.left.x
            right = 2 * probes.right.x
            direct_path = [(x, 0) for x in left:right]
            detour_path = vcat(
                [(left, y) for y in 0:-1:-2],
                [(x, -2) for x in (left + 1):right],
                [(right, y) for y in -1:0],
            )
            generic_direct = SG.dimer_height_increment_along_path(
                matching, direct_path)
            generic_detour = SG.dimer_height_increment_along_path(
                matching, detour_path)
            @test generic_direct == increment
            @test generic_detour == increment
            @test SG.dimer_height_increment_along_path(
                matching, reverse(direct_path)) == -increment
        end
    end
    matching = SG.build_temperley_matching(
        SG.sample_full_tree(SG.BaselineEnvironment(), 8, Xoshiro(94_100)))
    @test_throws ArgumentError SG.symmetric_probe_vertices(4, 8)
    @test_throws ArgumentError SG.spatial_height_increment(matching, 0)
    @test_throws ArgumentError SG.dimer_height_increment_along_path(matching, [(0, 0)])
end

@testset "paired square-grid shared environment" begin
    separations = [1, 2, 4]
    first = SG.sample_spatial_increment_pair(
        UInt64(95_001), 8, separations;
        environment_model=:directed_site_iid,
        distribution=:gamma,
        parameter=0.5,
    )
    second = SG.sample_spatial_increment_pair(
        UInt64(95_001), 8, separations;
        environment_model=:directed_site_iid,
        distribution=:gamma,
        parameter=0.5,
    )
    @test first == second
    @test first.diagnostics.replica_1_seed != first.diagnostics.replica_2_seed
    @test [row.separation for row in first.rows] == separations
    @test all(row.difference == row.increment_1 - row.increment_2 for row in first.rows)

    other = SG.sample_spatial_increment_pair(
        UInt64(95_002), 8, separations;
        environment_model=:directed_site_iid,
        distribution=:gamma,
        parameter=0.5,
    )
    @test first.diagnostics.environment_seed != other.diagnostics.environment_seed
end

@testset "square-grid finite-sample replica identity" begin
    first = [2, -1, 3, 0, 4, 1]
    second = [1, 2, -2, 0, 3, -1]
    difference = first .- second
    marginal = (var(first) + var(second)) / 2
    conditional = var(difference) / 2
    disorder = cov(first, second)
    @test isapprox(marginal - disorder, conditional; atol=1e-12)
end
