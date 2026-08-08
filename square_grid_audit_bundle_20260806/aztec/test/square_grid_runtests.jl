const SG = AztecDiamond.SquareGrid

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

@testset "square-grid central dimer height" begin
    tree = SG.sample_full_tree(SG.BaselineEnvironment(), 8, Xoshiro(94_001))
    matching = SG.build_temperley_matching(tree)
    for separation in (1, 2, 4, 8)
        probes = SG.symmetric_probe_vertices(8, separation)
        @test probes.right.x - probes.left.x == separation
        increment = SG.spatial_height_increment(matching, separation)
        @test -separation <= increment <= separation
    end
    @test_throws ArgumentError SG.symmetric_probe_vertices(4, 8)
    @test_throws ArgumentError SG.spatial_height_increment(matching, 0)
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
