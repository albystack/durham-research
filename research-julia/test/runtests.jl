using LERWResearch
using StableRNGs
using Test

@testset "configuration" begin
    @test distribution_spec("baseline", Dict{String,Float64}()) == (:baseline, nothing)
    @test distribution_spec("gamma", Dict("shape" => 0.5)) == (:gamma, 0.5)
    @test_throws ArgumentError distribution_spec("pareto", Dict("alpha" => 1.0))
    @test stable_seed("same", 1) == stable_seed("same", 1)
    @test stable_seed("same", 1) != stable_seed("same", 2)
end

@testset "four-weight site environment" begin
    environment = SiteIIDEnvironment(123, "gamma", Dict("shape" => 1.0))
    origin = LERWResearch.pack_point(Int32(0), Int32(0))
    first_weights = LERWResearch.site_weights!(environment, origin)
    second_weights = LERWResearch.site_weights!(environment, origin)
    @test first_weights == second_weights
    @test length(first_weights) == 4
    @test all(>(0), first_weights)
    @test length(unique(first_weights)) > 1

    replica = SiteIIDEnvironment(123, "gamma", Dict("shape" => 1.0))
    @test LERWResearch.site_weights!(replica, origin) == first_weights

    tiny_cache = SiteIIDEnvironment(123, "gamma", Dict("shape" => 1.0);
                                    cache_capacity=1)
    original = LERWResearch.site_weights!(tiny_cache, origin)
    neighbour = LERWResearch.pack_point(Int32(1), Int32(0))
    LERWResearch.site_weights!(tiny_cache, neighbour)
    @test LERWResearch.site_weights!(tiny_cache, origin) == original

    specifications = [
        ("baseline", Dict{String,Float64}()),
        ("gamma", Dict("shape" => 0.5)),
        ("exponential", Dict{String,Float64}()),
        ("lognormal", Dict("sigma" => 1.0)),
        ("pareto", Dict("alpha" => 2.0)),
        ("uniform", Dict("a" => 0.7)),
        ("beta", Dict("a" => 0.5)),
        ("weibull", Dict("shape" => 0.7)),
        ("inverse_gamma", Dict("alpha" => 2.2)),
        ("bernoulli", Dict("a" => 0.8)),
        ("triangular", Dict("a" => 0.8)),
    ]
    for (distribution, params) in specifications
        candidate = SiteIIDEnvironment(321, distribution, params)
        weights = LERWResearch.site_weights!(candidate, origin)
        @test all(isfinite, weights)
        @test all(>(0), weights)
    end
end

@testset "batch reproducibility" begin
    task = BatchConfig(0, "gamma", Dict("shape" => 1.0), 16, 0, 2, 2,
                       UInt64(12345), "site_iid")
    first = run_batch(task; max_steps=2_000_000)
    second = run_batch(task; max_steps=2_000_000)
    scientific_columns(row) = (row.environment_seed, row.walk_seed, row.winding,
                               row.loop_erased_path_length, row.raw_walk_length,
                               row.exit_x, row.exit_y, row.sampled_site_count, row.status)
    @test scientific_columns.(first) == scientific_columns.(second)
end

@testset "loop erasure and winding" begin
    environment = SiteIIDEnvironment(456, "baseline", Dict{String,Float64}())
    path, raw_steps = loop_erased_walk(8, environment, StableRNG(789); max_steps=1_000_000)
    @test raw_steps > 0
    @test length(path) == length(unique(path))
    x, y = LERWResearch.unpack_point(path[end])
    @test is_boundary(x, y, 8)

    points = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)]
    packed = [LERWResearch.pack_point(Int32(x), Int32(y)) for (x, y) in points]
    @test winding(packed) == 2
end

@testset "batch and analysis pipeline" begin
    mktempdir() do directory
        config = joinpath(directory, "tasks.csv")
        results = joinpath(directory, "results")
        analysis = joinpath(directory, "analysis")
        rows = generate_config(config, ["baseline"], [8, 16, 32];
                               batches=1, num_environments=3,
                               walks_per_environment=3, base_seed=100)
        @test length(rows) == 3
        for task_id in 0:2
            task = load_config_row(config, task_id)
            output = result_path(results, task)
            LERWResearch.write_csv_atomic(output, run_batch(task; max_steps=2_000_000))
            @test completed_result(output, 9)
        end
        report = analyze_results(config, results, analysis; bootstrap_reps=10)
        @test report.observations == 27
        @test isfile(joinpath(analysis, "summary.csv"))
        @test isfile(joinpath(analysis, "loglog_fits.csv"))
        @test isfile(joinpath(analysis, "scaling_model_comparison.csv"))
    end
end
