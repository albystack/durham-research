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

@testset "strict annealed reproduction" begin
    mktempdir() do directory
        campaign = joinpath(directory, "campaign.csv")
        rows = generate_strict_annealed_reproduction_config(campaign)
        @test length(rows) == 3_790
        @test all(row -> row.walks_per_environment == 1, rows)
        @test sum(row.num_environments for row in rows) == 379_000
        @test count(row -> row.L == 8192, rows) == 10

        config = joinpath(directory, "small.csv")
        results = joinpath(directory, "results")
        analysis = joinpath(directory, "analysis")
        generate_config(config, ["gamma:1.0"], [8, 16, 32];
                        batches=1, num_environments=4,
                        walks_per_environment=1, base_seed=200)
        campaign_report = run_campaign(config, results; strict_annealed=true,
                                       max_steps=2_000_000)
        @test campaign_report.new_tasks == 3
        resumed_report = run_campaign(config, results; strict_annealed=true,
                                      max_steps=2_000_000)
        @test resumed_report.already_complete == 3
        observed_seeds = String[]
        for task_id in 0:2
            task = load_config_row(config, task_id)
            batch_rows = NamedTuple.(LERWResearch.CSV.File(result_path(results, task)))
            @test all(row -> row.walk_id == 0, batch_rows)
            append!(observed_seeds, string.(getproperty.(batch_rows, :environment_seed)))
        end
        @test length(observed_seeds) == length(unique(observed_seeds))
        report = analyze_results(config, results, analysis;
                                 bootstrap_reps=10, strict_annealed=true)
        @test report.observations == 12
        @test report.fits == 1
        fits = collect(LERWResearch.CSV.File(joinpath(analysis, "loglog_fits.csv")))
        @test all(row -> row.variance_kind == "annealed", fits)
    end

    invalid = BatchConfig(99, "gamma", Dict("shape" => 1.0), 16, 0, 2, 2,
                          UInt64(123), "site_iid")
    @test_throws ArgumentError LERWResearch.require_strict_annealed([invalid])
end

@testset "double-dimer pairs" begin
    mktempdir() do directory
        frozen = joinpath(directory, "double_dimer_campaign.csv")
        frozen_rows = generate_double_dimer_reproduction_config(frozen)
        @test length(frozen_rows) == 3_790
        @test all(row -> row.walks_per_environment == 2, frozen_rows)
        @test sum(row.num_environments for row in frozen_rows) == 379_000
        @test sum(row.num_environments * row.walks_per_environment for row in frozen_rows) == 758_000
        pilot = generate_double_dimer_pilot_config(joinpath(directory, "pilot.csv");
                                                    environments_per_size=10)
        @test length(pilot) == 137
        @test sum(row.num_environments for row in pilot) == 1_370

        config = joinpath(directory, "small.csv")
        results = joinpath(directory, "results")
        analysis = joinpath(directory, "analysis")
        generate_config(config, ["gamma:1.0"], [8, 16, 32];
                        batches=1, num_environments=8,
                        walks_per_environment=2, base_seed=300)
        campaign = run_campaign(config, results; double_dimer=true,
                                max_steps=2_000_000)
        @test campaign.new_tasks == 3

        all_rows = NamedTuple[]
        for task_id in 0:2
            task = load_config_row(config, task_id)
            rows = NamedTuple.(LERWResearch.CSV.File(result_path(results, task); types=String))
            append!(all_rows, rows)
            by_environment = Dict{String,Vector{Any}}()
            for row in rows
                push!(get!(by_environment, row.environment_id, Any[]), row)
            end
            @test all(pair -> sort(parse.(Int, getproperty.(pair, :walk_id))) == [0, 1],
                      values(by_environment))
            @test all(pair -> length(unique(getproperty.(pair, :environment_seed))) == 1,
                      values(by_environment))
            @test all(pair -> length(unique(getproperty.(pair, :walk_seed))) == 2,
                      values(by_environment))
        end

        pairs = LERWResearch.make_double_dimer_pairs(all_rows)
        @test length(pairs) == 24
        @test all(row -> row.winding_difference == row.winding_1 - row.winding_2, pairs)

        report = analyze_results(config, results, analysis;
                                 bootstrap_reps=10, double_dimer=true)
        @test report.observations == 24
        @test report.raw_walks == 48
        @test report.fits == 1
        @test isfile(joinpath(analysis, "double_dimer_pairs.csv"))
        summary = collect(LERWResearch.CSV.File(joinpath(analysis, "summary.csv")))
        @test all(row -> row.pairs == 8, summary)
        @test all(row -> isapprox(row.double_dimer_variance,
                                  row.variance_identity_rhs; atol=1e-12), summary)
    end

    invalid = BatchConfig(100, "baseline", Dict{String,Float64}(), 16, 0, 2, 1,
                          UInt64(123), "site_iid")
    @test_throws ArgumentError LERWResearch.require_double_dimer([invalid])
end

@testset "heterogeneous table schemas" begin
    mktempdir() do directory
        path = joinpath(directory, "mixed.csv")
        rows = NamedTuple[(common=1, python_version="3.14"),
                          (common=2, julia_version="1.12")]
        LERWResearch.write_table(path, rows)
        text = read(path, String)
        @test occursin("python_version", text)
        @test occursin("julia_version", text)
        @test length(collect(LERWResearch.CSV.File(path))) == 2
    end
end
