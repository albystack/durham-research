using CSV
using LERWResearch
using StableRNGs
using Test

const LR = LERWResearch

function task(; task_id=0, distribution="gamma",
              params=Dict("shape" => 0.5), L=16, batch_id=0,
              num_environments=8, walks_per_environment=1,
              base_seed=UInt64(20260726), environment_model="site_iid")
    return BatchConfig(task_id, distribution, Dict{String,Float64}(params), L,
                       batch_id, num_environments, walks_per_environment,
                       base_seed, environment_model)
end

function deterministic_fields(row)
    ignored = Set((:runtime, :started_at_utc, :finished_at_utc))
    names = Tuple(name for name in keys(row) if name ∉ ignored)
    return NamedTuple{names}(Tuple(getproperty(row, name) for name in names))
end

@testset "LERWResearch" begin
    @testset "distribution parsing and validation" begin
        valid = [
            "baseline", "symmetric", "exponential", "exp", "gamma:0.5",
            "lognormal:1.0", "pareto:2.0", "uniform:0.7", "beta:0.5",
            "weibull:0.7", "inverse_gamma:2.2", "bernoulli:0.8",
            "triangular:0.8",
        ]
        for specification in valid
            name, params = LR.parse_distribution_argument(specification)
            model, parameter = LR.distribution_spec(name, params)
            @test model isa Symbol
            @test parameter === nothing || parameter > 0
            rng = StableRNG(123)
            @test all(isfinite(weight) && weight > 0
                      for weight in (LR.sample_weight(rng, model, parameter)
                                     for _ in 1:100))
        end
        @test LR.distribution_spec("gamma_edges", Dict("k" => 2.0)) == (:gamma, 2.0)
        @test LR.parse_params("") == Dict{String,Float64}()
        @test LR.parse_params("""{"sigma":1.0}""") == Dict("sigma" => 1.0)
        @test_throws ArgumentError LR.parse_params("[1,2]")
        @test_throws ArgumentError LR.distribution_spec("unknown", Dict())
        @test_throws ArgumentError LR.distribution_spec("baseline", Dict("a" => 1.0))
        @test_throws ArgumentError LR.distribution_spec("gamma", Dict())
        @test_throws ArgumentError LR.distribution_spec(
            "gamma", Dict("shape" => 1.0, "extra" => 2.0))
        @test_throws ArgumentError LR.distribution_spec("gamma", Dict("shape" => 0.0))
        @test_throws ArgumentError LR.distribution_spec("lognormal", Dict("sigma" => -1.0))
        @test_throws ArgumentError LR.distribution_spec("pareto", Dict("alpha" => 1.0))
        @test_throws ArgumentError LR.distribution_spec("inverse_gamma", Dict("alpha" => 0.8))
        @test_throws ArgumentError LR.distribution_spec("uniform", Dict("a" => 1.0))
        @test_throws ArgumentError LR.distribution_spec("bernoulli", Dict("a" => -0.1))
        @test_throws ArgumentError LR.parse_distribution_argument("gamma:not-a-number")
    end

    @testset "stable seeds" begin
        @test stable_seed("test", 123, "gamma") == UInt64(13875854726940453535)
        @test stable_seed(20260726, "temporal") == UInt64(17941603992211970450)
        @test stable_seed("a", "bc") != stable_seed("ab", "c")
        @test stable_seed("repeat", 42) == stable_seed("repeat", 42)
    end

    @testset "site-i.i.d. fixed weights" begin
        environment = SiteIIDEnvironment(123, "gamma", Dict("shape" => 0.5);
                                         cache_capacity=2)
        origin = LR.pack_point(Int32(0), Int32(0))
        first_weights = LR.site_weights!(environment, origin)
        @test LR.site_weights!(environment, origin) == first_weights
        # Even after direct-mapped eviction, the point seed regenerates exactly.
        for x in Int32(1):Int32(20)
            LR.site_weights!(environment, LR.pack_point(x, Int32(0)))
        end
        @test LR.site_weights!(environment, origin) == first_weights
    end

    @testset "temporal-i.i.d. resampling and reproducibility" begin
        point = LR.pack_point(Int32(0), Int32(0))
        environment = TemporalIIDEnvironment(123, "gamma", Dict("shape" => 0.5))
        first_weights = temporal_weights!(environment)
        second_weights = temporal_weights!(environment)
        @test first_weights != second_weights
        @test environment.sampled_weight_vectors == 2
        # The lattice point does not enter temporal weight generation.
        third_weights = LR.transition_weights!(environment, point)
        fourth_weights = LR.transition_weights!(environment, point)
        @test third_weights != fourth_weights

        function temporal_walk(environment_seed, direction_seed)
            env = TemporalIIDEnvironment(environment_seed, "lognormal",
                                         Dict("sigma" => 1.0))
            return loop_erased_walk_diagnostics(24, env, StableRNG(direction_seed))
        end
        path_1, raw_1, counts_1 = temporal_walk(1001, 2002)
        path_2, raw_2, counts_2 = temporal_walk(1001, 2002)
        @test path_1 == path_2
        @test raw_1 == raw_2
        @test counts_1 == counts_2
        path_3, raw_3, counts_3 = temporal_walk(1002, 2003)
        @test (path_3, raw_3, counts_3) != (path_1, raw_1, counts_1)
        @test sum(counts_1) == raw_1
    end

    @testset "temporal directional symmetry" begin
        environment = TemporalIIDEnvironment(555, "pareto", Dict("alpha" => 2.0))
        direction_rng = StableRNG(777)
        point = LR.pack_point(Int32(0), Int32(0))
        counts = zeros(Int, 4)
        for _ in 1:100_000
            _, direction = LR.next_point_and_direction(point, environment, direction_rng)
            counts[direction] += 1
        end
        frequencies = counts ./ sum(counts)
        @test all(abs.(frequencies .- 0.25) .< 0.01)
        @test environment.sampled_weight_vectors == 100_000
    end

    @testset "loop erasure" begin
        p(x, y) = LR.pack_point(Int32(x), Int32(y))
        raw = [p(0, 0), p(1, 0), p(1, 1), p(0, 1), p(0, 0), p(0, 1), p(0, 2)]
        @test loop_erase(raw) == [p(0, 0), p(0, 1), p(0, 2)]
        @test loop_erase(LR.PointKey[]) == LR.PointKey[]
        environment = SiteIIDEnvironment(1, "baseline", Dict(); cache_capacity=1)
        path, raw_steps = loop_erased_walk(1, environment, StableRNG(2))
        @test raw_steps == 1
        @test length(path) == 2
    end

    @testset "winding" begin
        p(x, y) = LR.pack_point(Int32(x), Int32(y))
        @test winding([p(0, 0)]) == 0
        @test winding([p(0, 0), p(1, 0), p(1, 1), p(0, 1)]) == 2
        @test winding([p(0, 0), p(1, 0), p(1, -1), p(0, -1)]) == -2
        @test_throws ArgumentError winding([p(0, 0), p(1, 0), p(0, 0)])
        @test_throws ArgumentError LR.direction_code(p(0, 0), p(2, 0))
    end

    @testset "batch reproducibility and schema" begin
        temporal_task = task(environment_model="temporal_iid", num_environments=12)
        rows_1 = run_batch(temporal_task)
        rows_2 = run_batch(temporal_task)
        @test isequal(deterministic_fields.(rows_1), deterministic_fields.(rows_2))
        @test length(rows_1) == 12
        required = (
            :schema_version, :environment_model, :environment_seed, :walk_seed,
            :weight_seed, :direction_seed, :winding, :loop_erased_path_length,
            :raw_walk_length, :north_steps, :east_steps, :south_steps, :west_steps,
            :sampled_site_count, :sampled_weight_vector_count, :status,
        )
        @test all(name -> hasproperty(first(rows_1), name), required)
        @test all(row -> row.schema_version == "batch_v7_julia", rows_1)
        @test all(row -> row.environment_model == "temporal_iid", rows_1)
        @test all(row -> row.walk_id == 0, rows_1)
        @test all(row -> row.sampled_site_count === missing, rows_1)
        @test all(row -> row.sampled_weight_vector_count == row.raw_walk_length, rows_1)
        @test all(row -> row.north_steps + row.east_steps + row.south_steps +
                         row.west_steps == row.raw_walk_length, rows_1)

        paired_task = task(environment_model="site_iid", walks_per_environment=2,
                           num_environments=3)
        paired_rows = run_batch(paired_task)
        @test length(paired_rows) == 6
        @test paired_rows[1].environment_seed == paired_rows[2].environment_seed
        @test paired_rows[1].direction_seed != paired_rows[2].direction_seed
        @test all(row -> row.sampled_site_count !== missing, paired_rows)
    end

    @testset "campaign validation" begin
        @test LR.require_strict_annealed([task(walks_per_environment=1)]) === nothing
        @test_throws ArgumentError LR.require_strict_annealed(
            [task(walks_per_environment=2)])
        @test LR.require_double_dimer([task(walks_per_environment=2)]) === nothing
        @test_throws ArgumentError LR.require_double_dimer(
            [task(walks_per_environment=1)])
        @test_throws ArgumentError LR.require_double_dimer(
            [task(environment_model="temporal_iid", walks_per_environment=2)])
        @test LR.require_temporal_iid(
            [task(environment_model="temporal_iid")]) === nothing
        @test_throws ArgumentError LR.require_temporal_iid([task()])
        @test_throws ArgumentError LR.require_temporal_iid(
            [task(environment_model="temporal_iid", walks_per_environment=2)])
        @test_throws ArgumentError run_batch(
            task(environment_model="temporal_iid", walks_per_environment=2))
        @test_throws ArgumentError run_batch(task(environment_model="unknown"))
        @test_throws ArgumentError LR.validate_environment_model("space_iid")

        mktempdir() do directory
            unknown = joinpath(directory, "unknown.csv")
            CSV.write(unknown, [(
                task_id=0, environment_model="unknown", distribution="baseline",
                distribution_params="{}", L=8, batch_id=0, num_environments=1,
                walks_per_environment=1, base_seed=1,
            )])
            @test_throws ArgumentError LR.load_config_row(unknown, 0)

            invalid_temporal = joinpath(directory, "invalid_temporal.csv")
            CSV.write(invalid_temporal, [(
                task_id=0, environment_model="temporal_iid", distribution="baseline",
                distribution_params="{}", L=8, batch_id=0, num_environments=1,
                walks_per_environment=2, base_seed=1,
            )])
            @test_throws ArgumentError LR.load_config_row(invalid_temporal, 0)
            @test_throws ArgumentError generate_config(
                joinpath(directory, "generated.csv"), ("baseline",), (8,);
                environment_model="temporal_iid", walks_per_environment=2)
        end
        @test occursin("temporal_iid", result_path(
            "results", task(environment_model="temporal_iid")))
        @test occursin("site_iid", result_path("results", task()))
    end

    @testset "Gamma length campaign and path-length scaling" begin
        mktempdir() do directory
            config = joinpath(directory, "gamma_length.csv")
            rows = generate_temporal_iid_gamma_length_config(
                config; extension_walks=10)
            @test length(rows) == 119
            @test sum(row.num_environments for row in rows) == 11_540
            @test sort(unique(row.L for row in rows)) ==
                [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 5000, 8192]
            @test all(row -> row.distribution == "gamma", rows)
            @test all(row -> row.environment_model == "temporal_iid", rows)
        end

        synthetic = [(
            distribution="gamma", distribution_params="""{"shape":0.5}""",
            L=L, log_L=log(L), observations=100,
            mean_loop_erased_path_length=3.0 * L^1.25,
        ) for L in (16, 32, 64, 128)]
        fit = only(LR.fit_path_length_exponents(synthetic))
        @test fit.length_exponent ≈ 1.25 atol=1e-12
        pointwise = LR.path_length_pointwise_rows(synthetic)
        @test all(row -> row.mean_length_over_L_to_5_over_4 ≈ 3.0,
                  pointwise)
    end

    @testset "atomic writing and restart-safe detection" begin
        mktempdir() do directory
            path = joinpath(directory, "nested", "batch.csv")
            rows = [(status="ok", value=1), (status="ok", value=2)]
            LR.write_csv_atomic(path, rows)
            @test isfile(path)
            @test LR.completed_result(path, 2)
            @test !LR.completed_result(path, 3)
            @test isempty(filter(name -> occursin(".tmp.", name), readdir(dirname(path))))
            LR.write_csv_atomic(path, [(status="failed", value=1)])
            @test !LR.completed_result(path, 1)
        end
    end

    @testset "end-to-end temporal campaign and analysis" begin
        mktempdir() do directory
            config = joinpath(directory, "config.csv")
            results = joinpath(directory, "results")
            analysis = joinpath(directory, "analysis")
            rows = generate_config(config, ("baseline", "gamma:0.5"), (4, 8, 16);
                batches=1, num_environments=20, walks_per_environment=1,
                base_seed=20260726, environment_model="temporal_iid")
            @test length(rows) == 6
            first_report = LR.run_campaign(config, results)
            @test first_report.new_tasks == 6
            tasks = LR.read_config(config)
            first_path = result_path(results, first(tasks))
            @test completed_result(first_path, first(tasks))
            stale_task = BatchConfig(999, first(tasks).distribution,
                first(tasks).distribution_params, first(tasks).L,
                first(tasks).batch_id, first(tasks).num_environments,
                first(tasks).walks_per_environment, first(tasks).base_seed,
                first(tasks).environment_model)
            @test !completed_result(first_path, stale_task)
            second_report = LR.run_campaign(config, results)
            @test second_report.already_complete == 6
            @test second_report.new_tasks == 0

            report = analyze_results(config, results, analysis;
                bootstrap_reps=20, bootstrap_seed=20260726)
            @test report.raw_walks == 120
            @test report.summary_rows == 6
            summary = collect(CSV.File(joinpath(analysis, "summary.csv")))
            @test length(summary) == 6
            @test :quenched_variance ∉ propertynames(first(summary))
            @test all(row -> row.environment_model == "temporal_iid", summary)
            @test all(row -> row.total_raw_steps ==
                             row.total_north_steps + row.total_east_steps +
                             row.total_south_steps + row.total_west_steps, summary)
            fits = collect(CSV.File(joinpath(analysis, "loglog_fits.csv")))
            @test length(fits) == 2
            @test all(row -> row.variance_kind == "annealed", fits)
            comparisons = collect(CSV.File(
                joinpath(analysis, "temporal_baseline_comparisons.csv")))
            @test length(comparisons) == 3
            @test all(row -> isfinite(row.variance_ratio), comparisons)
            diagnostics = collect(CSV.File(
                joinpath(analysis, "temporal_direction_diagnostics.csv")))
            @test length(diagnostics) == 6
            @test all(row -> isfinite(row.max_abs_frequency_deviation), diagnostics)
            length_fits = collect(CSV.File(
                joinpath(analysis, "path_length_fits.csv")))
            @test length(length_fits) == 2
            @test all(row -> isfinite(row.length_exponent), length_fits)
            @test isfile(joinpath(analysis, "path_length_pointwise.csv"))
            @test isfile(joinpath(
                analysis, "path_length_local_effective_exponents.csv"))
        end
    end
end
