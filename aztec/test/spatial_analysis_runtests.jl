module SpatialAnalyzerHarness
include(joinpath(@__DIR__, "..", "scripts", "analyze_spatial_campaign.jl"))
end
const SGA = SpatialAnalyzerHarness

@testset "spatial analysis validation and block GLS" begin
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
