const GSG = AztecDiamond.GlauberSquareGrid

include(joinpath(@__DIR__, "..", "scripts", "analyze_glauber_square_grid_scaling.jl"))
const GPSA = GlauberProductionScaling
include(joinpath(@__DIR__, "..", "scripts", "analyze_glauber_kasteleyn_campaign.jl"))
const GKA = GlauberKasteleynAnalysis

@testset "Glauber production environment-blocked scaling" begin
    rows = [
        (replica_1_mean=value, replica_2_mean=2value,
         conditional_variance=3.0) for value in 1.0:3.0
    ]
    estimate = GPSA.component_statistics(rows)
    @test estimate.conditional == 3.0
    @test estimate.disorder == 2.0
    @test estimate.total == 5.0

    sizes = [2, 4, 6, 8, 12, 16]
    x = log.(Float64.(sizes))
    exactly_quadratic = 1 .+ 0.4 .* x .+ 0.75 .* x .^ 2
    comparison = GPSA.fit_comparison(x, exactly_quadratic)
    @test isapprox(comparison.quadratic.coefficients[3], 0.75; atol=1e-11)
    @test comparison.delta_bic > 0
    @test comparison.quadratic_loocv_rmse < 1e-10

    exactly_log = 2 .+ 1.25 .* x
    null_comparison = GPSA.fit_comparison(x, exactly_log)
    @test abs(null_comparison.quadratic.coefficients[3]) < 1e-11
    @test null_comparison.delta_bic < 0
    @test null_comparison.log_loocv_rmse < 1e-10

    paired = [
        (replica_1_mean=Float64(index), replica_2_mean=Float64(index),
         conditional_variance=Float64(index)) for index in 1:8
    ]
    draws = GPSA.bootstrap_size_draws(Xoshiro(1234), Dict(4 => paired), 20)[4]
    @test length(draws.total) == 20
    @test all(draws.total .== draws.conditional .+ draws.disorder)
end

@testset "Kasteleyn/MCMC comparison joins exact environment identities" begin
    exact = [
        (L=4, environment_id=1, environment_seed=UInt64(11),
         conditional_mean=1.0, conditional_variance=2.0),
        (L=4, environment_id=2, environment_seed=UInt64(12),
         conditional_mean=2.0, conditional_variance=3.0),
    ]
    mcmc = [
        (L=4, environment_id=1, environment_seed=UInt64(11),
         replica_1_mean=0.8, replica_2_mean=1.0,
         replica_1_variance=2.0, replica_2_variance=2.0,
         replica_1_ess=100.0, replica_2_ess=100.0,
         conditional_variance=2.2, start_gap=-0.2,
         standardized_start_gap=1.0, minimum_pair_swap_acceptance=0.2,
         target_exchange_acceptance=0.8),
        (L=4, environment_id=2, environment_seed=UInt64(12),
         replica_1_mean=2.0, replica_2_mean=2.0,
         replica_1_variance=0.0, replica_2_variance=0.0,
         replica_1_ess=100.0, replica_2_ess=100.0,
         conditional_variance=2.8, start_gap=0.0,
         standardized_start_gap=0.0, minimum_pair_swap_acceptance=0.3,
         target_exchange_acceptance=0.9),
    ]
    comparison = GKA.comparison_rows(exact, mcmc)
    @test length(comparison) == 2
    @test isapprox(comparison[1].mcmc_minus_exact_mean, -0.1; atol=1e-12)
    @test comparison[1].estimated_mcmc_mean_standard_error > 0
    @test !comparison[1].zero_standard_error_inconsistent
    @test comparison[2].standardized_mean_error == 0
    @test !comparison[2].zero_standard_error_inconsistent
    summary = only(GKA.comparison_summary(comparison))
    @test summary.L == 4
    @test summary.environments == 2
end

@testset "square-grid Glauber height reference" begin
    @testset "extremal height boundaries encode valid dimer matchings" begin
        for L in 1:4
            maximum = GSG.max_height_configuration(L)
            minimum = GSG.min_height_configuration(L)
            max_validation = GSG.validate_height_configuration(maximum, L)
            min_validation = GSG.validate_height_configuration(minimum, L)
            @test max_validation.valid
            @test min_validation.valid
            @test max_validation.dimer_edges == 2L^2
            @test min_validation.dimer_edges == 2L^2
            @test maximum[1, :] == minimum[1, :]
            @test maximum[end, :] == minimum[end, :]
            @test maximum[:, 1] == minimum[:, 1]
            @test maximum[:, end] == minimum[:, end]
        end
    end

    @testset "tiny state spaces match domino-tiling counts" begin
        @test length(GSG.enumerate_height_configurations(1)) == 2
        @test length(GSG.enumerate_height_configurations(2)) == 36
    end

    @testset "Kasteleyn moments match exhaustive weighted enumeration" begin
        for L in 1:2, seed in 100_020:100_023
            weights = GSG.random_edge_weights(
                Xoshiro(seed + L), L; distribution=:gamma, parameter=0.7)
            states = GSG.enumerate_height_configurations(L)
            partition = sum(GSG.matching_weight(state, weights) for state in states)
            exact = GSG.exact_center_distribution(L, weights)
            determinantal = GSG.center_height_moments_kasteleyn(weights)

            @test isapprox(exp(determinantal.log_partition), partition;
                           rtol=2e-12, atol=1e-14)
            @test isapprox(determinantal.mean, exact.mean; rtol=2e-12, atol=2e-12)
            @test isapprox(determinantal.variance, exact.variance;
                           rtol=2e-12, atol=2e-12)
            @test determinantal.crossed_edges == L
            @test determinantal.matrix_order == 2L^2
            @test all(probability -> 0 <= probability <= 1,
                      determinantal.edge_probabilities)
            @test determinantal.relative_solve_residual < 1e-12
            @test determinantal.maximum_probability_imaginary_residual < 1e-12
            @test determinantal.maximum_covariance_imaginary_residual < 1e-12

            reverse_path = reverse(determinantal.path)
            reverse_difference = GSG.height_difference_moments_kasteleyn(
                weights, reverse_path)
            forward_mean = determinantal.mean - determinantal.boundary_height
            @test isapprox(reverse_difference.mean, -forward_mean; atol=2e-12)
            @test isapprox(reverse_difference.variance, determinantal.variance;
                           atol=2e-12)
        end

        uniform = GSG.center_height_moments_kasteleyn(GSG.constant_edge_weights(2))
        @test isapprox(exp(uniform.log_partition), 36.0; atol=1e-11)
        high_precision = setprecision(128) do
            GSG.center_height_moments_kasteleyn(
                GSG.constant_edge_weights(2); number_type=BigFloat)
        end
        @test isapprox(high_precision.mean, 0; atol=big"1e-35")
        @test isapprox(high_precision.variance, big(8) / 9; atol=big"1e-35")
        @test high_precision.relative_solve_residual < big"1e-35"
        @test_throws ArgumentError GSG.height_difference_moments_kasteleyn(
            GSG.constant_edge_weights(2), [(1, 2), (2, 3)])
        @test_throws ArgumentError GSG.height_difference_moments_kasteleyn(
            GSG.constant_edge_weights(2), [(1, 1), (2, 1)])

        moderate = GSG.center_height_moments_kasteleyn(
            GSG.random_edge_weights(Xoshiro(100_024), 6;
                                    distribution=:gamma, parameter=0.5))
        @test isfinite(moderate.mean)
        @test moderate.variance >= 0
        @test moderate.relative_solve_residual < 1e-11
    end

    @testset "local heat bath has the specified ac versus bd probability" begin
        height = GSG.max_height_configuration(1)
        weights = GSG.EdgeWeights(ones(2, 3), ones(3, 2))
        # At L=1, the sole selected face is flippable and uniform weights give
        # an exact half/half local heat bath.
        samples = 20_000
        ac_count = 0
        rng = Xoshiro(100_001)
        for _ in 1:samples
            height .= GSG.max_height_configuration(1)
            GSG.heatbath_update!(rng, height, weights)
            configurations = GSG.local_configurations(height, 2, 2)
            @test !isnothing(configurations)
            ac_count += height[2, 2] == configurations.ac
        end
        @test abs(ac_count / samples - 0.5) < 0.02

        # With a*c=6 and b*d=1 the first configuration has probability 6/7.
        biased = GSG.EdgeWeights([1.0 6.0 1.0; 1.0 1.0 1.0],
                                 [1.0 1.0; 1.0 1.0; 1.0 1.0])
        # Top and bottom are vertical[1,2] and vertical[2,2].
        biased.vertical[1, 2] = 2.0
        biased.vertical[2, 2] = 3.0
        biased.horizontal[2, 1] = 1.0
        biased.horizontal[2, 2] = 1.0
        ac_count = 0
        rng = Xoshiro(100_002)
        for _ in 1:samples
            height .= GSG.max_height_configuration(1)
            GSG.heatbath_update!(rng, height, biased)
            configurations = GSG.local_configurations(height, 2, 2)
            ac_count += height[2, 2] == configurations.ac
        end
        @test abs(ac_count / samples - 6 / 7) < 0.015
    end

    @testset "weighted Gibbs detailed balance is exact on L=2" begin
        weights = GSG.random_edge_weights(Xoshiro(100_003), 2; distribution=:gamma, parameter=0.7)
        @test GSG.detailed_balance_residual(2, weights) < 1e-12
        exact = GSG.exact_center_distribution(2, weights)
        @test isapprox(sum(values(exact.masses)), 1.0; atol=1e-12)
        @test exact.variance >= 0
    end

    @testset "replica exchange preserves the tempered weighted law" begin
        weights = GSG.random_edge_weights(Xoshiro(100_004), 2; distribution=:gamma, parameter=0.7)
        @test GSG.tempered_edge_weights(weights, 0.0).vertical == ones(4, 5)
        @test GSG.tempered_edge_weights(weights, 0.0).horizontal == ones(5, 4)
        states = GSG.enumerate_height_configurations(2)
        beta_low, beta_high = 0.25, 1.0
        for left in states, right in states
            log_ratio = GSG.replica_exchange_log_ratio(
                left, beta_low, right, beta_high, weights)
            forward = exp(beta_low * GSG.matching_log_weight(left, weights) +
                          beta_high * GSG.matching_log_weight(right, weights)) *
                      min(1.0, exp(log_ratio))
            reverse = exp(beta_low * GSG.matching_log_weight(right, weights) +
                          beta_high * GSG.matching_log_weight(left, weights)) *
                      min(1.0, exp(-log_ratio))
            @test isapprox(forward, reverse; rtol=1e-12, atol=1e-12)
        end

        start_left = GSG.max_height_configuration(2)
        start_right = GSG.min_height_configuration(2)
        beta = 0.5
        left_chain = GSG.accelerated_chain(copy(start_left), GSG.tempered_edge_weights(weights, beta))
        right_chain = GSG.accelerated_chain(copy(start_right), GSG.tempered_edge_weights(weights, beta))
        report = GSG.parallel_tempering_swap!(
            Xoshiro(100_005), left_chain, beta, right_chain, beta, weights)
        @test report.accepted
        @test left_chain.height == start_right
        @test right_chain.height == start_left
        @test GSG.validate_height_configuration(left_chain.height, 2).valid
        @test GSG.validate_height_configuration(right_chain.height, 2).valid

        # The exchange clock and alternating-pair parity must survive calls that
        # are shorter than the exchange interval.  Otherwise the beta=1 pair is
        # never attempted when each retained-sample call resets the schedule.
        betas = [0.0, 0.5, 1.0]
        chains = [GSG.accelerated_chain(GSG.max_height_configuration(2),
                                        GSG.tempered_edge_weights(weights, value))
                  for value in betas]
        before_boundary = GSG.run_parallel_tempering_updates!(
            Xoshiro(100_005), chains, betas, weights, 2;
            swap_interval_attempts=5)
        @test before_boundary.attempted_swaps_by_pair == [0, 0]
        @test before_boundary.attempts_since_swap == 2
        first_parity = GSG.run_parallel_tempering_updates!(
            Xoshiro(100_005), chains, betas, weights, 3;
            swap_interval_attempts=5, swap_round_offset=before_boundary.swap_round,
            attempts_since_swap=before_boundary.attempts_since_swap)
        @test first_parity.attempted_swaps_by_pair == [1, 0]
        second_parity = GSG.run_parallel_tempering_updates!(
            Xoshiro(100_005), chains, betas, weights, 5;
            swap_interval_attempts=5, swap_round_offset=first_parity.swap_round,
            attempts_since_swap=first_parity.attempts_since_swap)
        @test second_parity.attempted_swaps_by_pair == [0, 1]
        @test second_parity.attempts_since_swap == 0

        tempered_first = GSG.sample_center_height_chain_parallel_tempering(
            Xoshiro(100_006), 2, weights; start=:max, betas=[0.0, 0.5, 1.0],
            burn_in_attempts=100, thin_attempts=10, samples=20,
            swap_interval_attempts=5)
        tempered_second = GSG.sample_center_height_chain_parallel_tempering(
            Xoshiro(100_006), 2, weights; start=:max, betas=[0.0, 0.5, 1.0],
            burn_in_attempts=100, thin_attempts=10, samples=20,
            swap_interval_attempts=5)
        @test tempered_first.heights == tempered_second.heights
        @test 0 <= tempered_first.diagnostics.swap_acceptance_rate <= 1
        @test all(>(0), tempered_first.diagnostics.attempted_swaps_by_pair)
        @test tempered_first.diagnostics.target_exchange_attempts > 0
        @test GSG.validate_height_configuration(tempered_first.final_height, 2).valid
    end

    @testset "rejection-free acceleration preserves literal update time" begin
        weights = GSG.constant_edge_weights(1)
        changed = 0
        rng = Xoshiro(100_050)
        for _ in 1:20_000
            chain = GSG.accelerated_chain(GSG.max_height_configuration(1), weights)
            GSG.run_accelerated_updates!(rng, chain, 1)
            changed += chain.height[2, 2] != GSG.max_height_configuration(1)[2, 2]
        end
        # One literal L=1 update changes the configuration with probability 1/2.
        @test abs(changed / 20_000 - 0.5) < 0.02

        chain = GSG.accelerated_chain(GSG.max_height_configuration(2),
                                      GSG.random_edge_weights(Xoshiro(100_051), 2))
        GSG.run_accelerated_updates!(Xoshiro(100_052), chain, 10_000)
        @test GSG.validate_height_configuration(chain.height, 2).valid

        order_two_weights = GSG.constant_edge_weights(2)
        first = GSG.sample_center_height_chain_accelerated(Xoshiro(100_053), 2, order_two_weights;
            burn_in_attempts=100, thin_attempts=10, samples=12)
        second = GSG.sample_center_height_chain_accelerated(Xoshiro(100_053), 2, order_two_weights;
            burn_in_attempts=100, thin_attempts=10, samples=12)
        @test first.heights == second.heights
    end

    @testset "short trace diagnostics and extremal starts are explicit" begin
        weights = GSG.constant_edge_weights(2)
        comparison = GSG.compare_extremal_starts(Xoshiro(100_101), Xoshiro(100_102), 2, weights;
            burn_in_attempts=500, thin_attempts=20, samples=40)
        @test isfinite(comparison.maximum.diagnostics.integrated_autocorrelation_time)
        @test 1 <= comparison.maximum.diagnostics.effective_sample_size <= 40
        @test 0 <= comparison.maximum.diagnostics.changed_rate <= 1
        @test GSG.validate_height_configuration(comparison.maximum.final_height, 2).valid
        @test GSG.validate_height_configuration(comparison.minimum.final_height, 2).valid
    end

    @testset "seeded paired chains are reproducible and preserve the environment" begin
        first = GSG.sample_center_height_pair(100_007, 2;
            distribution=:gamma, parameter=0.5,
            burn_in_attempts=100, thin_attempts=10, samples=12)
        second = GSG.sample_center_height_pair(100_007, 2;
            distribution=:gamma, parameter=0.5,
            burn_in_attempts=100, thin_attempts=10, samples=12)
        @test first.environment_seed == second.environment_seed
        @test first.replica_1_seed != first.replica_2_seed
        @test first.replica_1.heights == second.replica_1.heights
        @test first.replica_2.heights == second.replica_2.heights
        @test GSG.validate_height_configuration(first.replica_1.final_height, 2).valid
        @test GSG.validate_height_configuration(first.replica_2.final_height, 2).valid
    end
end
