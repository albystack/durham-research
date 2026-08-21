const GSG = AztecDiamond.GlauberSquareGrid

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
