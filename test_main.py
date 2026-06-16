import random
import unittest

from main import (
    RandomEnvironment,
    grouped_stats,
    loop_erase,
    parse_model_config,
    run_trial,
    winding,
)


class LerwProjectTests(unittest.TestCase):
    """Small tests that protect the main modelling assumptions."""

    def test_parse_model_config_accepts_plain_and_parameterised_models(self) -> None:
        self.assertEqual(parse_model_config("symmetric"), ("symmetric", None))
        self.assertEqual(parse_model_config("gamma_balanced:1"), ("gamma_balanced", 1.0))

    def test_loop_erase_removes_chronological_loop(self) -> None:
        path = [(0, 0), (1, 0), (1, 1), (1, 0), (2, 0)]

        self.assertEqual(loop_erase(path), [(0, 0), (1, 0), (2, 0)])

    def test_winding_counts_left_minus_right_turns(self) -> None:
        # Directions: E, N, E, N, W.
        # Turns: left, right, left, left, so winding is 2.
        path = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)]

        self.assertEqual(winding(path), 2)

    def test_symmetric_probabilities_are_equal(self) -> None:
        env = RandomEnvironment(random.Random(1), model="symmetric")

        self.assertEqual(env.transition_probabilities((0, 0)), (0.25, 0.25, 0.25, 0.25))

    def test_random_environment_is_fixed_after_first_sample(self) -> None:
        env = RandomEnvironment(random.Random(1), model="gamma4", k=1.0)

        first = env.weights_at((2, 3))
        second = env.weights_at((2, 3))

        self.assertEqual(first, second)
        self.assertEqual(env.sampled_site_count, 1)

    def test_different_sites_receive_different_samples(self) -> None:
        env = RandomEnvironment(random.Random(1), model="gamma_balanced", k=1.0)

        first_site = env.weights_at((0, 0))
        second_site = env.weights_at((1, 0))

        self.assertNotEqual(first_site, second_site)
        self.assertEqual(env.sampled_site_count, 2)

    def test_opposite_pair_models_are_locally_balanced(self) -> None:
        for model, k in [
            ("gamma_balanced", 1.0),
            ("lognormal_balanced", 1.0),
            ("pareto_balanced", 2.0),
            ("uniform_balanced", 0.5),
        ]:
            with self.subTest(model=model):
                env = RandomEnvironment(random.Random(1), model=model, k=k)

                p_n, p_e, p_s, p_w = env.transition_probabilities((0, 0))

                self.assertAlmostEqual(p_n, p_s)
                self.assertAlmostEqual(p_e, p_w)

    def test_small_symmetric_trial_hits_boundary(self) -> None:
        result = run_trial(
            L=4,
            model="symmetric",
            k=None,
            sample_index=0,
            seed=123,
            max_steps=10_000,
        )

        self.assertGreater(result.raw_steps, 0)
        self.assertGreater(result.lerw_steps, 0)
        self.assertTrue(abs(result.hit_x) == 4 or abs(result.hit_y) == 4)

    def test_grouped_stats_summarises_rows(self) -> None:
        rows = [
            {
                "model": "symmetric",
                "k": "",
                "L": "4",
                "winding": "1",
                "raw_steps": "10",
                "lerw_steps": "5",
                "hit_x": "4",
                "hit_y": "0",
            },
            {
                "model": "symmetric",
                "k": "",
                "L": "4",
                "winding": "3",
                "raw_steps": "12",
                "lerw_steps": "6",
                "hit_x": "0",
                "hit_y": "4",
            },
        ]

        stats = grouped_stats(rows)

        self.assertEqual(len(stats), 1)
        self.assertEqual(stats[0].samples, 2)
        self.assertEqual(stats[0].mean_winding, 2)


if __name__ == "__main__":
    unittest.main()
