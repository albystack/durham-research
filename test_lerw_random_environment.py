import random
import unittest

from lerw_random_environment import (
    RandomEnvironment,
    loop_erase,
    run_trial,
    winding,
)


class LerwRandomEnvironmentTests(unittest.TestCase):
    def test_loop_erase_removes_chronological_loop(self) -> None:
        path = [(0, 0), (1, 0), (1, 1), (1, 0), (2, 0)]

        self.assertEqual(loop_erase(path), [(0, 0), (1, 0), (2, 0)])

    def test_winding_counts_left_minus_right_turns(self) -> None:
        # Directions: E, N, E, N, W
        # Turns: left, right, left, left -> winding = 2.
        path = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)]

        self.assertEqual(winding(path), 2)

    def test_symmetric_probabilities_are_equal(self) -> None:
        env = RandomEnvironment(random.Random(1), model="symmetric")

        self.assertEqual(env.transition_probabilities((0, 0)), (0.25, 0.25, 0.25, 0.25))

    def test_gamma_probabilities_are_positive_and_normalised(self) -> None:
        env = RandomEnvironment(random.Random(1), model="gamma", k=1.0)

        probabilities = env.transition_probabilities((0, 0))

        self.assertTrue(all(probability > 0 for probability in probabilities))
        self.assertAlmostEqual(sum(probabilities), 1.0)

    def test_gamma4_probabilities_are_positive_and_normalised(self) -> None:
        env = RandomEnvironment(random.Random(1), model="gamma4", k=1.0)

        probabilities = env.transition_probabilities((0, 0))

        self.assertTrue(all(probability > 0 for probability in probabilities))
        self.assertAlmostEqual(sum(probabilities), 1.0)

    def test_environment_is_fixed_after_first_sample(self) -> None:
        env = RandomEnvironment(random.Random(1), model="gamma", k=1.0)

        first = env.weights_at((2, 3))
        second = env.weights_at((2, 3))

        self.assertEqual(first, second)
        self.assertEqual(env.sampled_site_count, 1)

    def test_gamma4_environment_is_fixed_after_first_sample(self) -> None:
        env = RandomEnvironment(random.Random(1), model="gamma4", k=1.0)

        first = env.weights_at((2, 3))
        second = env.weights_at((2, 3))

        self.assertEqual(first, second)
        self.assertEqual(len(first), 4)
        self.assertEqual(env.sampled_site_count, 1)

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


if __name__ == "__main__":
    unittest.main()
