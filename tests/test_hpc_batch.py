from __future__ import annotations

import csv
import math
import random
import tempfile
import unittest
from pathlib import Path

from analyze_hpc_results import fit_loglog_exponent, make_local_exponents
from generate_hpc_config import generate_rows, write_config
from run_hpc_batch import completed_result, load_config_row, main as run_hpc_main, result_path
from simulation import (
    SiteDirectionalEnvironment,
    SiteIIDEnvironment,
    is_boundary,
    loop_erase,
    loop_erased_walk_until_boundary,
    winding,
)


class SiteDirectionalEnvironmentTests(unittest.TestCase):
    def test_site_weights_are_positive_and_probabilities_sum_to_one(self) -> None:
        environment = SiteDirectionalEnvironment(123, "gamma_edges", 1.0)

        weights = environment.weights_at((3, -2))
        probabilities = environment.transition_probabilities((3, -2))

        self.assertTrue(all(weight > 0.0 for weight in weights))
        self.assertAlmostEqual(sum(probabilities), 1.0)

    def test_baseline_probabilities_are_quarters(self) -> None:
        environment = SiteDirectionalEnvironment(123, "symmetric", None)

        self.assertEqual(environment.weights_at((0, 0)), (1.0, 1.0, 1.0, 1.0))
        self.assertEqual(environment.transition_probabilities((0, 0)), (0.25, 0.25, 0.25, 0.25))

    def test_site_iid_environment_samples_all_four_directions(self) -> None:
        environment = SiteIIDEnvironment(123, "gamma_edges", 1.0)

        weights = environment.weights_at((0, 0))

        self.assertEqual(len(weights), 4)
        self.assertTrue(all(weight > 0.0 for weight in weights))
        self.assertNotEqual(weights[0], 1.0)
        self.assertNotEqual(weights[1], 1.0)

    def test_fixed_environment_remains_fixed_during_walk(self) -> None:
        environment = SiteDirectionalEnvironment(456, "gamma_edges", 1.0)
        before = environment.weights_at((0, 0))

        loop_erased_walk_until_boundary(4, environment, random.Random(789), max_steps=10_000)
        after = environment.weights_at((0, 0))

        self.assertEqual(before, after)

    def test_different_environment_seeds_give_different_sites(self) -> None:
        first = SiteDirectionalEnvironment(1, "gamma_edges", 1.0).weights_at((0, 0))
        second = SiteDirectionalEnvironment(2, "gamma_edges", 1.0).weights_at((0, 0))

        self.assertNotEqual(first, second)

    def test_walk_seed_reproducibility_is_separate_from_environment_seed(self) -> None:
        env_seed = 77
        walk_seed = 88
        first_path, first_raw = loop_erased_walk_until_boundary(
            5,
            SiteDirectionalEnvironment(env_seed, "symmetric", None),
            random.Random(walk_seed),
            max_steps=10_000,
        )
        second_path, second_raw = loop_erased_walk_until_boundary(
            5,
            SiteDirectionalEnvironment(env_seed, "symmetric", None),
            random.Random(walk_seed),
            max_steps=10_000,
        )

        self.assertEqual(first_path, second_path)
        self.assertEqual(first_raw, second_raw)


class LerwPrimitiveTests(unittest.TestCase):
    def test_loop_erasure_is_self_avoiding(self) -> None:
        path = [(0, 0), (1, 0), (1, 1), (1, 0), (2, 0), (2, 1)]

        erased = loop_erase(path)

        self.assertEqual(len(erased), len(set(erased)))

    def test_winding_counts_hand_constructed_turns(self) -> None:
        path = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)]

        self.assertEqual(winding(path), 2)

    def test_lerw_stops_on_boundary(self) -> None:
        path, raw_steps = loop_erased_walk_until_boundary(
            4,
            SiteDirectionalEnvironment(12, "symmetric", None),
            random.Random(34),
            max_steps=10_000,
        )

        self.assertGreater(raw_steps, 0)
        self.assertTrue(is_boundary(path[-1], 4))


class HpcBatchTests(unittest.TestCase):
    def test_task_id_maps_to_configuration_row(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            config_path = Path(directory) / "tasks.csv"
            rows = generate_rows(
                distributions=["baseline", "gamma:1.0"],
                sizes=[64, 128],
                batches=1,
                num_environments=2,
                walks_per_environment=3,
                base_seed=20260623,
                environment_model="site_iid",
            )
            write_config(config_path, rows)

            task = load_config_row(config_path, 2)

            self.assertEqual(task.task_id, 2)
            self.assertEqual(task.distribution, "gamma")
            self.assertEqual(task.distribution_params, {"shape": 1.0})
            self.assertEqual(task.environment_model, "site_iid")
            self.assertEqual(task.L, 64)

    def test_completed_batch_is_not_overwritten_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config_path = root / "tasks.csv"
            output_dir = root / "results"
            rows = generate_rows(
                distributions=["baseline"],
                sizes=[4],
                batches=1,
                num_environments=1,
                walks_per_environment=2,
                base_seed=100,
                environment_model="site_iid",
            )
            write_config(config_path, rows)

            first_code = run_hpc_main(
                [
                    "--task-id",
                    "0",
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(output_dir),
                ]
            )
            task = load_config_row(config_path, 0)
            path = result_path(output_dir, task)
            first_contents = path.read_text(encoding="utf-8")
            second_code = run_hpc_main(
                [
                    "--task-id",
                    "0",
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(output_dir),
                ]
            )
            second_contents = path.read_text(encoding="utf-8")

            self.assertEqual(first_code, 0)
            self.assertEqual(second_code, 0)
            self.assertEqual(first_contents, second_contents)
            self.assertTrue(completed_result(path, expected_rows=2))

            with path.open(newline="", encoding="utf-8") as handle:
                result_rows = list(csv.DictReader(handle))
            self.assertEqual({row["status"] for row in result_rows}, {"ok"})
            self.assertEqual({row["environment_id"] for row in result_rows}, {"0"})
            self.assertEqual({row["walk_id"] for row in result_rows}, {"0", "1"})


class HpcAnalysisTests(unittest.TestCase):
    def _synthetic_summary(self, power: float, constant: float) -> list[dict[str, float | str | int]]:
        rows: list[dict[str, float | str | int]] = []
        for L in [32, 64, 128, 256, 512, 1024]:
            log_l = math.log(L)
            rows.append(
                {
                    "distribution": "synthetic",
                    "distribution_params": "{}",
                    "L": L,
                    "log_L": log_l,
                    "log_log_L": math.log(log_l),
                    "annealed_variance": constant * (log_l**power),
                    "quenched_variance": constant * (log_l**power),
                }
            )
        return rows

    def test_synthetic_log_l_variance_has_exponent_one(self) -> None:
        fits = fit_loglog_exponent(self._synthetic_summary(power=1.0, constant=3.0), "annealed_variance", None)

        self.assertEqual(len(fits), 1)
        self.assertAlmostEqual(float(fits[0]["p"]), 1.0, places=12)

    def test_synthetic_log_l_squared_variance_has_exponent_two(self) -> None:
        fits = fit_loglog_exponent(self._synthetic_summary(power=2.0, constant=2.0), "annealed_variance", None)

        self.assertEqual(len(fits), 1)
        self.assertAlmostEqual(float(fits[0]["p"]), 2.0, places=12)

    def test_local_effective_exponent_uses_log_log_denominator(self) -> None:
        local_rows = make_local_exponents(self._synthetic_summary(power=2.0, constant=2.0), "annealed_variance")

        self.assertTrue(local_rows)
        for row in local_rows:
            self.assertAlmostEqual(float(row["p_local"]), 2.0, places=12)


if __name__ == "__main__":
    unittest.main()
