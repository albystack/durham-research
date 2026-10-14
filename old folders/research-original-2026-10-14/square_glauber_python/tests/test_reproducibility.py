import numpy as np
import pytest

from square_glauber.environment import sample_environment
from square_glauber.experiment import derive_seeds, run_paired_environment
from square_glauber.glauber import run_chain
from square_glauber.matching import DimerState


def test_environment_and_chain_reproducibility() -> None:
    first_environment = sample_environment(6, np.random.default_rng(301), "gamma", gamma_shape=0.8)
    second_environment = sample_environment(6, np.random.default_rng(301), "gamma", gamma_shape=0.8)
    np.testing.assert_array_equal(
        first_environment.horizontal_weights, second_environment.horizontal_weights
    )
    np.testing.assert_array_equal(first_environment.vertical_weights, second_environment.vertical_weights)
    kwargs = dict(
        environment=first_environment,
        initial_matching=DimerState.all_horizontal(6),
        burnin_sweeps=20,
        measurement_gap_sweeps=3,
        num_measurements=5,
    )
    first = run_chain(rng=np.random.default_rng(302), **kwargs)
    second = run_chain(rng=np.random.default_rng(302), **kwargs)
    assert first.center_heights == second.center_heights
    np.testing.assert_array_equal(
        first.final_state.horizontal_occupied, second.final_state.horizontal_occupied
    )
    np.testing.assert_array_equal(first.final_state.vertical_occupied, second.final_state.vertical_occupied)
    assert first.diagnostics() | {"elapsed_seconds": 0} == second.diagnostics() | {
        "elapsed_seconds": 0
    }


def test_environment_arrays_are_immutable() -> None:
    environment = sample_environment(4, np.random.default_rng(303), "gamma")
    with pytest.raises(ValueError):
        environment.horizontal_weights[0, 0] = 2.0
    with pytest.raises(ValueError):
        environment.vertical_weights[0, 0] = 2.0


def test_chain_sweep_accounting_and_optional_trace() -> None:
    environment = sample_environment(4, np.random.default_rng(304), "uniform")
    result = run_chain(
        environment,
        DimerState.all_horizontal(4),
        np.random.default_rng(305),
        burnin_sweeps=3,
        measurement_gap_sweeps=2,
        num_measurements=3,
        record_trace=True,
    )
    expected_sweeps = 3 + 2 * 2
    assert result.counts.attempted_updates == expected_sweeps * (4 - 1) ** 2
    assert result.counts.null_moves == result.counts.attempted_updates - result.counts.actual_changes
    assert len(result.trace) == expected_sweeps + 1
    assert result.trace[0].phase == "initial"
    assert len(result.center_heights) == 3


def test_paired_chains_share_environment_but_have_independent_seeds() -> None:
    first = run_paired_environment(
        L=4,
        model="gamma",
        gamma_shape=1.0,
        environment_index=7,
        master_seed=20260907,
        burnin_sweeps=10,
        measurement_gap_sweeps=2,
        num_measurements=3,
    )
    second = run_paired_environment(
        L=4,
        model="gamma",
        gamma_shape=1.0,
        environment_index=7,
        master_seed=20260907,
        burnin_sweeps=10,
        measurement_gap_sweeps=2,
        num_measurements=3,
    )
    assert first.seeds == derive_seeds(20260907, 4, "gamma", 7)
    assert first.seeds.chain1_seed != first.seeds.chain2_seed
    assert first.chain1.center_heights == second.chain1.center_heights
    assert first.chain2.center_heights == second.chain2.center_heights
    np.testing.assert_array_equal(
        first.environment.horizontal_weights, second.environment.horizontal_weights
    )
