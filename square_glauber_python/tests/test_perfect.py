import numpy as np
import pytest

from square_glauber.accelerated import NUMBA_AVAILABLE
from square_glauber.environment import sample_environment
from square_glauber.exact import enumerate_perfect_matchings, rectangular_face_heights
from square_glauber.extrema import (
    construct_extremal_state,
    face_flip_height_delta,
    true_extremal_states,
)
from square_glauber.height import face_heights
from square_glauber.matching import DimerState
from square_glauber.order_audit import local_pattern_monotonicity_truth_table
from square_glauber.perfect import CFTPNotCoalesced, MonotoneCFTPSampler


def test_flip_height_direction_matches_reconstructed_height() -> None:
    state = DimerState.all_horizontal(6)
    for i in range(5):
        for j in range(5):
            orientation = state.face_orientation(i, j)
            if orientation is None:
                continue
            before = face_heights(state)
            flipped = state.copy()
            flipped.flip_face(i, j)
            observed = int(face_heights(flipped)[i, j] - before[i, j])
            assert observed == face_flip_height_delta(i, j, orientation)


@pytest.mark.parametrize("L", [2, 4])
def test_constructed_extrema_equal_exact_pointwise_extrema(L: int) -> None:
    matchings = enumerate_perfect_matchings(L, L)
    exact_heights = np.stack(
        [rectangular_face_heights(item, L, L) for item in matchings]
    )
    minimum, maximum = true_extremal_states(L)
    np.testing.assert_array_equal(face_heights(minimum.state), exact_heights.min(axis=0))
    np.testing.assert_array_equal(face_heights(maximum.state), exact_heights.max(axis=0))


def test_extremal_construction_is_independent_of_valid_start() -> None:
    for direction in ("minimum", "maximum"):
        horizontal = construct_extremal_state(
            6, direction, initial_state=DimerState.all_horizontal(6)
        )
        vertical = construct_extremal_state(
            6, direction, initial_state=DimerState.all_vertical(6)
        )
        np.testing.assert_array_equal(
            horizontal.state.horizontal_occupied,
            vertical.state.horizontal_occupied,
        )
        np.testing.assert_array_equal(
            horizontal.state.vertical_occupied,
            vertical.state.vertical_occupied,
        )


def test_complete_local_pattern_truth_table_preserves_order() -> None:
    table = local_pattern_monotonicity_truth_table()
    assert len(table) == 496
    assert table["order_preserved"].all()
    assert set(table["center_gap_before"]) == {0, 4, 8}
    assert set(table["heatbath_choice"]) == {"horizontal", "vertical"}


def test_reference_cftp_is_reproducible_and_refuses_uncertified_state() -> None:
    environment = sample_environment(4, np.random.default_rng(6001), "gamma")
    sampler = MonotoneCFTPSampler(
        environment, backend="reference", chunk_attempts=23
    )
    first = sampler.sample(6002, max_horizon_sweeps=2048)
    second = sampler.sample(6002, max_horizon_sweeps=2048)
    assert first.certified and second.certified
    assert first.horizon_sweeps == second.horizon_sweeps
    assert first.blocks == second.blocks
    np.testing.assert_array_equal(
        first.state.horizontal_occupied, second.state.horizontal_occupied
    )
    np.testing.assert_array_equal(
        first.state.vertical_occupied, second.state.vertical_occupied
    )
    with pytest.raises(CFTPNotCoalesced) as failure:
        sampler.sample(6002, max_horizon_sweeps=1)
    assert second.iterations[: len(failure.value.iterations)] == failure.value.iterations


def test_perfect_pair_requires_independent_seeds() -> None:
    environment = sample_environment(2, np.random.default_rng(6101), "uniform")
    sampler = MonotoneCFTPSampler(environment, backend="reference")
    with pytest.raises(ValueError, match="distinct"):
        sampler.sample_pair(7, 7)
    pair = sampler.sample_pair(7, 8, max_horizon_sweeps=128)
    assert pair.first.certified and pair.second.certified


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="optional Numba backend is not installed")
def test_reference_and_numba_cftp_are_identical() -> None:
    # Keep the ordinary suite fast and horizon-insensitive.  Larger grids are
    # already covered checkpoint-by-checkpoint in test_accelerated.py; this
    # test exercises the CFTP replay/doubling wrapper itself.
    environment = sample_environment(4, np.random.default_rng(6201), "gamma")
    reference = MonotoneCFTPSampler(
        environment, backend="reference", chunk_attempts=197
    ).sample(6202, max_horizon_sweeps=4096)
    accelerated = MonotoneCFTPSampler(
        environment, backend="numba", chunk_attempts=113
    ).sample(6202, max_horizon_sweeps=4096)
    assert reference.horizon_sweeps == accelerated.horizon_sweeps
    assert reference.blocks == accelerated.blocks
    assert reference.lower_update_counts == accelerated.lower_update_counts
    assert reference.upper_update_counts == accelerated.upper_update_counts
    np.testing.assert_array_equal(
        reference.state.horizontal_occupied,
        accelerated.state.horizontal_occupied,
    )
    np.testing.assert_array_equal(
        reference.state.vertical_occupied,
        accelerated.state.vertical_occupied,
    )
