import numpy as np
import pytest

from square_glauber.accelerated import (
    AcceleratedRandomStreams,
    NUMBA_AVAILABLE,
    horizontal_probability_grid,
    run_attempts_accelerated,
    run_supplied_attempts_accelerated,
)
from square_glauber.environment import sample_environment
from square_glauber.glauber import UpdateCounts, run_supplied_attempts
from square_glauber.height import face_heights
from square_glauber.matching import DimerState


def test_reference_supplied_stream_validation_and_reproducibility() -> None:
    environment = sample_environment(4, np.random.default_rng(100), "gamma")
    faces = np.asarray([0, 1, 4, 8, 3, 2], dtype=np.int64)
    uniforms = np.asarray([0.1, 0.9, 0.5, 0.25, 0.75, 0.4])
    first = DimerState.all_horizontal(4)
    second = first.copy()
    first_counts = run_supplied_attempts(first, environment, faces, uniforms)
    second_counts = run_supplied_attempts(second, environment, faces, uniforms)
    np.testing.assert_array_equal(first.horizontal_occupied, second.horizontal_occupied)
    np.testing.assert_array_equal(first.vertical_occupied, second.vertical_occupied)
    assert first_counts == second_counts
    with pytest.raises(ValueError, match="equal length"):
        run_supplied_attempts(first, environment, faces, uniforms[:-1])
    with pytest.raises(ValueError, match="must lie"):
        run_supplied_attempts(first, environment, np.asarray([9]), np.asarray([0.5]))
    with pytest.raises(ValueError, match=r"\[0,1\)"):
        run_supplied_attempts(first, environment, np.asarray([0]), np.asarray([1.0]))


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="optional Numba backend is not installed")
@pytest.mark.parametrize("L", [4, 6, 8, 12])
@pytest.mark.parametrize("model", ["uniform", "gamma"])
@pytest.mark.parametrize("start_name", ["all_horizontal", "all_vertical"])
def test_numba_matches_reference_at_every_checkpoint(
    L: int, model: str, start_name: str
) -> None:
    environment = sample_environment(L, np.random.default_rng(1000 + L), model)
    stream = np.random.default_rng(2000 + L + int(model == "gamma"))
    attempts = 5_003
    faces = stream.integers(0, (L - 1) ** 2, size=attempts, dtype=np.int64)
    uniforms = stream.random(attempts)
    constructor = getattr(DimerState, start_name)
    reference = constructor(L)
    accelerated = constructor(L)
    reference_counts = UpdateCounts()
    accelerated_counts = UpdateCounts()
    probabilities = horizontal_probability_grid(environment)
    checkpoints = [1, 8, 39, 296, 997, 2048, attempts]
    begin = 0
    for end in checkpoints:
        run_supplied_attempts(
            reference,
            environment,
            faces[begin:end],
            uniforms[begin:end],
            reference_counts,
        )
        run_supplied_attempts_accelerated(
            accelerated,
            environment,
            faces[begin:end],
            uniforms[begin:end],
            accelerated_counts,
            horizontal_probabilities=probabilities,
        )
        np.testing.assert_array_equal(
            reference.horizontal_occupied, accelerated.horizontal_occupied
        )
        np.testing.assert_array_equal(
            reference.vertical_occupied, accelerated.vertical_occupied
        )
        assert reference_counts == accelerated_counts
        begin = end
    reference.validate()
    accelerated.validate()


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="optional Numba backend is not installed")
def test_accelerated_rng_api_is_reproducible_and_preserves_matching() -> None:
    environment = sample_environment(8, np.random.default_rng(3001), "gamma")
    first = DimerState.all_vertical(8)
    second = first.copy()
    first_counts = run_attempts_accelerated(
        first, environment, np.random.default_rng(3002), 25_000
    )
    second_counts = run_attempts_accelerated(
        second, environment, np.random.default_rng(3002), 25_000
    )
    np.testing.assert_array_equal(first.horizontal_occupied, second.horizontal_occupied)
    np.testing.assert_array_equal(first.vertical_occupied, second.vertical_occupied)
    assert first_counts == second_counts
    first.validate()


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="optional Numba backend is not installed")
def test_accelerated_rng_trajectory_is_invariant_to_chunking() -> None:
    environment = sample_environment(8, np.random.default_rng(3101), "gamma")
    first = DimerState.all_horizontal(8)
    second = first.copy()
    first_streams = AcceleratedRandomStreams.from_generator(np.random.default_rng(3102))
    second_streams = AcceleratedRandomStreams.from_generator(np.random.default_rng(3102))
    first_counts = run_attempts_accelerated(
        first, environment, first_streams, 25_000, chunk_attempts=25_000
    )
    second_counts = UpdateCounts()
    for block in (1, 7, 101, 997, 4_003, 19_891):
        run_attempts_accelerated(
            second,
            environment,
            second_streams,
            block,
            second_counts,
            chunk_attempts=113,
        )
    np.testing.assert_array_equal(first.horizontal_occupied, second.horizontal_occupied)
    np.testing.assert_array_equal(first.vertical_occupied, second.vertical_occupied)
    assert first_counts == second_counts


def test_common_randomness_preserves_extremal_height_order_in_reference_cases() -> None:
    for L in (4, 6, 8):
        environment = sample_environment(L, np.random.default_rng(4000 + L), "gamma")
        stream = np.random.default_rng(5000 + L)
        lower = DimerState.all_horizontal(L)
        upper = DimerState.all_vertical(L)
        for _ in range(2_000):
            face = np.asarray(
                [stream.integers(0, (L - 1) ** 2)], dtype=np.int64
            )
            uniform = np.asarray([stream.random()])
            run_supplied_attempts(lower, environment, face, uniform)
            run_supplied_attempts(upper, environment, face, uniform)
            assert np.all(face_heights(lower) <= face_heights(upper))
