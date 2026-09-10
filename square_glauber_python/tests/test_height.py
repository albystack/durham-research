import numpy as np
import pytest

from square_glauber.environment import sample_environment
from square_glauber.geometry import center_face_index, spatial_separations
from square_glauber.glauber import heatbath_update
from square_glauber.height import (
    face_heights,
    height_constraints,
    spatial_increments,
    validate_height_rules,
)
from square_glauber.matching import DimerState


@pytest.mark.parametrize("L", [2, 4, 6, 8])
@pytest.mark.parametrize("start", [DimerState.all_horizontal, DimerState.all_vertical])
def test_height_is_path_consistent_and_obeys_occupancy_rule(L: int, start) -> None:
    state = start(L)
    heights = face_heights(state)
    assert heights.shape == (L - 1, L - 1)
    assert validate_height_rules(state, heights)
    # All constrained crossings have exactly the supervisor-specified magnitude.
    for item in height_constraints(state):
        if item.left_face is not None and item.right_face is not None:
            difference = int(heights[item.right_face]) - int(heights[item.left_face])
            assert abs(difference) == (3 if item.occupied else 1)


def test_fixed_boundary_cut_retains_additive_fluctuation() -> None:
    horizontal = DimerState.all_horizontal(2)
    vertical = DimerState.all_vertical(2)
    assert int(face_heights(horizontal)[0, 0]) == -3
    assert int(face_heights(vertical)[0, 0]) == 1
    assert abs(int(face_heights(horizontal)[0, 0] - face_heights(vertical)[0, 0])) == 4


def test_extremal_L6_matrices_are_frozen_for_convention_audit() -> None:
    expected_horizontal = np.array(
        [
            [-3, -2, -3, -2, -3],
            [0, -1, 0, -1, 0],
            [-3, -2, -3, -2, -3],
            [0, -1, 0, -1, 0],
            [-3, -2, -3, -2, -3],
        ]
    )
    expected_vertical = np.array(
        [
            [1, -2, 1, -2, 1],
            [0, -1, 0, -1, 0],
            [1, -2, 1, -2, 1],
            [0, -1, 0, -1, 0],
            [1, -2, 1, -2, 1],
        ]
    )
    np.testing.assert_array_equal(face_heights(DimerState.all_horizontal(6)), expected_horizontal)
    np.testing.assert_array_equal(face_heights(DimerState.all_vertical(6)), expected_vertical)


def test_every_legal_local_flip_changes_only_selected_face_by_four() -> None:
    L = 6
    state = DimerState.all_horizontal(L)
    environment = sample_environment(L, np.random.default_rng(201), "gamma", gamma_shape=0.5)
    rng = np.random.default_rng(202)
    for _ in range(100):
        before = face_heights(state)
        for i in range(L - 1):
            for j in range(L - 1):
                if state.is_face_flippable(i, j):
                    flipped = state.copy()
                    flipped.flip_face(i, j)
                    difference = face_heights(flipped) - before
                    changed = np.argwhere(difference != 0)
                    assert changed.tolist() == [[i, j]]
                    assert abs(int(difference[i, j])) == 4
                    assert flipped.validate()
        heatbath_update(state, environment, rng)


def test_center_and_spatial_increments() -> None:
    L = 16
    state = DimerState.all_horizontal(L)
    heights = face_heights(state)
    assert center_face_index(L) == (7, 7)
    separations = spatial_separations(L)
    assert [item.r for item in separations] == [1, 2, 4]
    increments = spatial_increments(heights, separations)
    assert set(increments) == {1, 2, 4}
    center_row, center_column = center_face_index(L)
    for r, value in increments.items():
        assert value == int(heights[center_row, center_column + r] - heights[center_row, center_column])
