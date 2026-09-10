import numpy as np
import pytest

from square_glauber.environment import sample_environment
from square_glauber.glauber import heatbath_update
from square_glauber.matching import DimerState


@pytest.mark.parametrize("L", [2, 4, 6, 12])
def test_extremal_matchings_validate(L: int) -> None:
    assert DimerState.all_horizontal(L).validate()
    assert DimerState.all_vertical(L).validate()


def test_odd_L_is_refused() -> None:
    with pytest.raises(ValueError, match="even"):
        DimerState.all_horizontal(5)


def test_validator_checks_shapes_dtype_count_and_degree() -> None:
    valid = DimerState.all_horizontal(4)
    wrong_shape = DimerState(4, np.zeros((4, 4), dtype=bool), valid.vertical_occupied.copy())
    with pytest.raises(ValueError, match="shape"):
        wrong_shape.validate()
    wrong_dtype = DimerState(
        4, valid.horizontal_occupied.astype(np.int8), valid.vertical_occupied.copy()
    )
    with pytest.raises(TypeError, match="boolean"):
        wrong_dtype.validate()
    missing = valid.copy()
    missing.horizontal_occupied[0, 0] = False
    with pytest.raises(ValueError, match="dimers"):
        missing.validate()
    bad_degree = valid.copy()
    bad_degree.horizontal_occupied[0, 0] = False
    bad_degree.horizontal_occupied[0, 1] = True
    with pytest.raises(ValueError, match="vertex"):
        bad_degree.validate()


def test_face_orientation_and_flip() -> None:
    state = DimerState.all_horizontal(4)
    assert state.face_orientation(0, 0) == "horizontal"
    assert state.face_orientation(0, 1) is None
    changed = state.set_face_orientation(0, 0, "vertical")
    assert changed
    assert state.face_orientation(0, 0) == "vertical"
    assert state.validate()
    assert not state.set_face_orientation(0, 0, "vertical")
    with pytest.raises(ValueError, match="not flippable"):
        state.flip_face(0, 1)


def test_random_updates_preserve_perfect_matching() -> None:
    L = 8
    state = DimerState.all_horizontal(L)
    environment = sample_environment(L, np.random.default_rng(101), "gamma", gamma_shape=0.7)
    rng = np.random.default_rng(102)
    for _ in range(2_000):
        heatbath_update(state, environment, rng)
        assert state.validate()
