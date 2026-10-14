import numpy as np
import pytest

from square_glauber.environment import DimerEnvironment
from square_glauber.exact import (
    deterministic_test_weights,
    enumerate_perfect_matchings,
    exact_center_height_distribution,
    matching_occupancy,
    rectangular_face_heights,
    standard_exact_validations,
    transition_matrix,
)
from square_glauber.height import face_heights
from square_glauber.glauber import heatbath_update
from square_glauber.matching import DimerState


@pytest.mark.parametrize(
    ("rows", "columns", "expected"), [(2, 2, 2), (2, 4, 5), (4, 4, 36)]
)
def test_tiny_grid_matching_counts(rows: int, columns: int, expected: int) -> None:
    assert len(enumerate_perfect_matchings(rows, columns)) == expected


def test_weighted_2x2_analytic_heatbath_probability() -> None:
    horizontal = np.array([[2.0], [3.0]])
    vertical = np.array([[1.0, 1.0]])
    environment = DimerEnvironment(2, horizontal, vertical)
    assert environment.horizontal_probability(0, 0) == pytest.approx(6 / 7, abs=1e-15)
    states, matrix = transition_matrix(2, 2, horizontal, vertical)
    horizontal_index = next(
        index
        for index, matching in enumerate(states)
        if ((0, 0), (0, 1)) in matching
    )
    assert matrix[horizontal_index, horizontal_index] == pytest.approx(6 / 7, abs=1e-15)
    assert matrix[1 - horizontal_index, horizontal_index] == pytest.approx(6 / 7, abs=1e-15)


def test_update_reports_exact_face_probability() -> None:
    environment = DimerEnvironment(
        2,
        np.array([[2.0], [3.0]]),
        np.array([[1.0, 1.0]]),
    )
    result = heatbath_update(DimerState.all_horizontal(2), environment, np.random.default_rng(7))
    assert result.flippable
    assert result.probability_horizontal == pytest.approx(6 / 7, abs=1e-15)


def test_exact_random_scan_kernels() -> None:
    expected = {(2, 2): 2, (2, 4): 5, (4, 4): 36}
    for report in standard_exact_validations(20260907):
        assert report.matching_count == expected[(report.rows, report.columns)]
        assert report.row_sum_residual < 1e-12
        assert report.stationarity_residual < 1e-12
        assert report.detailed_balance_residual < 1e-12
        assert report.illegal_transition_count == 0


@pytest.mark.parametrize(("rows", "columns"), [(2, 2), (2, 4), (4, 4)])
def test_exact_center_height_distribution(rows: int, columns: int) -> None:
    horizontal, vertical = deterministic_test_weights(rows, columns)
    result = exact_center_height_distribution(
        rows, columns, horizontal, vertical
    )
    assert result.matching_count == {(2, 2): 2, (2, 4): 5, (4, 4): 36}[
        (rows, columns)
    ]
    assert sum(result.pmf.values()) == pytest.approx(1.0, abs=1e-14)
    reconstructed_mean = sum(height * probability for height, probability in result.pmf.items())
    reconstructed_variance = sum(
        (height - reconstructed_mean) ** 2 * probability
        for height, probability in result.pmf.items()
    )
    assert result.mean == pytest.approx(reconstructed_mean, abs=1e-14)
    assert result.variance == pytest.approx(reconstructed_variance, abs=1e-14)


def test_rectangular_height_integrator_matches_square_implementation() -> None:
    for L in (2, 4):
        for matching in enumerate_perfect_matchings(L, L):
            horizontal, vertical = matching_occupancy(matching, L, L)
            square = DimerState(L, horizontal, vertical)
            np.testing.assert_array_equal(
                rectangular_face_heights(matching, L, L), face_heights(square)
            )


def test_deterministic_weight_exact_pmfs_are_frozen() -> None:
    expected = {
        (2, 2): {-3: 0.5992683311985683, 1: 0.40073166880143174},
        (2, 4): {-2: 0.8168578915645399, 2: 0.18314210843546005},
        (4, 4): {
            -5: 0.03380410895876429,
            -1: 0.9368436119205726,
            3: 0.029352279120663263,
        },
    }
    for dimensions, expected_pmf in expected.items():
        horizontal, vertical = deterministic_test_weights(*dimensions)
        result = exact_center_height_distribution(
            *dimensions, horizontal, vertical
        )
        assert result.pmf == pytest.approx(expected_pmf, abs=2e-15)
