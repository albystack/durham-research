import numpy as np

from square_glauber.exact import deterministic_test_weights
from square_glauber.order_audit import (
    exact_extremality_case,
    exhaustive_monotonicity_case,
    first_order_violation,
)


def test_exact_tiny_grid_extrema_identify_nonextremal_standard_starts() -> None:
    two_by_two = exact_extremality_case(2, 2)
    assert two_by_two["state_count"] == 2
    assert two_by_two["all_horizontal_is_global_minimum"]
    assert two_by_two["all_vertical_is_global_maximum"]

    two_by_four = exact_extremality_case(2, 4)
    assert two_by_four["state_count"] == 5
    assert two_by_four["unique_pointwise_minimum_exists"]
    assert two_by_four["unique_pointwise_maximum_exists"]
    assert two_by_four["all_horizontal_is_global_minimum"]
    assert not two_by_four["all_vertical_is_global_maximum"]
    assert two_by_four["componentwise_maximum_height"] == [[1, 2, 1]]

    four_by_four = exact_extremality_case(4, 4)
    assert four_by_four["state_count"] == 36
    assert four_by_four["unique_pointwise_minimum_exists"]
    assert four_by_four["unique_pointwise_maximum_exists"]
    assert not four_by_four["all_horizontal_is_global_minimum"]
    assert not four_by_four["all_vertical_is_global_maximum"]
    assert four_by_four["componentwise_minimum_height"] == [
        [-3, -2, -3],
        [-4, -5, -4],
        [-3, -2, -3],
    ]
    assert four_by_four["componentwise_maximum_height"] == [
        [1, 2, 1],
        [0, 3, 0],
        [1, 2, 1],
    ]


def test_exhaustive_common_update_monotonicity_on_weighted_4x4() -> None:
    horizontal, vertical = deterministic_test_weights(4, 4)
    report, counterexample = exhaustive_monotonicity_case(
        4,
        4,
        "unit_test_deterministic",
        horizontal,
        vertical,
    )
    assert counterexample is None
    assert report["state_count"] == 36
    assert report["order_preserved"]
    assert report["regime_checks_completed"] > 0
    assert report["only_lower_flippable_cases"] > 0
    assert report["only_upper_flippable_cases"] > 0
    assert report["both_flippable_equal_threshold_checks"] > 0
    assert report["different_threshold_cases"] == 0
    assert report["flip_graph_connected"]


def test_artificial_order_violation_is_reported() -> None:
    lower = np.asarray([[0, 5], [1, 2]], dtype=np.int64)
    upper = np.asarray([[0, 4], [1, 3]], dtype=np.int64)
    assert first_order_violation(lower, upper) == (0, 1)
    assert first_order_violation(upper, upper) is None
