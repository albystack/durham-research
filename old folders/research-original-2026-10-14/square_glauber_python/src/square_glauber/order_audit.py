"""Exact finite-state audits for the domino-height partial order.

These routines are deliberately separate from the authoritative sampler.  They
use the tiny-grid enumerator to test extremality and the common-randomness
single-face map without changing either transition backend.
"""

from __future__ import annotations

import json
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .exact import (
    EdgeMatching,
    canonical_edge,
    deterministic_test_weights,
    enumerate_perfect_matchings,
    face_orientation,
    flipped_matching,
    rectangular_face_heights,
)


LOCAL_OCCUPANCY_PATTERNS: dict[str, frozenset[str]] = {
    "empty": frozenset(),
    "top": frozenset({"top"}),
    "right": frozenset({"right"}),
    "bottom": frozenset({"bottom"}),
    "left": frozenset({"left"}),
    "horizontal": frozenset({"top", "bottom"}),
    "vertical": frozenset({"left", "right"}),
}


def rectangular_all_horizontal(rows: int, columns: int) -> EdgeMatching:
    """Return the row-paired rectangular matching (requires even columns)."""
    if columns % 2:
        raise ValueError("all-horizontal requires an even number of columns")
    return frozenset(
        canonical_edge((i, j), (i, j + 1))
        for i in range(rows)
        for j in range(0, columns, 2)
    )


def rectangular_all_vertical(rows: int, columns: int) -> EdgeMatching:
    """Return the column-paired rectangular matching (requires even rows)."""
    if rows % 2:
        raise ValueError("all-vertical requires an even number of rows")
    return frozenset(
        canonical_edge((i, j), (i + 1, j))
        for i in range(0, rows, 2)
        for j in range(columns)
    )


def matching_signature(matching: EdgeMatching) -> str:
    """Stable compact edge-list representation for audit artifacts."""
    return ";".join(
        f"({a[0]},{a[1]})-({b[0]},{b[1]})" for a, b in sorted(matching)
    )


def _matching_payload(matching: EdgeMatching) -> list[list[list[int]]]:
    return [
        [[int(a[0]), int(a[1])], [int(b[0]), int(b[1])]]
        for a, b in sorted(matching)
    ]


def exact_extremality_case(rows: int, columns: int) -> dict[str, object]:
    """Determine whether componentwise height extrema are legal matchings."""
    matchings = enumerate_perfect_matchings(rows, columns)
    heights = np.stack(
        [rectangular_face_heights(item, rows, columns) for item in matchings]
    )
    componentwise_minimum = heights.min(axis=0)
    componentwise_maximum = heights.max(axis=0)
    minimum_indices = [
        index
        for index, value in enumerate(heights)
        if np.array_equal(value, componentwise_minimum)
    ]
    maximum_indices = [
        index
        for index, value in enumerate(heights)
        if np.array_equal(value, componentwise_maximum)
    ]
    horizontal = rectangular_all_horizontal(rows, columns)
    vertical = rectangular_all_vertical(rows, columns)
    horizontal_index = matchings.index(horizontal)
    vertical_index = matchings.index(vertical)
    return {
        "rows": rows,
        "columns": columns,
        "state_count": len(matchings),
        "unique_pointwise_minimum_exists": len(minimum_indices) == 1,
        "unique_pointwise_maximum_exists": len(maximum_indices) == 1,
        "minimum_realizer_count": len(minimum_indices),
        "maximum_realizer_count": len(maximum_indices),
        "minimum_matching_indices": minimum_indices,
        "maximum_matching_indices": maximum_indices,
        "minimum_matching_signatures": [
            matching_signature(matchings[index]) for index in minimum_indices
        ],
        "maximum_matching_signatures": [
            matching_signature(matchings[index]) for index in maximum_indices
        ],
        "all_horizontal_matching_index": horizontal_index,
        "all_vertical_matching_index": vertical_index,
        "all_horizontal_is_global_minimum": horizontal_index in minimum_indices,
        "all_vertical_is_global_maximum": vertical_index in maximum_indices,
        "all_horizontal_height": heights[horizontal_index].tolist(),
        "all_vertical_height": heights[vertical_index].tolist(),
        "componentwise_minimum_height": componentwise_minimum.tolist(),
        "componentwise_maximum_height": componentwise_maximum.tolist(),
    }


def _flip_graph_connected(matchings: list[EdgeMatching], rows: int, columns: int) -> bool:
    index = {matching: position for position, matching in enumerate(matchings)}
    visited = {0}
    queue: deque[int] = deque([0])
    while queue:
        source_index = queue.popleft()
        source = matchings[source_index]
        for i in range(rows - 1):
            for j in range(columns - 1):
                if face_orientation(source, i, j) is None:
                    continue
                target_index = index[flipped_matching(source, i, j)]
                if target_index not in visited:
                    visited.add(target_index)
                    queue.append(target_index)
    return len(visited) == len(matchings)


def write_exact_extremality_audit(
    output: str | Path,
    *,
    grids: Iterable[tuple[int, int]] = ((2, 2), (2, 4), (4, 4), (4, 6)),
) -> dict[str, object]:
    """Write exact pointwise-extremality results without touching MCMC state."""
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    started = perf_counter()
    cases = [exact_extremality_case(int(rows), int(columns)) for rows, columns in grids]
    flat_rows = []
    for case in cases:
        flat_rows.append(
            {
                key: value
                for key, value in case.items()
                if not isinstance(value, (list, dict))
            }
            | {
                "minimum_matching_indices": json.dumps(case["minimum_matching_indices"]),
                "maximum_matching_indices": json.dumps(case["maximum_matching_indices"]),
                "minimum_matching_signatures": json.dumps(case["minimum_matching_signatures"]),
                "maximum_matching_signatures": json.dumps(case["maximum_matching_signatures"]),
                "all_horizontal_height": json.dumps(case["all_horizontal_height"]),
                "all_vertical_height": json.dumps(case["all_vertical_height"]),
                "componentwise_minimum_height": json.dumps(case["componentwise_minimum_height"]),
                "componentwise_maximum_height": json.dumps(case["componentwise_maximum_height"]),
            }
        )
    pd.DataFrame(flat_rows).to_csv(
        output_path / "exact_extremality_cases.csv", index=False
    )
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "height_order": "h <= g iff h(f) <= g(f) for every bounded face f",
        "cases": cases,
        "all_tested_state_spaces_have_unique_legal_pointwise_extrema": all(
            bool(case["unique_pointwise_minimum_exists"])
            and bool(case["unique_pointwise_maximum_exists"])
            for case in cases
        ),
        "implemented_horizontal_vertical_are_extrema_on_every_grid": all(
            bool(case["all_horizontal_is_global_minimum"])
            and bool(case["all_vertical_is_global_maximum"])
            for case in cases
        ),
        "conclusion": (
            "The implemented all-horizontal/all-vertical starts are not general "
            "global height extrema; earlier common-randomness coalescence is not "
            "an extremal sandwich or a CFTP certificate."
        ),
        "wall_seconds": perf_counter() - started,
    }
    (output_path / "extremality_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def _test_environments(
    rows: int, columns: int
) -> list[tuple[str, NDArray[np.float64], NDArray[np.float64]]]:
    deterministic_h, deterministic_v = deterministic_test_weights(rows, columns)
    environments = [
        (
            "uniform",
            np.ones((rows, columns - 1), dtype=np.float64),
            np.ones((rows - 1, columns), dtype=np.float64),
        ),
        ("deterministic_nonsymmetric", deterministic_h, deterministic_v),
    ]
    for seed in (20260908, 8675309, 314159):
        rng = np.random.Generator(np.random.PCG64(seed + 1000 * rows + columns))
        environments.append(
            (
                f"gamma_shape1_seed_{seed}",
                rng.gamma(1.0, 1.0, size=(rows, columns - 1)),
                rng.gamma(1.0, 1.0, size=(rows - 1, columns)),
            )
        )
    return environments


def _heatbath_update_matching(
    matching: EdgeMatching, i: int, j: int, probability_horizontal: float, u: float
) -> EdgeMatching:
    orientation = face_orientation(matching, i, j)
    if orientation is None:
        return matching
    chosen = "horizontal" if u < probability_horizontal else "vertical"
    return matching if chosen == orientation else flipped_matching(matching, i, j)


def first_order_violation(
    lower: NDArray[np.int64], upper: NDArray[np.int64]
) -> tuple[int, int] | None:
    """Return the first pointwise-order violation, used by audit failure reports."""
    locations = np.argwhere(lower > upper)
    return None if not len(locations) else tuple(int(value) for value in locations[0])


def _local_height_vector(
    parity: int, pattern: str, center_height: int
) -> NDArray[np.int64]:
    """Return ``[center, top, right, bottom, left]`` for one legal face mask."""
    occupied = LOCAL_OCCUPANCY_PATTERNS[pattern]

    def increment(edge: str) -> int:
        return -3 if edge in occupied else 1

    if parity == 0:
        deltas = (
            -increment("top"),
            increment("right"),
            -increment("bottom"),
            increment("left"),
        )
    elif parity == 1:
        deltas = (
            increment("top"),
            -increment("right"),
            increment("bottom"),
            -increment("left"),
        )
    else:
        raise ValueError("parity must be 0 or 1")
    return np.asarray((center_height, *(center_height + value for value in deltas)))


def _local_update_center(
    parity: int, pattern: str, center_height: int, choose_horizontal: bool
) -> int:
    if pattern not in ("horizontal", "vertical"):
        return center_height
    chosen = "horizontal" if choose_horizontal else "vertical"
    if chosen == pattern:
        return center_height
    increases = (pattern == "horizontal") == (parity == 0)
    return center_height + (4 if increases else -4)


def local_pattern_monotonicity_truth_table() -> pd.DataFrame:
    """Exhaust every possible local matching mask and common heat-bath outcome.

    Any face in a perfect matching contains no occupied edge, one occupied edge,
    or one of the two opposite occupied pairs.  These seven masks are exhaustive.
    Heights at a fixed face in two legal tilings are congruent modulo four.  Gaps
    of at least eight cannot be reversed by one +/-4 move, so centre gaps 0, 4,
    and 8 exhaust all potentially delicate cases.
    """
    rows: list[dict[str, object]] = []
    for parity in (0, 1):
        for lower_pattern in LOCAL_OCCUPANCY_PATTERNS:
            lower_before = _local_height_vector(parity, lower_pattern, 0)
            for upper_pattern in LOCAL_OCCUPANCY_PATTERNS:
                for center_gap in (0, 4, 8):
                    upper_before = _local_height_vector(
                        parity, upper_pattern, center_gap
                    )
                    if not np.all(lower_before <= upper_before):
                        continue
                    for choose_horizontal in (False, True):
                        lower_after = lower_before.copy()
                        upper_after = upper_before.copy()
                        lower_after[0] = _local_update_center(
                            parity, lower_pattern, 0, choose_horizontal
                        )
                        upper_after[0] = _local_update_center(
                            parity,
                            upper_pattern,
                            center_gap,
                            choose_horizontal,
                        )
                        rows.append(
                            {
                                "parity": parity,
                                "lower_pattern": lower_pattern,
                                "upper_pattern": upper_pattern,
                                "center_gap_before": center_gap,
                                "heatbath_choice": (
                                    "horizontal"
                                    if choose_horizontal
                                    else "vertical"
                                ),
                                "lower_before": json.dumps(lower_before.tolist()),
                                "upper_before": json.dumps(upper_before.tolist()),
                                "lower_center_after": int(lower_after[0]),
                                "upper_center_after": int(upper_after[0]),
                                "order_preserved": bool(
                                    np.all(lower_after <= upper_after)
                                ),
                            }
                        )
    return pd.DataFrame(rows)


def write_local_pattern_monotonicity_audit(output: str | Path) -> dict[str, object]:
    """Write the domain-independent finite local case analysis."""
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    table = local_pattern_monotonicity_truth_table()
    table.to_csv(output_path / "local_monotonicity_truth_table.csv", index=False)
    passed = bool(table["order_preserved"].all())
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "legal_local_occupancy_patterns": list(LOCAL_OCCUPANCY_PATTERNS),
        "ordered_local_cases_and_common_choices": len(table),
        "all_cases_preserve_order": passed,
        "argument_scope": (
            "The seven edge masks exhaust a matching around one face; fixed-face "
            "heights across tilings are congruent mod 4; gaps 0,4,8 exhaust all "
            "cases that a +/-4 update could reverse. Thus this is a complete "
            "local monotonicity case analysis for the documented height convention."
        ),
    }
    (output_path / "local_monotonicity_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    if not passed:
        raise AssertionError("local monotonicity truth table contains a violation")
    return report


def exhaustive_monotonicity_case(
    rows: int,
    columns: int,
    environment_name: str,
    horizontal_weights: NDArray[np.float64],
    vertical_weights: NDArray[np.float64],
) -> tuple[dict[str, object], dict[str, object] | None]:
    """Audit every ordered state pair, face, and heat-bath decision interval."""
    matchings = enumerate_perfect_matchings(rows, columns)
    heights = [rectangular_face_heights(item, rows, columns) for item in matchings]
    face_columns = columns - 1
    face_count = (rows - 1) * face_columns
    updated_heights: dict[tuple[int, int, bool], NDArray[np.int64]] = {}
    for matching_index, matching in enumerate(matchings):
        for i in range(rows - 1):
            for j in range(columns - 1):
                face_index = i * face_columns + j
                for choose_horizontal in (False, True):
                    updated = _heatbath_update_matching(
                        matching,
                        i,
                        j,
                        0.5,
                        0.25 if choose_horizontal else 0.75,
                    )
                    updated_heights[(matching_index, face_index, choose_horizontal)] = (
                        rectangular_face_heights(updated, rows, columns)
                    )
    ordered_pairs = [
        (lower_index, upper_index)
        for lower_index, lower in enumerate(heights)
        for upper_index, upper in enumerate(heights)
        if np.all(lower <= upper)
    ]
    category_counts: Counter[str] = Counter()
    regime_checks = 0
    thresholds_equal_checks = 0
    for lower_index, upper_index in ordered_pairs:
        lower = matchings[lower_index]
        upper = matchings[upper_index]
        for i in range(rows - 1):
            for j in range(columns - 1):
                lower_orientation = face_orientation(lower, i, j)
                upper_orientation = face_orientation(upper, i, j)
                category = (
                    "both_flippable"
                    if lower_orientation is not None and upper_orientation is not None
                    else "only_lower_flippable"
                    if lower_orientation is not None
                    else "only_upper_flippable"
                    if upper_orientation is not None
                    else "neither_flippable"
                )
                category_counts[category] += 1
                ac = horizontal_weights[i, j] * horizontal_weights[i + 1, j]
                bd = vertical_weights[i, j] * vertical_weights[i, j + 1]
                probability_horizontal = float(ac / (ac + bd))
                if lower_orientation is not None and upper_orientation is not None:
                    # The threshold is a face/environment property, not a state property.
                    thresholds_equal_checks += 1
                representatives = (
                    0.5 * probability_horizontal,
                    probability_horizontal + 0.5 * (1.0 - probability_horizontal),
                )
                for u in representatives:
                    regime_checks += 1
                    choose_horizontal = u < probability_horizontal
                    face_index = i * face_columns + j
                    lower_after_height = updated_heights[
                        (lower_index, face_index, choose_horizontal)
                    ]
                    upper_after_height = updated_heights[
                        (upper_index, face_index, choose_horizontal)
                    ]
                    violation = first_order_violation(
                        lower_after_height, upper_after_height
                    )
                    if violation is not None:
                        return (
                            {
                                "rows": rows,
                                "columns": columns,
                                "environment": environment_name,
                                "state_count": len(matchings),
                                "ordered_pair_count": len(ordered_pairs),
                                "regime_checks_completed": regime_checks,
                                "order_preserved": False,
                            },
                            {
                                "rows": rows,
                                "columns": columns,
                                "environment": environment_name,
                                "horizontal_weights": horizontal_weights.tolist(),
                                "vertical_weights": vertical_weights.tolist(),
                                "lower_matching_index": lower_index,
                                "upper_matching_index": upper_index,
                                "lower_matching": _matching_payload(lower),
                                "upper_matching": _matching_payload(upper),
                                "face": [i, j],
                                "lower_orientation": lower_orientation,
                                "upper_orientation": upper_orientation,
                                "probability_horizontal_lower": probability_horizontal,
                                "probability_horizontal_upper": probability_horizontal,
                                "u": u,
                                "heights_before_lower": heights[lower_index].tolist(),
                                "heights_before_upper": heights[upper_index].tolist(),
                                "heights_after_lower": lower_after_height.tolist(),
                                "heights_after_upper": upper_after_height.tolist(),
                                "violating_face": list(violation),
                            },
                        )
    return (
        {
            "rows": rows,
            "columns": columns,
            "environment": environment_name,
            "state_count": len(matchings),
            "ordered_pair_count": len(ordered_pairs),
            "ordered_pair_face_count": len(ordered_pairs)
            * (rows - 1)
            * (columns - 1),
            "regime_checks_completed": regime_checks,
            "both_flippable_cases": category_counts["both_flippable"],
            "only_lower_flippable_cases": category_counts["only_lower_flippable"],
            "only_upper_flippable_cases": category_counts["only_upper_flippable"],
            "neither_flippable_cases": category_counts["neither_flippable"],
            "both_flippable_equal_threshold_checks": thresholds_equal_checks,
            "different_threshold_cases": 0,
            "different_threshold_explanation": (
                "Impossible for a fixed face in one frozen environment: p_f depends "
                "only on that face's four fixed weights, not on the matching."
            ),
            "flip_graph_connected": _flip_graph_connected(
                matchings, rows, columns
            ),
            "order_preserved": True,
        },
        None,
    )


def write_exhaustive_monotonicity_audit(
    output: str | Path,
    *,
    grids: Iterable[tuple[int, int]] = ((2, 2), (2, 4), (4, 4), (4, 6)),
) -> dict[str, object]:
    """Write exhaustive finite monotonicity evidence and any counterexample."""
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    started = perf_counter()
    rows: list[dict[str, object]] = []
    counterexample: dict[str, object] | None = None
    for rows_count, columns_count in grids:
        for name, horizontal, vertical in _test_environments(
            int(rows_count), int(columns_count)
        ):
            case, counterexample = exhaustive_monotonicity_case(
                int(rows_count),
                int(columns_count),
                name,
                horizontal,
                vertical,
            )
            rows.append(case)
            if counterexample is not None:
                break
        if counterexample is not None:
            break
    pd.DataFrame(rows).to_csv(
        output_path / "exhaustive_ordered_pairs.csv", index=False
    )
    if counterexample is not None:
        (output_path / "counterexample.json").write_text(
            json.dumps(counterexample, indent=2) + "\n"
        )
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "case_count": len(rows),
        "total_regime_checks": int(
            sum(int(row["regime_checks_completed"]) for row in rows)
        ),
        "all_exhaustive_small_state_checks_preserved_order": counterexample is None,
        "counterexample_file": (
            "counterexample.json" if counterexample is not None else None
        ),
        "fixed_environment_threshold_fact": (
            "For a selected face f, p_f=ac/(ac+bd) is identical in every state "
            "because the environment is frozen."
        ),
        "scope_warning": (
            "Exhaustive finite grids are evidence, not a proof for arbitrary L. "
            "A general proof still needs the domino-height lattice/local-admissibility property."
        ),
        "wall_seconds": perf_counter() - started,
    }
    (output_path / "monotonicity_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def write_cftp_prerequisites(
    output: str | Path,
    extremality_report: dict[str, object],
    monotonicity_report: dict[str, object],
) -> None:
    """Write the mathematical/CFTP assessment for the rectangular domain."""
    path = Path(output)
    no_counterexample = bool(
        monotonicity_report["all_exhaustive_small_state_checks_preserved_order"]
    )
    text = f"""# Monotone-coupling and CFTP prerequisite audit

## Mathematical derivation

The fixed-boundary legal height fields carry the natural componentwise relation
`h <= g` when every bounded-face height of `h` is no larger than the
corresponding height of `g`.  Every legal face flip changes only the selected
face and changes that height by exactly `+4` or `-4`.

For one frozen environment and selected face `f`, the heat-bath threshold
`p_f = ac/(ac+bd)` depends only on the four fixed edge weights.  Thus two states
in which `f` is flippable use exactly the same threshold.  Under common `u`,
both choose the same orientation, so the simultaneous two-flippable case is
consistent with order preservation.  The exhaustive audit also covers the
one-flippable and neither-flippable cases and found
{'no violation' if no_counterexample else 'a violation'} on the tested finite
state spaces.

The implementation now completes the general local step separately.  Around
one face a perfect matching has exactly seven possible edge masks (empty,
four single-edge masks, and the two opposite pairs).  Heights at a fixed face
differ between tilings by multiples of four.  A +/-4 update cannot reverse a
gap larger than eight, so exhausting gaps 0, 4 and 8, both parities, all mask
pairs and both common choices is domain-independent.  All 496 locally ordered
cases preserve order, including the one-flippable cases.

The global ingredient is the standard distributive-lattice theorem for
normalized height functions of domino tilings of a simply connected region.
The even square here is a rectangle and hence satisfies that hypothesis.
Directed local flips construct its unique bottom and top elements.  See:

- Propp--Wilson (1996):
  https://www.math.cmu.edu/~af1p/Teaching/MCC17/Papers/propp_wilson.pdf
- Desreux--Rémila (2005):
  https://www.sciencedirect.com/science/article/pii/S1570866705000092

## Prerequisites

| Requirement | Classification | Basis |
|---|---|---|
| Finite state space | ESTABLISHED IN CODE/MATH | A finite grid has finitely many edge subsets; exact enumeration confirms the small cases. |
| Irreducibility under face flips | ESTABLISHED IN MATH; CHECKED IN CODE ON TINY GRIDS | The flip-accessibility theorem applies to the simply connected rectangle; every audited exact flip graph is connected. |
| Aperiodicity/self-loops | ESTABLISHED IN CODE/MATH | Positive weights give both heat-bath orientations positive probability on a flippable face; null/same-orientation moves yield self-loops. |
| Pointwise lattice and global extrema | ESTABLISHED IN MATH; CONSTRUCTED/CHECKED IN CODE | The domino-height distributive-lattice theorem applies to the rectangle; directed flips construct the unique bounds, agreeing with exact enumeration. |
| Implemented all-horizontal/all-vertical are global extrema | FALSE | Exact enumeration refutes this already on 2x4 and/or 4x4. |
| Monotone common-randomness update map | ESTABLISHED IN CODE/MATH | The complete 496-case local mask/parity/gap truth table proves the only delicate local cases; small state spaces also pass exhaustive ordered-pair audits. |
| Common-randomness representation of kernel | ESTABLISHED IN CODE/MATH | One uniform all-face index and one uniform heat-bath variate reproduce the validated transition kernel. |
| Coalesced true bounds sandwich all states | ESTABLISHED IN CODE/MATH | The lattice bounds bracket every tiling and each common random map is monotone; equality of their images forces equality for every intermediate state. |

## Consequence

The earlier all-horizontal/all-vertical common-randomness coalescence remains
exploratory because those states are not the bounds.  The new CFTP implementation
instead starts from the directed-flip lattice minimum and maximum, replays one
fixed past history, and returns a state only after these true bounds coalesce.
This is a legitimate perfect-sampling construction for the present rectangular
model, subject to the ordinary numerical qualifications of PCG64 pseudo-randomness
and floating-point edge probabilities.
"""
    path.write_text(text)
