"""Transparent tiny-grid enumeration and exact heat-bath kernel checks."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np
from numpy.typing import NDArray

Vertex = tuple[int, int]
Edge = tuple[Vertex, Vertex]
EdgeMatching = frozenset[Edge]


def canonical_edge(first: Vertex, second: Vertex) -> Edge:
    return (first, second) if first < second else (second, first)


def _neighbours(vertex: Vertex, rows: int, columns: int) -> Iterable[Vertex]:
    i, j = vertex
    if i > 0:
        yield i - 1, j
    if i + 1 < rows:
        yield i + 1, j
    if j > 0:
        yield i, j - 1
    if j + 1 < columns:
        yield i, j + 1


def enumerate_perfect_matchings(rows: int, columns: int | None = None) -> list[EdgeMatching]:
    """Recursively enumerate tiny rectangular-grid perfect matchings."""
    columns = rows if columns is None else columns
    if rows <= 0 or columns <= 0:
        raise ValueError("grid dimensions must be positive")
    if rows * columns % 2:
        return []
    unmatched = frozenset((i, j) for i in range(rows) for j in range(columns))
    result: list[EdgeMatching] = []

    def visit(remaining: frozenset[Vertex], selected: tuple[Edge, ...]) -> None:
        if not remaining:
            result.append(frozenset(selected))
            return
        first = min(remaining)
        for neighbour in _neighbours(first, rows, columns):
            if neighbour in remaining:
                visit(
                    remaining.difference((first, neighbour)),
                    selected + (canonical_edge(first, neighbour),),
                )

    visit(unmatched, ())
    return result


def matching_weight(
    matching: EdgeMatching,
    horizontal_weights: NDArray[np.float64],
    vertical_weights: NDArray[np.float64],
) -> float:
    """Product of occupied edge weights for one enumerated matching."""
    log_weight = 0.0
    for first, second in matching:
        if first[0] == second[0]:
            i, j = first
            log_weight += float(np.log(horizontal_weights[i, j]))
        else:
            i, j = first
            log_weight += float(np.log(vertical_weights[i, j]))
    return float(np.exp(log_weight))


def gibbs_probabilities(
    matchings: list[EdgeMatching],
    horizontal_weights: NDArray[np.float64],
    vertical_weights: NDArray[np.float64],
) -> NDArray[np.float64]:
    logs = []
    for matching in matchings:
        value = 0.0
        for first, second in matching:
            if first[0] == second[0]:
                value += float(np.log(horizontal_weights[first[0], first[1]]))
            else:
                value += float(np.log(vertical_weights[first[0], first[1]]))
        logs.append(value)
    log_array = np.asarray(logs, dtype=np.float64)
    shifted = np.exp(log_array - np.max(log_array))
    return shifted / shifted.sum()


def matching_occupancy(
    matching: EdgeMatching, rows: int, columns: int
) -> tuple[NDArray[np.bool_], NDArray[np.bool_]]:
    """Convert an enumerated matching to rectangular occupancy arrays."""
    horizontal = np.zeros((rows, columns - 1), dtype=np.bool_)
    vertical = np.zeros((rows - 1, columns), dtype=np.bool_)
    for first, second in matching:
        if first[0] == second[0]:
            horizontal[first[0], first[1]] = True
        else:
            vertical[first[0], first[1]] = True
    degree = np.zeros((rows, columns), dtype=np.int8)
    degree[:, :-1] += horizontal
    degree[:, 1:] += horizontal
    degree[:-1, :] += vertical
    degree[1:, :] += vertical
    if not np.all(degree == 1):
        raise ValueError("enumerated edge set is not a rectangular-grid perfect matching")
    return horizontal, vertical


def rectangular_face_heights(
    matching: EdgeMatching, rows: int, columns: int
) -> NDArray[np.int64]:
    """Integrate the fixed-cut height convention on a tiny rectangle.

    The exterior reference crosses the top-left horizontal edge into face
    ``(0,0)``. Interior dual-edge increments use the same white-even,
    white-to-black, left-to-right ``+1/-3`` convention as ``height.face_heights``.
    """
    if rows < 2 or columns < 2:
        raise ValueError("height requires at least one bounded face")
    horizontal, vertical = matching_occupancy(matching, rows, columns)
    face_rows = rows - 1
    face_columns = columns - 1
    face_count = face_rows * face_columns
    outside = face_count
    adjacency: list[list[tuple[int, int]]] = [[] for _ in range(face_count + 1)]

    def face_id(i: int, j: int) -> int:
        return i * face_columns + j

    def connect(first: int, second: int, difference: int) -> None:
        """Record ``h(second)-h(first)=difference`` in both directions."""
        adjacency[first].append((second, difference))
        adjacency[second].append((first, -difference))

    anchor_increment = -3 if horizontal[0, 0] else 1
    connect(outside, face_id(0, 0), anchor_increment)

    # Vertically adjacent bounded faces cross a horizontal primal edge.
    for face_i in range(face_rows - 1):
        edge_i = face_i + 1
        for j in range(face_columns):
            increment = -3 if horizontal[edge_i, j] else 1
            above_to_below = increment if (edge_i + j) % 2 == 0 else -increment
            connect(face_id(face_i, j), face_id(face_i + 1, j), above_to_below)

    # Horizontally adjacent bounded faces cross a vertical primal edge.
    for i in range(face_rows):
        for face_j in range(face_columns - 1):
            edge_j = face_j + 1
            increment = -3 if vertical[i, edge_j] else 1
            west_to_east = -increment if (i + edge_j) % 2 == 0 else increment
            connect(face_id(i, face_j), face_id(i, face_j + 1), west_to_east)

    assigned: list[int | None] = [None] * (face_count + 1)
    assigned[outside] = 0
    queue: deque[int] = deque([outside])
    while queue:
        source = queue.popleft()
        source_height = assigned[source]
        assert source_height is not None
        for target, increment in adjacency[source]:
            implied = source_height + increment
            if assigned[target] is None:
                assigned[target] = implied
                queue.append(target)
            elif assigned[target] != implied:
                raise ValueError("rectangular height is path-inconsistent")
    if any(value is None for value in assigned):
        raise ValueError("rectangular bounded-face dual graph is disconnected")
    return np.asarray(assigned[:-1], dtype=np.int64).reshape((face_rows, face_columns))


@dataclass(frozen=True)
class ExactCenterHeightDistribution:
    rows: int
    columns: int
    matching_count: int
    center_index: tuple[int, int]
    mean: float
    variance: float
    pmf: dict[int, float]


def exact_center_height_distribution(
    rows: int,
    columns: int,
    horizontal_weights: NDArray[np.float64],
    vertical_weights: NDArray[np.float64],
) -> ExactCenterHeightDistribution:
    """Enumerate the exact Gibbs PMF and moments of the central bounded face."""
    matchings = enumerate_perfect_matchings(rows, columns)
    probabilities = gibbs_probabilities(matchings, horizontal_weights, vertical_weights)
    center = ((rows - 2) // 2, (columns - 2) // 2)
    masses: dict[int, float] = {}
    for matching, probability in zip(matchings, probabilities, strict=True):
        height = int(rectangular_face_heights(matching, rows, columns)[center])
        masses[height] = masses.get(height, 0.0) + float(probability)
    pmf = dict(sorted(masses.items()))
    mean = sum(height * probability for height, probability in pmf.items())
    variance = sum((height - mean) ** 2 * probability for height, probability in pmf.items())
    return ExactCenterHeightDistribution(
        rows,
        columns,
        len(matchings),
        center,
        float(mean),
        float(variance),
        pmf,
    )


def deterministic_test_weights(
    rows: int, columns: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Positive, non-symmetric, deterministic weights used by exact audits."""
    horizontal = np.fromfunction(
        lambda i, j: 0.55 + ((7 * i + 11 * j + 3) % 13) / 7.0,
        (rows, columns - 1),
        dtype=int,
    )
    vertical = np.fromfunction(
        lambda i, j: 0.65 + ((5 * i + 9 * j + 4) % 11) / 6.0,
        (rows - 1, columns),
        dtype=int,
    )
    return np.asarray(horizontal, dtype=np.float64), np.asarray(vertical, dtype=np.float64)


def face_orientation(
    matching: EdgeMatching, i: int, j: int
) -> Literal["horizontal", "vertical"] | None:
    top = canonical_edge((i, j), (i, j + 1))
    right = canonical_edge((i, j + 1), (i + 1, j + 1))
    bottom = canonical_edge((i + 1, j), (i + 1, j + 1))
    left = canonical_edge((i, j), (i + 1, j))
    if top in matching and bottom in matching and right not in matching and left not in matching:
        return "horizontal"
    if left in matching and right in matching and top not in matching and bottom not in matching:
        return "vertical"
    return None


def flipped_matching(matching: EdgeMatching, i: int, j: int) -> EdgeMatching:
    orientation = face_orientation(matching, i, j)
    if orientation is None:
        raise ValueError("face is not flippable")
    top = canonical_edge((i, j), (i, j + 1))
    right = canonical_edge((i, j + 1), (i + 1, j + 1))
    bottom = canonical_edge((i + 1, j), (i + 1, j + 1))
    left = canonical_edge((i, j), (i + 1, j))
    removed = (top, bottom) if orientation == "horizontal" else (left, right)
    added = (left, right) if orientation == "horizontal" else (top, bottom)
    return frozenset(matching.difference(removed).union(added))


def transition_matrix(
    rows: int,
    columns: int,
    horizontal_weights: NDArray[np.float64],
    vertical_weights: NDArray[np.float64],
    *,
    matchings: list[EdgeMatching] | None = None,
) -> tuple[list[EdgeMatching], NDArray[np.float64]]:
    """Exact literal random-all-faces one-step heat-bath matrix."""
    if horizontal_weights.shape != (rows, columns - 1):
        raise ValueError("horizontal weight shape does not match the grid")
    if vertical_weights.shape != (rows - 1, columns):
        raise ValueError("vertical weight shape does not match the grid")
    if not np.all(horizontal_weights > 0) or not np.all(vertical_weights > 0):
        raise ValueError("all exact-kernel weights must be positive")
    states = enumerate_perfect_matchings(rows, columns) if matchings is None else matchings
    index = {matching: position for position, matching in enumerate(states)}
    face_count = (rows - 1) * (columns - 1)
    matrix = np.zeros((len(states), len(states)), dtype=np.float64)
    for source_index, matching in enumerate(states):
        for i in range(rows - 1):
            for j in range(columns - 1):
                orientation = face_orientation(matching, i, j)
                if orientation is None:
                    matrix[source_index, source_index] += 1.0 / face_count
                    continue
                ac = horizontal_weights[i, j] * horizontal_weights[i + 1, j]
                bd = vertical_weights[i, j] * vertical_weights[i, j + 1]
                probability_horizontal = float(ac / (ac + bd))
                target = flipped_matching(matching, i, j)
                if target not in index:
                    raise AssertionError("a local face flip produced an unenumerated matching")
                target_index = index[target]
                probability_same = (
                    probability_horizontal
                    if orientation == "horizontal"
                    else 1.0 - probability_horizontal
                )
                matrix[source_index, source_index] += probability_same / face_count
                matrix[source_index, target_index] += (1.0 - probability_same) / face_count
    return states, matrix


@dataclass(frozen=True)
class ExactValidation:
    rows: int
    columns: int
    matching_count: int
    row_sum_residual: float
    stationarity_residual: float
    detailed_balance_residual: float
    illegal_transition_count: int


def validate_exact_kernel(
    rows: int,
    columns: int,
    horizontal_weights: NDArray[np.float64],
    vertical_weights: NDArray[np.float64],
) -> ExactValidation:
    """Check stochasticity, legal support, stationarity, and detailed balance."""
    states, matrix = transition_matrix(rows, columns, horizontal_weights, vertical_weights)
    pi = gibbs_probabilities(states, horizontal_weights, vertical_weights)
    row_residual = float(np.max(np.abs(matrix.sum(axis=1) - 1.0)))
    stationarity = float(np.max(np.abs(pi @ matrix - pi)))
    flux = pi[:, None] * matrix
    detailed_balance = float(np.max(np.abs(flux - flux.T)))

    illegal = 0
    for i, source in enumerate(states):
        for j, target in enumerate(states):
            if i == j or matrix[i, j] == 0:
                continue
            difference = source.symmetric_difference(target)
            if len(difference) != 4:
                illegal += 1
                continue
            legal_flip = any(
                face_orientation(source, face_i, face_j) is not None
                and flipped_matching(source, face_i, face_j) == target
                for face_i in range(rows - 1)
                for face_j in range(columns - 1)
            )
            illegal += int(not legal_flip)
    return ExactValidation(
        rows,
        columns,
        len(states),
        row_residual,
        stationarity,
        detailed_balance,
        illegal,
    )


def standard_exact_validations(seed: int = 20260907) -> tuple[ExactValidation, ...]:
    """Run deterministic arbitrary-weight checks for 2x2, 2x4, and 4x4."""
    seed_sequence = np.random.SeedSequence(seed)
    reports: list[ExactValidation] = []
    for (rows, columns), child in zip(((2, 2), (2, 4), (4, 4)), seed_sequence.spawn(3)):
        rng = np.random.Generator(np.random.PCG64(child))
        horizontal = rng.gamma(0.7, 1.0 / 0.7, size=(rows, columns - 1))
        vertical = rng.gamma(0.7, 1.0 / 0.7, size=(rows - 1, columns))
        reports.append(validate_exact_kernel(rows, columns, horizontal, vertical))
    return tuple(reports)
