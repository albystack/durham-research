"""Dimer height function on bounded faces with a fixed exterior reference.

Vertices use array coordinates ``(row, column)`` and are coloured white when
``(row + column)`` is even. Orient every primal edge from its white endpoint
to its black endpoint. The dual crossing from the geometrical left side of
that oriented edge to its right side is clockwise around the white endpoint,
so its height increment is ``+1`` for an unoccupied edge and ``-3`` for an
occupied dimer. Reversing the crossing reverses the sign. This is precisely
the supervisor-specified white-clockwise convention; black is the opposite.

The induced finite graph has boundary vertices of degree two or three. Thus
one cannot impose the degree-four ``+1/-3`` rule across *every* boundary edge
while identifying all boundary sectors with one combinatorial exterior face.
We instead fix one configuration-independent dual boundary cut: exterior
height zero crosses the top-left horizontal edge into bounded face ``(0,0)``.
All other faces are integrated using interior primal edges. This preserves the
matching-dependent additive component (the anchor crossing is not reset per
configuration) and is the standard fixed-cut equivalent of an exterior
reference for the bounded-face height field.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .geometry import SpatialSeparation, center_face_index
from .matching import DimerState

IntArray = NDArray[np.int64]
EdgeKind = Literal["horizontal", "vertical"]
Face = tuple[int, int] | None  # None denotes the single exterior face.


@dataclass(frozen=True)
class HeightConstraint:
    """One oriented dual-edge equation ``h[right] - h[left] = increment``."""

    kind: EdgeKind
    edge_index: tuple[int, int]
    left_face: Face
    right_face: Face
    occupied: bool
    increment: int


def _bounded_face(L: int, i: int, j: int) -> Face:
    if 0 <= i < L - 1 and 0 <= j < L - 1:
        return i, j
    return None


def height_constraints(state: DimerState) -> tuple[HeightConstraint, ...]:
    """Return every primal-edge constraint in the documented orientation."""
    L = state.L
    constraints: list[HeightConstraint] = []

    # Horizontal edge (i,j)--(i,j+1). If its left endpoint is white, the
    # white-to-black direction is east: north/above is geometrical left and
    # south/below is right. If the right endpoint is white, reverse both.
    for i in range(L):
        for j in range(L - 1):
            above = _bounded_face(L, i - 1, j)
            below = _bounded_face(L, i, j)
            if (i + j) % 2 == 0:
                left_face, right_face = above, below
            else:
                left_face, right_face = below, above
            occupied = bool(state.horizontal_occupied[i, j])
            constraints.append(
                HeightConstraint(
                    "horizontal",
                    (i, j),
                    left_face,
                    right_face,
                    occupied,
                    -3 if occupied else 1,
                )
            )

    # Vertical edge (i,j)--(i+1,j). With a white top endpoint the direction is
    # south: east/right is geometrical left and west/left is right.
    for i in range(L - 1):
        for j in range(L):
            west = _bounded_face(L, i, j - 1)
            east = _bounded_face(L, i, j)
            if (i + j) % 2 == 0:
                left_face, right_face = east, west
            else:
                left_face, right_face = west, east
            occupied = bool(state.vertical_occupied[i, j])
            constraints.append(
                HeightConstraint(
                    "vertical",
                    (i, j),
                    left_face,
                    right_face,
                    occupied,
                    -3 if occupied else 1,
                )
            )
    return tuple(constraints)


def _integration_constraints(state: DimerState) -> tuple[HeightConstraint, ...]:
    """Interior dual edges plus the single fixed top-left boundary cut."""
    result = []
    for item in height_constraints(state):
        both_bounded = item.left_face is not None and item.right_face is not None
        anchor = item.kind == "horizontal" and item.edge_index == (0, 0)
        if both_bounded or anchor:
            result.append(item)
    return tuple(result)


def face_heights(state: DimerState, *, validate_matching: bool = True) -> IntArray:
    """Integrate from exterior height zero through the fixed boundary cut.

    The dual graph is traversed by breadth-first search. Every already-visited
    face is checked again, so multiple dual paths must imply exactly the same
    integer height or a ``ValueError`` is raised.
    """
    if validate_matching:
        state.validate()
    L = state.L
    outside = (L - 1) ** 2

    def face_id(face: Face) -> int:
        return outside if face is None else face[0] * (L - 1) + face[1]

    adjacency: list[list[tuple[int, int]]] = [[] for _ in range(outside + 1)]
    constraints = _integration_constraints(state)
    for item in constraints:
        left = face_id(item.left_face)
        right = face_id(item.right_face)
        adjacency[left].append((right, item.increment))
        adjacency[right].append((left, -item.increment))

    assigned: list[int | None] = [None] * (outside + 1)
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
                raise ValueError(
                    "height is path-inconsistent: one dual face has implied "
                    f"heights {assigned[target]} and {implied}"
                )
    if any(value is None for value in assigned):
        raise ValueError("dual face graph is unexpectedly disconnected")
    if assigned[outside] != 0:
        raise ValueError("exterior height reference changed")

    heights = np.asarray(assigned[:-1], dtype=np.int64).reshape((L - 1, L - 1))
    validate_height_rules(state, heights, constraints=constraints)
    return heights


def validate_height_rules(
    state: DimerState,
    heights: IntArray,
    *,
    constraints: tuple[HeightConstraint, ...] | None = None,
) -> bool:
    """Check shape and every interior/anchor occupancy equation."""
    expected = (state.L - 1, state.L - 1)
    if heights.shape != expected:
        raise ValueError(f"height array must have shape {expected}, got {heights.shape}")
    if not np.issubdtype(heights.dtype, np.integer):
        raise TypeError("height array must have integer dtype")

    def value(face: Face) -> int:
        return 0 if face is None else int(heights[face])

    for item in constraints if constraints is not None else _integration_constraints(state):
        difference = value(item.right_face) - value(item.left_face)
        if difference != item.increment:
            raise ValueError(
                f"height rule fails at {item.kind} edge {item.edge_index}: "
                f"got {difference}, expected {item.increment}"
            )
        if abs(difference) != (3 if item.occupied else 1):
            raise ValueError("height magnitude does not agree with edge occupancy")
    return True


def center_height(state: DimerState) -> int:
    heights = face_heights(state)
    return int(heights[center_face_index(state.L)])


def spatial_increments(
    heights: IntArray,
    separations: tuple[SpatialSeparation, ...],
) -> dict[int, int]:
    """Return rightward central increments keyed by actual integer offset."""
    L = heights.shape[0] + 1
    row, column = center_face_index(L)
    reference = int(heights[row, column])
    return {item.r: int(heights[row, column + item.r]) - reference for item in separations}
