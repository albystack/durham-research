"""True pointwise minimum and maximum square-grid dimer height states.

For the fixed height convention, a legal face flip changes only that face by
``+4`` or ``-4``.  Repeated directed flips terminate because the sum of bounded
face heights changes strictly.  For a simply connected domino-tileable region,
the domino-height lattice theorem identifies the resulting no-down-flip and
no-up-flip states as the unique global minimum and maximum.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Literal

import numpy as np

from .height import face_heights
from .matching import DimerState, FaceOrientation

HeightDirection = Literal["minimum", "maximum"]


@dataclass(frozen=True)
class ExtremalConstruction:
    direction: HeightDirection
    state: DimerState
    directed_flips: int
    processed_faces: int
    starting_state: str


def face_flip_height_delta(
    i: int, j: int, orientation: FaceOrientation
) -> int:
    """Return the selected-face height change caused by flipping orientation."""
    if orientation not in ("horizontal", "vertical"):
        raise ValueError("orientation must be horizontal or vertical")
    parity_even = (i + j) % 2 == 0
    increases = (orientation == "horizontal") == parity_even
    return 4 if increases else -4


def _has_directed_flip(
    state: DimerState, i: int, j: int, direction: HeightDirection
) -> bool:
    orientation = state.face_orientation(i, j)
    if orientation is None:
        return False
    delta = face_flip_height_delta(i, j, orientation)
    return delta < 0 if direction == "minimum" else delta > 0


def construct_extremal_state(
    L: int,
    direction: HeightDirection,
    *,
    initial_state: DimerState | None = None,
) -> ExtremalConstruction:
    """Construct the unique directed-flip sink from any supplied legal tiling."""
    if direction not in ("minimum", "maximum"):
        raise ValueError("direction must be 'minimum' or 'maximum'")
    if initial_state is None:
        state = DimerState.all_horizontal(L)
        starting_state = "all_horizontal"
    else:
        if initial_state.L != L:
            raise ValueError("initial_state has the wrong L")
        initial_state.validate()
        state = initial_state.copy()
        starting_state = "supplied"

    side = L - 1
    queue = deque((i, j) for i in range(side) for j in range(side))
    queued = np.ones((side, side), dtype=np.bool_)
    directed_flips = 0
    processed = 0
    while queue:
        i, j = queue.popleft()
        queued[i, j] = False
        processed += 1
        if not _has_directed_flip(state, i, j, direction):
            continue
        state.flip_face(i, j)
        directed_flips += 1
        # Only this face and edge-neighbouring faces can change flippability.
        for ni, nj in ((i, j), (i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
            if 0 <= ni < side and 0 <= nj < side and not queued[ni, nj]:
                queue.append((ni, nj))
                queued[ni, nj] = True

    state.validate()
    for i in range(side):
        for j in range(side):
            if _has_directed_flip(state, i, j, direction):
                raise AssertionError("directed-flip construction stopped before a sink")
    return ExtremalConstruction(
        direction, state, directed_flips, processed, starting_state
    )


def true_extremal_states(L: int) -> tuple[ExtremalConstruction, ExtremalConstruction]:
    """Return the true minimum and maximum, verified to be pointwise ordered."""
    minimum = construct_extremal_state(L, "minimum")
    maximum = construct_extremal_state(L, "maximum")
    minimum_heights = face_heights(minimum.state)
    maximum_heights = face_heights(maximum.state)
    if not np.all(minimum_heights <= maximum_heights):
        raise AssertionError("constructed height extrema are not pointwise ordered")
    return minimum, maximum
