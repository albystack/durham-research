"""Authoritative perfect-matching state and local face operations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .geometry import validate_L

BoolArray = NDArray[np.bool_]
FaceOrientation = Literal["horizontal", "vertical"]


@dataclass
class DimerState:
    """Occupied undirected edges of an even square-grid perfect matching."""

    L: int
    horizontal_occupied: BoolArray
    vertical_occupied: BoolArray

    def __post_init__(self) -> None:
        validate_L(self.L)
        if not isinstance(self.horizontal_occupied, np.ndarray):
            raise TypeError("horizontal_occupied must be a NumPy ndarray")
        if not isinstance(self.vertical_occupied, np.ndarray):
            raise TypeError("vertical_occupied must be a NumPy ndarray")

    @classmethod
    def all_horizontal(cls, L: int) -> DimerState:
        """Pair columns ``(0,1), (2,3), ...`` independently in every row."""
        validate_L(L)
        horizontal = np.zeros((L, L - 1), dtype=np.bool_)
        horizontal[:, 0::2] = True
        state = cls(L, horizontal, np.zeros((L - 1, L), dtype=np.bool_))
        state.validate()
        return state

    @classmethod
    def all_vertical(cls, L: int) -> DimerState:
        """Pair rows ``(0,1), (2,3), ...`` independently in every column."""
        validate_L(L)
        vertical = np.zeros((L - 1, L), dtype=np.bool_)
        vertical[0::2, :] = True
        state = cls(L, np.zeros((L, L - 1), dtype=np.bool_), vertical)
        state.validate()
        return state

    def copy(self) -> DimerState:
        return DimerState(
            self.L,
            self.horizontal_occupied.copy(),
            self.vertical_occupied.copy(),
        )

    def validate(self) -> bool:
        """Raise on any shape, dtype, dimer-count, or vertex-degree violation."""
        expected_horizontal = (self.L, self.L - 1)
        expected_vertical = (self.L - 1, self.L)
        if self.horizontal_occupied.shape != expected_horizontal:
            raise ValueError(
                f"horizontal_occupied must have shape {expected_horizontal}, "
                f"got {self.horizontal_occupied.shape}"
            )
        if self.vertical_occupied.shape != expected_vertical:
            raise ValueError(
                f"vertical_occupied must have shape {expected_vertical}, "
                f"got {self.vertical_occupied.shape}"
            )
        if self.horizontal_occupied.dtype != np.bool_:
            raise TypeError("horizontal_occupied must have boolean dtype")
        if self.vertical_occupied.dtype != np.bool_:
            raise TypeError("vertical_occupied must have boolean dtype")
        total = int(self.horizontal_occupied.sum() + self.vertical_occupied.sum())
        expected_total = self.L * self.L // 2
        if total != expected_total:
            raise ValueError(f"matching has {total} dimers; expected {expected_total}")
        degree = np.zeros((self.L, self.L), dtype=np.int8)
        degree[:, :-1] += self.horizontal_occupied
        degree[:, 1:] += self.horizontal_occupied
        degree[:-1, :] += self.vertical_occupied
        degree[1:, :] += self.vertical_occupied
        bad = np.argwhere(degree != 1)
        if bad.size:
            i, j = (int(value) for value in bad[0])
            raise ValueError(f"vertex {(i, j)} has occupied degree {int(degree[i, j])}, not 1")
        return True

    def face_orientation(self, i: int, j: int) -> FaceOrientation | None:
        """Return the occupied opposite pair, or ``None`` if not flippable."""
        if not (0 <= i < self.L - 1 and 0 <= j < self.L - 1):
            raise IndexError(f"face {(i, j)} is outside the bounded-face array")
        top = bool(self.horizontal_occupied[i, j])
        right = bool(self.vertical_occupied[i, j + 1])
        bottom = bool(self.horizontal_occupied[i + 1, j])
        left = bool(self.vertical_occupied[i, j])
        if top and bottom and not right and not left:
            return "horizontal"
        if left and right and not top and not bottom:
            return "vertical"
        return None

    def is_face_flippable(self, i: int, j: int) -> bool:
        return self.face_orientation(i, j) is not None

    def set_face_orientation(self, i: int, j: int, orientation: FaceOrientation) -> bool:
        """Set a flippable face to one of its two orientations in O(1).

        Returns whether the state changed. The operation is refused unless the
        current face is flippable, which prevents accidental nonlocal damage.
        """
        if orientation not in ("horizontal", "vertical"):
            raise ValueError("orientation must be 'horizontal' or 'vertical'")
        current = self.face_orientation(i, j)
        if current is None:
            raise ValueError(f"face {(i, j)} is not flippable")
        if current == orientation:
            return False
        horizontal = orientation == "horizontal"
        self.horizontal_occupied[i, j] = horizontal
        self.horizontal_occupied[i + 1, j] = horizontal
        self.vertical_occupied[i, j] = not horizontal
        self.vertical_occupied[i, j + 1] = not horizontal
        return True

    def flip_face(self, i: int, j: int) -> FaceOrientation:
        """Rotate the two dimers and return the new orientation."""
        current = self.face_orientation(i, j)
        if current is None:
            raise ValueError(f"face {(i, j)} is not flippable")
        new: FaceOrientation = "vertical" if current == "horizontal" else "horizontal"
        self.set_face_orientation(i, j, new)
        return new
