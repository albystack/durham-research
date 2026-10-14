"""Geometry helpers for an even ``L x L`` square grid of vertices."""

from __future__ import annotations

from dataclasses import dataclass
from math import floor
from typing import Iterable


DEFAULT_RELATIVE_SEPARATIONS = (1 / 32, 1 / 16, 1 / 8, 1 / 4)


def validate_L(L: int) -> int:
    """Return ``L`` as an int after enforcing the v1 square-grid contract."""
    if isinstance(L, bool) or not isinstance(L, int):
        raise TypeError("L must be an integer")
    if L < 2:
        raise ValueError("L must be at least 2")
    if L % 2:
        raise ValueError("L must be even in the v1 square-grid implementation")
    return L


def center_face_index(L: int) -> tuple[int, int]:
    """Index of the unique central bounded face for even ``L``."""
    validate_L(L)
    center = (L - 2) // 2
    return center, center


@dataclass(frozen=True)
class SpatialSeparation:
    """A requested relative separation and its actual integer face offset."""

    requested_rho: float
    r: int


def spatial_separations(
    L: int,
    relative_separations: Iterable[float] = DEFAULT_RELATIVE_SEPARATIONS,
) -> tuple[SpatialSeparation, ...]:
    """Convert fractions of the vertex side length to valid rightward offsets.

    Half-integers are rounded upward. Zero, duplicate, and out-of-domain
    offsets are omitted. When requested fractions collide, the first one is
    retained and the actual integer ``r`` remains the operative label.
    """
    validate_L(L)
    _, center_col = center_face_index(L)
    maximum = (L - 2) - center_col
    result: list[SpatialSeparation] = []
    seen: set[int] = set()
    for value in relative_separations:
        rho = float(value)
        if not 0 < rho <= 1:
            raise ValueError("relative separations must lie in (0, 1]")
        r = floor(L * rho + 0.5)
        if r <= 0 or r > maximum or r in seen:
            continue
        result.append(SpatialSeparation(rho, r))
        seen.add(r)
    return tuple(result)
