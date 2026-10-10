"""Geometry and disorder for the prospective Aztec/square comparison.

This module deliberately leaves the September model, seeds and rightward probes
unchanged. New campaigns use an explicit specification and a new seed namespace.

Figure 1 of Duits--Van Peski (arXiv:2512.03033) has random red (a) and yellow (b)
edges. In the unrotated lattice used by AztecDiamond.covered_cells, these point
north and east from black vertices. Here black vertices have odd row+column;
rows increase downwards. Other incident directions have weight one.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from math import floor, hypot

import numpy as np

from .environment import DimerEnvironment
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, center_face_index, validate_L


@dataclass(frozen=True)
class ComparisonModel:
    """Full law, including deterministic controls; Gamma uses shape/scale."""

    name: str = "aztec_gamma"
    alpha: float = 0.2
    beta: float = 0.25
    scale: float = 1.0

    def __post_init__(self) -> None:
        if self.name not in ("aztec_gamma", "aztec_mean", "uniform"):
            raise ValueError("model must be aztec_gamma, aztec_mean or uniform")
        if not all(
            np.isfinite(x) and x > 0 for x in (self.alpha, self.beta, self.scale)
        ):
            raise ValueError("Gamma shapes and scale must be finite and positive")

    def identity(self) -> dict:
        return {"schema": 1, "placement": "black-odd:a-north,b-east:v1", **asdict(self)}


def family_masks(L: int) -> tuple[np.ndarray, np.ndarray]:
    """Return b-horizontal and a-vertical masks (random edges only).

    Horizontal left endpoint black: (i+j) odd. Vertical bottom endpoint
    black: (i+j+1) odd, hence top coordinate (i+j) even.
    """
    validate_L(L)
    h_i, h_j = np.indices((L, L - 1))
    v_i, v_j = np.indices((L - 1, L))
    return (h_i + h_j) % 2 == 1, (v_i + v_j) % 2 == 0


def comparison_environment(
    L: int, rng: np.random.Generator, model: ComparisonModel
) -> DimerEnvironment:
    """Restrict the stated independent two-family law to an even square."""
    b_mask, a_mask = family_masks(L)
    horizontal = np.ones((L, L - 1))
    vertical = np.ones((L - 1, L))
    if model.name == "aztec_gamma":
        # Draw each undirected random edge once; do not normalise each family.
        vertical[a_mask] = rng.gamma(model.alpha, model.scale, a_mask.sum())
        horizontal[b_mask] = rng.gamma(model.beta, model.scale, b_mask.sum())
    elif model.name == "aztec_mean":
        vertical[a_mask] = model.alpha * model.scale
        horizontal[b_mask] = model.beta * model.scale
    return DimerEnvironment(L, horizontal, vertical, model=model.name)


def comparison_seeds(
    master_seed: int, L: int, model: ComparisonModel, environment_index: int
) -> tuple[int, int, int]:
    """Parameter-aware PCG64 seed namespace, independent of legacy campaigns."""
    validate_L(L)
    if master_seed < 0 or environment_index < 0:
        raise ValueError("seeds and environment indices must be nonnegative")
    encoded = json.dumps(model.identity(), sort_keys=True).encode()
    words = np.frombuffer(hashlib.sha256(encoded).digest(), dtype="<u4").tolist()
    root = np.random.SeedSequence([0xA27EC, master_seed, L, environment_index, *words])
    seeds = tuple(
        int(child.generate_state(1, dtype=np.uint64)[0]) for child in root.spawn(3)
    )
    if len(set(seeds)) != 3:
        raise RuntimeError("seed collision")
    return seeds


@dataclass(frozen=True)
class SymmetricProbe:
    """Endpoints are bounded-face array indices, not vertex indices.

    r is total displacement in each active coordinate. Diagonal length is
    sqrt(2)*r in units of one primal edge. Odd separations have midpoint offset
    -1/2 in each active coordinate (a fixed, documented tie-break).
    """

    requested_rho: float
    r: int
    axis: str
    start: tuple[int, int]
    end: tuple[int, int]

    @property
    def distance(self) -> float:
        return hypot(self.end[0] - self.start[0], self.end[1] - self.start[1])

    @property
    def label(self) -> str:
        return f"{self.axis}:rho={self.requested_rho:.12g}"


def symmetric_probes(
    L: int, fractions=DEFAULT_RELATIVE_SEPARATIONS, axes=("diagonal", "horizontal")
) -> tuple[SymmetricProbe, ...]:
    """Use round-half-up and omit zero/duplicate offsets within each axis.

    Aztec's even-order central staggered row has physical displacement (r,-r)
    in (x,y), which is (r,r) in row/column arrays. Thus diagonal is primary.
    """
    center, _ = center_face_index(L)
    probes = []
    for axis in axes:
        if axis not in ("diagonal", "horizontal"):
            raise ValueError("axis must be diagonal or horizontal")
        seen = set()
        for fraction in fractions:
            rho = float(fraction)
            if not np.isfinite(rho) or not 0 < rho <= 1:
                raise ValueError("fractions must be finite and in (0,1]")
            r = floor(L * rho + 0.5)
            if r == 0 or r in seen:
                continue
            lo = center - (r + 1) // 2
            hi = lo + r
            if lo < 0 or hi >= L - 1:
                continue
            seen.add(r)
            start = (lo, lo) if axis == "diagonal" else (center, lo)
            end = (hi, hi) if axis == "diagonal" else (center, hi)
            probes.append(SymmetricProbe(rho, r, axis, start, end))
    return tuple(probes)


def measure_probes(
    heights: np.ndarray, probes: tuple[SymmetricProbe, ...]
) -> dict[str, int]:
    """Height differences cancel the common boundary reference."""
    if heights.ndim != 2 or heights.shape[0] != heights.shape[1]:
        raise ValueError("expected square bounded-face height array")
    return {p.label: int(heights[p.end]) - int(heights[p.start]) for p in probes}


def aztec_face_xy(order: int, row: int, column: int) -> tuple[int, int]:
    """Physical Aztec face position in units of one primal edge.

    Julia indices are one-based. Every step in a height-integration column
    crosses alternately a horizontal and a vertical edge. Starting at
    (column-1, order+1-column), move south, west, south, west, ... .
    """
    if order < 1 or not 1 <= row <= 2 * order + 1 or not 1 <= column <= order:
        raise ValueError("invalid staggered face index")
    steps = row - 1
    return column - 1 - steps // 2, order + 1 - column - (steps + 1) // 2
