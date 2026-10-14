#!/usr/bin/env python3
"""Weighted domino tilings of an even L-by-L square of vertices

For each fixed random environment, sample two tilings independently. Measure
height differences at offsets L/8 and L/4 along each diagonal coordinate;
physical distances are sqrt(2) times these offsets. The random edge families
have Gamma shapes 0.2 and 0.25, scale 1; other weights are one.

Run: python square.py --size 8 --environments 2 --output square-example
Requires Python 3.11+ and NumPy. The output folder must be NEW

"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import random
import resource
import sys
import time
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import floor, hypot
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

BoolArray = NDArray[np.bool_]
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
EdgeKind = Literal["horizontal", "vertical"]
Face = tuple[int, int] | None  # None is the exterior height reference.
DEFAULT_RELATIVE_SEPARATIONS = (1 / 32, 1 / 16, 1 / 8, 1 / 4)
VERSION = "embedded-square-shuffle-v1"
NAMESPACE = "matched-shuffle-20261004"


# 1. The square and its random edge weights


def validate_L(L: int) -> int:
    """Require an even number of vertices along each side."""
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


@dataclass
class DimerState:
    """A matching: True marks an occupied horizontal or vertical edge."""

    L: int
    horizontal_occupied: BoolArray
    vertical_occupied: BoolArray

    def __post_init__(self) -> None:
        validate_L(self.L)
        if not isinstance(self.horizontal_occupied, np.ndarray):
            raise TypeError("horizontal_occupied must be a NumPy ndarray")
        if not isinstance(self.vertical_occupied, np.ndarray):
            raise TypeError("vertical_occupied must be a NumPy ndarray")

    def validate(self) -> bool:
        """Check that every vertex belongs to exactly one dimer."""
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
            raise ValueError(
                f"vertex {(i, j)} has occupied degree {int(degree[i, j])}, not 1"
            )
        return True


@dataclass(frozen=True)
class DimerEnvironment:
    """Positive edge weights, copied and held fixed while tilings are sampled."""

    L: int
    horizontal_weights: FloatArray
    vertical_weights: FloatArray
    model: str = "custom"
    gamma_shape: float | None = None

    def __post_init__(self) -> None:
        L = validate_L(self.L)
        horizontal = np.array(self.horizontal_weights, dtype=np.float64, copy=True)
        vertical = np.array(self.vertical_weights, dtype=np.float64, copy=True)
        if horizontal.shape != (L, L - 1):
            raise ValueError(
                f"horizontal_weights must have shape {(L, L - 1)}, got {horizontal.shape}"
            )
        if vertical.shape != (L - 1, L):
            raise ValueError(
                f"vertical_weights must have shape {(L - 1, L)}, got {vertical.shape}"
            )
        if not np.all(np.isfinite(horizontal)) or not np.all(horizontal > 0):
            raise ValueError("horizontal edge weights must be finite and positive")
        if not np.all(np.isfinite(vertical)) or not np.all(vertical > 0):
            raise ValueError("vertical edge weights must be finite and positive")
        horizontal.setflags(write=False)
        vertical.setflags(write=False)
        object.__setattr__(self, "horizontal_weights", horizontal)
        object.__setattr__(self, "vertical_weights", vertical)


@dataclass(frozen=True)
class ComparisonModel:
    """Gamma parameters, or the fixed-mean and uniform control models."""

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


def family_masks(L: int) -> tuple[np.ndarray, np.ndarray]:
    """Locate the north (a) and east (b) edges of each black vertex."""
    validate_L(L)
    h_i, h_j = np.indices((L, L - 1))
    v_i, v_j = np.indices((L - 1, L))
    return (h_i + h_j) % 2 == 1, (v_i + v_j) % 2 == 0


def comparison_environment(
    L: int, rng: np.random.Generator, model: ComparisonModel
) -> DimerEnvironment:
    """Draw the two Gamma families once; all other edge weights equal one."""
    # North edges use a; east edges use b, relative to black vertices.
    b_mask, a_mask = family_masks(L)
    horizontal = np.ones((L, L - 1))
    vertical = np.ones((L - 1, L))
    if model.name == "aztec_gamma":
        vertical[a_mask] = rng.gamma(model.alpha, model.scale, a_mask.sum())
        horizontal[b_mask] = rng.gamma(model.beta, model.scale, b_mask.sum())
    elif model.name == "aztec_mean":
        vertical[a_mask] = model.alpha * model.scale
        horizontal[b_mask] = model.beta * model.scale
    return DimerEnvironment(L, horizontal, vertical, model=model.name)


# 2. Weighted shuffling inside an Aztec diamond

# Zero crossing weights isolate the square. The exterior can be tiled,
# so summing over its tilings leaves the desired weighted law on the square.


def randbelow_bits(rng, bound: int) -> int:
    """Draw a uniform integer from 0 to bound - 1 using fair random bits."""
    if not isinstance(bound, int) or bound <= 0:
        raise ValueError("positive integer bound required")
    if bound == 1:
        return 0
    bits = (bound - 1).bit_length()
    while True:
        value = rng.getrandbits(bits)
        if value < bound:
            return value


def bernoulli_ratio(rng, numerator: int, denominator: int) -> bool:
    """Draw a Bernoulli coin without rounding its probability to a float."""
    if (
        not isinstance(numerator, int)
        or not isinstance(denominator, int)
        or not 0 <= numerator <= denominator
        or denominator <= 0
    ):
        raise ValueError("require integer 0 <= numerator <= denominator > 0")
    return randbelow_bits(rng, denominator) < numerator


@dataclass(frozen=True)
class Diamond:
    """Vertices, edges and active 2-by-2 cells of the surrounding diamond."""

    order: int
    vertices: frozenset
    edges: tuple
    faces: tuple
    face_edges: tuple
    edge_ends: tuple
    incident: dict
    edge_index: dict


@lru_cache(maxsize=128)
def diamond(n):
    """Build the host graph. Coordinates are (row, column); 0=right, 1=down."""
    if not 1 <= n <= 127:
        raise ValueError("prototype supports diamond orders 1..127")
    vertices = frozenset(
        (r, c)
        for r in range(1 - n, n + 1)
        for c in range(1 - n, n + 1)
        if abs(2 * r - 1) + abs(2 * c - 1) <= 2 * n
    )
    edges, ends = [], []
    incident = {v: [] for v in vertices}
    for r, c in sorted(vertices):
        for direction, end in [(0, (r, c + 1)), (1, (r + 1, c))]:
            if end in vertices:
                index = len(edges)
                edges.append((r, c, direction))
                ends.append(((r, c), end))
                incident[(r, c)].append(index)
                incident[end].append(index)
    indices = {edge: i for i, edge in enumerate(edges)}
    faces, face_edges = [], []
    for r, c in sorted(vertices):
        corners = [(r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c)]
        if (r + c + n) % 2 == 1 and all(v in vertices for v in corners):
            faces.append((r, c))
            face_edges.append(
                tuple(
                    indices[e]
                    for e in [(r, c, 0), (r, c + 1, 1), (r + 1, c, 0), (r, c, 1)]
                )
            )
    assert len(vertices) == 2 * n * (n + 1) and len(faces) == n * n
    assert sorted(i for f in face_edges for i in f) == list(range(len(edges)))
    return Diamond(
        n,
        vertices,
        tuple(edges),
        tuple(faces),
        tuple(face_edges),
        tuple(ends),
        incident,
        indices,
    )


def square_embedding(environment, cast=Fraction):
    """Isolate the square with zero crossing weights.

    Also construct a tiling of the exterior, to check it can be covered."""
    L = environment.L
    g = diamond(L - 1)
    lo, hi = 1 - L // 2, L // 2
    inside = lambda v: lo <= v[0] <= hi and lo <= v[1] <= hi
    weights = []
    for (r, c, direction), (a, b) in zip(g.edges, g.edge_ends):
        if inside(a) and inside(b):
            array = (
                environment.horizontal_weights
                if direction == 0
                else environment.vertical_weights
            )
            weights.append(cast(float(array[r - lo, c - lo])))
        else:
            weights.append(cast(int(inside(a) == inside(b))))
    outside = set(g.vertices) - {v for v in g.vertices if inside(v)}
    witness = []
    for r in sorted({v[0] for v in outside if v[0] < lo or v[0] > hi}):
        columns = sorted(c for rr, c in outside if rr == r)
        assert len(columns) % 2 == 0
        for j in range(0, len(columns), 2):
            assert columns[j + 1] == columns[j] + 1
            witness.append(g.edge_index[(r, columns[j], 0)])
    covered = {v for i in witness for v in g.edge_ends[i]}
    remaining = outside - covered
    for c in sorted({v[1] for v in remaining}):
        rows = sorted(r for r, cc in remaining if cc == c)
        assert len(rows) % 2 == 0
        for j in range(0, len(rows), 2):
            assert rows[j + 1] == rows[j] + 1
            witness.append(g.edge_index[(rows[j], c, 1)])
    all_ends = [v for i in witness for v in g.edge_ends[i]]
    assert len(all_ends) == len(set(all_ends)) and set(all_ends) == outside
    return weights, tuple(witness)


def reduce_level(n, weights):
    """Replace each edge weight by its opposite weight divided by ac + bd.

    Restrict to the smaller diamond; reject singular cells and underflow."""
    if n <= 1:
        raise ValueError("cannot reduce order one")
    g, smaller = diamond(n), diamond(n - 1)
    number = type(weights[0])
    precomputed = [number(0)] * len(g.edges)
    for ids in g.face_edges:
        w = [weights[i] for i in ids]
        has_pair = (w[0] > 0 and w[2] > 0) or (w[1] > 0 and w[3] > 0)
        if not has_pair:
            raise ValueError("singular zero cell: outside supported embedding")
        products = [w[0] * w[2], w[1] * w[3]]
        if any(w[i] > 0 and w[i + 2] > 0 and products[i] == 0 for i in (0, 1)):
            raise ArithmeticError("positive pair product underflow")
        dp = sum(products)
        if not dp > 0 or (number is float and not math.isfinite(dp)):
            raise ArithmeticError("nonzero double product underflow/overflow")
        new = [w[(i + 2) % 4] / dp for i in range(4)]
        if any(w[(i + 2) % 4] > 0 and new[i] == 0 for i in range(4)):
            raise ArithmeticError("positive reduced edge underflow")
        if number is float and any(not math.isfinite(x) or x < 0 for x in new):
            raise ArithmeticError("invalid reduced weight")
        for i, value in zip(ids, new):
            precomputed[i] = value
    return [precomputed[g.edge_index[e]] for e in smaller.edges]


def coin_probability(weights, ids):
    """Horizontal-pair probability ac / (ac + bd), using exact stored weights."""
    a, b, c, d = [Fraction(weights[i]) for i in ids]
    horizontal, vertical = a * c, b * d
    if horizontal + vertical == 0:
        raise ValueError("reachable empty cell has zero double product")
    return horizontal / (horizontal + vertical)


def move_and_holes(n, previous):
    """Delete pairs, move single edges across their cell, and find empty cells."""
    g = diamond(n)
    occupied = (
        set() if n == 1 else {g.edge_index[diamond(n - 1).edges[i]] for i in previous}
    )
    moved, holes = set(), []
    for ids in g.face_edges:
        present = [j for j, i in enumerate(ids) if i in occupied]
        if len(present) == 1:
            moved.add(ids[(present[0] + 2) % 4])
        elif len(present) == 0:
            holes.append(ids)
        elif len(present) > 2:
            raise ValueError("invalid old matching")

    # Only cells empty before the move are filled; deleted pairs stay deleted.
    covered = [v for i in moved for v in g.edge_ends[i]]
    hole_vertices = [
        v for ids in holes for i in (ids[0], ids[2]) for v in g.edge_ends[i]
    ]
    if len(covered + hole_vertices) != len(set(covered + hole_vertices)):
        raise ValueError("overlap between moved edges and creation cells")
    return moved, holes


def validate_host(n, matching, weights):
    """Require a full matching with no occupied zero-weight edge."""
    g = diamond(n)
    vertices = [v for i in matching for v in g.edge_ends[i]]
    if len(vertices) != len(set(vertices)) or set(vertices) != set(g.vertices):
        raise ValueError("host matching does not cover vertices exactly once")
    if any(weights[i] <= 0 for i in matching):
        raise ValueError("host matching occupies a forbidden edge")


class EmbeddedSquareSampler:
    """Reduce weights once, then shuffle a fresh tiling for each seed.

    Fraction arithmetic is the small-grid reference. The campaign uses
    float reductions; its creation coins are exact for those stored weights."""

    def __init__(self, environment, cast=Fraction):
        self.L = environment.L
        self.order = self.L - 1
        self.levels = {}
        weights, self.exterior_witness = square_embedding(environment, cast)
        self.levels[self.order] = weights
        for n in range(self.order, 1, -1):
            weights = reduce_level(n, weights)
            self.levels[n - 1] = weights

    def extract(self, matching):
        """Keep the dimers inside the square and check that none cross its boundary."""
        lo, L = 1 - self.L // 2, self.L
        h, v = np.zeros((L, L - 1), bool), np.zeros((L - 1, L), bool)
        g = diamond(self.order)
        for i in matching:
            r, c, direction = g.edges[i]
            ends = g.edge_ends[i]
            inside = [lo <= a < lo + L and lo <= b < lo + L for a, b in ends]
            if any(inside) and not all(inside):
                raise ValueError("matching crosses target square boundary")
            if all(inside):
                (h if direction == 0 else v)[r - lo, c - lo] = True
        result = DimerState(L, h, v)
        result.validate()
        return result

    def sample(self, seed):
        """Grow a matching from order 1, filling each empty cell with a random pair."""
        rng = random.Random(int(seed))
        matching = set()
        for n in range(1, self.order + 1):
            matching, holes = move_and_holes(n, matching)
            weights = self.levels[n]
            for ids in holes:
                p = coin_probability(weights, ids)
                horizontal = bernoulli_ratio(rng, p.numerator, p.denominator)
                matching.update([ids[0], ids[2]] if horizontal else [ids[1], ids[3]])
            validate_host(n, matching, weights)
        return self.extract(matching)


# 3. Heights and the two diagonal measurements

# Vertices with even row + column are white. Clockwise around white,
# an empty edge contributes +1 and an occupied edge -3. Anchor the height
# through the top-left boundary edge, then integrate across interior edges.


@dataclass(frozen=True)
class HeightConstraint:
    """One height equation: h[right] - h[left] = increment."""

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
    """Orient edges white to black; crossing left to right adds 1 or -3."""
    L = state.L
    constraints: list[HeightConstraint] = []

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
    """Use interior crossings and one fixed boundary crossing as the anchor."""
    result = []
    for item in height_constraints(state):
        both_bounded = item.left_face is not None and item.right_face is not None
        anchor = item.kind == "horizontal" and item.edge_index == (0, 0)
        if both_bounded or anchor:
            result.append(item)
    return tuple(result)


def face_heights(state: DimerState, *, validate_matching: bool = True) -> IntArray:
    """Integrate height differences from the fixed boundary anchor.

    Check every alternative path gives the same integer height."""
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
    """Check every interior crossing and the boundary anchor."""
    expected = (state.L - 1, state.L - 1)
    if heights.shape != expected:
        raise ValueError(
            f"height array must have shape {expected}, got {heights.shape}"
        )
    if not np.issubdtype(heights.dtype, np.integer):
        raise TypeError("height array must have integer dtype")

    def value(face: Face) -> int:
        return 0 if face is None else int(heights[face])

    for item in (
        constraints if constraints is not None else _integration_constraints(state)
    ):
        difference = value(item.right_face) - value(item.left_face)
        if difference != item.increment:
            raise ValueError(
                f"height rule fails at {item.kind} edge {item.edge_index}: "
                f"got {difference}, expected {item.increment}"
            )
        if abs(difference) != (3 if item.occupied else 1):
            raise ValueError("height magnitude does not agree with edge occupancy")
    return True


@dataclass(frozen=True)
class SymmetricProbe:
    """Two face indices, centred on the grid; r is displacement per coordinate."""

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
    """Choose centred face pairs, rounding L * fraction half upward."""
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
    """Return h(end) - h(start); the additive height reference cancels."""
    if heights.ndim != 2 or heights.shape[0] != heights.shape[1]:
        raise ValueError("expected square bounded-face height array")
    return {p.label: int(heights[p.end]) - int(heights[p.start]) for p in probes}


# 4. Reproducible batches and saved results

FIELDS = [
    "domain",
    "size",
    "model",
    "law_id",
    "environment_id",
    "probe",
    "axis",
    "r",
    "distance",
    "height_scale",
    "value1",
    "value2",
    "environment_seed",
    "replica1_seed",
    "replica2_seed",
]


def peak_rss_bytes():
    """Peak memory use in bytes (the operating systems report different units)."""
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seed(master, domain, model, size, index, role):
    """Derive a repeatable seed for each environment or replica."""
    text = f"{NAMESPACE}:{master}:{domain}:{model}:{size}:{index}:{role}"
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "big")


def write_json(path, data):
    """Replace the output only after the complete JSON has been written."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)


def run_square(plan_path, cell, directory, validation=False):
    """For each environment, draw two tilings and save their height increments."""
    plan = json.loads(plan_path.read_text())
    count = plan["validation_pairs"] if validation else cell["pairs"]
    master = plan["validation_seed"] if validation else plan["master_seed"]
    directory.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    status = {
        "state": "STARTED",
        "completed": 0,
        "cell": cell,
        "pairs": count,
        "validation_only": validation,
        "plan_sha256": digest(plan_path),
        "sampler_version": VERSION,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "master_seed": master,
    }
    write_json(directory / "manifest.json", status)
    size, model = cell["size"], cell["model"]
    probes = symmetric_probes(size, fractions=(0.125, 0.25), axes=("diagonal",))
    law = f"{model}:alpha=0.2:beta=0.25:scale=1:a-north,b-east"
    temporary = directory / "measurements.csv.tmp"
    try:
        with (
            temporary.open("w", newline="") as out,
            (directory / "certificates.jsonl").open("w") as certificates,
        ):
            writer = csv.DictWriter(out, fieldnames=FIELDS)
            writer.writeheader()
            for index in range(1, count + 1):
                seeds = [
                    seed(master, "square", model, size, index, role)
                    for role in ("environment", "replica1", "replica2")
                ]
                env = comparison_environment(
                    size, np.random.default_rng(seeds[0]), ComparisonModel(name=model)
                )
                case = directory / "cases" / f"env-{index:05d}"
                case.mkdir(parents=True)
                np.savez_compressed(
                    case / "weights.npz",
                    horizontal=env.horizontal_weights,
                    vertical=env.vertical_weights,
                )

                # Share the weights, but draw each replica from its own seed.
                sampler = EmbeddedSquareSampler(env, cast=float)
                states = [sampler.sample(s) for s in seeds[1:]]
                values = [measure_probes(face_heights(s), probes) for s in states]
                for replica, state in enumerate(states, 1):
                    state.validate()
                    np.savez_compressed(
                        case / f"sample{replica}.npz",
                        horizontal=state.horizontal_occupied,
                        vertical=state.vertical_occupied,
                    )
                rows = []
                for probe in probes:
                    rows.append(
                        {
                            "domain": "square",
                            "size": size,
                            "model": model,
                            "law_id": law,
                            "environment_id": f"square-{model}-n{size}-env{index}",
                            "probe": probe.label,
                            "axis": "diagonal",
                            "r": probe.r,
                            "distance": probe.distance,
                            "height_scale": 1,
                            "value1": values[0][probe.label],
                            "value2": values[1][probe.label],
                            "environment_seed": str(seeds[0]),
                            "replica1_seed": str(seeds[1]),
                            "replica2_seed": str(seeds[2]),
                        }
                    )
                write_json(case / "rows.json", rows)
                certificate = {
                    "index": index,
                    "files": {
                        name: digest(case / name)
                        for name in (
                            "weights.npz",
                            "sample1.npz",
                            "sample2.npz",
                            "rows.json",
                        )
                    },
                }
                certificates.write(json.dumps(certificate) + "\n")
                writer.writerows(rows)
                status["completed"] = index
                if index % 64 == 0:
                    out.flush()
                    certificates.flush()
                    write_json(directory / "status.json", status)
                    print(f"square {model} L={size}: {index}/{count}", flush=True)
                del sampler, states
        temporary.replace(directory / "measurements.csv")
        status.update(
            state="COMPLETE",
            measurements_sha256=digest(directory / "measurements.csv"),
            certificates_sha256=digest(directory / "certificates.jsonl"),
            arithmetic="Float64 reductions; exact rational creation coins for represented weights; finite numerical audit, not universal exactness",
        )
    except Exception as error:
        status.update(state="ERROR", reason=f"{type(error).__name__}: {error}")
        raise
    finally:
        status.update(
            seconds=time.perf_counter() - started, peak_rss_bytes=peak_rss_bytes()
        )
        write_json(directory / "status.json", status)
    print(json.dumps(status), flush=True)
    return status


def main():
    """Read options and run a small batch; no arguments means two pairs at L=8."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--size", type=int, default=8, help="Even vertex side length, 8..128"
    )
    parser.add_argument(
        "--environments", type=int, default=2, help="Pairs in this batch, 1..512"
    )
    parser.add_argument(
        "--model",
        choices=("aztec_gamma", "aztec_mean", "uniform"),
        default="aztec_gamma",
    )
    parser.add_argument("--seed", type=int, default=2026101301)
    parser.add_argument("--output", type=Path, default=Path("square-example"))
    args = parser.parse_args()
    if args.size not in (8, 16, 24, 32, 48, 64, 96, 128):
        parser.error("Use a supported size: 8, 16, 24, 32, 48, 64, 96 or 128")
    if not 1 <= args.environments <= 512 or args.seed < 0:
        parser.error("Need 1..512 environments and a nonnegative seed")
    args.output.mkdir(parents=True, exist_ok=False)
    cell = dict(
        index=0,
        domain="square",
        model=args.model,
        size=args.size,
        pairs=args.environments,
    )
    plan = dict(
        namespace=NAMESPACE,
        master_seed=args.seed,
        validation_seed=args.seed,
        validation_pairs=args.environments,
        cells=[cell],
        purpose="Standalone demonstration",
    )
    plan_path = args.output / "plan.json"
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    run_square(plan_path, cell, args.output / "samples", validation=True)


if __name__ == "__main__":
    main()
