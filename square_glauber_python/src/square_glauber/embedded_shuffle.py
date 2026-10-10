"""Experimental square sampling by embedded weighted domino shuffling.

Reference: Janvresse, de la Rue & Velenik, A Note on Domino Shuffling, Sec. 2.
A square of side L is isolated inside diamond order L-1 by zero crossing
weights. Its exterior is tileable and independent, so marginalising it leaves
the square Gibbs law. No square edges or boundary conditions are replaced.

Rational arithmetic supplies a small-grid oracle. Float arithmetic is a
separately audited approximation for the reduction; creation coins use exact
products of represented weights. Structural zeros are never epsilon-clipped. Singular zero cells are rejected.
"""

import math
import random
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

import numpy as np

from .integer_sampling import bernoulli_ratio
from .matching import DimerState

VERSION = "embedded-square-shuffle-v1"


@dataclass(frozen=True)
class Diamond:
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
    """Coordinates are (row,column); edge direction 0=right, 1=down."""
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
    """Return host weights plus an explicit exterior matching witness.

    North/south exterior strips are matched horizontally; remaining west/east
    strips vertically. This independently checks that the complement is tileable.
    """
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
    """Complement cells and restrict to the smaller diamond.

    This prototype supports zero edges only while every cell double product
    stays positive. For this regime the usual positive-weight distribution
    identity extends by continuity to zero edges. Singular zero-cell rules
    are deliberately unsupported: fail rather than invent a probability.
    The centred even-square embedding has passed this support check at every
    tested level. Arithmetic underflow is also a failure, never a new zero.
    """
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
    """Exact creation probability for the represented weights, including tails."""
    a, b, c, d = [Fraction(weights[i]) for i in ids]
    horizontal, vertical = a * c, b * d
    if horizontal + vertical == 0:
        raise ValueError("reachable empty cell has zero double product")
    return horizontal / (horizontal + vertical)


def move_and_holes(n, previous):
    """Delete pairs and move single edges across active cells, then find holes."""
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
    # Creation is determined by OLD empty active cells. Looking only for
    # unmatched vertices after pair deletion would also refill deleted cells
    # and create overlapping candidate holes.
    covered = [v for i in moved for v in g.edge_ends[i]]
    hole_vertices = [
        v for ids in holes for i in (ids[0], ids[2]) for v in g.edge_ends[i]
    ]
    if len(covered + hole_vertices) != len(set(covered + hole_vertices)):
        raise ValueError("overlap between moved edges and creation cells")
    return moved, holes


def validate_host(n, matching, weights):
    g = diamond(n)
    vertices = [v for i in matching for v in g.edge_ends[i]]
    if len(vertices) != len(set(vertices)) or set(vertices) != set(g.vertices):
        raise ValueError("host matching does not cover vertices exactly once")
    if any(weights[i] <= 0 for i in matching):
        raise ValueError("host matching occupies a forbidden edge")


class EmbeddedSquareSampler:
    """Rational reference or audited floating reduction; no CFTP burn-in.

    Construction is cubic in side length in arithmetic operations, rather
    than exponential in frontier width. Rational bit complexity is not bounded
    by that operation count. Use only within a documented validation scope;
    floating reductions do not inherit the rational oracle's exactness.
    """

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
