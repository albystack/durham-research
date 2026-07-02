"""Simulation code for quenched random edge environments and LERW winding.

The project uses the following computational proxy for the dimer height at the
origin.  In a box [-L, L]^2, run a random walk from the origin until it first
hits the boundary.  Loop-erase that path.  The winding of the loop-erased path
is the tree-side observable motivated by Temperley's bijection.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from config import DIRECTIONS, STEP_VECTOR, parse_model, sample_edge_weight, stable_seed


Point = tuple[int, int]
Edge = tuple[Point, Point]


# Direction codes are arranged so that a +1 change is a left turn and a -1
# change is a right turn.  The winding is counted in quarter-turns.
DIRECTION_CODE = {
    (1, 0): 0,   # east
    (0, 1): 1,   # north
    (-1, 0): 2,  # west
    (0, -1): 3,  # south
}


@dataclass(frozen=True)
class TrialTask:
    """Everything needed to run one independent sample."""

    L: int
    model: str
    parameter: float | None
    sample_number: int
    seed: int
    max_steps: int | None


class RandomEdgeEnvironment:
    """Lazy quenched nearest-neighbour conductance environment.

    Each undirected edge receives one positive weight the first time the walk
    needs that edge.  The same value is reused later, even if the walk crosses
    the edge in the opposite direction.  This matches the idea of fixed random
    edge weights rather than fresh random transition probabilities at every
    visit.
    """

    def __init__(self, rng: random.Random, model: str, parameter: float | None):
        self.rng = rng
        self.model = model
        self.parameter = parameter
        self.edge_weights: dict[Edge, float] = {}

    @staticmethod
    def edge_key(a: Point, b: Point) -> Edge:
        """Canonical key for an undirected edge."""

        return (a, b) if a <= b else (b, a)

    def weight_between(self, a: Point, b: Point) -> float:
        """Return the fixed random weight on edge {a,b}, sampling if needed."""

        key = self.edge_key(a, b)
        if key not in self.edge_weights:
            self.edge_weights[key] = sample_edge_weight(self.rng, self.model, self.parameter)
        return self.edge_weights[key]

    def weights_at(self, point: Point) -> tuple[float, float, float, float]:
        """Return incident edge weights in N,E,S,W order."""

        weights: list[float] = []
        for direction in DIRECTIONS:
            dx, dy = STEP_VECTOR[direction]
            neighbour = (point[0] + dx, point[1] + dy)
            weights.append(self.weight_between(point, neighbour))
        return weights[0], weights[1], weights[2], weights[3]

    @property
    def sampled_edge_count(self) -> int:
        """Number of distinct edge weights exposed during the trial."""

        return len(self.edge_weights)


class SiteDirectionalEnvironment:
    """Lazy fixed site environment with directional weights N=E=1, S=u, W=v.

    This is the older drift-prone model from the June 8 email.  It is kept for
    reproducibility, but the HPC runner defaults to `SiteIIDEnvironment`.

    The value at each site is generated from a seed that depends only on the
    environment seed, model, parameter, and lattice coordinate.  This matters
    for quenched batches: adding more walks or changing walk order must not
    change already-defined site weights in the sampled environment.
    """

    def __init__(self, environment_seed: int, model: str, parameter: float | None):
        self.environment_seed = environment_seed
        self.model = model
        self.parameter = parameter
        self.site_weights: dict[Point, tuple[float, float]] = {}

    def _site_rng(self, point: Point) -> random.Random:
        parameter_key = "" if self.parameter is None else f"{self.parameter:.17g}"
        seed = stable_seed(self.environment_seed, self.model, parameter_key, point[0], point[1])
        return random.Random(seed)

    def uv_at(self, point: Point) -> tuple[float, float]:
        """Return the fixed positive (u, v) pair at a site."""

        if point not in self.site_weights:
            rng = self._site_rng(point)
            u = sample_edge_weight(rng, self.model, self.parameter)
            v = sample_edge_weight(rng, self.model, self.parameter)
            if u <= 0.0 or v <= 0.0:
                raise ValueError(f"non-positive site weights at {point}: u={u}, v={v}")
            self.site_weights[point] = (u, v)
        return self.site_weights[point]

    def weights_at(self, point: Point) -> tuple[float, float, float, float]:
        """Return directional weights in N,E,S,W order."""

        u, v = self.uv_at(point)
        return 1.0, 1.0, u, v

    def transition_probabilities(self, point: Point) -> tuple[float, float, float, float]:
        """Return normalised transition probabilities in N,E,S,W order."""

        weights = self.weights_at(point)
        total = sum(weights)
        return tuple(weight / total for weight in weights)

    @property
    def sampled_site_count(self) -> int:
        """Number of distinct site weights exposed so far."""

        return len(self.site_weights)


class SiteIIDEnvironment:
    """Lazy fixed site environment with IID N,E,S,W positive weights.

    Each site receives four independent samples from the same positive
    distribution, one for each outgoing direction, and transition probabilities
    are obtained by normalising these four weights.  This is the locally
    balanced model requested after the north-east drift diagnostic: no compass
    direction is distinguished in the sampling law, while each site remains
    fixed for all walks in a quenched batch.
    """

    def __init__(self, environment_seed: int, model: str, parameter: float | None):
        self.environment_seed = environment_seed
        self.model = model
        self.parameter = parameter
        self.site_weights: dict[Point, tuple[float, float, float, float]] = {}

    def _site_rng(self, point: Point) -> random.Random:
        parameter_key = "" if self.parameter is None else f"{self.parameter:.17g}"
        seed = stable_seed("site_iid", self.environment_seed, self.model, parameter_key, point[0], point[1])
        return random.Random(seed)

    def weights_at(self, point: Point) -> tuple[float, float, float, float]:
        """Return fixed IID directional weights in N,E,S,W order."""

        if point not in self.site_weights:
            rng = self._site_rng(point)
            weights = tuple(sample_edge_weight(rng, self.model, self.parameter) for _ in DIRECTIONS)
            if any(weight <= 0.0 for weight in weights):
                raise ValueError(f"non-positive site weights at {point}: {weights}")
            self.site_weights[point] = weights
        return self.site_weights[point]

    def transition_probabilities(self, point: Point) -> tuple[float, float, float, float]:
        """Return normalised transition probabilities in N,E,S,W order."""

        weights = self.weights_at(point)
        total = sum(weights)
        return tuple(weight / total for weight in weights)

    @property
    def sampled_site_count(self) -> int:
        """Number of distinct site weights exposed so far."""

        return len(self.site_weights)


def is_boundary(point: Point, L: int) -> bool:
    """True exactly on the boundary of the square box [-L,L]^2."""

    x, y = point
    return abs(x) == L or abs(y) == L


def choose_next_point(point: Point, environment: RandomEdgeEnvironment, rng: random.Random) -> Point:
    """Choose the next walk step using the four incident edge weights."""

    weights = environment.weights_at(point)
    total_weight = sum(weights)

    # Draw from the categorical distribution without building a probability
    # list.  This is both clear and avoids tiny rounding errors in normalising.
    threshold = rng.random() * total_weight
    running_total = 0.0
    for direction, weight in zip(DIRECTIONS, weights):
        running_total += weight
        if threshold <= running_total:
            dx, dy = STEP_VECTOR[direction]
            return (point[0] + dx, point[1] + dy)

    # Floating point roundoff can only leave us infinitesimally above the final
    # cumulative weight.  Falling back to west keeps the function total.
    dx, dy = STEP_VECTOR["W"]
    return (point[0] + dx, point[1] + dy)


def random_walk_until_boundary(
    L: int,
    environment: RandomEdgeEnvironment,
    rng: random.Random,
    max_steps: int | None,
) -> list[Point]:
    """Run a nearest-neighbour walk from the origin to the box boundary."""

    if L < 1:
        raise ValueError("L must be at least 1")
    if max_steps is None:
        # The safety cap should almost never matter in ordinary runs.  It is
        # here to catch pathological parameter choices or accidental bugs.
        max_steps = max(100_000, 2_000 * L * L)

    path = [(0, 0)]
    current = (0, 0)
    for _ in range(max_steps):
        if is_boundary(current, L):
            return path
        current = choose_next_point(current, environment, rng)
        path.append(current)

    raise RuntimeError(f"walk did not hit the boundary within {max_steps} steps")


def loop_erased_walk_until_boundary(
    L: int,
    environment: RandomEdgeEnvironment,
    rng: random.Random,
    max_steps: int | None,
) -> tuple[list[Point], int]:
    """Run a walk and maintain its loop erasure online.

    This is equivalent to storing the full random walk and calling
    `loop_erase` afterwards, but it avoids keeping enormous raw paths in memory
    and removes a second pass over the path.  For the publication probe at
    L=512, this matters.
    """

    if L < 1:
        raise ValueError("L must be at least 1")
    if max_steps is None:
        max_steps = max(100_000, 2_000 * L * L)

    current = (0, 0)
    raw_steps = 0
    erased: list[Point] = [current]
    position: dict[Point, int] = {current: 0}

    for _ in range(max_steps):
        if is_boundary(current, L):
            return erased, raw_steps

        current = choose_next_point(current, environment, rng)
        raw_steps += 1

        if current in position:
            keep_until = position[current]
            for removed in erased[keep_until + 1 :]:
                del position[removed]
            erased = erased[: keep_until + 1]
        else:
            position[current] = len(erased)
            erased.append(current)

    raise RuntimeError(f"walk did not hit the boundary within {max_steps} steps")


def loop_erase(path: list[Point]) -> list[Point]:
    """Chronological loop erasure of a lattice path.

    When the path returns to a point already in the erased path, everything
    after the previous visit is deleted.  The result is a self-avoiding path
    with the same start and final boundary hit.
    """

    erased: list[Point] = []
    position: dict[Point, int] = {}

    for point in path:
        if point in position:
            keep_until = position[point]
            for removed in erased[keep_until + 1 :]:
                del position[removed]
            erased = erased[: keep_until + 1]
            continue

        position[point] = len(erased)
        erased.append(point)

    return erased


def winding(path: list[Point]) -> int:
    """Left-turn count minus right-turn count for a self-avoiding path."""

    if len(path) < 3:
        return 0

    directions: list[int] = []
    for start, end in zip(path, path[1:]):
        step = (end[0] - start[0], end[1] - start[1])
        directions.append(DIRECTION_CODE[step])

    turns = 0
    for before, after in zip(directions, directions[1:]):
        turn = (after - before) % 4
        if turn == 1:
            turns += 1
        elif turn == 3:
            turns -= 1
        elif turn == 2:
            raise ValueError("loop-erased path contains an immediate U-turn")
    return turns


def winding_angle(path: list[Point]) -> float:
    """The same winding measured in radians."""

    return 0.5 * math.pi * winding(path)


def make_trial_tasks(
    sizes: list[int],
    samples_per_size: int,
    model_strings: list[str],
    base_seed: int,
    max_steps: int | None,
) -> list[TrialTask]:
    """Expand a preset into independent trial tasks."""

    tasks: list[TrialTask] = []
    for model_string in model_strings:
        model, parameter = parse_model(model_string)
        for L in sizes:
            for sample_number in range(samples_per_size):
                seed = stable_seed(base_seed, model, parameter, L, sample_number)
                tasks.append(
                    TrialTask(
                        L=L,
                        model=model,
                        parameter=parameter,
                        sample_number=sample_number,
                        seed=seed,
                        max_steps=max_steps,
                    )
                )
    return tasks


def run_trial(task: TrialTask) -> dict[str, int | float | str]:
    """Run one random environment, one walk, and one loop erasure."""

    # Use separate random streams so the quenched environment samples are not
    # entangled with the walk's transition draws in an avoidable way.
    edge_rng = random.Random(stable_seed(task.seed, "edges"))
    walk_rng = random.Random(stable_seed(task.seed, "walk"))

    environment = RandomEdgeEnvironment(edge_rng, task.model, task.parameter)
    erased_path, raw_steps = loop_erased_walk_until_boundary(
        task.L,
        environment,
        walk_rng,
        task.max_steps,
    )

    hit_x, hit_y = erased_path[-1]
    winding_count = winding(erased_path)
    angle = 0.5 * math.pi * winding_count

    return {
        "model": task.model,
        "parameter": "" if task.parameter is None else task.parameter,
        "L": task.L,
        "sample": task.sample_number,
        "seed": task.seed,
        "raw_steps": raw_steps,
        "lerw_steps": len(erased_path) - 1,
        "winding": winding_count,
        "winding_angle": angle,
        "height_proxy": winding_count,
        "height_proxy_angle": angle,
        "hit_x": hit_x,
        "hit_y": hit_y,
        "environment_edges": environment.sampled_edge_count,
    }
