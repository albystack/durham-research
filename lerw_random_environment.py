#!/usr/bin/env python3
"""Simulate LERW winding in symmetric and fixed random environments.

The main model is the one in project_brief.md:

    w_N = 1, w_E = 1, w_S = u_x, w_W = v_x

with u_x and v_x positive, mean-one random variables sampled once per site and
then kept fixed for the whole walk.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, variance
from typing import Iterable

Point = tuple[int, int]

# Direction codes are cyclic counterclockwise, which makes winding simple.
# E -> N is a left turn (+1); E -> S is a right turn (-1).
DIRECTION_CODES: dict[Point, int] = {
    (1, 0): 0,   # E
    (0, 1): 1,   # N
    (-1, 0): 2,  # W
    (0, -1): 3,  # S
}

# Sampling order follows the email: N, E, S, W.
STEP_VECTORS: tuple[tuple[str, Point], ...] = (
    ("N", (0, 1)),
    ("E", (1, 0)),
    ("S", (0, -1)),
    ("W", (-1, 0)),
)


@dataclass(frozen=True)
class TrialResult:
    L: int
    model: str
    k: float | None
    sample_index: int
    seed: int
    raw_steps: int
    lerw_steps: int
    winding: int
    hit_x: int
    hit_y: int
    environment_sites: int


class RandomEnvironment:
    """Fixed site-dependent environment, sampled lazily.

    Laziness is equivalent to pre-sampling the whole box because the site
    weights are independent. It avoids spending time on sites the walk never
    reaches.
    """

    def __init__(self, rng: random.Random, model: str, k: float | None = None):
        if model not in {"symmetric", "gamma", "exponential", "gamma4"}:
            raise ValueError(f"unknown model: {model}")
        if model in {"gamma", "gamma4"} and (k is None or k <= 0):
            raise ValueError(f"{model} model requires k > 0")

        self.rng = rng
        self.model = model
        self.k = k
        self._weights: dict[Point, tuple[float, ...]] = {}

    @property
    def sampled_site_count(self) -> int:
        return len(self._weights)

    def weights_at(self, point: Point) -> tuple[float, ...]:
        """Return the fixed random weights for one site."""

        if self.model == "symmetric":
            return (1.0, 1.0)

        if point not in self._weights:
            if self.model == "exponential":
                u = self.rng.expovariate(1.0)
                v = self.rng.expovariate(1.0)
                weights = (u, v)
            elif self.model == "gamma4":
                assert self.k is not None
                weights = tuple(
                    self.rng.gammavariate(self.k, 1.0 / self.k)
                    for _ in STEP_VECTORS
                )
            else:
                assert self.k is not None
                u = self.rng.gammavariate(self.k, 1.0 / self.k)
                v = self.rng.gammavariate(self.k, 1.0 / self.k)
                weights = (u, v)
            self._weights[point] = weights

        return self._weights[point]

    def transition_weights(self, point: Point) -> tuple[float, float, float, float]:
        """Return unnormalised weights in N, E, S, W order."""

        if self.model == "symmetric":
            return (1.0, 1.0, 1.0, 1.0)
        if self.model == "gamma4":
            weights = self.weights_at(point)
            return (weights[0], weights[1], weights[2], weights[3])
        u, v = self.weights_at(point)
        return (1.0, 1.0, u, v)

    def transition_probabilities(self, point: Point) -> tuple[float, float, float, float]:
        weights = self.transition_weights(point)
        total = sum(weights)
        return tuple(weight / total for weight in weights)


def is_boundary(point: Point, L: int) -> bool:
    x, y = point
    return abs(x) == L or abs(y) == L


def validate_L(L: int) -> None:
    if L < 1:
        raise ValueError("L must be at least 1")


def choose_step(point: Point, env: RandomEnvironment) -> Point:
    """Sample the next step from the transition weights at `point`."""

    weights = env.transition_weights(point)
    total = sum(weights)
    threshold = env.rng.random() * total

    cumulative = 0.0
    for (_, step), weight in zip(STEP_VECTORS, weights):
        cumulative += weight
        if threshold <= cumulative:
            return (point[0] + step[0], point[1] + step[1])

    # Floating-point fallback for threshold extremely close to total.
    step = STEP_VECTORS[-1][1]
    return (point[0] + step[0], point[1] + step[1])


def random_walk_to_boundary(
    L: int,
    env: RandomEnvironment,
    max_steps: int | None = None,
) -> list[Point]:
    """Run the weighted walk from the origin until it first hits the boundary."""

    validate_L(L)
    if max_steps is None:
        max_steps = max(100_000, 1_000 * L * L)

    path: list[Point] = [(0, 0)]
    current = (0, 0)

    for _ in range(max_steps):
        if is_boundary(current, L):
            return path
        current = choose_step(current, env)
        path.append(current)

    raise RuntimeError(f"walk failed to hit boundary within {max_steps} steps")


def loop_erase(path: Iterable[Point]) -> list[Point]:
    """Chronologically erase loops from a nearest-neighbour path."""

    erased: list[Point] = []
    index: dict[Point, int] = {}

    for point in path:
        if point in index:
            keep_until = index[point]
            for removed in erased[keep_until + 1 :]:
                del index[removed]
            erased = erased[: keep_until + 1]
        else:
            index[point] = len(erased)
            erased.append(point)

    return erased


def direction_code(start: Point, end: Point) -> int:
    step = (end[0] - start[0], end[1] - start[1])
    try:
        return DIRECTION_CODES[step]
    except KeyError as exc:
        raise ValueError(f"points are not nearest neighbours: {start} -> {end}") from exc


def winding(path: list[Point]) -> int:
    """Return number of left turns minus number of right turns."""

    if len(path) < 3:
        return 0

    directions = [
        direction_code(path[i], path[i + 1])
        for i in range(len(path) - 1)
    ]

    total = 0
    for before, after in zip(directions, directions[1:]):
        turn = (after - before) % 4
        if turn == 1:
            total += 1
        elif turn == 3:
            total -= 1
        elif turn == 2:
            raise ValueError("loop-erased path contains an immediate U-turn")

    return total


def run_trial(
    L: int,
    model: str,
    k: float | None,
    sample_index: int,
    seed: int,
    max_steps: int | None = None,
) -> TrialResult:
    rng = random.Random(seed)
    env = RandomEnvironment(rng=rng, model=model, k=k)
    raw_path = random_walk_to_boundary(L=L, env=env, max_steps=max_steps)
    erased_path = loop_erase(raw_path)
    hit_x, hit_y = erased_path[-1]

    return TrialResult(
        L=L,
        model=model,
        k=k,
        sample_index=sample_index,
        seed=seed,
        raw_steps=len(raw_path) - 1,
        lerw_steps=len(erased_path) - 1,
        winding=winding(erased_path),
        hit_x=hit_x,
        hit_y=hit_y,
        environment_sites=env.sampled_site_count,
    )


def sample_seeds(base_seed: int, count: int) -> list[int]:
    rng = random.Random(base_seed)
    return [rng.randrange(0, 2**63) for _ in range(count)]


def run_experiment(
    sizes: list[int],
    samples: int,
    model: str,
    k: float | None,
    seed: int,
    max_steps: int | None,
) -> list[TrialResult]:
    total_trials = len(sizes) * samples
    seeds = sample_seeds(seed, total_trials)
    results: list[TrialResult] = []
    seed_index = 0

    for L in sizes:
        for sample_index in range(samples):
            result = run_trial(
                L=L,
                model=model,
                k=k,
                sample_index=sample_index,
                seed=seeds[seed_index],
                max_steps=max_steps,
            )
            results.append(result)
            seed_index += 1

    return results


def write_csv(path: Path, results: list[TrialResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(TrialResult.__dataclass_fields__.keys())

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)


def print_summary(results: list[TrialResult]) -> None:
    by_size: dict[int, list[TrialResult]] = {}
    for result in results:
        by_size.setdefault(result.L, []).append(result)

    print("L,samples,mean_winding,var_winding,mean_raw_steps,mean_lerw_steps")
    for L in sorted(by_size):
        group = by_size[L]
        windings = [item.winding for item in group]
        var_winding = variance(windings) if len(windings) > 1 else 0.0
        print(
            f"{L},{len(group)},"
            f"{mean(windings):.6g},{var_winding:.6g},"
            f"{mean(item.raw_steps for item in group):.6g},"
            f"{mean(item.lerw_steps for item in group):.6g}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Simulate LERW winding in symmetric or fixed random environments."
    )
    parser.add_argument("--sizes", type=int, nargs="+", default=[16, 32, 64])
    parser.add_argument("--samples", type=int, default=50)
    parser.add_argument(
        "--model",
        choices=["symmetric", "gamma", "exponential", "gamma4"],
        default="gamma",
    )
    parser.add_argument(
        "--k",
        type=float,
        default=1.0,
        help="Gamma shape parameter. Scale is set to 1/k so the mean is 1.",
    )
    parser.add_argument("--seed", type=int, default=20260612)
    parser.add_argument("--max-steps", type=int, default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/debug_lerw.csv"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    k = None if args.model in {"symmetric", "exponential"} else args.k

    if args.samples < 1:
        raise ValueError("--samples must be at least 1")
    for L in args.sizes:
        validate_L(L)

    results = run_experiment(
        sizes=args.sizes,
        samples=args.samples,
        model=args.model,
        k=k,
        seed=args.seed,
        max_steps=args.max_steps,
    )
    write_csv(args.output, results)
    print_summary(results)
    print(f"\nwrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()
