#!/usr/bin/env python3
"""One-file runner for the LERW random-environment project.

The goal of this file is to be readable first and clever second. It contains
the simulation, analysis, and simple SVG plotting code so the project can be
run from one entry point:

    python3 main.py list
    python3 main.py run smoke
    python3 main.py analyse results/smoke_balanced_distributions.csv

Experiment setup lives in `experiments.py`. Tests live in `test_main.py`.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, variance
from typing import Iterable

from experiments import EXPERIMENTS


# A lattice point is always represented as (x, y).
Point = tuple[int, int]


# Direction codes go counterclockwise. This makes winding easy:
# E -> N is +1, E -> S is -1.
DIRECTION_CODES: dict[Point, int] = {
    (1, 0): 0,   # E
    (0, 1): 1,   # N
    (-1, 0): 2,  # W
    (0, -1): 3,  # S
}


# The transition weights are stored in this order everywhere in the code.
STEP_VECTORS: tuple[tuple[str, Point], ...] = (
    ("N", (0, 1)),
    ("E", (1, 0)),
    ("S", (0, -1)),
    ("W", (-1, 0)),
)


# Model groups. Keeping these sets explicit makes validation easier to read.
TWO_WEIGHT_MODELS = {"gamma", "exponential"}
FOUR_WEIGHT_MODELS = {"gamma4", "lognormal4", "pareto4", "uniform4"}
PAIR_BALANCED_MODELS = {
    "gamma_balanced",
    "lognormal_balanced",
    "pareto_balanced",
    "uniform_balanced",
}
PARAMETERISED_MODELS = {
    "gamma",
    *FOUR_WEIGHT_MODELS,
    *PAIR_BALANCED_MODELS,
}
SUPPORTED_MODELS = {
    "symmetric",
    *TWO_WEIGHT_MODELS,
    *FOUR_WEIGHT_MODELS,
    *PAIR_BALANCED_MODELS,
}


# Labels used in tables and plot legends.
MODEL_PARAMETER_LABELS = {
    "gamma": "k",
    "gamma4": "k",
    "lognormal4": "sigma",
    "pareto4": "alpha",
    "uniform4": "a",
    "gamma_balanced": "k",
    "lognormal_balanced": "sigma",
    "pareto_balanced": "alpha",
    "uniform_balanced": "a",
}


# Stable plot colors. Add a model here if you add a new model family.
COLORS = {
    "symmetric": "#111827",
    "gamma": "#dc2626",
    "gamma4": "#2563eb",
    "lognormal4": "#9333ea",
    "pareto4": "#ea580c",
    "uniform4": "#059669",
    "gamma_balanced": "#0f766e",
    "lognormal_balanced": "#7c3aed",
    "pareto_balanced": "#c2410c",
    "uniform_balanced": "#047857",
}


@dataclass(frozen=True)
class TrialResult:
    """One completed random-walk trial, exactly as it is written to CSV."""

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


@dataclass(frozen=True)
class GroupStats:
    """Summary statistics for one model/parameter/size group."""

    model: str
    k: str
    L: int
    samples: int
    mean_winding: float
    var_winding: float
    mean_raw_steps: float
    mean_lerw_steps: float
    mean_hit_x: float
    mean_hit_y: float


@dataclass(frozen=True)
class Series:
    """One line on a plot."""

    label: str
    points: list[tuple[float, float]]
    color: str


def model_requires_parameter(model: str) -> bool:
    """Return whether a model needs a distribution parameter."""

    return model in PARAMETERISED_MODELS


def validate_model_parameter(model: str, k: float | None) -> None:
    """Reject unknown models and invalid distribution parameters early."""

    if model not in SUPPORTED_MODELS:
        raise ValueError(f"unknown model: {model}")
    if model in PARAMETERISED_MODELS and (k is None or k <= 0):
        raise ValueError(f"{model} requires a positive parameter")
    if model in {"pareto4", "pareto_balanced"} and k is not None and k <= 1:
        raise ValueError(f"{model} requires alpha > 1 for finite mean")
    if model in {"uniform4", "uniform_balanced"} and k is not None and k >= 1:
        raise ValueError(f"{model} requires 0 < a < 1")


def parse_model_config(config: str) -> tuple[str, float | None]:
    """Parse strings like 'symmetric' or 'gamma_balanced:1'."""

    if ":" in config:
        model, raw_k = config.split(":", 1)
        k = float(raw_k)
    else:
        model = config
        k = None

    if model_requires_parameter(model) and k is None:
        raise ValueError(f"{model!r} needs a parameter, for example {model}:1")
    if not model_requires_parameter(model) and k is not None:
        raise ValueError(f"{model!r} should not include a parameter")

    validate_model_parameter(model, k)
    return model, k


def fmt_k(k: str | float | None) -> str:
    """Format a parameter without noisy trailing decimals."""

    if k in {None, "", "None"}:
        return "None"
    return f"{float(k):g}"


def model_label(model: str, k: str | float | None) -> str:
    """Human-readable model label for reports and figures."""

    if model == "symmetric":
        return "symmetric"
    parameter = MODEL_PARAMETER_LABELS.get(model, "parameter")
    return f"{model}, {parameter}={fmt_k(k)}"


def config_label(model: str, k: float | None) -> str:
    """Filesystem-friendly label for checkpoint CSVs."""

    if k is None:
        return model
    return f"{model}_{MODEL_PARAMETER_LABELS.get(model, 'p')}{fmt_k(k)}"


def stable_config_seed(base_seed: int, model: str, k: float | None) -> int:
    """Make one repeatable seed per model config from the experiment seed."""

    label = f"{base_seed}:{model}:{k}"
    digest = hashlib.sha256(label.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")


class RandomEnvironment:
    """Fixed site-dependent transition environment.

    The environment is sampled lazily: a site's weights are generated the first
    time the walk visits that site, then reused forever. Because site weights
    are independent, this is equivalent to sampling the whole box in advance.
    """

    def __init__(self, rng: random.Random, model: str, k: float | None = None):
        validate_model_parameter(model, k)
        self.rng = rng
        self.model = model
        self.k = k
        self._weights: dict[Point, tuple[float, ...]] = {}

    @property
    def sampled_site_count(self) -> int:
        """Number of environment sites sampled so far."""

        return len(self._weights)

    def sample_two_weights(self) -> tuple[float, float]:
        """Sample the u, v weights in the original drift-diagnostic model."""

        if self.model == "exponential":
            return (self.rng.expovariate(1.0), self.rng.expovariate(1.0))

        assert self.k is not None
        return (
            self.rng.gammavariate(self.k, 1.0 / self.k),
            self.rng.gammavariate(self.k, 1.0 / self.k),
        )

    def sample_four_weights(self) -> tuple[float, float, float, float]:
        """Sample independent N, E, S, W weights for comparison models."""

        assert self.k is not None

        if self.model == "gamma4":
            return tuple(
                self.rng.gammavariate(self.k, 1.0 / self.k)
                for _ in STEP_VECTORS
            )
        if self.model == "lognormal4":
            sigma = self.k
            mu = -0.5 * sigma * sigma
            return tuple(self.rng.lognormvariate(mu, sigma) for _ in STEP_VECTORS)
        if self.model == "pareto4":
            alpha = self.k
            scale = (alpha - 1.0) / alpha
            return tuple(scale * self.rng.paretovariate(alpha) for _ in STEP_VECTORS)
        if self.model == "uniform4":
            half_width = self.k
            return tuple(
                self.rng.uniform(1.0 - half_width, 1.0 + half_width)
                for _ in STEP_VECTORS
            )

        raise ValueError(f"{self.model} is not a four-weight model")

    def sample_balanced_pair_weights(self) -> tuple[float, float]:
        """Sample vertical and horizontal weights for exact local balance."""

        assert self.k is not None

        if self.model == "gamma_balanced":
            return (
                self.rng.gammavariate(self.k, 1.0 / self.k),
                self.rng.gammavariate(self.k, 1.0 / self.k),
            )
        if self.model == "lognormal_balanced":
            sigma = self.k
            mu = -0.5 * sigma * sigma
            return (
                self.rng.lognormvariate(mu, sigma),
                self.rng.lognormvariate(mu, sigma),
            )
        if self.model == "pareto_balanced":
            alpha = self.k
            scale = (alpha - 1.0) / alpha
            return (
                scale * self.rng.paretovariate(alpha),
                scale * self.rng.paretovariate(alpha),
            )
        if self.model == "uniform_balanced":
            half_width = self.k
            return (
                self.rng.uniform(1.0 - half_width, 1.0 + half_width),
                self.rng.uniform(1.0 - half_width, 1.0 + half_width),
            )

        raise ValueError(f"{self.model} is not an opposite-pair model")

    def weights_at(self, point: Point) -> tuple[float, ...]:
        """Return the fixed random weights for one site."""

        if self.model == "symmetric":
            return (1.0, 1.0)

        if point not in self._weights:
            if self.model in FOUR_WEIGHT_MODELS:
                weights = self.sample_four_weights()
            elif self.model in PAIR_BALANCED_MODELS:
                weights = self.sample_balanced_pair_weights()
            else:
                weights = self.sample_two_weights()
            self._weights[point] = weights

        return self._weights[point]

    def transition_weights(self, point: Point) -> tuple[float, float, float, float]:
        """Return unnormalised weights in N, E, S, W order."""

        if self.model == "symmetric":
            return (1.0, 1.0, 1.0, 1.0)
        if self.model in FOUR_WEIGHT_MODELS:
            w_n, w_e, w_s, w_w = self.weights_at(point)
            return (w_n, w_e, w_s, w_w)
        if self.model in PAIR_BALANCED_MODELS:
            vertical, horizontal = self.weights_at(point)
            return (vertical, horizontal, vertical, horizontal)

        # Original two-weight model: north/east are fixed at one.
        u, v = self.weights_at(point)
        return (1.0, 1.0, u, v)

    def transition_probabilities(self, point: Point) -> tuple[float, float, float, float]:
        """Return normalised transition probabilities in N, E, S, W order."""

        weights = self.transition_weights(point)
        total = sum(weights)
        return tuple(weight / total for weight in weights)


def validate_L(L: int) -> None:
    """Validate the square-box half-width L."""

    if L < 1:
        raise ValueError("L must be at least 1")


def is_boundary(point: Point, L: int) -> bool:
    """Return whether a point lies on the boundary of [-L, L]^2."""

    x, y = point
    return abs(x) == L or abs(y) == L


def choose_step(point: Point, env: RandomEnvironment) -> Point:
    """Sample the next lattice point from the local transition weights."""

    weights = env.transition_weights(point)
    total = sum(weights)
    threshold = env.rng.random() * total

    cumulative = 0.0
    for (_, step), weight in zip(STEP_VECTORS, weights):
        cumulative += weight
        if threshold <= cumulative:
            return (point[0] + step[0], point[1] + step[1])

    # Floating-point fallback for thresholds extremely close to `total`.
    step = STEP_VECTORS[-1][1]
    return (point[0] + step[0], point[1] + step[1])


def random_walk_to_boundary(
    L: int,
    env: RandomEnvironment,
    max_steps: int | None = None,
) -> list[Point]:
    """Run a weighted walk from the origin until it first hits the boundary."""

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
            # We returned to an old point, so everything after it is erased.
            keep_until = index[point]
            for removed in erased[keep_until + 1 :]:
                del index[removed]
            erased = erased[: keep_until + 1]
        else:
            index[point] = len(erased)
            erased.append(point)

    return erased


def direction_code(start: Point, end: Point) -> int:
    """Encode a nearest-neighbour step as E=0, N=1, W=2, S=3."""

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
    """Run one seeded trial and return all measured quantities."""

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
    """Expand one base seed into independent trial seeds."""

    rng = random.Random(base_seed)
    return [rng.randrange(0, 2**63) for _ in range(count)]


def run_model_config(
    sizes: list[int],
    samples: int,
    model: str,
    k: float | None,
    seed: int,
    max_steps: int | None,
) -> list[TrialResult]:
    """Run every requested size/sample pair for one model configuration."""

    total_trials = len(sizes) * samples
    seeds = sample_seeds(seed, total_trials)
    results: list[TrialResult] = []
    seed_index = 0

    for L in sizes:
        validate_L(L)
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


def write_results_csv(path: Path, results: list[TrialResult]) -> None:
    """Write trial-level simulation rows."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(TrialResult.__dataclass_fields__.keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    """Read one or more result CSVs."""

    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def grouped_stats(rows: list[dict[str, str]]) -> list[GroupStats]:
    """Aggregate trial rows by model, parameter, and box size."""

    groups: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row["model"], row["k"] or "None", int(row["L"]))
        groups[key].append(row)

    stats: list[GroupStats] = []
    for (model, k, L), group in sorted(groups.items(), key=lambda item: item[0]):
        windings = [int(row["winding"]) for row in group]
        stats.append(
            GroupStats(
                model=model,
                k=k,
                L=L,
                samples=len(group),
                mean_winding=mean(windings),
                var_winding=variance(windings) if len(windings) > 1 else 0.0,
                mean_raw_steps=mean(int(row["raw_steps"]) for row in group),
                mean_lerw_steps=mean(int(row["lerw_steps"]) for row in group),
                mean_hit_x=mean(int(row["hit_x"]) for row in group),
                mean_hit_y=mean(int(row["hit_y"]) for row in group),
            )
        )
    return stats


def linear_fit(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Return slope, intercept, and R^2 for y = slope*x + intercept."""

    if len(xs) != len(ys) or len(xs) < 2:
        raise ValueError("need at least two points for a linear fit")

    x_bar = mean(xs)
    y_bar = mean(ys)
    ss_xx = sum((x - x_bar) ** 2 for x in xs)
    if ss_xx == 0:
        raise ValueError("x values are all identical")

    slope = sum((x - x_bar) * (y - y_bar) for x, y in zip(xs, ys)) / ss_xx
    intercept = y_bar - slope * x_bar

    fitted = [slope * x + intercept for x in xs]
    ss_res = sum((y - y_hat) ** 2 for y, y_hat in zip(ys, fitted))
    ss_tot = sum((y - y_bar) ** 2 for y in ys)
    r_squared = 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot
    return slope, intercept, r_squared


def write_group_stats_csv(path: Path, stats: list[GroupStats]) -> None:
    """Write the grouped statistics table used by the report and plots."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "model",
                "k",
                "L",
                "samples",
                "mean_winding",
                "var_winding",
                "mean_raw_steps",
                "mean_lerw_steps",
                "mean_hit_x",
                "mean_hit_y",
                "mean_hit_x_over_L",
                "mean_hit_y_over_L",
            ]
        )
        for item in stats:
            writer.writerow(
                [
                    item.model,
                    item.k,
                    item.L,
                    item.samples,
                    item.mean_winding,
                    item.var_winding,
                    item.mean_raw_steps,
                    item.mean_lerw_steps,
                    item.mean_hit_x,
                    item.mean_hit_y,
                    item.mean_hit_x / item.L,
                    item.mean_hit_y / item.L,
                ]
            )


def write_fit_csv(path: Path, stats: list[GroupStats]) -> list[dict[str, str | float]]:
    """Write descriptive fits of winding variance against log scales."""

    by_model: dict[tuple[str, str], list[GroupStats]] = defaultdict(list)
    for item in stats:
        by_model[(item.model, item.k)].append(item)

    rows: list[dict[str, str | float]] = []
    for (model, k), items in sorted(by_model.items()):
        items = sorted(items, key=lambda item: item.L)
        if len(items) < 3:
            continue
        ys = [item.var_winding for item in items]
        for fit_name, xs in (
            ("logL", [math.log(item.L) for item in items]),
            ("logL_squared", [math.log(item.L) ** 2 for item in items]),
        ):
            slope, intercept, r2 = linear_fit(xs, ys)
            fitted = [slope * x + intercept for x in xs]
            rss = sum((y - y_hat) ** 2 for y, y_hat in zip(ys, fitted))
            rows.append(
                {
                    "model": model,
                    "k": k,
                    "fit": fit_name,
                    "slope": slope,
                    "intercept": intercept,
                    "r_squared": r2,
                    "rss": rss,
                }
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["model", "k", "fit", "slope", "intercept", "r_squared", "rss"],
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def key_color(model: str, k: str | float | None) -> str:
    """Pick a repeatable color for one model line."""

    if model in COLORS:
        return COLORS[model]
    return "#2563eb"


def stat_lookup(stats: list[GroupStats]) -> dict[tuple[str, str, int], GroupStats]:
    """Index stats for quick plot lookups."""

    return {(item.model, item.k, item.L): item for item in stats}


def preferred_selected(stats: list[GroupStats], limit: int = 6) -> list[tuple[str, str]]:
    """Choose the most useful model lines for multi-model figures."""

    available = {(item.model, item.k) for item in stats}
    preferred = [
        ("symmetric", "None"),
        ("gamma_balanced", "1.0"),
        ("gamma_balanced", "0.5"),
        ("lognormal_balanced", "0.5"),
        ("lognormal_balanced", "1.0"),
        ("pareto_balanced", "2.0"),
        ("pareto_balanced", "3.0"),
        ("uniform_balanced", "0.5"),
        ("gamma4", "1.0"),
        ("gamma4", "0.5"),
        ("gamma", "1.0"),
    ]

    selected = [item for item in preferred if item in available]
    for item in sorted(available):
        if len(selected) >= limit:
            break
        if item not in selected:
            selected.append(item)
    return selected[:limit]


def selected_series(
    stats: list[GroupStats],
    value_name: str,
    transform_x,
    selected: list[tuple[str, str]],
) -> list[Series]:
    """Build plot lines from grouped statistics."""

    lookup = stat_lookup(stats)
    sizes = sorted({item.L for item in stats})
    series_list: list[Series] = []

    for model, k in selected:
        points = []
        for L in sizes:
            item = lookup.get((model, k, L))
            if item is not None:
                points.append((transform_x(L), getattr(item, value_name)))
        if points:
            series_list.append(
                Series(
                    label=model_label(model, k),
                    points=points,
                    color=key_color(model, k),
                )
            )
    return series_list


def nice_ticks(v_min: float, v_max: float, count: int = 6) -> list[float]:
    """Choose simple numeric tick positions for SVG plots."""

    if not math.isfinite(v_min) or not math.isfinite(v_max):
        return []
    if v_min == v_max:
        return [v_min]

    span = v_max - v_min
    rough = span / max(count - 1, 1)
    magnitude = 10 ** math.floor(math.log10(abs(rough)))
    fraction = rough / magnitude
    if fraction <= 1:
        nice_fraction = 1
    elif fraction <= 2:
        nice_fraction = 2
    elif fraction <= 5:
        nice_fraction = 5
    else:
        nice_fraction = 10

    step = nice_fraction * magnitude
    start = math.floor(v_min / step) * step
    end = math.ceil(v_max / step) * step
    ticks = []
    value = start
    while value <= end + step * 0.5:
        ticks.append(0.0 if abs(value) < 1e-12 else value)
        value += step
    return ticks


def format_tick(value: float) -> str:
    """Format plot tick labels compactly."""

    if abs(value) >= 1000:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2f}".rstrip("0").rstrip(".")


def svg_text(
    x: float,
    y: float,
    text: str,
    size: int = 13,
    anchor: str = "middle",
    weight: str = "400",
    fill: str = "#111827",
) -> str:
    """Create escaped SVG text."""

    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
        f'fill="{fill}">{html.escape(text)}</text>'
    )


def svg_doc(width: int, height: int, body: list[str]) -> str:
    """Wrap SVG body fragments in a complete document."""

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        + "\n".join(body)
        + "\n</svg>\n"
    )


def write_svg(path: Path, width: int, height: int, body: list[str]) -> None:
    """Write one SVG file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(svg_doc(width, height, body), encoding="utf-8")


def line_plot(
    path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    series_list: list[Series],
    x_ticks: list[tuple[float, str]] | None = None,
    width: int = 1050,
    height: int = 760,
) -> None:
    """Write a simple multi-series SVG line plot."""

    if not series_list:
        return

    left, right, top, bottom = 82, 285, 58, 78
    plot_w = width - left - right
    plot_h = height - top - bottom

    xs = [x for series in series_list for x, _ in series.points]
    ys = [y for series in series_list for _, y in series.points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if x_min == x_max:
        x_min -= 1
        x_max += 1
    y_pad = (y_max - y_min) * 0.08 if y_max != y_min else 1.0
    y_min -= y_pad
    y_max += y_pad

    def px(x: float) -> float:
        return left + (x - x_min) / (x_max - x_min) * plot_w

    def py(y: float) -> float:
        return top + (y_max - y) / (y_max - y_min) * plot_h

    body: list[str] = []
    body.append(svg_text(width / 2, 30, title, size=19, weight="700"))
    body.append(
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" '
        'fill="#f9fafb" stroke="#d1d5db"/>'
    )

    for tick in nice_ticks(y_min, y_max):
        y = py(tick)
        body.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" '
            f'y2="{y:.2f}" stroke="#e5e7eb"/>'
        )
        body.append(svg_text(left - 10, y + 4, format_tick(tick), size=11, anchor="end"))

    if x_ticks is None:
        x_ticks = [(tick, format_tick(tick)) for tick in nice_ticks(x_min, x_max)]
    for tick, label in x_ticks:
        x = px(tick)
        body.append(
            f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" '
            f'y2="{top + plot_h}" stroke="#eef2f7"/>'
        )
        body.append(svg_text(x, top + plot_h + 21, label, size=11))

    body.append(
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" '
        f'y2="{top + plot_h}" stroke="#111827"/>'
    )
    body.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" '
        f'y2="{top + plot_h}" stroke="#111827"/>'
    )
    body.append(svg_text(left + plot_w / 2, height - 25, xlabel))
    body.append(
        f'<text x="22" y="{top + plot_h / 2:.2f}" font-family="Arial, sans-serif" '
        f'font-size="13" text-anchor="middle" fill="#111827" '
        f'transform="rotate(-90 22 {top + plot_h / 2:.2f})">{html.escape(ylabel)}</text>'
    )

    for series in series_list:
        points = " ".join(f"{px(x):.2f},{py(y):.2f}" for x, y in series.points)
        body.append(
            f'<polyline points="{points}" fill="none" stroke="{series.color}" '
            'stroke-width="2.4"/>'
        )
        for x, y in series.points:
            body.append(
                f'<circle cx="{px(x):.2f}" cy="{py(y):.2f}" r="4" '
                f'fill="{series.color}" stroke="#ffffff" stroke-width="1"/>'
            )

    legend_x = width - right + 28
    legend_y = top + 12
    for index, series in enumerate(series_list):
        y = legend_y + index * 24
        body.append(
            f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 22}" '
            f'y2="{y}" stroke="{series.color}" stroke-width="3"/>'
        )
        body.append(svg_text(legend_x + 30, y + 4, series.label, size=12, anchor="start"))

    write_svg(path, width, height, body)


def boundary_hit_plot(
    path: Path,
    rows: list[dict[str, str]],
    selected: list[tuple[str, str]],
    L: int,
    max_points: int = 650,
) -> None:
    """Write a small-multiple scatter plot of boundary-hit locations."""

    width, height = 960, 760
    margin = 42
    panel_gap = 24
    panel_w = (width - 2 * margin - panel_gap * (len(selected) - 1)) / len(selected)
    panel_size = min(panel_w, height - 118)
    top = 62
    rng = random.Random(123)
    body: list[str] = [svg_text(width / 2, 31, f"Boundary Hit Locations, L={L}", size=19, weight="700")]

    for idx, (model, k) in enumerate(selected):
        x0 = margin + idx * (panel_w + panel_gap)
        y0 = top
        group = [
            row
            for row in rows
            if row["model"] == model and (row["k"] or "None") == k and int(row["L"]) == L
        ]
        if len(group) > max_points:
            group = rng.sample(group, max_points)

        body.append(
            f'<rect x="{x0:.2f}" y="{y0}" width="{panel_size:.2f}" '
            f'height="{panel_size:.2f}" fill="#f9fafb" stroke="#d1d5db"/>'
        )
        body.append(svg_text(x0 + panel_size / 2, y0 - 14, model_label(model, k), size=13, weight="700"))

        def px(value: float) -> float:
            return x0 + (value / L + 1) / 2 * panel_size

        def py(value: float) -> float:
            return y0 + (1 - (value / L + 1) / 2) * panel_size

        body.append(f'<line x1="{px(0):.2f}" y1="{y0}" x2="{px(0):.2f}" y2="{y0 + panel_size}" stroke="#e5e7eb"/>')
        body.append(f'<line x1="{x0:.2f}" y1="{py(0):.2f}" x2="{x0 + panel_size:.2f}" y2="{py(0):.2f}" stroke="#e5e7eb"/>')

        for row in group:
            x = int(row["hit_x"])
            y = int(row["hit_y"])
            body.append(f'<circle cx="{px(x):.2f}" cy="{py(y):.2f}" r="2.2" fill="{key_color(model, k)}" opacity="0.55"/>')

        if group:
            mean_x = mean(int(row["hit_x"]) for row in group)
            mean_y = mean(int(row["hit_y"]) for row in group)
            body.append(f'<circle cx="{px(mean_x):.2f}" cy="{py(mean_y):.2f}" r="5.2" fill="#111827" stroke="#ffffff" stroke-width="1.2"/>')
            body.append(svg_text(x0 + panel_size / 2, y0 + panel_size + 24, f"mean=({mean_x / L:.2f}, {mean_y / L:.2f})", size=12))

    write_svg(path, width, height, body)


def write_analysis_report(
    path: Path,
    source_csv: Path,
    stats: list[GroupStats],
    fit_rows: list[dict[str, str | float]],
    figure_names: list[str],
) -> None:
    """Write a readable markdown report next to the generated figures."""

    largest_L = max(item.L for item in stats)
    stats_by_key = stat_lookup(stats)
    selected = preferred_selected(stats)

    lines: list[str] = [
        "# Numerical Analysis Dossier",
        "",
        f"Source data: `{source_csv}`",
        "",
        "## Dataset",
        "",
        f"- Models: {', '.join(sorted({item.model for item in stats}))}",
        f"- Box sizes: {', '.join(str(size) for size in sorted({item.L for item in stats}))}",
        f"- Total grouped cells: {len(stats)}",
        "",
        "## Main Diagnostics",
        "",
    ]

    for model, k in selected:
        item = stats_by_key.get((model, k, largest_L))
        if item is None:
            continue
        lines.append(
            f"- `{model_label(model, k)}` at L={largest_L}: "
            f"Var(W)={item.var_winding:.3f}, mean W={item.mean_winding:.3f}, "
            f"mean hit/L=({item.mean_hit_x / item.L:.3f}, {item.mean_hit_y / item.L:.3f}), "
            f"mean raw steps={item.mean_raw_steps:.1f}."
        )

    lines.extend(
        [
            "",
            "Boundary-hit displacement is the main drift diagnostic. A strong",
            "displacement means the model should be treated cautiously for the",
            "centred random-environment question.",
            "",
            "## Scaling Fits",
            "",
            "These fits are descriptive finite-size diagnostics, not asymptotic proofs.",
            "",
            "| model | k | fit | slope | intercept | R^2 | RSS |",
            "|---|---:|---|---:|---:|---:|---:|",
        ]
    )

    for row in fit_rows:
        lines.append(
            f"| {row['model']} | {row['k']} | {row['fit']} | "
            f"{float(row['slope']):.4g} | {float(row['intercept']):.4g} | "
            f"{float(row['r_squared']):.4f} | {float(row['rss']):.4g} |"
        )

    lines.extend(["", "## Figures", ""])
    for name in figure_names:
        lines.append(f"![{name}](figures/{name})")
        lines.append("")

    lines.extend(
        [
            "## Reading The Results",
            "",
            "- Symmetric is the calibration case.",
            "- Two-weight gamma is mainly a drift diagnostic.",
            "- Opposite-pair balanced models are the main no-local-drift candidates.",
            "- Larger sizes and more samples are needed before making asymptotic claims.",
        ]
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def analyse_results(csv_path: Path, output_dir: Path) -> None:
    """Create tables, figures, and a markdown dossier from one result CSV."""

    rows = read_rows([csv_path])
    stats = grouped_stats(rows)
    if not stats:
        raise ValueError(f"no rows found in {csv_path}")

    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    write_group_stats_csv(tables_dir / "group_stats.csv", stats)
    fit_rows = write_fit_csv(tables_dir / "scaling_fits.csv", stats)

    sizes = sorted({item.L for item in stats})
    x_ticks_log = [(math.log(L), str(L)) for L in sizes]
    x_ticks_log2 = [(math.log(L) ** 2, str(L)) for L in sizes]
    selected = preferred_selected(stats)
    figure_names: list[str] = []

    plots = [
        (
            "variance_vs_logL.svg",
            "Winding Variance vs log L",
            "box size L (shown at log L positions)",
            "sample variance of winding",
            selected_series(stats, "var_winding", math.log, selected),
            x_ticks_log,
        ),
        (
            "variance_vs_logL2.svg",
            "Winding Variance vs (log L)^2",
            "box size L (shown at (log L)^2 positions)",
            "sample variance of winding",
            selected_series(stats, "var_winding", lambda L: math.log(L) ** 2, selected),
            x_ticks_log2,
        ),
        (
            "mean_hit_x_over_L.svg",
            "Drift Diagnostic: Mean Boundary Hit x/L",
            "box size L (shown at log L positions)",
            "mean hit x / L",
            [
                Series(series.label, [(x, y / math.exp(x)) for x, y in series.points], series.color)
                for series in selected_series(stats, "mean_hit_x", math.log, selected)
            ],
            x_ticks_log,
        ),
        (
            "mean_hit_y_over_L.svg",
            "Drift Diagnostic: Mean Boundary Hit y/L",
            "box size L (shown at log L positions)",
            "mean hit y / L",
            [
                Series(series.label, [(x, y / math.exp(x)) for x, y in series.points], series.color)
                for series in selected_series(stats, "mean_hit_y", math.log, selected)
            ],
            x_ticks_log,
        ),
        (
            "raw_walk_length.svg",
            "Mean Raw Walk Length",
            "box size L (shown at log L positions)",
            "mean raw steps",
            selected_series(stats, "mean_raw_steps", math.log, selected),
            x_ticks_log,
        ),
        (
            "lerw_path_length.svg",
            "Mean Loop-Erased Path Length",
            "box size L (shown at log L positions)",
            "mean LERW steps",
            selected_series(stats, "mean_lerw_steps", math.log, selected),
            x_ticks_log,
        ),
    ]

    for filename, title, xlabel, ylabel, series_list, ticks in plots:
        if series_list:
            line_plot(figures_dir / filename, title, xlabel, ylabel, series_list, ticks)
            figure_names.append(filename)

    mid_L = 128 if 128 in sizes else max(sizes)
    boundary_selected = selected[: min(3, len(selected))]
    boundary_hit_plot(figures_dir / f"boundary_hits_L{mid_L}.svg", rows, boundary_selected, mid_L)
    figure_names.append(f"boundary_hits_L{mid_L}.svg")

    write_analysis_report(
        output_dir / "analysis_dossier.md",
        csv_path,
        stats,
        fit_rows,
        figure_names,
    )

    print(f"wrote tables to {tables_dir}")
    print(f"wrote {len(figure_names)} figures to {figures_dir}")
    print(f"wrote report to {output_dir / 'analysis_dossier.md'}")


def run_experiment(name: str, analyse: bool) -> Path:
    """Run one preset from experiments.py and optionally analyse it."""

    if name not in EXPERIMENTS:
        available = ", ".join(sorted(EXPERIMENTS))
        raise ValueError(f"unknown experiment {name!r}; available: {available}")

    experiment = EXPERIMENTS[name]
    sizes = [int(size) for size in experiment["sizes"]]
    samples = int(experiment["samples"])
    seed = int(experiment["seed"])
    output = Path(str(experiment["output"]))
    max_steps = experiment.get("max_steps")
    max_steps = None if max_steps is None else int(max_steps)
    checkpoint_per_config = bool(experiment.get("checkpoint_per_config", False))
    model_configs = [parse_model_config(str(item)) for item in experiment["model_configs"]]

    if samples < 1:
        raise ValueError("samples must be at least 1")

    print(f"running experiment: {name}")
    print(f"sizes={sizes}, samples_per_size={samples}")
    all_results: list[TrialResult] = []

    for model, k in model_configs:
        label = config_label(model, k)
        print(f"\n== {label} ==", flush=True)
        results = run_model_config(
            sizes=sizes,
            samples=samples,
            model=model,
            k=k,
            seed=stable_config_seed(seed, model, k),
            max_steps=max_steps,
        )
        all_results.extend(results)
        print(f"completed {len(results)} trials", flush=True)

        if checkpoint_per_config:
            checkpoint = output.with_name(f"{output.stem}_{label}{output.suffix}")
            write_results_csv(checkpoint, results)
            print(f"checkpoint wrote {checkpoint}", flush=True)

    write_results_csv(output, all_results)
    print(f"\nwrote {len(all_results)} total rows to {output}")

    if analyse:
        default_analysis_dir = Path(f"analysis_{name}")
        analyse_results(output, default_analysis_dir)

    return output


def list_experiments() -> None:
    """Print the available experiment presets."""

    for name in sorted(EXPERIMENTS):
        experiment = EXPERIMENTS[name]
        print(f"{name}: {experiment['description']}")
        print(f"  output: {experiment['output']}")
        print(f"  sizes: {experiment['sizes']}")
        print(f"  samples: {experiment['samples']}")


def parse_args() -> argparse.Namespace:
    """Parse the small command-line interface."""

    parser = argparse.ArgumentParser(description="Run and analyse LERW experiments.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("list", help="show experiment presets from experiments.py")

    run_parser = subparsers.add_parser("run", help="run a named experiment")
    run_parser.add_argument("experiment", choices=sorted(EXPERIMENTS))
    run_parser.add_argument(
        "--analyse",
        action="store_true",
        help="also generate tables, figures, and a dossier after the run",
    )

    analyse_parser = subparsers.add_parser("analyse", help="analyse an existing CSV")
    analyse_parser.add_argument("csv_file", type=Path)
    analyse_parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_output"),
        help="where generated figures, tables, and report should go",
    )

    return parser.parse_args()


def main() -> None:
    """Command-line entry point."""

    args = parse_args()
    if args.command == "list":
        list_experiments()
    elif args.command == "run":
        run_experiment(args.experiment, analyse=args.analyse)
    elif args.command == "analyse":
        analyse_results(args.csv_file, args.output_dir)
    else:
        raise ValueError(f"unknown command: {args.command}")


if __name__ == "__main__":
    main()
