#!/usr/bin/env python3
"""Generate CSV task tables for Slurm LERW batches."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from config import stable_seed


FIELDNAMES = [
    "task_id",
    "environment_model",
    "distribution",
    "distribution_params",
    "L",
    "batch_id",
    "num_environments",
    "walks_per_environment",
    "base_seed",
]


PARAMETER_NAME_BY_DISTRIBUTION = {
    "gamma": "shape",
    "gamma_edges": "shape",
    "lognormal": "sigma",
    "lognormal_edges": "sigma",
    "pareto": "alpha",
    "pareto_edges": "alpha",
    "uniform": "a",
    "uniform_edges": "a",
    "beta": "a",
    "beta_edges": "a",
    "weibull": "shape",
    "weibull_edges": "shape",
    "inverse_gamma": "alpha",
    "inverse_gamma_edges": "alpha",
    "bernoulli": "a",
    "bernoulli_edges": "a",
    "triangular": "a",
    "triangular_edges": "a",
}


@dataclass(frozen=True)
class ConfigPreset:
    """Named distribution/size grids for common Slurm array tables."""

    sizes: list[int]
    distributions: list[str]


CONFIG_PRESETS = {
    "smoke": ConfigPreset(
        sizes=[64, 128],
        distributions=["baseline", "gamma:1.0"],
    ),
    "hpc_full": ConfigPreset(
        sizes=[32, 64, 128, 256, 512, 1024],
        distributions=[
            "baseline",
            "gamma:0.5",
            "gamma:1.0",
            "gamma:2.0",
            "exponential",
            "lognormal:0.5",
            "lognormal:1.0",
            "pareto:2.0",
            "pareto:3.0",
            "uniform:0.7",
            "beta:0.5",
            "weibull:0.7",
            "inverse_gamma:2.2",
            "bernoulli:0.8",
            "triangular:0.8",
        ],
    ),
}


def canonical_json(value: dict[str, Any]) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def parse_distribution_spec(spec: str) -> tuple[str, dict[str, Any]]:
    """Parse distribution specs such as baseline, gamma:1, lognormal:0.6."""

    if ":" not in spec:
        return spec, {}
    distribution, raw_value = spec.split(":", 1)
    if distribution not in PARAMETER_NAME_BY_DISTRIBUTION:
        raise ValueError(f"do not know the parameter name for distribution {distribution!r}")
    return distribution, {PARAMETER_NAME_BY_DISTRIBUTION[distribution]: float(raw_value)}


def generate_rows(
    distributions: list[str],
    sizes: list[int],
    batches: int,
    num_environments: int,
    walks_per_environment: int,
    base_seed: int,
    environment_model: str,
) -> list[dict[str, Any]]:
    """Expand CLI arguments into one independent task row per array task."""

    parsed_distributions = [parse_distribution_spec(spec) for spec in distributions]
    rows: list[dict[str, Any]] = []
    task_id = 0
    for distribution, params in parsed_distributions:
        params_json = canonical_json(params)
        for L in sizes:
            for batch_id in range(batches):
                rows.append(
                    {
                        "task_id": task_id,
                        "environment_model": environment_model,
                        "distribution": distribution,
                        "distribution_params": params_json,
                        "L": L,
                        "batch_id": batch_id,
                        "num_environments": num_environments,
                        "walks_per_environment": walks_per_environment,
                        "base_seed": stable_seed(base_seed, task_id),
                    }
                )
                task_id += 1
    return rows


def write_config(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate an HPC LERW task CSV.")
    parser.add_argument(
        "--preset",
        choices=sorted(CONFIG_PRESETS),
        default=None,
        help="use a named distribution/size grid; counts still come from the other flags",
    )
    parser.add_argument("--output", type=Path, default=Path("configs/hpc_test.csv"))
    parser.add_argument("--sizes", type=int, nargs="+", default=[64, 128])
    parser.add_argument("--distributions", nargs="+", default=["baseline", "gamma:1.0"])
    parser.add_argument("--batches", type=int, default=1)
    parser.add_argument("--num-environments", type=int, default=2)
    parser.add_argument("--walks-per-environment", type=int, default=3)
    parser.add_argument("--base-seed", type=int, default=20260623)
    parser.add_argument(
        "--environment-model",
        choices=["site_iid", "site_drift_ne_sw"],
        default="site_iid",
        help="site_iid is the balanced four-weight model; site_drift_ne_sw is the older N=E=1 model",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.preset is not None:
        preset = CONFIG_PRESETS[args.preset]
        args.sizes = preset.sizes
        args.distributions = preset.distributions
    rows = generate_rows(
        distributions=args.distributions,
        sizes=args.sizes,
        batches=args.batches,
        num_environments=args.num_environments,
        walks_per_environment=args.walks_per_environment,
        base_seed=args.base_seed,
        environment_model=args.environment_model,
    )
    write_config(args.output, rows)
    print(f"Wrote {len(rows)} tasks to {args.output}")
    print(f"Slurm array range: 0-{len(rows) - 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
