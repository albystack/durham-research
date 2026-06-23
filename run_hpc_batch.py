#!/usr/bin/env python3
"""Run one Slurm-array LERW batch from a CSV configuration row."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import random
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

from config import stable_seed, validate_model
from simulation import SiteDirectionalEnvironment, SiteIIDEnvironment, loop_erased_walk_until_boundary, winding


SCHEMA_VERSION = "hpc_batch_v2"
DEFAULT_ENVIRONMENT_MODEL = "site_iid"


@dataclass(frozen=True)
class BatchConfig:
    """One row from the HPC experiment configuration table."""

    task_id: int
    distribution: str
    distribution_params: dict[str, Any]
    L: int
    batch_id: int
    num_environments: int
    walks_per_environment: int
    base_seed: int
    environment_model: str = DEFAULT_ENVIRONMENT_MODEL

    @property
    def params_json(self) -> str:
        return canonical_json(self.distribution_params)


PARAMETER_ALIASES: dict[str, tuple[str, str | None, tuple[str, ...]]] = {
    "baseline": ("symmetric", None, ()),
    "symmetric": ("symmetric", None, ()),
    "gamma": ("gamma_edges", "shape", ("shape", "k", "parameter")),
    "gamma_edges": ("gamma_edges", "shape", ("shape", "k", "parameter")),
    "exponential": ("exp_edges", None, ()),
    "exp": ("exp_edges", None, ()),
    "exp_edges": ("exp_edges", None, ()),
    "lognormal": ("lognormal_edges", "sigma", ("sigma", "parameter")),
    "lognormal_edges": ("lognormal_edges", "sigma", ("sigma", "parameter")),
    "pareto": ("pareto_edges", "alpha", ("alpha", "parameter")),
    "pareto_edges": ("pareto_edges", "alpha", ("alpha", "parameter")),
    "uniform": ("uniform_edges", "a", ("a", "parameter")),
    "uniform_edges": ("uniform_edges", "a", ("a", "parameter")),
    "beta": ("beta_edges", "a", ("a", "parameter")),
    "beta_edges": ("beta_edges", "a", ("a", "parameter")),
    "weibull": ("weibull_edges", "shape", ("shape", "k", "parameter")),
    "weibull_edges": ("weibull_edges", "shape", ("shape", "k", "parameter")),
    "inverse_gamma": ("inverse_gamma_edges", "alpha", ("alpha", "parameter")),
    "inverse_gamma_edges": ("inverse_gamma_edges", "alpha", ("alpha", "parameter")),
    "bernoulli": ("bernoulli_edges", "a", ("a", "parameter")),
    "bernoulli_edges": ("bernoulli_edges", "a", ("a", "parameter")),
    "triangular": ("triangular_edges", "a", ("a", "parameter")),
    "triangular_edges": ("triangular_edges", "a", ("a", "parameter")),
}


def canonical_json(value: dict[str, Any]) -> str:
    """Stable compact JSON for seeds, output metadata, and CSV cells."""

    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def parse_params(text: str | None) -> dict[str, Any]:
    """Parse the JSON object stored in the distribution_params CSV column."""

    if text is None or text.strip() == "":
        return {}
    value = json.loads(text)
    if not isinstance(value, dict):
        raise ValueError("distribution_params must be a JSON object")
    return value


def _extract_parameter(
    distribution: str,
    params: dict[str, Any],
    canonical_name: str | None,
    aliases: tuple[str, ...],
) -> float | None:
    if canonical_name is None:
        if params:
            raise ValueError(f"{distribution} does not take distribution_params")
        return None
    for name in aliases:
        if name in params:
            return float(params[name])
    raise ValueError(f"{distribution} needs parameter `{canonical_name}` in distribution_params")


def distribution_to_model(distribution: str, params: dict[str, Any]) -> tuple[str, float | None]:
    """Map public HPC distribution names to the existing model catalogue."""

    key = distribution.strip()
    if key not in PARAMETER_ALIASES:
        available = ", ".join(sorted(PARAMETER_ALIASES))
        raise ValueError(f"unknown distribution {distribution!r}; available: {available}")
    model, canonical_name, aliases = PARAMETER_ALIASES[key]
    parameter = _extract_parameter(key, params, canonical_name, aliases)
    validate_model(model, parameter)
    return model, parameter


def load_config_row(config_path: Path, task_id: int) -> BatchConfig:
    """Read the requested task row from a CSV configuration table."""

    with config_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if int(row["task_id"]) != task_id:
                continue
            return BatchConfig(
                task_id=int(row["task_id"]),
                distribution=str(row["distribution"]),
                distribution_params=parse_params(row.get("distribution_params")),
                L=int(row["L"]),
                batch_id=int(row["batch_id"]),
                num_environments=int(row["num_environments"]),
                walks_per_environment=int(row["walks_per_environment"]),
                base_seed=int(row["base_seed"]),
                environment_model=row.get("environment_model") or DEFAULT_ENVIRONMENT_MODEL,
            )
    raise ValueError(f"task_id {task_id} not found in {config_path}")


@lru_cache(maxsize=1)
def code_version() -> str:
    """Best-effort git commit identifier for result metadata."""

    try:
        completed = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return "unknown"
    return completed.stdout.strip() or "unknown"


def result_path(output_dir: Path, task: BatchConfig) -> Path:
    """Return the per-task output path."""

    distribution_to_model(task.distribution, task.distribution_params)
    distribution_dir = safe_distribution_name(task.distribution, task.distribution_params)
    return output_dir / safe_environment_name(task.environment_model) / distribution_dir / f"L_{task.L:04d}" / f"batch_{task.batch_id:04d}.csv"


def _safe_token(value: object) -> str:
    if isinstance(value, float):
        text = f"{value:g}"
    else:
        text = str(value)
    text = text.replace(".", "p").replace("-", "m")
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_")


def safe_distribution_name(distribution: str, params: dict[str, Any]) -> str:
    """Directory-safe public distribution label."""

    base = _safe_token(distribution)
    if not params:
        return base
    suffix = "_".join(f"{_safe_token(key)}_{_safe_token(params[key])}" for key in sorted(params))
    return f"{base}_{suffix}"


def safe_environment_name(environment_model: str) -> str:
    """Directory-safe environment model label."""

    return _safe_token(environment_model)


def expected_observations(task: BatchConfig) -> int:
    return task.num_environments * task.walks_per_environment


def completed_result(path: Path, expected_rows: int) -> bool:
    """True when an existing batch output has all expected successful rows."""

    if not path.exists():
        return False
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != expected_rows:
        return False
    return all(row.get("status") == "ok" for row in rows)


def write_csv_atomic(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a CSV through a sibling temporary file and atomically replace."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    fieldnames = list(rows[0].keys()) if rows else []
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _seed_for(task: BatchConfig, *parts: object) -> int:
    return stable_seed(
        task.base_seed,
        task.environment_model,
        task.distribution,
        task.params_json,
        task.L,
        task.batch_id,
        *parts,
    )


def _empty_row(
    task: BatchConfig,
    model: str,
    parameter: float | None,
    environment_id: int,
    walk_id: int,
    environment_seed: int,
    walk_seed: int,
    started_at: str,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "code_version": code_version(),
        "python_version": platform.python_version(),
        "task_id": task.task_id,
        "batch_id": task.batch_id,
        "distribution": task.distribution,
        "distribution_params": task.params_json,
        "environment_model": task.environment_model,
        "model": model,
        "model_parameter": "" if parameter is None else parameter,
        "L": task.L,
        "environment_id": environment_id,
        "walk_id": walk_id,
        "environment_seed": environment_seed,
        "walk_seed": walk_seed,
        "winding": "",
        "loop_erased_path_length": "",
        "raw_walk_length": "",
        "exit_location": "",
        "exit_x": "",
        "exit_y": "",
        "runtime": "",
        "sampled_site_count": "",
        "status": "",
        "error_type": "",
        "error_message": "",
        "started_at_utc": started_at,
        "finished_at_utc": "",
    }


def run_batch(task: BatchConfig, max_steps: int | None = None) -> list[dict[str, Any]]:
    """Run all observations in one task row."""

    model, parameter = distribution_to_model(task.distribution, task.distribution_params)
    rows: list[dict[str, Any]] = []

    for environment_offset in range(task.num_environments):
        environment_id = task.batch_id * task.num_environments + environment_offset
        environment_seed = _seed_for(task, "environment", environment_id)
        if task.environment_model == "site_iid":
            environment = SiteIIDEnvironment(environment_seed, model, parameter)
        elif task.environment_model == "site_drift_ne_sw":
            environment = SiteDirectionalEnvironment(environment_seed, model, parameter)
        else:
            raise ValueError(
                f"unknown environment_model {task.environment_model!r}; "
                "use site_iid or site_drift_ne_sw"
            )

        for walk_id in range(task.walks_per_environment):
            walk_seed = _seed_for(task, "walk", environment_id, walk_id)
            walk_rng = random.Random(walk_seed)
            started_at = datetime.now(timezone.utc).isoformat()
            row = _empty_row(
                task,
                model,
                parameter,
                environment_id,
                walk_id,
                environment_seed,
                walk_seed,
                started_at,
            )
            start = time.perf_counter()
            try:
                erased_path, raw_steps = loop_erased_walk_until_boundary(
                    task.L,
                    environment,
                    walk_rng,
                    max_steps,
                )
                exit_x, exit_y = erased_path[-1]
                row.update(
                    {
                        "winding": winding(erased_path),
                        "loop_erased_path_length": len(erased_path) - 1,
                        "raw_walk_length": raw_steps,
                        "exit_location": f"{exit_x},{exit_y}",
                        "exit_x": exit_x,
                        "exit_y": exit_y,
                        "sampled_site_count": environment.sampled_site_count,
                        "status": "ok",
                    }
                )
            except Exception as exc:  # noqa: BLE001 - failures are recorded per observation.
                row.update(
                    {
                        "status": "failed",
                        "error_type": type(exc).__name__,
                        "error_message": str(exc),
                        "sampled_site_count": environment.sampled_site_count,
                    }
                )
            finally:
                row["runtime"] = time.perf_counter() - start
                row["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
                rows.append(row)

    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run one HPC LERW batch task.")
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--max-steps", type=int, default=None)
    parser.add_argument("--force", action="store_true", help="overwrite even completed result files")
    parser.add_argument(
        "--rerun-failed",
        action="store_true",
        help="overwrite an existing incomplete or failed result file",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    task = load_config_row(args.config, args.task_id)
    path = result_path(args.output_dir, task)
    expected = expected_observations(task)

    if path.exists() and not args.force:
        if completed_result(path, expected):
            print(f"Task {task.task_id} already complete: {path}")
            return 0
        if not args.rerun_failed:
            print(
                f"Existing output is not complete: {path}\n"
                "Use --rerun-failed to replace it, or --force to overwrite any result.",
                file=sys.stderr,
            )
            return 3

    start = time.perf_counter()
    rows = run_batch(task, max_steps=args.max_steps)
    write_csv_atomic(path, rows)
    seconds = time.perf_counter() - start
    failures = sum(1 for row in rows if row["status"] != "ok")
    print(
        f"Task {task.task_id} wrote {len(rows)} rows to {path} "
        f"in {seconds:.2f}s; failures={failures}",
        flush=True,
    )
    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
