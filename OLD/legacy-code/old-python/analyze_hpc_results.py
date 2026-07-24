#!/usr/bin/env python3
"""Combine and analyse Slurm LERW batch outputs."""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev, variance
from typing import Any

from run_hpc_batch import (
    BatchConfig,
    completed_result,
    expected_observations,
    load_config_row,
    result_path,
)


def parameter_key(value: object) -> str:
    if value in ("", None, "None"):
        return "{}"
    return str(value)


def variance_standard_error(sample_variance: float, sample_count: int) -> float:
    if sample_count <= 1:
        return float("nan")
    return math.sqrt(2.0 / (sample_count - 1)) * sample_variance


def read_config_table(path: Path) -> list[BatchConfig]:
    """Read all task rows using the same parser as the batch runner."""

    with path.open(newline="", encoding="utf-8") as handle:
        task_ids = [int(row["task_id"]) for row in csv.DictReader(handle)]
    return [load_config_row(path, task_id) for task_id in task_ids]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def validate_and_collect(
    config_rows: list[BatchConfig],
    results_dir: Path,
    allow_incomplete: bool,
) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    """Validate expected batch files and return available successful rows."""

    validation_rows: list[dict[str, Any]] = []
    raw_rows: list[dict[str, str]] = []
    for task in config_rows:
        path = result_path(results_dir, task)
        expected_rows = expected_observations(task)
        exists = path.exists()
        row_count = 0
        ok_count = 0
        failed_count = 0
        complete = False
        if exists:
            rows = read_csv_rows(path)
            row_count = len(rows)
            ok_count = sum(1 for row in rows if row.get("status") == "ok")
            failed_count = row_count - ok_count
            complete = completed_result(path, expected_rows)
            raw_rows.extend(row for row in rows if row.get("status") == "ok")
        validation_rows.append(
            {
                "task_id": task.task_id,
                "distribution": task.distribution,
                "distribution_params": task.params_json,
                "L": task.L,
                "batch_id": task.batch_id,
                "path": str(path),
                "exists": exists,
                "expected_rows": expected_rows,
                "row_count": row_count,
                "ok_count": ok_count,
                "failed_count": failed_count,
                "complete": complete,
            }
        )

    incomplete = [row for row in validation_rows if not row["complete"]]
    if incomplete and not allow_incomplete:
        missing = ", ".join(str(row["task_id"]) for row in incomplete[:10])
        extra = "" if len(incomplete) <= 10 else f" and {len(incomplete) - 10} more"
        raise RuntimeError(f"incomplete or missing tasks: {missing}{extra}; use --allow-incomplete to analyse available rows")
    return raw_rows, validation_rows


def duplicate_rows(raw_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Find duplicate observation keys."""

    counts: dict[tuple[str, ...], int] = defaultdict(int)
    for row in raw_rows:
        key = (
            row["distribution"],
            row["distribution_params"],
            row["L"],
            row["batch_id"],
            row["environment_id"],
            row["walk_id"],
        )
        counts[key] += 1
    duplicates = []
    for key, count in sorted(counts.items()):
        if count > 1:
            duplicates.append(
                {
                    "distribution": key[0],
                    "distribution_params": key[1],
                    "L": key[2],
                    "batch_id": key[3],
                    "environment_id": key[4],
                    "walk_id": key[5],
                    "count": count,
                }
            )
    return duplicates


def group_raw_rows(raw_rows: list[dict[str, str]]) -> dict[tuple[str, str, int], list[dict[str, str]]]:
    groups: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        groups[(row["distribution"], parameter_key(row["distribution_params"]), int(row["L"]))].append(row)
    return groups


def make_summary(raw_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Compute annealed and quenched variance summaries by distribution and L."""

    summary_rows: list[dict[str, Any]] = []
    for (distribution, params, L), rows in sorted(group_raw_rows(raw_rows).items()):
        windings = [float(row["winding"]) for row in rows]
        raw_lengths = [float(row["raw_walk_length"]) for row in rows if row["raw_walk_length"] != ""]
        lerw_lengths = [
            float(row["loop_erased_path_length"]) for row in rows if row["loop_erased_path_length"] != ""
        ]
        runtime_values = [float(row["runtime"]) for row in rows if row["runtime"] != ""]
        exit_x_values = [float(row["exit_x"]) / L for row in rows if row.get("exit_x", "") != ""]
        exit_y_values = [float(row["exit_y"]) / L for row in rows if row.get("exit_y", "") != ""]
        by_environment: dict[str, list[float]] = defaultdict(list)
        for row in rows:
            by_environment[row["environment_id"]].append(float(row["winding"]))

        env_variances = [variance(values) for values in by_environment.values() if len(values) >= 2]
        env_means = [mean(values) for values in by_environment.values()]
        annealed_variance = variance(windings) if len(windings) >= 2 else float("nan")
        quenched_variance = mean(env_variances) if env_variances else float("nan")
        between_env_variance = variance(env_means) if len(env_means) >= 2 else float("nan")
        if len(env_variances) >= 2:
            quenched_se = stdev(env_variances) / math.sqrt(len(env_variances))
        else:
            quenched_se = float("nan")

        summary_rows.append(
            {
                "distribution": distribution,
                "distribution_params": params,
                "L": L,
                "log_L": math.log(L),
                "log_log_L": math.log(math.log(L)),
                "observations": len(windings),
                "environments": len(by_environment),
                "walks_per_environment_min": min(len(values) for values in by_environment.values()),
                "walks_per_environment_max": max(len(values) for values in by_environment.values()),
                "mean_winding": mean(windings),
                "annealed_variance": annealed_variance,
                "annealed_variance_se": variance_standard_error(annealed_variance, len(windings)),
                "quenched_variance": quenched_variance,
                "quenched_variance_se": quenched_se,
                "environment_mean_variance": between_env_variance,
                "mean_exit_x_over_L": mean(exit_x_values) if exit_x_values else float("nan"),
                "mean_exit_y_over_L": mean(exit_y_values) if exit_y_values else float("nan"),
                "mean_raw_walk_length": mean(raw_lengths) if raw_lengths else float("nan"),
                "mean_loop_erased_path_length": mean(lerw_lengths) if lerw_lengths else float("nan"),
                "mean_runtime": mean(runtime_values) if runtime_values else float("nan"),
            }
        )
    return summary_rows


def _linear_fit(xs: list[float], ys: list[float]) -> dict[str, float]:
    n = len(xs)
    x_bar = mean(xs)
    y_bar = mean(ys)
    sxx = sum((x - x_bar) ** 2 for x in xs)
    if n < 2 or sxx == 0:
        return {
            "intercept": float("nan"),
            "slope": float("nan"),
            "slope_se": float("nan"),
            "sse": float("nan"),
            "r2": float("nan"),
            "aic": float("nan"),
            "bic": float("nan"),
            "residual_mean": float("nan"),
            "residual_std": float("nan"),
            "max_abs_residual": float("nan"),
            "durbin_watson": float("nan"),
        }
    slope = sum((x - x_bar) * (y - y_bar) for x, y in zip(xs, ys)) / sxx
    intercept = y_bar - slope * x_bar
    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    sse = sum(residual * residual for residual in residuals)
    total = sum((y - y_bar) ** 2 for y in ys)
    r2 = float("nan") if total == 0 else 1.0 - sse / total
    degrees = n - 2
    residual_variance = sse / degrees if degrees > 0 else float("nan")
    slope_se = math.sqrt(residual_variance / sxx) if degrees > 0 else float("nan")
    safe_sse = max(sse, 1e-300)
    parameter_count = 2
    aic = n * math.log(safe_sse / n) + 2 * parameter_count
    bic = n * math.log(safe_sse / n) + parameter_count * math.log(n)
    residual_std = stdev(residuals) if len(residuals) >= 2 else 0.0
    denominator = sum(residual * residual for residual in residuals)
    if denominator == 0 or n < 2:
        durbin_watson = float("nan")
    else:
        durbin_watson = sum((residuals[i] - residuals[i - 1]) ** 2 for i in range(1, n)) / denominator
    return {
        "intercept": intercept,
        "slope": slope,
        "slope_se": slope_se,
        "sse": sse,
        "r2": r2,
        "aic": aic,
        "bic": bic,
        "residual_mean": mean(residuals),
        "residual_std": residual_std,
        "max_abs_residual": max(abs(residual) for residual in residuals),
        "durbin_watson": durbin_watson,
    }


def fit_loglog_exponent(
    summary_rows: list[dict[str, Any]],
    variance_column: str,
    min_L: int | None,
) -> list[dict[str, Any]]:
    """Fit log V(L) = a + p log(log L) for each distribution."""

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        if min_L is not None and int(row["L"]) < min_L:
            continue
        value = float(row[variance_column])
        if value > 0.0 and math.isfinite(value):
            grouped[(str(row["distribution"]), str(row["distribution_params"]))].append(row)

    fit_rows: list[dict[str, Any]] = []
    for (distribution, params), rows in sorted(grouped.items()):
        rows = sorted(rows, key=lambda row: int(row["L"]))
        xs = [float(row["log_log_L"]) for row in rows]
        ys = [math.log(float(row[variance_column])) for row in rows]
        fit = _linear_fit(xs, ys)
        slope = fit["slope"]
        slope_se = fit["slope_se"]
        fit_rows.append(
            {
                "distribution": distribution,
                "distribution_params": params,
                "variance_kind": variance_column.replace("_variance", ""),
                "fit_min_L": "" if min_L is None else min_L,
                "n_sizes": len(rows),
                "L_values": " ".join(str(row["L"]) for row in rows),
                "intercept_log_C": fit["intercept"],
                "p": slope,
                "p_se": slope_se,
                "p_ci_low": slope - 1.96 * slope_se if math.isfinite(slope_se) else float("nan"),
                "p_ci_high": slope + 1.96 * slope_se if math.isfinite(slope_se) else float("nan"),
                "r2": fit["r2"],
                "sse": fit["sse"],
                "aic": fit["aic"],
                "bic": fit["bic"],
                "residual_mean": fit["residual_mean"],
                "residual_std": fit["residual_std"],
                "max_abs_residual": fit["max_abs_residual"],
                "durbin_watson": fit["durbin_watson"],
            }
        )
    return fit_rows


def bootstrap_exponent_ci(
    raw_rows: list[dict[str, str]],
    summary_rows: list[dict[str, Any]],
    variance_kind: str,
    min_L: int | None,
    reps: int,
    seed: int,
) -> dict[tuple[str, str], tuple[float, float]]:
    """Bootstrap exponent intervals by resampling observations or environments."""

    if reps <= 1:
        return {}
    rng = random.Random(seed)
    summaries_by_model: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        if min_L is None or int(row["L"]) >= min_L:
            summaries_by_model[(str(row["distribution"]), str(row["distribution_params"]))].append(row)

    raw_by_group: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        raw_by_group[(row["distribution"], parameter_key(row["distribution_params"]), int(row["L"]))].append(row)

    intervals: dict[tuple[str, str], tuple[float, float]] = {}
    for key, rows in summaries_by_model.items():
        rows = sorted(rows, key=lambda row: int(row["L"]))
        if len(rows) < 3:
            continue
        slopes: list[float] = []
        for _ in range(reps):
            pseudo_summary: list[dict[str, Any]] = []
            for row in rows:
                L = int(row["L"])
                sample_rows = raw_by_group[(key[0], key[1], L)]
                if variance_kind == "annealed":
                    windings = [float(sample_rows[rng.randrange(len(sample_rows))]["winding"]) for _ in sample_rows]
                    value = variance(windings) if len(windings) >= 2 else float("nan")
                else:
                    by_env: dict[str, list[float]] = defaultdict(list)
                    for sample_row in sample_rows:
                        by_env[sample_row["environment_id"]].append(float(sample_row["winding"]))
                    env_items = list(by_env.items())
                    sampled_env_vars = []
                    for _env_id, values in (env_items[rng.randrange(len(env_items))] for _ in env_items):
                        if len(values) >= 2:
                            sampled_walks = [values[rng.randrange(len(values))] for _ in values]
                            sampled_env_vars.append(variance(sampled_walks))
                    value = mean(sampled_env_vars) if sampled_env_vars else float("nan")
                pseudo = dict(row)
                pseudo[f"{variance_kind}_variance"] = value
                pseudo_summary.append(pseudo)
            fit_rows = fit_loglog_exponent(pseudo_summary, f"{variance_kind}_variance", min_L=None)
            if fit_rows and math.isfinite(float(fit_rows[0]["p"])):
                slopes.append(float(fit_rows[0]["p"]))
        if slopes:
            slopes.sort()
            low_index = int(0.025 * (len(slopes) - 1))
            high_index = int(0.975 * (len(slopes) - 1))
            intervals[key] = (slopes[low_index], slopes[high_index])
    return intervals


def add_bootstrap_columns(
    fit_rows: list[dict[str, Any]],
    intervals: dict[tuple[str, str], tuple[float, float]],
) -> None:
    for row in fit_rows:
        interval = intervals.get((str(row["distribution"]), str(row["distribution_params"])))
        if interval is None:
            row["p_bootstrap_ci_low"] = float("nan")
            row["p_bootstrap_ci_high"] = float("nan")
        else:
            row["p_bootstrap_ci_low"] = interval[0]
            row["p_bootstrap_ci_high"] = interval[1]


def make_pointwise_ratios(summary_rows: list[dict[str, Any]], variance_column: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in sorted(summary_rows, key=lambda item: (item["distribution"], item["distribution_params"], int(item["L"]))):
        value = float(row[variance_column])
        log_log_l = float(row["log_log_L"])
        rows.append(
            {
                "distribution": row["distribution"],
                "distribution_params": row["distribution_params"],
                "variance_kind": variance_column.replace("_variance", ""),
                "L": row["L"],
                "variance": value,
                "log_variance": math.log(value) if value > 0.0 else float("nan"),
                "log_log_L": log_log_l,
                "log_variance_over_log_log_L": math.log(value) / log_log_l if value > 0.0 else float("nan"),
            }
        )
    return rows


def make_local_exponents(summary_rows: list[dict[str, Any]], variance_column: str) -> list[dict[str, Any]]:
    """Compute local p using log(log L_{i+1}) - log(log L_i) in the denominator."""

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        value = float(row[variance_column])
        if value > 0.0 and math.isfinite(value):
            grouped[(str(row["distribution"]), str(row["distribution_params"]))].append(row)

    rows: list[dict[str, Any]] = []
    for (distribution, params), group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["L"]))
        for before, after in zip(ordered, ordered[1:]):
            log_var_before = math.log(float(before[variance_column]))
            log_var_after = math.log(float(after[variance_column]))
            denom = float(after["log_log_L"]) - float(before["log_log_L"])
            rows.append(
                {
                    "distribution": distribution,
                    "distribution_params": params,
                    "variance_kind": variance_column.replace("_variance", ""),
                    "L_left": before["L"],
                    "L_right": after["L"],
                    "L_mid_geometric": math.sqrt(int(before["L"]) * int(after["L"])),
                    "p_local": (log_var_after - log_var_before) / denom if denom != 0.0 else float("nan"),
                }
            )
    return rows


def write_plots(output_dir: Path, summary_rows: list[dict[str, Any]], fit_rows: list[dict[str, Any]]) -> str:
    """Write the required HPC analysis plots when matplotlib is installed."""

    plot_cache = output_dir / "plot_cache"
    plot_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(plot_cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(plot_cache / "xdg"))

    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except Exception as exc:  # noqa: BLE001 - optional plotting stack.
        message = f"plots skipped: {exc}"
        (output_dir / "plots_skipped.txt").write_text(message + "\n", encoding="utf-8")
        return message

    if not summary_rows:
        return "plots skipped: no summary rows"

    def label(row: dict[str, Any]) -> str:
        params = row["distribution_params"]
        return str(row["distribution"]) if params == "{}" else f"{row['distribution']} {params}"

    def grouped_summary() -> dict[str, list[dict[str, Any]]]:
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in summary_rows:
            groups[label(row)].append(row)
        return groups

    groups = grouped_summary()
    plot_specs = [
        ("variance_vs_L.png", "L", lambda row: float(row["L"]), "L"),
        ("variance_vs_log_L.png", "log L", lambda row: float(row["log_L"]), "log L"),
        (
            "variance_vs_log_L_squared.png",
            "(log L)^2",
            lambda row: float(row["log_L"]) ** 2,
            "(log L)^2",
        ),
    ]
    for filename, title_x, x_fn, xlabel in plot_specs:
        _, ax = plt.subplots(figsize=(9, 6))
        for group_label, rows in sorted(groups.items()):
            ordered = sorted(rows, key=lambda row: int(row["L"]))
            ax.plot([x_fn(row) for row in ordered], [float(row["annealed_variance"]) for row in ordered], marker="o", label=group_label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Annealed Var(winding)")
        ax.set_title(f"Variance vs {title_x}")
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(output_dir / filename, dpi=180)
        plt.close()

    _, ax = plt.subplots(figsize=(9, 6))
    for group_label, rows in sorted(groups.items()):
        ordered = [row for row in sorted(rows, key=lambda row: int(row["L"])) if float(row["annealed_variance"]) > 0.0]
        ax.plot(
            [float(row["log_log_L"]) for row in ordered],
            [math.log(float(row["annealed_variance"])) for row in ordered],
            marker="o",
            label=group_label,
        )
    ax.set_xlabel("log(log L)")
    ax.set_ylabel("log Var(winding)")
    ax.set_title("Log variance vs log log L")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_dir / "log_variance_vs_log_log_L.png", dpi=180)
    plt.close()

    local_rows = make_local_exponents(summary_rows, "annealed_variance")
    local_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in local_rows:
        local_groups[str(row["distribution"]) if row["distribution_params"] == "{}" else f"{row['distribution']} {row['distribution_params']}"].append(row)
    _, ax = plt.subplots(figsize=(9, 6))
    for group_label, rows in sorted(local_groups.items()):
        ordered = sorted(rows, key=lambda row: float(row["L_mid_geometric"]))
        ax.plot([float(row["L_mid_geometric"]) for row in ordered], [float(row["p_local"]) for row in ordered], marker="o", label=group_label)
    ax.set_xscale("log")
    ax.set_xlabel("geometric mean of adjacent L values")
    ax.set_ylabel("local effective exponent")
    ax.set_title("Local effective exponent vs L")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_dir / "local_effective_exponent_vs_L.png", dpi=180)
    plt.close()

    exponent_rows = [row for row in fit_rows if row["variance_kind"] == "annealed"]
    if exponent_rows:
        _, ax = plt.subplots(figsize=(10, 6))
        labels = [
            str(row["distribution"]) if row["distribution_params"] == "{}" else f"{row['distribution']} {row['distribution_params']}"
            for row in exponent_rows
        ]
        estimates = [float(row["p"]) for row in exponent_rows]
        lower = [max(0.0, estimate - float(row["p_ci_low"])) if math.isfinite(float(row["p_ci_low"])) else 0.0 for estimate, row in zip(estimates, exponent_rows)]
        upper = [max(0.0, float(row["p_ci_high"]) - estimate) if math.isfinite(float(row["p_ci_high"])) else 0.0 for estimate, row in zip(estimates, exponent_rows)]
        x_values = list(range(len(exponent_rows)))
        ax.errorbar(x_values, estimates, yerr=[lower, upper], fmt="o", capsize=4)
        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        ax.axhline(2.0, linestyle=":", linewidth=1.0)
        ax.set_xticks(x_values)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel("estimated p")
        ax.set_title("Estimated exponent by distribution")
        plt.tight_layout()
        plt.savefig(output_dir / "estimated_exponent_by_distribution.png", dpi=180)
        plt.close()

    skipped_path = output_dir / "plots_skipped.txt"
    if skipped_path.exists():
        skipped_path.unlink()
    return "plots written"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Combine and analyse HPC LERW batch outputs.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, default=Path("results_hpc"))
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_hpc"))
    parser.add_argument("--fit-min-L", type=int, default=None)
    parser.add_argument("--bootstrap-reps", type=int, default=0)
    parser.add_argument("--bootstrap-seed", type=int, default=20260623)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    config_rows = read_config_table(args.config)
    raw_rows, validation_rows = validate_and_collect(config_rows, args.results_dir, args.allow_incomplete)
    duplicates = duplicate_rows(raw_rows)
    if duplicates:
        write_csv(args.output_dir / "duplicates.csv", duplicates)
        raise RuntimeError(f"found {len(duplicates)} duplicate observation keys")

    summary_rows = make_summary(raw_rows)
    annealed_fits = fit_loglog_exponent(summary_rows, "annealed_variance", args.fit_min_L)
    quenched_fits = fit_loglog_exponent(summary_rows, "quenched_variance", args.fit_min_L)
    if args.bootstrap_reps > 1:
        add_bootstrap_columns(
            annealed_fits,
            bootstrap_exponent_ci(
                raw_rows,
                summary_rows,
                variance_kind="annealed",
                min_L=args.fit_min_L,
                reps=args.bootstrap_reps,
                seed=args.bootstrap_seed,
            ),
        )
        add_bootstrap_columns(
            quenched_fits,
            bootstrap_exponent_ci(
                raw_rows,
                summary_rows,
                variance_kind="quenched",
                min_L=args.fit_min_L,
                reps=args.bootstrap_reps,
                seed=args.bootstrap_seed + 7919,
            ),
        )
    else:
        add_bootstrap_columns(annealed_fits, {})
        add_bootstrap_columns(quenched_fits, {})

    fit_rows = annealed_fits + quenched_fits
    pointwise_rows = make_pointwise_ratios(summary_rows, "annealed_variance") + make_pointwise_ratios(
        summary_rows,
        "quenched_variance",
    )
    local_rows = make_local_exponents(summary_rows, "annealed_variance") + make_local_exponents(
        summary_rows,
        "quenched_variance",
    )

    write_csv(args.output_dir / "validation.csv", validation_rows)
    write_csv(args.output_dir / "combined_raw.csv", raw_rows)
    write_csv(args.output_dir / "summary.csv", summary_rows)
    write_csv(args.output_dir / "loglog_fits.csv", fit_rows)
    write_csv(args.output_dir / "pointwise_ratios.csv", pointwise_rows)
    write_csv(args.output_dir / "local_effective_exponents.csv", local_rows)

    plot_status = "plots skipped by --no-plots" if args.no_plots else write_plots(args.output_dir, summary_rows, fit_rows)
    print(f"Analysed {len(raw_rows)} observations from {len(config_rows)} tasks")
    print(f"Summary rows: {len(summary_rows)}")
    print(f"Fit rows: {len(fit_rows)}")
    print(f"Output: {args.output_dir}")
    print(plot_status)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
