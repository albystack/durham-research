#!/usr/bin/env python3
"""Summarise LERW winding CSVs and compare simple scaling fits."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, variance


@dataclass(frozen=True)
class GroupStats:
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


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def grouped_stats(rows: list[dict[str, str]]) -> list[GroupStats]:
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


def print_group_table(stats: list[GroupStats]) -> None:
    print(
        "model,k,L,samples,mean_winding,var_winding,"
        "mean_raw_steps,mean_lerw_steps,mean_hit_x,mean_hit_y"
    )
    for item in stats:
        print(
            f"{item.model},{item.k},{item.L},{item.samples},"
            f"{item.mean_winding:.6g},{item.var_winding:.6g},"
            f"{item.mean_raw_steps:.6g},{item.mean_lerw_steps:.6g},"
            f"{item.mean_hit_x:.6g},{item.mean_hit_y:.6g}"
        )


def print_scaling_fits(stats: list[GroupStats]) -> None:
    by_model: dict[tuple[str, str], list[GroupStats]] = defaultdict(list)
    for item in stats:
        by_model[(item.model, item.k)].append(item)

    print("\nScaling fits")
    print("model,k,fit,slope,intercept,r_squared")
    for (model, k), items in sorted(by_model.items()):
        items = sorted(items, key=lambda item: item.L)
        if len(items) < 2:
            continue

        ys = [item.var_winding for item in items]
        log_xs = [math.log(item.L) for item in items]
        log_sq_xs = [math.log(item.L) ** 2 for item in items]

        log_fit = linear_fit(log_xs, ys)
        log_sq_fit = linear_fit(log_sq_xs, ys)

        print(
            f"{model},{k},logL,"
            f"{log_fit[0]:.6g},{log_fit[1]:.6g},{log_fit[2]:.6g}"
        )
        print(
            f"{model},{k},logL_squared,"
            f"{log_sq_fit[0]:.6g},{log_sq_fit[1]:.6g},{log_sq_fit[2]:.6g}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarise LERW winding simulation CSV files."
    )
    parser.add_argument("csv_files", type=Path, nargs="+")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_rows(args.csv_files)
    stats = grouped_stats(rows)
    print_group_table(stats)
    print_scaling_fits(stats)


if __name__ == "__main__":
    main()
