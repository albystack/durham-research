"""Height-convention figures and exact centre-height distribution audits."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .environment import DimerEnvironment
from .exact import (
    deterministic_test_weights,
    exact_center_height_distribution,
)
from .glauber import run_chain
from .height import face_heights
from .matching import DimerState


def _draw_matching_height(
    axis: plt.Axes,
    state: DimerState,
    heights: np.ndarray,
    title: str,
) -> None:
    """Draw grid, vertex colours, occupied dimers, face heights, and anchor."""
    L = state.L
    for row in range(L):
        axis.plot([0, L - 1], [row, row], color="0.82", linewidth=0.7, zorder=0)
    for column in range(L):
        axis.plot([column, column], [0, L - 1], color="0.82", linewidth=0.7, zorder=0)

    for row in range(L):
        for column in range(L - 1):
            if state.horizontal_occupied[row, column]:
                axis.plot(
                    [column, column + 1],
                    [row, row],
                    color="#1769aa",
                    linewidth=4.0,
                    solid_capstyle="round",
                    zorder=2,
                )
    for row in range(L - 1):
        for column in range(L):
            if state.vertical_occupied[row, column]:
                axis.plot(
                    [column, column],
                    [row, row + 1],
                    color="#1769aa",
                    linewidth=4.0,
                    solid_capstyle="round",
                    zorder=2,
                )

    white_x, white_y, black_x, black_y = [], [], [], []
    for row in range(L):
        for column in range(L):
            if (row + column) % 2 == 0:
                white_x.append(column)
                white_y.append(row)
            else:
                black_x.append(column)
                black_y.append(row)
    axis.scatter(white_x, white_y, s=20, facecolor="white", edgecolor="black", zorder=3)
    axis.scatter(black_x, black_y, s=20, facecolor="black", edgecolor="black", zorder=3)

    for row in range(L - 1):
        for column in range(L - 1):
            axis.text(
                column + 0.5,
                row + 0.5,
                str(int(heights[row, column])),
                color="#9b1c31",
                fontsize=9 if L <= 6 else 7.5,
                ha="center",
                va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.5},
                zorder=4,
            )
    axis.annotate(
        "fixed cut\n$h_{ext}=0$",
        xy=(0.5, 0.02),
        xytext=(0.5, -0.9),
        ha="center",
        va="top",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "color": "#9b1c31"},
        color="#9b1c31",
    )
    axis.set_title(title)
    axis.set_aspect("equal")
    axis.set_xlim(-0.35, L - 0.65)
    axis.set_ylim(L - 0.65, -1.25)
    axis.set_xticks(range(L))
    axis.set_yticks(range(L))
    axis.set_xlabel("vertex column")
    axis.set_ylabel("vertex row")


def write_height_convention_audit(
    output: str | Path, sizes: Iterable[int] = (6, 8)
) -> dict[str, object]:
    """Write full extremal height matrices and labelled dimer figures."""
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    matrix_rows: list[dict[str, object]] = []
    report_lines = [
        "# Height-convention audit",
        "",
        "White vertices satisfy `(row + column) % 2 == 0`. Primal edges are",
        "oriented white-to-black. A left-to-right dual crossing is clockwise",
        "around white and contributes +1 if unoccupied or -3 if occupied.",
        "The exterior reference h=0 crosses the fixed top-left horizontal edge.",
        "",
        "Swapping white and black reverses every oriented crossing, so with the",
        "same zero reference the entire height matrix is multiplied by -1.",
        "Changing only the numerical additive reference from 0 to C adds the",
        "same constant C to every bounded-face height. The boundary cut itself",
        "is part of the convention and is not changed in this comparison.",
        "",
    ]
    audited_sizes = []
    for raw_L in sizes:
        L = int(raw_L)
        audited_sizes.append(L)
        states = {
            "all_horizontal": DimerState.all_horizontal(L),
            "all_vertical": DimerState.all_vertical(L),
        }
        figure, axes = plt.subplots(1, 2, figsize=(2 * max(5.2, L * 0.8), max(5.0, L * 0.8)))
        report_lines.extend((f"## L={L}", ""))
        for axis, (orientation, state) in zip(axes, states.items(), strict=True):
            heights = face_heights(state)
            np.savetxt(
                output_path / f"L{L}_{orientation}_height_matrix.csv",
                heights,
                fmt="%d",
                delimiter=",",
            )
            (output_path / f"L{L}_{orientation}_height_matrix.txt").write_text(
                np.array2string(heights, separator=" ") + "\n"
            )
            for row in range(L - 1):
                for column in range(L - 1):
                    matrix_rows.append(
                        {
                            "L": L,
                            "orientation": orientation,
                            "face_row": row,
                            "face_column": column,
                            "height": int(heights[row, column]),
                        }
                    )
            _draw_matching_height(axis, state, heights, orientation.replace("_", " "))
            report_lines.extend(
                (
                    f"### {orientation.replace('_', ' ')}",
                    "",
                    "```text",
                    np.array2string(heights, separator=" "),
                    "```",
                    "",
                )
            )
        figure.suptitle(
            f"L={L}: occupied dimers (blue) and bounded-face heights (red)", fontsize=12
        )
        figure.tight_layout()
        figure.savefig(output_path / f"L{L}_extremal_height_audit.png", dpi=180)
        plt.close(figure)
    pd.DataFrame(matrix_rows).to_csv(output_path / "height_matrices_long.csv", index=False)
    (output_path / "height_convention_audit.md").write_text("\n".join(report_lines) + "\n")
    metadata = {
        "sizes": audited_sizes,
        "white_vertex_rule": "(row + column) % 2 == 0",
        "crossing_rule": "white-to-black left-to-right: +1 unoccupied, -3 occupied",
        "additive_reference": "exterior h=0 through fixed top-left horizontal boundary cut",
        "colour_swap_effect": "global sign reversal",
        "reference_shift_effect": "global additive constant",
    }
    (output_path / "height_convention_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    return metadata


def write_exact_center_height_audit(
    output: str | Path,
    *,
    seed: int = 20260908,
    burnin_sweeps: int = 10_000,
    measurement_gap_sweeps: int = 2,
    num_measurements: int = 100_000,
) -> dict[str, object]:
    """Write exact tiny-grid PMFs and a long diagnostic 4x4 chain comparison."""
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    exact_rows = []
    pmf_rows = []
    weight_rows = []
    distributions = {}
    for rows, columns in ((2, 2), (2, 4), (4, 4)):
        horizontal, vertical = deterministic_test_weights(rows, columns)
        distribution = exact_center_height_distribution(
            rows, columns, horizontal, vertical
        )
        distributions[(rows, columns)] = distribution
        exact_rows.append(
            {
                "rows": rows,
                "columns": columns,
                "matching_count": distribution.matching_count,
                "center_row": distribution.center_index[0],
                "center_column": distribution.center_index[1],
                "exact_mean": distribution.mean,
                "exact_variance": distribution.variance,
            }
        )
        for height, probability in distribution.pmf.items():
            pmf_rows.append(
                {
                    "rows": rows,
                    "columns": columns,
                    "center_height": height,
                    "exact_probability": probability,
                    "empirical_probability": np.nan,
                }
            )
        for i in range(rows):
            for j in range(columns - 1):
                weight_rows.append(
                    {"rows": rows, "columns": columns, "edge": "horizontal", "i": i, "j": j, "weight": horizontal[i, j]}
                )
        for i in range(rows - 1):
            for j in range(columns):
                weight_rows.append(
                    {"rows": rows, "columns": columns, "edge": "vertical", "i": i, "j": j, "weight": vertical[i, j]}
                )

    horizontal, vertical = deterministic_test_weights(4, 4)
    environment = DimerEnvironment(
        4, horizontal, vertical, model="deterministic_exact_audit"
    )
    started = perf_counter()
    chain = run_chain(
        environment,
        DimerState.all_horizontal(4),
        np.random.Generator(np.random.PCG64(seed)),
        burnin_sweeps=burnin_sweeps,
        measurement_gap_sweeps=measurement_gap_sweeps,
        num_measurements=num_measurements,
    )
    wall_seconds = perf_counter() - started
    empirical_values, empirical_counts = np.unique(chain.center_heights, return_counts=True)
    empirical_pmf = {
        int(value): float(count / num_measurements)
        for value, count in zip(empirical_values, empirical_counts, strict=True)
    }
    exact_4 = distributions[(4, 4)]
    union = sorted(set(exact_4.pmf) | set(empirical_pmf))
    for row in pmf_rows:
        if row["rows"] == 4 and row["columns"] == 4:
            row["empirical_probability"] = empirical_pmf.get(int(row["center_height"]), 0.0)
    empirical_array = np.asarray(chain.center_heights, dtype=np.float64)
    empirical_mean = float(empirical_array.mean())
    empirical_variance = float(empirical_array.var(ddof=0))
    total_variation = 0.5 * sum(
        abs(exact_4.pmf.get(value, 0.0) - empirical_pmf.get(value, 0.0)) for value in union
    )
    maximum_probability_error = max(
        abs(exact_4.pmf.get(value, 0.0) - empirical_pmf.get(value, 0.0)) for value in union
    )
    comparison = {
        "seed": seed,
        "burnin_sweeps": burnin_sweeps,
        "measurement_gap_sweeps": measurement_gap_sweeps,
        "num_measurements": num_measurements,
        "exact_mean": exact_4.mean,
        "empirical_mean": empirical_mean,
        "mean_error": empirical_mean - exact_4.mean,
        "exact_variance": exact_4.variance,
        "empirical_variance": empirical_variance,
        "variance_error": empirical_variance - exact_4.variance,
        "total_variation_distance": total_variation,
        "maximum_probability_error": maximum_probability_error,
        "runtime_seconds": wall_seconds,
        **chain.diagnostics(),
    }
    pd.DataFrame(exact_rows).to_csv(output_path / "exact_center_moments.csv", index=False)
    pd.DataFrame(pmf_rows).to_csv(output_path / "exact_center_pmf.csv", index=False)
    pd.DataFrame(weight_rows).to_csv(output_path / "deterministic_test_weights.csv", index=False)
    pd.DataFrame([comparison]).to_csv(output_path / "four_by_four_chain_comparison.csv", index=False)
    (output_path / "four_by_four_chain_comparison.json").write_text(
        json.dumps(comparison, indent=2) + "\n"
    )

    figure, axis = plt.subplots(figsize=(6.5, 4.2))
    positions = np.arange(len(union))
    width = 0.36
    axis.bar(
        positions - width / 2,
        [exact_4.pmf.get(value, 0.0) for value in union],
        width,
        label="exact Gibbs",
    )
    axis.bar(
        positions + width / 2,
        [empirical_pmf.get(value, 0.0) for value in union],
        width,
        label="long Glauber diagnostic",
    )
    axis.set_xticks(positions, [str(value) for value in union])
    axis.set(xlabel="centre height", ylabel="probability", title="4x4 exact vs empirical centre height")
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_path / "four_by_four_exact_vs_empirical.png", dpi=180)
    plt.close(figure)

    report = {
        "weights": "deterministic positive non-symmetric test arrays saved in CSV",
        "exact_distributions": [asdict(distributions[key]) for key in ((2, 2), (2, 4), (4, 4))],
        "four_by_four_monte_carlo": comparison,
        "warning": "The long-chain comparison is a diagnostic, not a normal unit test or proof of mixing.",
    }
    (output_path / "exact_height_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
