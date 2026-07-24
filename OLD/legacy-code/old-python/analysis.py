"""Analysis, scaling fits, and plotting for simulation output.

The code uses `pandas`, `numpy`, `matplotlib`, and `seaborn` when they are
installed.  If they are missing, the command-line runner still writes CSV
tables with standard-library code and records that plots were skipped.
"""

from __future__ import annotations

import csv
import math
import os
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean, variance
from typing import Any

from config import model_color, model_label, safe_model_name


def _optional_imports() -> tuple[Any | None, Any | None]:
    """Return pandas/numpy modules when available."""

    try:
        import pandas as pd  # type: ignore
        import numpy as np  # type: ignore
    except Exception:
        return None, None
    return pd, np


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write dictionaries as a CSV file with stable column order."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parameter_key(value: object) -> str:
    """Normalise empty and numeric parameter values for grouping."""

    if value in ("", None, "None"):
        return ""
    return f"{float(value):g}"


def make_summary_table(raw_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Aggregate trial rows by model, parameter, and L.

    Pandas handles the group-by path when available.  The fallback below keeps
    the repository runnable in a bare Python installation.
    """

    pd, _ = _optional_imports()
    if pd is not None and raw_rows:
        frame = pd.DataFrame(raw_rows)
        for column in ["L", "raw_steps", "lerw_steps", "winding", "winding_angle", "hit_x", "hit_y"]:
            frame[column] = pd.to_numeric(frame[column])
        frame["parameter"] = frame["parameter"].apply(parameter_key)
        grouped = frame.groupby(["model", "parameter", "L"], as_index=False)
        summary = grouped.agg(
            samples=("winding", "size"),
            mean_winding=("winding", "mean"),
            variance_winding=("winding", "var"),
            mean_winding_angle=("winding_angle", "mean"),
            variance_winding_angle=("winding_angle", "var"),
            mean_raw_steps=("raw_steps", "mean"),
            mean_lerw_steps=("lerw_steps", "mean"),
            mean_hit_x_over_L=("hit_x", lambda values: float(values.mean())),
            mean_hit_y_over_L=("hit_y", lambda values: float(values.mean())),
        )
        summary["logL"] = summary["L"].apply(math.log)
        summary["mean_hit_x_over_L"] = summary["mean_hit_x_over_L"] / summary["L"]
        summary["mean_hit_y_over_L"] = summary["mean_hit_y_over_L"] / summary["L"]
        summary["variance_winding_se"] = summary.apply(
            lambda row: variance_standard_error(float(row["variance_winding"]), int(row["samples"])),
            axis=1,
        )
        summary["variance_winding_angle_se"] = summary.apply(
            lambda row: variance_standard_error(float(row["variance_winding_angle"]), int(row["samples"])),
            axis=1,
        )
        summary["model_label"] = summary.apply(
            lambda row: model_label(str(row["model"]), str(row["parameter"])),
            axis=1,
        )
        summary = summary.fillna(0.0)
        columns = [
            "model",
            "parameter",
            "model_label",
            "L",
            "logL",
            "samples",
            "mean_winding",
            "variance_winding",
            "variance_winding_se",
            "mean_winding_angle",
            "variance_winding_angle",
            "variance_winding_angle_se",
            "mean_raw_steps",
            "mean_lerw_steps",
            "mean_hit_x_over_L",
            "mean_hit_y_over_L",
        ]
        return summary[columns].sort_values(["model", "parameter", "L"]).to_dict("records")

    groups: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        key = (str(row["model"]), parameter_key(row["parameter"]), int(row["L"]))
        groups[key].append(row)

    summary_rows: list[dict[str, Any]] = []
    for (model, parameter, L), rows in sorted(groups.items()):
        windings = [float(row["winding"]) for row in rows]
        angles = [float(row["winding_angle"]) for row in rows]
        winding_variance = variance(windings) if len(windings) > 1 else 0.0
        angle_variance = variance(angles) if len(angles) > 1 else 0.0
        summary_rows.append(
            {
                "model": model,
                "parameter": parameter,
                "model_label": model_label(model, parameter),
                "L": L,
                "logL": math.log(L),
                "samples": len(rows),
                "mean_winding": mean(windings),
                "variance_winding": winding_variance,
                "variance_winding_se": variance_standard_error(winding_variance, len(rows)),
                "mean_winding_angle": mean(angles),
                "variance_winding_angle": angle_variance,
                "variance_winding_angle_se": variance_standard_error(angle_variance, len(rows)),
                "mean_raw_steps": mean(float(row["raw_steps"]) for row in rows),
                "mean_lerw_steps": mean(float(row["lerw_steps"]) for row in rows),
                "mean_hit_x_over_L": mean(float(row["hit_x"]) / L for row in rows),
                "mean_hit_y_over_L": mean(float(row["hit_y"]) / L for row in rows),
            }
        )
    return summary_rows


def variance_standard_error(sample_variance: float, sample_count: int) -> float:
    """Normal-approximation standard error for a sample variance.

    This is a readable finite-size error bar for the plots, not a claim that
    the winding distribution is exactly Gaussian.  The scaling constants still
    use bootstrap intervals computed from raw winding samples.
    """

    if sample_count <= 1:
        return 0.0
    return math.sqrt(2.0 / (sample_count - 1)) * sample_variance


def _fit(xs: list[float], ys: list[float], intercept: bool) -> tuple[float, float, float]:
    """Fit y = C x or y = A + C x and return intercept, slope, SSE."""

    _, np = _optional_imports()
    if np is not None and xs:
        x_array = np.asarray(xs, dtype=float)
        y_array = np.asarray(ys, dtype=float)
        design = np.column_stack([np.ones_like(x_array), x_array]) if intercept else x_array[:, None]
        coefficients, *_ = np.linalg.lstsq(design, y_array, rcond=None)
        if intercept:
            a = float(coefficients[0])
            c = float(coefficients[1])
            fitted = a + c * x_array
        else:
            a = 0.0
            c = float(coefficients[0])
            fitted = c * x_array
        sse = float(np.sum((y_array - fitted) ** 2))
        return a, c, sse

    if not xs:
        return 0.0, 0.0, float("inf")
    if intercept:
        x_bar = mean(xs)
        y_bar = mean(ys)
        denominator = sum((x - x_bar) ** 2 for x in xs)
        if denominator == 0:
            return 0.0, 0.0, float("inf")
        c = sum((x - x_bar) * (y - y_bar) for x, y in zip(xs, ys)) / denominator
        a = y_bar - c * x_bar
        sse = sum((y - (a + c * x)) ** 2 for x, y in zip(xs, ys))
        return a, c, sse

    denominator = sum(x * x for x in xs)
    if denominator == 0:
        return 0.0, 0.0, float("inf")
    c = sum(x * y for x, y in zip(xs, ys)) / denominator
    sse = sum((y - c * x) ** 2 for x, y in zip(xs, ys))
    return 0.0, c, sse


def _r_squared(ys: list[float], sse: float) -> float:
    """Coefficient of determination for descriptive comparison."""

    if len(ys) < 2:
        return 0.0
    y_bar = mean(ys)
    total = sum((y - y_bar) ** 2 for y in ys)
    return 0.0 if total == 0 else 1.0 - sse / total


def _aic_bic(n: int, sse: float, parameter_count: int) -> tuple[float, float]:
    """Small diagnostic scores; lower is better."""

    safe_sse = max(sse, 1e-12)
    aic = n * math.log(safe_sse / n) + 2 * parameter_count
    bic = n * math.log(safe_sse / n) + parameter_count * math.log(n)
    return aic, bic


def _bootstrap_slope_ci(
    grouped_windings: dict[int, list[float]],
    transform,
    reps: int,
    seed: int,
) -> tuple[float, float]:
    """Bootstrap a through-origin slope confidence interval.

    This resamples windings within each L, recomputes the variance points, and
    refits the slope.  It is a finite-sample error bar, not an asymptotic proof.
    """

    if reps <= 1:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    slopes: list[float] = []
    sizes = sorted(grouped_windings)

    for _ in range(reps):
        xs: list[float] = []
        ys: list[float] = []
        for L in sizes:
            values = grouped_windings[L]
            if len(values) < 2:
                continue
            sample = [values[rng.randrange(len(values))] for _ in range(len(values))]
            xs.append(transform(math.log(L)))
            ys.append(variance(sample))
        if len(xs) >= 2:
            _, slope, _ = _fit(xs, ys, intercept=False)
            slopes.append(slope)

    if not slopes:
        return float("nan"), float("nan")
    slopes.sort()
    low = slopes[int(0.025 * (len(slopes) - 1))]
    high = slopes[int(0.975 * (len(slopes) - 1))]
    return low, high


def make_fit_table(
    summary_rows: list[dict[str, Any]],
    raw_rows: list[dict[str, Any]],
    bootstrap_reps: int,
    bootstrap_seed: int,
) -> list[dict[str, Any]]:
    """Compare C log L and C (log L)^2 scaling for each model."""

    summaries: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in summary_rows:
        summaries[(str(row["model"]), parameter_key(row["parameter"]))].append(row)

    raw_groups: dict[tuple[str, str, int], list[float]] = defaultdict(list)
    for row in raw_rows:
        raw_groups[(str(row["model"]), parameter_key(row["parameter"]), int(row["L"]))].append(
            float(row["winding"])
        )

    fit_rows: list[dict[str, Any]] = []
    for (model, parameter), rows in sorted(summaries.items()):
        rows = sorted(rows, key=lambda row: int(row["L"]))
        if len(rows) < 2:
            continue

        x_log = [float(row["logL"]) for row in rows]
        x_log2 = [x * x for x in x_log]
        y = [float(row["variance_winding"]) for row in rows]
        n = len(rows)

        _, c_log, sse_log = _fit(x_log, y, intercept=False)
        _, c_log2, sse_log2 = _fit(x_log2, y, intercept=False)
        a_log, slope_log, sse_log_i = _fit(x_log, y, intercept=True)
        a_log2, slope_log2, sse_log2_i = _fit(x_log2, y, intercept=True)

        grouped_windings = {int(row["L"]): raw_groups[(model, parameter, int(row["L"]))] for row in rows}
        ci_log = _bootstrap_slope_ci(
            grouped_windings,
            transform=lambda log_l: log_l,
            reps=bootstrap_reps,
            seed=bootstrap_seed + 17,
        )
        ci_log2 = _bootstrap_slope_ci(
            grouped_windings,
            transform=lambda log_l: log_l * log_l,
            reps=bootstrap_reps,
            seed=bootstrap_seed + 71,
        )

        aic_log, bic_log = _aic_bic(n, sse_log, 1)
        aic_log2, bic_log2 = _aic_bic(n, sse_log2, 1)
        aic_log_i, bic_log_i = _aic_bic(n, sse_log_i, 2)
        aic_log2_i, bic_log2_i = _aic_bic(n, sse_log2_i, 2)

        fit_rows.append(
            {
                "model": model,
                "parameter": parameter,
                "n_sizes": n,
                "C_logL": c_log,
                "C_logL_ci_low": ci_log[0],
                "C_logL_ci_high": ci_log[1],
                "SSE_logL_origin": sse_log,
                "R2_logL_origin": _r_squared(y, sse_log),
                "AIC_logL_origin": aic_log,
                "BIC_logL_origin": bic_log,
                "C_logL2": c_log2,
                "C_logL2_ci_low": ci_log2[0],
                "C_logL2_ci_high": ci_log2[1],
                "SSE_logL2_origin": sse_log2,
                "R2_logL2_origin": _r_squared(y, sse_log2),
                "AIC_logL2_origin": aic_log2,
                "BIC_logL2_origin": bic_log2,
                "intercept_logL": a_log,
                "slope_logL_intercept": slope_log,
                "SSE_logL_intercept": sse_log_i,
                "AIC_logL_intercept": aic_log_i,
                "BIC_logL_intercept": bic_log_i,
                "intercept_logL2": a_log2,
                "slope_logL2_intercept": slope_log2,
                "SSE_logL2_intercept": sse_log2_i,
                "AIC_logL2_intercept": aic_log2_i,
                "BIC_logL2_intercept": bic_log2_i,
                "better_origin_fit_by_sse": "C log L" if sse_log <= sse_log2 else "C (log L)^2",
                "better_origin_fit_by_aic": "C log L" if aic_log <= aic_log2 else "C (log L)^2",
                "better_origin_fit_by_bic": "C log L" if bic_log <= bic_log2 else "C (log L)^2",
            }
        )

    return fit_rows


def write_plots(output_dir: Path, summary_rows: list[dict[str, Any]], fit_rows: list[dict[str, Any]]) -> str:
    """Write PNG plots with seaborn/matplotlib when available."""

    plot_cache = output_dir / "plot_cache"
    plot_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(plot_cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(plot_cache / "xdg"))

    try:
        import matplotlib  # type: ignore

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt  # type: ignore
        import pandas as pd  # type: ignore
        import seaborn as sns  # type: ignore
    except Exception:
        message = "plots skipped; install requirements.txt for matplotlib/seaborn figures"
        (output_dir / "plots_skipped.txt").write_text(message + "\n", encoding="utf-8")
        return message

    if not summary_rows:
        return "plots skipped; no summary rows"

    sns.set_theme(style="whitegrid", context="talk")
    summary = pd.DataFrame(summary_rows)
    summary["parameter"] = summary["parameter"].astype(str)
    summary["label"] = summary.apply(lambda row: model_label(row["model"], row["parameter"]), axis=1)

    # One compact all-model view.  The x-axis is log L because the two
    # theoretical competitors are log L and (log L)^2.
    _, ax = plt.subplots(figsize=(12, 7))
    for label, rows in summary.groupby("label", sort=True):
        rows = rows.sort_values("L")
        color = model_color(str(rows.iloc[0]["model"]))
        ax.errorbar(
            rows["logL"],
            rows["variance_winding"],
            yerr=rows["variance_winding_se"],
            marker="o",
            linewidth=2.0,
            capsize=3,
            label=label,
            color=color,
        )
    ax.set_title("LERW winding variance in random edge environments")
    ax.set_xlabel("log L")
    ax.set_ylabel("Var(winding)")
    ax.legend(title="model", fontsize=8, title_fontsize=9, ncol=1, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    all_path = output_dir / "variance_scaling_all_models.png"
    plt.savefig(all_path, dpi=180)
    plt.close()

    # Per-model figures include the fitted C log L and C (log L)^2 curves.
    fit_map = {(str(row["model"]), str(row["parameter"])): row for row in fit_rows}
    for (model, parameter), rows in summary.groupby(["model", "parameter"], sort=True):
        rows = rows.sort_values("L")
        fit = fit_map.get((str(model), str(parameter)))
        if fit is None:
            continue
        x_min = float(rows["logL"].min())
        x_max = float(rows["logL"].max())
        dense_x = [x_min + (x_max - x_min) * i / 150 for i in range(151)]
        _, ax = plt.subplots(figsize=(9, 6))
        ax.errorbar(
            rows["logL"],
            rows["variance_winding"],
            yerr=rows["variance_winding_se"],
            marker="o",
            linewidth=2.0,
            capsize=4,
            color=model_color(str(model)),
            label="data with SE",
        )
        ax.plot(dense_x, [float(fit["C_logL"]) * x for x in dense_x], label="C log L", color="#dc2626")
        ax.plot(
            dense_x,
            [float(fit["C_logL2"]) * x * x for x in dense_x],
            label="C (log L)^2",
            color="#2563eb",
            linestyle="--",
        )
        ax.set_title(model_label(str(model), str(parameter)))
        ax.set_xlabel("log L")
        ax.set_ylabel("Var(winding)")
        ax.legend()
        plt.tight_layout()
        path = output_dir / f"variance_fit_{safe_model_name(str(model), str(parameter))}.png"
        plt.savefig(path, dpi=180)
        plt.close()

    return "plots written"


def write_run_report(
    path: Path,
    preset_name: str,
    raw_rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    fit_rows: list[dict[str, Any]],
    seconds: float,
    plot_status: str,
) -> None:
    """Small Markdown summary saved next to the result tables."""

    symmetric_fit = next((row for row in fit_rows if row["model"] == "symmetric"), None)
    log2_candidates = [
        row
        for row in fit_rows
        if row["model"] != "symmetric" and row["better_origin_fit_by_bic"] == "C (log L)^2"
    ]

    lines = [
        f"# Run report: `{preset_name}`",
        "",
        f"Trials: {len(raw_rows)}",
        f"Summary rows: {len(summary_rows)}",
        f"Runtime seconds: {seconds:.2f}",
        f"Plot status: {plot_status}",
        "",
        "Height convention:",
        "",
        "The reported height proxy is the quarter-turn winding of the loop-erased path. "
        "A physical dimer-height normalization may multiply this observable by a fixed "
        "orientation-dependent constant; that changes fitted constants but not the "
        "`log L` versus `(log L)^2` exponent comparison.",
        "",
        "Annealed/quenched convention:",
        "",
        "Each trial samples one fixed edge-weight environment and one walk in that environment. "
        "The current variance is annealed over both environment and walk randomness. "
        "A quenched-only variance would require several walks in each fixed environment.",
        "",
    ]
    if symmetric_fit is not None:
        lines.extend(
            [
                "Symmetric baseline:",
                "",
                f"- `C_logL = {float(symmetric_fit['C_logL']):.6g}` "
                f"(bootstrap 95% CI {float(symmetric_fit['C_logL_ci_low']):.6g}, "
                f"{float(symmetric_fit['C_logL_ci_high']):.6g})",
                f"- BIC prefers `{symmetric_fit['better_origin_fit_by_bic']}` for the baseline.",
                "",
            ]
        )
    lines.extend(
        [
            "Potential super-rough candidates:",
            "",
        ]
    )
    if log2_candidates:
        for row in log2_candidates:
            lines.append(
                f"- {model_label(str(row['model']), str(row['parameter']))}: "
                f"BIC prefers `{row['better_origin_fit_by_bic']}`, "
                f"`C_logL2 = {float(row['C_logL2']):.6g}` "
                f"(bootstrap 95% CI {float(row['C_logL2_ci_low']):.6g}, "
                f"{float(row['C_logL2_ci_high']):.6g})"
            )
    else:
        lines.append("- No non-baseline model preferred `C (log L)^2` by BIC in this run.")
    lines.extend(
        [
            "",
        "Fit comparison:",
        "",
            "| model | parameter | better by SSE | better by AIC | better by BIC | C log L 95% CI | C (log L)^2 95% CI |",
            "|---|---:|---|---|---|---:|---:|",
        ]
    )
    for row in fit_rows:
        lines.append(
            f"| {row['model']} | {row['parameter'] or '-'} | "
            f"{row['better_origin_fit_by_sse']} | {row['better_origin_fit_by_aic']} | "
            f"{row['better_origin_fit_by_bic']} | "
            f"[{float(row['C_logL_ci_low']):.4g}, {float(row['C_logL_ci_high']):.4g}] | "
            f"[{float(row['C_logL2_ci_low']):.4g}, {float(row['C_logL2_ci_high']):.4g}] |"
        )
    lines.extend(
        [
            "",
            "Interpretation note: winding is used as the tree-side height proxy. "
            "The finite-size fits are diagnostics for log L versus (log L)^2 growth, not proofs.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
