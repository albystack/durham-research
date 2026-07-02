#!/usr/bin/env python3
"""Create scaling figures from the Julia aggregate-analysis CSVs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D


LABELS = {
    ("baseline", "{}"): "Baseline",
    ("gamma", '{"shape":0.5}'): "Gamma(shape=0.5)",
    ("gamma", '{"shape":1.0}'): "Gamma(shape=1)",
    ("gamma", '{"shape":2.0}'): "Gamma(shape=2)",
    ("exponential", "{}"): "Exponential",
    ("lognormal", '{"sigma":0.5}'): "Lognormal(σ=0.5)",
    ("lognormal", '{"sigma":1.0}'): "Lognormal(σ=1)",
    ("pareto", '{"alpha":2.0}'): "Pareto(α=2)",
    ("pareto", '{"alpha":3.0}'): "Pareto(α=3)",
    ("uniform", '{"a":0.7}'): "Uniform(a=0.7)",
    ("beta", '{"a":0.5}'): "Beta(a=0.5)",
    ("weibull", '{"shape":0.7}'): "Weibull(shape=0.7)",
    ("inverse_gamma", '{"alpha":2.2}'): "Inverse gamma(α=2.2)",
    ("bernoulli", '{"a":0.8}'): "Bernoulli(a=0.8)",
    ("triangular", '{"a":0.8}'): "Triangular(a=0.8)",
}

ORDER = list(LABELS)
DATA_COLOR = "#202124"
LOG_COLOR = "#1769aa"
LOG2_COLOR = "#d97706"
GRID_COLOR = "#d8dee7"


def canonical_params(value: str) -> str:
    parsed = json.loads(value)
    return json.dumps(parsed, separators=(",", ":"), sort_keys=True)


def key_columns(frame: pd.DataFrame) -> pd.Series:
    return list(zip(frame["distribution"], frame["distribution_params"], strict=True))


def prepare(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(path / "summary.csv")
    fits = pd.read_csv(path / "loglog_fits.csv")
    models = pd.read_csv(path / "scaling_model_comparison.csv")
    for frame in (summary, fits, models):
        frame["distribution_params"] = frame["distribution_params"].map(canonical_params)
        frame["key"] = key_columns(frame)
    return summary, fits, models


def model_row(models: pd.DataFrame, key: tuple[str, str], kind: str, model: str):
    selected = models[
        (models["key"] == key)
        & (models["variance_kind"] == kind)
        & (models["model"] == model)
    ]
    if len(selected) != 1:
        raise ValueError(f"expected one model row for {key}, {kind}, {model}")
    return selected.iloc[0]


def fit_row(fits: pd.DataFrame, key: tuple[str, str], kind: str):
    selected = fits[(fits["key"] == key) & (fits["variance_kind"] == kind)]
    if len(selected) != 1:
        raise ValueError(f"expected one exponent row for {key}, {kind}")
    return selected.iloc[0]


def plot_panel(ax, summary, fits, models, key, kind, *, annotate=True):
    rows = summary[summary["key"] == key].sort_values("L")
    variance = f"{kind}_variance"
    standard_error = f"{kind}_variance_se"
    x = np.log2(rows["L"].to_numpy(dtype=float))
    y = rows[variance].to_numpy(dtype=float)
    yerr = 1.96 * rows[standard_error].to_numpy(dtype=float)

    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt="o",
        color=DATA_COLOR,
        ecolor="#697386",
        elinewidth=1.0,
        capsize=2.5,
        markersize=4.2,
        zorder=4,
    )

    smooth_l = np.geomspace(rows["L"].min(), rows["L"].max(), 240)
    smooth_x = np.log2(smooth_l)
    log_model = model_row(models, key, kind, "a_plus_b_log_L")
    log2_model = model_row(models, key, kind, "a_plus_b_log_L_squared")
    ax.plot(
        smooth_x,
        log_model["intercept"] + log_model["coefficient"] * np.log(smooth_l),
        color=LOG_COLOR,
        linewidth=1.8,
    )
    ax.plot(
        smooth_x,
        log2_model["intercept"]
        + log2_model["coefficient"] * np.log(smooth_l) ** 2,
        color=LOG2_COLOR,
        linewidth=1.8,
        linestyle="--",
    )

    ticks = np.log2(rows["L"].to_numpy(dtype=float))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(int(value)) for value in rows["L"]], rotation=45, ha="right")
    ax.grid(True, color=GRID_COLOR, linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.set_title(LABELS[key], fontsize=10.5, fontweight="semibold")
    ax.set_xlabel("Box half-size L")
    ax.set_ylabel("Winding variance")

    if annotate:
        fit = fit_row(fits, key, kind)
        delta_bic = log2_model["bic"] - log_model["bic"]
        ax.text(
            0.03,
            0.96,
            f"p={fit['p']:.3f}\nΔBIC={delta_bic:+.2f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.86, "edgecolor": "none"},
        )


def legend_handles():
    return [
        Line2D([0], [0], marker="o", color=DATA_COLOR, linewidth=0, label="Estimate ± 1.96 SE", markersize=5),
        Line2D([0], [0], color=LOG_COLOR, linewidth=2, label="Fitted a + b log L"),
        Line2D([0], [0], color=LOG2_COLOR, linestyle="--", linewidth=2, label="Fitted a + b(log L)²"),
    ]


def scaling_grid(summary, fits, models, kind):
    fig, axes = plt.subplots(4, 4, figsize=(16, 14), sharey=True)
    axes = axes.ravel()
    for ax, key in zip(axes, ORDER, strict=False):
        plot_panel(ax, summary, fits, models, key, kind)
    axes[-1].axis("off")
    fig.suptitle(
        f"{kind.capitalize()} winding variance: logarithmic versus squared-log fits",
        fontsize=18,
        fontweight="semibold",
        y=0.988,
    )
    fig.legend(handles=legend_handles(), loc="upper center", bbox_to_anchor=(0.5, 0.963), ncol=3, frameon=False)
    fig.text(
        0.5,
        0.008,
        "Error bars are ±1.96 clustered/jackknife SE. ΔBIC = BIC(log²) − BIC(log); positive values favor log L.",
        ha="center",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0.02, 0.035, 0.99, 0.94))
    return fig


def exponent_forest(fits):
    fig, axes = plt.subplots(1, 2, figsize=(13, 8.5), sharey=True)
    y = np.arange(len(ORDER))
    for ax, kind in zip(axes, ("annealed", "quenched"), strict=True):
        values = [fit_row(fits, key, kind) for key in ORDER]
        p = np.array([row["p"] for row in values])
        low = np.array([row["p_bootstrap_ci_low"] for row in values])
        high = np.array([row["p_bootstrap_ci_high"] for row in values])
        ax.errorbar(
            p,
            y,
            xerr=np.vstack((p - low, high - p)),
            fmt="o",
            color=LOG_COLOR,
            ecolor="#697386",
            capsize=3,
            markersize=5,
        )
        ax.axvline(1, color="#18864b", linewidth=1.8, label="p=1 (log L)")
        ax.axvline(2, color="#b42318", linewidth=1.8, linestyle="--", label="p=2 ((log L)²)")
        ax.set_title(kind.capitalize(), fontweight="semibold")
        ax.set_xlabel("Fitted exponent p with 95% clustered-bootstrap CI")
        ax.set_xlim(min(-0.1, np.nanmin(low) - 0.1), max(2.15, np.nanmax(high) + 0.1))
        ax.grid(True, axis="x", color=GRID_COLOR, linewidth=0.7)
        ax.invert_yaxis()
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([LABELS[key] for key in ORDER])
    axes[1].tick_params(labelleft=False)
    axes[1].legend(loc="lower right", frameon=False)
    fig.suptitle("Scaling exponent estimates across all distributions", fontsize=17, fontweight="semibold")
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.95))
    return fig


def bic_figure(models):
    fig, axes = plt.subplots(1, 2, figsize=(13, 8.5), sharey=True)
    y = np.arange(len(ORDER))
    for ax, kind in zip(axes, ("annealed", "quenched"), strict=True):
        deltas = []
        for key in ORDER:
            log_model = model_row(models, key, kind, "a_plus_b_log_L")
            log2_model = model_row(models, key, kind, "a_plus_b_log_L_squared")
            deltas.append(log2_model["bic"] - log_model["bic"])
        colors = [LOG_COLOR if value >= 0 else LOG2_COLOR for value in deltas]
        ax.barh(y, deltas, color=colors, alpha=0.9)
        ax.axvline(0, color=DATA_COLOR, linewidth=1)
        ax.set_title(kind.capitalize(), fontweight="semibold")
        ax.set_xlabel("ΔBIC = BIC(log²) − BIC(log)")
        ax.grid(True, axis="x", color=GRID_COLOR, linewidth=0.7)
        ax.invert_yaxis()
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([LABELS[key] for key in ORDER])
    axes[1].tick_params(labelleft=False)
    fig.suptitle("Additive scaling-model comparison", fontsize=17, fontweight="semibold")
    fig.text(0.5, 0.015, "Positive values favor a + b log L; negative values favor a + b(log L)².", ha="center")
    fig.tight_layout(rect=(0.02, 0.04, 0.99, 0.95))
    return fig


def save_figure(fig, stem: Path):
    fig.savefig(stem.with_suffix(".png"), dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    detail_dir = args.output_dir / "by_distribution"
    detail_dir.mkdir(exist_ok=True)

    summary, fits, models = prepare(args.analysis_dir)
    figures = {
        "annealed_scaling_all_distributions": scaling_grid(summary, fits, models, "annealed"),
        "quenched_scaling_all_distributions": scaling_grid(summary, fits, models, "quenched"),
        "scaling_exponent_forest": exponent_forest(fits),
        "bic_model_comparison": bic_figure(models),
    }
    for name, fig in figures.items():
        save_figure(fig, args.output_dir / name)

    with PdfPages(args.output_dir / "professor_scaling_overview.pdf") as pdf:
        for fig in figures.values():
            pdf.savefig(fig, bbox_inches="tight", facecolor="white")

    with PdfPages(args.output_dir / "all_distributions_detailed.pdf") as pdf:
        for key in ORDER:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
            for ax, kind in zip(axes, ("annealed", "quenched"), strict=True):
                plot_panel(ax, summary, fits, models, key, kind)
                ax.set_title(f"{LABELS[key]} — {kind}", fontweight="semibold")
            fig.legend(handles=legend_handles(), loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=3, frameon=False)
            fig.tight_layout(rect=(0.01, 0.01, 0.99, 0.90))
            pdf.savefig(fig, bbox_inches="tight", facecolor="white")
            safe_name = LABELS[key].lower().replace(" ", "_").replace("(", "_").replace(")", "").replace("=", "_").replace(".", "p").replace("σ", "sigma").replace("α", "alpha")
            fig.savefig(detail_dir / f"{safe_name}.png", dpi=220, bbox_inches="tight", facecolor="white")
            plt.close(fig)

    for fig in figures.values():
        plt.close(fig)
    print(f"Wrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
