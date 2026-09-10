"""Exploratory environment-blocked analysis for local pilot output."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GROUP_COLUMNS = ["model", "gamma_shape", "L"]
PAIR_STAT_COLUMNS = [
    "mean_1",
    "mean_2",
    "variance_1",
    "variance_2",
    "connected_component",
    "disorder_covariance",
    "correlation",
]


def _pair_statistics(first: Iterable[float], second: Iterable[float]) -> dict[str, float | int]:
    x = np.asarray(list(first), dtype=np.float64)
    y = np.asarray(list(second), dtype=np.float64)
    if x.shape != y.shape or x.ndim != 1 or x.size == 0:
        raise ValueError("paired statistics require equally sized nonempty vectors")
    if x.size == 1:
        variance_x = variance_y = connected = covariance = correlation = np.nan
    else:
        variance_x = float(np.var(x, ddof=1))
        variance_y = float(np.var(y, ddof=1))
        connected = float(0.5 * np.var(x - y, ddof=1))
        covariance = float(np.cov(x, y, ddof=1)[0, 1])
        correlation = (
            float(np.corrcoef(x, y)[0, 1])
            if variance_x > 0 and variance_y > 0
            else np.nan
        )
    return {
        "environments": int(x.size),
        "mean_1": float(x.mean()),
        "mean_2": float(y.mean()),
        "variance_1": variance_x,
        "variance_2": variance_y,
        "connected_component": connected,
        "disorder_covariance": covariance,
        "correlation": correlation,
    }


def aggregate_environment_measurements(
    center: pd.DataFrame, spatial: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reduce repeated MCMC records to one pair of means per environment."""
    identity = [
        "master_seed",
        "environment_seed",
        "chain1_seed",
        "chain2_seed",
        "L",
        "model",
        "gamma_shape",
        "environment_index",
        "environment_id",
    ]
    center_environment = (
        center.groupby(identity, dropna=False, as_index=False)[["H1_center", "H2_center"]]
        .mean()
        .rename(columns={"H1_center": "H1_environment_mean", "H2_center": "H2_environment_mean"})
    )
    if spatial.empty:
        spatial_environment = spatial.copy()
    else:
        spatial_identity = identity + ["requested_rho", "r"]
        spatial_environment = (
            spatial.groupby(spatial_identity, dropna=False, as_index=False)[["DeltaH1", "DeltaH2"]]
            .mean()
            .rename(
                columns={
                    "DeltaH1": "DeltaH1_environment_mean",
                    "DeltaH2": "DeltaH2_environment_mean",
                }
            )
        )
    return center_environment, spatial_environment


def summarize_center(center_environment: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in center_environment.groupby(GROUP_COLUMNS, dropna=False, sort=True):
        rows.append(
            dict(zip(GROUP_COLUMNS, keys, strict=True))
            | _pair_statistics(group["H1_environment_mean"], group["H2_environment_mean"])
        )
    return pd.DataFrame(rows)


def summarize_spatial(spatial_environment: pd.DataFrame) -> pd.DataFrame:
    if spatial_environment.empty:
        return pd.DataFrame()
    group_columns = GROUP_COLUMNS + ["requested_rho", "r"]
    rows = []
    for keys, group in spatial_environment.groupby(group_columns, dropna=False, sort=True):
        rows.append(
            dict(zip(group_columns, keys, strict=True))
            | _pair_statistics(
                group["DeltaH1_environment_mean"], group["DeltaH2_environment_mean"]
            )
        )
    return pd.DataFrame(rows)


def _fit(x: np.ndarray, y: np.ndarray, degree: int) -> dict[str, float | int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    parameter_count = degree + 1
    if x.size < parameter_count:
        return {
            "n": int(x.size),
            "a": np.nan,
            "b": np.nan,
            "c": np.nan,
            "sse": np.nan,
            "aic": np.nan,
            "bic": np.nan,
        }
    design = np.column_stack([x**power for power in range(parameter_count)])
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    sse = float(residual @ residual)
    safe_sse = max(sse, np.finfo(float).tiny)
    aic = float(x.size * np.log(safe_sse / x.size) + 2 * parameter_count)
    bic = float(x.size * np.log(safe_sse / x.size) + parameter_count * np.log(x.size))
    return {
        "n": int(x.size),
        "a": float(coefficients[0]),
        "b": float(coefficients[1]),
        "c": float(coefficients[2]) if degree == 2 else 0.0,
        "sse": sse,
        "aic": aic,
        "bic": bic,
    }


def fit_center_models(center_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in center_summary.groupby(["model", "gamma_shape"], dropna=False, sort=True):
        x = np.log(group["L"].to_numpy(dtype=np.float64))
        y = group["disorder_covariance"].to_numpy(dtype=np.float64)
        for degree, name in ((1, "a+b_logL"), (2, "a+b_logL+c_logL2")):
            rows.append(
                {"model": keys[0], "gamma_shape": keys[1], "fit": name}
                | _fit(x, y, degree)
            )
    return pd.DataFrame(rows)


def fit_pooled_spatial_models(spatial_summary: pd.DataFrame) -> pd.DataFrame:
    """Fit rho-specific intercept/slopes and optionally one shared quadratic term."""
    if spatial_summary.empty:
        return pd.DataFrame()
    rows: list[dict[str, object]] = []
    for keys, group in spatial_summary.groupby(["model", "gamma_shape"], dropna=False, sort=True):
        clean = group[np.isfinite(group["disorder_covariance"])].copy()
        rho_values = sorted(float(value) for value in clean["requested_rho"].unique())
        if not rho_values:
            continue
        x = np.log(clean["L"].to_numpy(dtype=np.float64))
        y = clean["disorder_covariance"].to_numpy(dtype=np.float64)
        rho = clean["requested_rho"].to_numpy(dtype=np.float64)
        for quadratic, fit_name in ((False, "rho_intercept+slope"), (True, "rho_intercept+slope+shared_logL2")):
            columns = []
            labels = []
            for value in rho_values:
                indicator = (rho == value).astype(np.float64)
                columns.extend((indicator, indicator * x))
                labels.extend((f"a_rho_{value:g}", f"b_rho_{value:g}"))
            if quadratic:
                columns.append(x**2)
                labels.append("shared_c")
            design = np.column_stack(columns)
            if len(y) < design.shape[1]:
                coefficients = np.full(design.shape[1], np.nan)
                sse = aic = bic = np.nan
            else:
                coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
                residual = y - design @ coefficients
                sse = float(residual @ residual)
                safe_sse = max(sse, np.finfo(float).tiny)
                aic = float(len(y) * np.log(safe_sse / len(y)) + 2 * design.shape[1])
                bic = float(
                    len(y) * np.log(safe_sse / len(y)) + design.shape[1] * np.log(len(y))
                )
            coefficient_map = {label: float(value) for label, value in zip(labels, coefficients, strict=True)}
            rows.append(
                {
                    "model": keys[0],
                    "gamma_shape": keys[1],
                    "fit": fit_name,
                    "n": len(y),
                    "parameters": design.shape[1],
                    "shared_c": coefficient_map.get("shared_c", 0.0),
                    "sse": sse,
                    "aic": aic,
                    "bic": bic,
                    "coefficients_json": json.dumps(coefficient_map, sort_keys=True),
                }
            )
    return pd.DataFrame(rows)


def _bootstrap(
    center_environment: pd.DataFrame,
    spatial_environment: pd.DataFrame,
    replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resample whole environments, keeping both chains and all observables joined."""
    if replicates <= 0:
        return pd.DataFrame(), pd.DataFrame()
    rng = np.random.Generator(np.random.PCG64(seed))
    center_draws: list[pd.DataFrame] = []
    spatial_draws: list[pd.DataFrame] = []
    grouped = list(center_environment.groupby(GROUP_COLUMNS, dropna=False, sort=True))
    for replicate in range(replicates):
        center_parts = []
        spatial_parts = []
        for keys, center_group in grouped:
            center_group = center_group.reset_index(drop=True)
            chosen = rng.integers(0, len(center_group), size=len(center_group))
            sampled = center_group.iloc[chosen].copy()
            # Give repeated selections unique bootstrap IDs while preserving the
            # joint chain pair in each selected source environment.
            sampled["environment_id"] = [f"bootstrap-{replicate}-{i}" for i in range(len(sampled))]
            center_parts.append(sampled)
            if not spatial_environment.empty:
                lookup = spatial_environment[
                    (spatial_environment["model"] == keys[0])
                    & (spatial_environment["L"] == keys[2])
                    & (
                        spatial_environment["gamma_shape"].isna()
                        if pd.isna(keys[1])
                        else spatial_environment["gamma_shape"] == keys[1]
                    )
                ]
                choices = center_group.iloc[chosen][["environment_id"]].reset_index(
                    drop=True
                )
                choices = choices.rename(columns={"environment_id": "source_environment_id"})
                choices["bootstrap_environment_id"] = [
                    f"bootstrap-{replicate}-{index}" for index in range(len(choices))
                ]
                selected_spatial = choices.merge(
                    lookup,
                    left_on="source_environment_id",
                    right_on="environment_id",
                    how="left",
                    validate="many_to_many",
                )
                selected_spatial = selected_spatial.drop(
                    columns=["source_environment_id", "environment_id"]
                ).rename(columns={"bootstrap_environment_id": "environment_id"})
                spatial_parts.append(selected_spatial)
        center_summary = summarize_center(pd.concat(center_parts, ignore_index=True))
        center_summary["replicate"] = replicate
        center_draws.append(center_summary)
        if spatial_parts:
            spatial_summary = summarize_spatial(pd.concat(spatial_parts, ignore_index=True))
            spatial_summary["replicate"] = replicate
            spatial_draws.append(spatial_summary)
    return (
        pd.concat(center_draws, ignore_index=True),
        pd.concat(spatial_draws, ignore_index=True) if spatial_draws else pd.DataFrame(),
    )


def _bootstrap_intervals(
    center_draws: pd.DataFrame, spatial_draws: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    component_rows = []
    for observable, frame, extra in (
        ("center", center_draws, []),
        ("spatial", spatial_draws, ["requested_rho", "r"]),
    ):
        if frame.empty:
            continue
        columns = GROUP_COLUMNS + extra
        for keys, group in frame.groupby(columns, dropna=False, sort=True):
            base = dict(zip(columns, keys if isinstance(keys, tuple) else (keys,), strict=True))
            for statistic in ("connected_component", "disorder_covariance"):
                values = group[statistic].to_numpy(dtype=np.float64)
                values = values[np.isfinite(values)]
                component_rows.append(
                    base
                    | {
                        "observable": observable,
                        "statistic": statistic,
                        "bootstrap_replicates": len(values),
                        "ci_lower": float(np.quantile(values, 0.025)) if len(values) else np.nan,
                        "ci_upper": float(np.quantile(values, 0.975)) if len(values) else np.nan,
                        "probability_positive": float(np.mean(values > 0)) if len(values) else np.nan,
                    }
                )

    fit_rows = []
    for replicate, frame in center_draws.groupby("replicate", sort=True):
        fitted = fit_center_models(frame)
        if not fitted.empty:
            fitted["replicate"] = replicate
            fit_rows.append(fitted)
    if spatial_draws.empty:
        pooled_draws = pd.DataFrame()
    else:
        pooled_parts = []
        for replicate, frame in spatial_draws.groupby("replicate", sort=True):
            fitted = fit_pooled_spatial_models(frame)
            if not fitted.empty:
                fitted["replicate"] = replicate
                pooled_parts.append(fitted)
        pooled_draws = pd.concat(pooled_parts, ignore_index=True) if pooled_parts else pd.DataFrame()

    fit_draws = pd.concat(fit_rows, ignore_index=True) if fit_rows else pd.DataFrame()
    interval_rows = []
    if not fit_draws.empty:
        for keys, group in fit_draws.groupby(["model", "gamma_shape", "fit"], dropna=False):
            for coefficient in ("a", "b", "c"):
                values = group[coefficient].to_numpy(dtype=np.float64)
                values = values[np.isfinite(values)]
                interval_rows.append(
                    {
                        "observable": "center",
                        "model": keys[0],
                        "gamma_shape": keys[1],
                        "fit": keys[2],
                        "coefficient": coefficient,
                        "bootstrap_replicates": len(values),
                        "ci_lower": float(np.quantile(values, 0.025)) if len(values) else np.nan,
                        "ci_upper": float(np.quantile(values, 0.975)) if len(values) else np.nan,
                        "probability_positive": float(np.mean(values > 0)) if len(values) else np.nan,
                    }
                )
    if not pooled_draws.empty:
        for keys, group in pooled_draws.groupby(["model", "gamma_shape", "fit"], dropna=False):
            values = group["shared_c"].to_numpy(dtype=np.float64)
            values = values[np.isfinite(values)]
            interval_rows.append(
                {
                    "observable": "spatial_pooled",
                    "model": keys[0],
                    "gamma_shape": keys[1],
                    "fit": keys[2],
                    "coefficient": "shared_c",
                    "bootstrap_replicates": len(values),
                    "ci_lower": float(np.quantile(values, 0.025)) if len(values) else np.nan,
                    "ci_upper": float(np.quantile(values, 0.975)) if len(values) else np.nan,
                    "probability_positive": float(np.mean(values > 0)) if len(values) else np.nan,
                }
            )
    return pd.DataFrame(component_rows), pd.DataFrame(interval_rows)


def _series_label(model: str, gamma_shape: float) -> str:
    return model if pd.isna(gamma_shape) else f"{model}(k={gamma_shape:g})"


def _plot_outputs(
    center_summary: pd.DataFrame,
    center_fits: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    output: Path,
) -> None:
    figure, axis = plt.subplots(figsize=(7, 4.5))
    for keys, group in center_summary.groupby(["model", "gamma_shape"], dropna=False):
        ordered = group.sort_values("L")
        label = _series_label(keys[0], keys[1])
        x = np.log(ordered["L"].to_numpy(dtype=float))
        axis.plot(x, ordered["disorder_covariance"], "o", label=f"{label} data")
        fit_group = center_fits[
            (center_fits["model"] == keys[0])
            & (center_fits["gamma_shape"].isna() if pd.isna(keys[1]) else center_fits["gamma_shape"] == keys[1])
        ]
        grid = np.linspace(x.min(), x.max(), 100)
        for _, fit in fit_group.iterrows():
            if np.isfinite(fit["a"]):
                prediction = fit["a"] + fit["b"] * grid + fit["c"] * grid**2
                axis.plot(grid, prediction, label=f"{label} {fit['fit']}")
    axis.axhline(0, color="black", linewidth=0.7)
    axis.set(xlabel="log L", ylabel="Cov(H1, H2)", title="Exploratory disorder covariance fits")
    axis.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(output / "center_covariance_fits.png", dpi=160)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7, 4.5))
    for keys, group in center_summary.groupby(["model", "gamma_shape"], dropna=False):
        ordered = group.sort_values("L")
        axis.plot(
            np.log(ordered["L"]),
            ordered["connected_component"],
            "o-",
            label=_series_label(keys[0], keys[1]),
        )
    axis.set(xlabel="log L", ylabel="0.5 Var(H1-H2)", title="Connected component")
    axis.legend()
    figure.tight_layout()
    figure.savefig(output / "connected_component.png", dpi=160)
    plt.close(figure)

    uniform = center_summary[center_summary["model"] == "uniform"].sort_values("L")
    figure, axis = plt.subplots(figsize=(7, 4.5))
    if not uniform.empty:
        axis.plot(np.log(uniform["L"]), uniform["disorder_covariance"], "o-")
    axis.axhline(0, color="black", linewidth=0.7)
    axis.set(xlabel="log L", ylabel="Cov(H1, H2)", title="Uniform negative control")
    figure.tight_layout()
    figure.savefig(output / "uniform_control_covariance.png", dpi=160)
    plt.close(figure)

    if not spatial_summary.empty:
        figure, axis = plt.subplots(figsize=(7, 4.5))
        for keys, group in spatial_summary.groupby(
            ["model", "gamma_shape", "requested_rho"], dropna=False
        ):
            ordered = group.sort_values("L")
            axis.plot(
                np.log(ordered["L"]),
                ordered["disorder_covariance"],
                "o-",
                label=f"{_series_label(keys[0], keys[1])}, rho={keys[2]:g}",
            )
        axis.axhline(0, color="black", linewidth=0.7)
        axis.set(xlabel="log L", ylabel="Cov(DeltaH1, DeltaH2)", title="Spatial disorder covariance")
        axis.legend(fontsize=6)
        figure.tight_layout()
        figure.savefig(output / "spatial_covariance.png", dpi=160)
        plt.close(figure)


def analyze_pilot(
    input_path: str | Path,
    output_path: str | Path,
    *,
    bootstrap_replicates: int = 500,
    bootstrap_seed: int = 20260908,
) -> dict[str, object]:
    """Analyze a pilot while treating ``environment_id`` as the sampling unit."""
    input_directory = Path(input_path)
    output_directory = Path(output_path)
    campaign_summary_path = input_directory / "campaign_summary.json"
    if campaign_summary_path.exists():
        campaign_summary = json.loads(campaign_summary_path.read_text())
        if not campaign_summary.get("all_pairs_cftp_certified", False):
            raise ValueError(
                "refusing scientific analysis of an incomplete perfect-sampling "
                "campaign; extend every missing/censored environment first"
            )
    output_directory.mkdir(parents=True, exist_ok=True)
    center = pd.read_csv(input_directory / "center_measurements.csv")
    spatial_file = input_directory / "spatial_measurements.csv"
    spatial = pd.read_csv(spatial_file) if spatial_file.exists() and spatial_file.stat().st_size else pd.DataFrame()
    diagnostics = pd.read_csv(input_directory / "chain_diagnostics.csv")
    center_environment, spatial_environment = aggregate_environment_measurements(center, spatial)
    center_summary = summarize_center(center_environment)
    spatial_summary = summarize_spatial(spatial_environment)
    center_fits = fit_center_models(center_summary)
    pooled_fits = fit_pooled_spatial_models(spatial_summary)
    center_cutoff_parts = []
    spatial_cutoff_parts = []
    for minimum_L in (8, 16, 20):
        center_cutoff = fit_center_models(
            center_summary.loc[center_summary["L"] >= minimum_L]
        )
        center_cutoff["minimum_L"] = minimum_L
        center_cutoff_parts.append(center_cutoff)
        spatial_cutoff = fit_pooled_spatial_models(
            spatial_summary.loc[spatial_summary["L"] >= minimum_L]
        )
        spatial_cutoff["minimum_L"] = minimum_L
        spatial_cutoff_parts.append(spatial_cutoff)
    center_cutoff_fits = pd.concat(center_cutoff_parts, ignore_index=True)
    spatial_cutoff_fits = pd.concat(spatial_cutoff_parts, ignore_index=True)

    diagnostic_numeric = [
        column
        for column in diagnostics.columns
        if (column.startswith("chain1_") or column.startswith("chain2_"))
        and pd.api.types.is_numeric_dtype(diagnostics[column])
    ]
    diagnostic_summary = (
        diagnostics.groupby(GROUP_COLUMNS, dropna=False, as_index=False)[diagnostic_numeric].mean()
    )
    center_draws, spatial_draws = _bootstrap(
        center_environment, spatial_environment, bootstrap_replicates, bootstrap_seed
    )
    component_intervals, fit_intervals = _bootstrap_intervals(center_draws, spatial_draws)

    center_environment.to_csv(output_directory / "environment_center_means.csv", index=False)
    spatial_environment.to_csv(output_directory / "environment_spatial_means.csv", index=False)
    center_summary.to_csv(output_directory / "center_summary.csv", index=False)
    spatial_summary.to_csv(output_directory / "spatial_summary.csv", index=False)
    diagnostic_summary.to_csv(output_directory / "diagnostic_summary.csv", index=False)
    center_fits.to_csv(output_directory / "center_fit_comparison.csv", index=False)
    pooled_fits.to_csv(output_directory / "spatial_pooled_fit_comparison.csv", index=False)
    center_cutoff_fits.to_csv(
        output_directory / "center_fit_cutoff_sensitivity.csv", index=False
    )
    spatial_cutoff_fits.to_csv(
        output_directory / "spatial_pooled_fit_cutoff_sensitivity.csv", index=False
    )
    component_intervals.to_csv(output_directory / "bootstrap_component_intervals.csv", index=False)
    fit_intervals.to_csv(output_directory / "bootstrap_fit_intervals.csv", index=False)
    _plot_outputs(center_summary, center_fits, spatial_summary, output_directory)

    report = {
        "exploratory": True,
        "warning": "Small sizes or agreement between short traces do not establish mixing or a (log L)^2 term.",
        "statistical_unit": "whole environment after within-environment chain averaging",
        "bootstrap": "whole environments resampled jointly across chain replicas and all observables",
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": bootstrap_seed,
        "center_environment_rows": len(center_environment),
        "spatial_environment_rows": len(spatial_environment),
    }
    (output_directory / "analysis_metadata.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
