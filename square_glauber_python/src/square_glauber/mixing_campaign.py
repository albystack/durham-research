"""Output and orchestration for local extremal-start mixing campaigns."""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .diagnostics import MixingComparison, MixingChainTrace, autocorrelation, compare_extremal_starts
from .environment import sample_environment
from .experiment import derive_seeds
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, spatial_separations, validate_L


def trace_frame(trace: MixingChainTrace) -> pd.DataFrame:
    data: dict[str, object] = {
        "start": trace.start,
        "sweep": trace.sweeps,
        "center_height": trace.center_heights,
        "height_sum": trace.height_sums,
    }
    for r, values in trace.spatial_increments.items():
        data[f"DeltaH_r{r}"] = values
    return pd.DataFrame(data)


def acf_frame(
    trace: MixingChainTrace, record_every_sweeps: int, max_lag_records: int
) -> pd.DataFrame:
    rows = []
    for observable, values in trace.observable_series().items():
        acf = autocorrelation(values, max_lag=max_lag_records)
        for lag, value in enumerate(acf):
            rows.append(
                {
                    "start": trace.start,
                    "observable": observable,
                    "lag_records": lag,
                    "lag_sweeps": lag * record_every_sweeps,
                    "acf": float(value),
                }
            )
    return pd.DataFrame(rows)


def plot_comparison(
    comparison: MixingComparison, output: str | Path, title: str
) -> None:
    frames = [trace_frame(comparison.horizontal_start), trace_frame(comparison.vertical_start)]
    traces = pd.concat(frames, ignore_index=True)
    observable_columns = [
        "center_height",
        "height_sum",
        *[column for column in traces if column.startswith("DeltaH_")],
    ]
    figure, axes = plt.subplots(
        len(observable_columns),
        1,
        figsize=(9, max(5.5, 2.25 * len(observable_columns))),
        sharex=True,
        squeeze=False,
    )
    for axis, column in zip(axes[:, 0], observable_columns, strict=True):
        for start, group in traces.groupby("start"):
            axis.plot(
                group["sweep"],
                group[column],
                linewidth=0.45,
                alpha=0.35,
                label=f"{start} raw",
            )
            rolling_points = max(5, len(group) // 200)
            rolling = group[column].rolling(rolling_points, center=True).mean()
            axis.plot(
                group["sweep"],
                rolling,
                linewidth=1.1,
                alpha=0.95,
                label=f"{start} rolling mean",
            )
        axis.axvline(
            float(traces["sweep"].max()) / 2.0,
            color="black",
            linestyle="--",
            linewidth=0.7,
            alpha=0.65,
            label="post-transient start",
        )
        axis.set_ylabel(column)
        axis.legend(fontsize=7, loc="upper right")
    axes[-1, 0].set_xlabel("attempted-update sweeps")
    figure.suptitle(title + "\nDiagnostic only — trace agreement is not proof of convergence")
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def _finite_or_none(value: object) -> object:
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _finite_or_none(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_finite_or_none(item) for item in value]
    return value


def write_mixing_case(
    comparison: MixingComparison,
    output: str | Path,
    *,
    metadata: dict[str, object],
    max_lag_records: int = 200,
    late_fraction: float = 0.5,
    window_count: int = 4,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Write one two-start case and return flat case/observable summaries."""
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    trace = pd.concat(
        [trace_frame(comparison.horizontal_start), trace_frame(comparison.vertical_start)],
        ignore_index=True,
    )
    acf = pd.concat(
        [
            acf_frame(
                comparison.horizontal_start,
                comparison.record_every_sweeps,
                max_lag_records,
            ),
            acf_frame(
                comparison.vertical_start,
                comparison.record_every_sweeps,
                max_lag_records,
            ),
        ],
        ignore_index=True,
    )
    trace.to_csv(output_path / "traces.csv.gz", index=False, compression="gzip")
    acf.to_csv(output_path / "acf.csv.gz", index=False, compression="gzip")
    plot_comparison(comparison, output_path / "traces.png", str(metadata["environment_id"]))
    summary = comparison.summary(
        late_fraction=late_fraction, window_count=window_count
    ) | metadata
    (output_path / "summary.json").write_text(
        json.dumps(_finite_or_none(summary), indent=2) + "\n"
    )

    window_rows = []
    for start_name in ("horizontal_start", "vertical_start"):
        for observable, details in summary[start_name]["observables"].items():
            for window in details["post_transient_windows"]:
                window_rows.append(
                    metadata
                    | {
                        "start": summary[start_name]["start"],
                        "observable": observable,
                        **window,
                    }
                )
    pd.DataFrame(window_rows).to_csv(output_path / "window_summary.csv", index=False)

    horizontal = summary["horizontal_start"]
    vertical = summary["vertical_start"]
    case_row = metadata | {
        "record_every_sweeps": comparison.record_every_sweeps,
        "trace_points": horizontal["trace_points"],
        "horizontal_elapsed_seconds": horizontal["elapsed_seconds"],
        "vertical_elapsed_seconds": vertical["elapsed_seconds"],
        "total_chain_seconds": horizontal["elapsed_seconds"] + vertical["elapsed_seconds"],
        "horizontal_flippable_proposal_rate": horizontal["flippable_proposal_rate"],
        "vertical_flippable_proposal_rate": vertical["flippable_proposal_rate"],
        "horizontal_actual_change_rate": horizontal["actual_change_rate"],
        "vertical_actual_change_rate": vertical["actual_change_rate"],
        "center_split_rhat_late_half": summary["split_rhat_center_late_half"],
        "height_sum_split_rhat_late_half": summary["split_rhat_height_sum_late_half"],
    }
    observable_rows = []
    for observable, details in summary["observables"].items():
        first = details["horizontal_start"]
        second = details["vertical_start"]
        numerical_screen = (
            np.isfinite(details["split_rhat_late_half"])
            and details["split_rhat_late_half"] <= 1.01
            and details["absolute_standardized_late_mean_gap"] <= 2.0
            and min(first["effective_sample_size_late"], second["effective_sample_size_late"])
            >= 100
        )
        observable_rows.append(
            metadata
            | {
                "observable": observable,
                "horizontal_late_mean": first["late_mean"],
                "vertical_late_mean": second["late_mean"],
                "late_mean_gap_horizontal_minus_vertical": details[
                    "late_mean_gap_horizontal_minus_vertical"
                ],
                "horizontal_late_variance": first["late_variance"],
                "vertical_late_variance": second["late_variance"],
                "late_variance_difference_horizontal_minus_vertical": details[
                    "late_variance_difference_horizontal_minus_vertical"
                ],
                "horizontal_lag1_acf": first["late_lag1_autocorrelation"],
                "vertical_lag1_acf": second["late_lag1_autocorrelation"],
                "horizontal_iact_records": first[
                    "integrated_autocorrelation_time_records"
                ],
                "vertical_iact_records": second[
                    "integrated_autocorrelation_time_records"
                ],
                "horizontal_iact_sweeps": first[
                    "integrated_autocorrelation_time_sweeps"
                ],
                "vertical_iact_sweeps": second[
                    "integrated_autocorrelation_time_sweeps"
                ],
                "horizontal_effective_sample_size_late": first[
                    "effective_sample_size_late"
                ],
                "vertical_effective_sample_size_late": second[
                    "effective_sample_size_late"
                ],
                "horizontal_late_mean_mcse": first["late_mean_mcse"],
                "vertical_late_mean_mcse": second["late_mean_mcse"],
                "z_gap_horizontal_minus_vertical": details[
                    "z_gap_horizontal_minus_vertical"
                ],
                "absolute_standardized_late_mean_gap": details[
                    "absolute_standardized_late_mean_gap"
                ],
                "split_rhat_late_half": details["split_rhat_late_half"],
                "conservative_numerical_screen": (
                    "pass_requires_visual_review" if numerical_screen else "review"
                ),
            }
        )
    return case_row, observable_rows


def run_mixing_campaign(
    *,
    sizes: Iterable[int],
    gamma_environments: int,
    uniform_replicates: int,
    gamma_shape: float,
    sweeps: int,
    record_every: int,
    master_seed: int,
    output: str | Path,
    relative_separations: Iterable[float] = DEFAULT_RELATIVE_SEPARATIONS,
    max_lag_records: int = 200,
) -> dict[str, object]:
    """Run the local uniform/Gamma extremal-start diagnostic campaign."""
    sizes = tuple(validate_L(int(L)) for L in sizes)
    relative_separations = tuple(float(value) for value in relative_separations)
    if gamma_environments <= 0 or uniform_replicates <= 0:
        raise ValueError("Gamma environments and uniform replicates must be positive")
    if sweeps <= 0 or record_every <= 0:
        raise ValueError("sweeps and record_every must be positive")
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    case_rows: list[dict[str, object]] = []
    observable_rows: list[dict[str, object]] = []
    started = perf_counter()
    for L in sizes:
        separations = spatial_separations(L, relative_separations)
        for model, case_count in (("uniform", uniform_replicates), ("gamma", gamma_environments)):
            for environment_index in range(case_count):
                seeds = derive_seeds(master_seed, L, model, environment_index)
                environment = sample_environment(
                    L,
                    np.random.Generator(np.random.PCG64(seeds.environment_seed)),
                    model,
                    gamma_shape=gamma_shape,
                )
                comparison = compare_extremal_starts(
                    environment,
                    np.random.Generator(np.random.PCG64(seeds.chain1_seed)),
                    np.random.Generator(np.random.PCG64(seeds.chain2_seed)),
                    sweeps=sweeps,
                    record_every=record_every,
                    separations=separations,
                )
                environment_id = f"{model}-L{L}-diagnostic-{environment_index:02d}"
                case_directory = output_path / environment_id
                case_directory.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(
                    case_directory / "environment_weights.npz",
                    horizontal_weights=environment.horizontal_weights,
                    vertical_weights=environment.vertical_weights,
                )
                metadata = {
                    "environment_id": environment_id,
                    "L": L,
                    "model": model,
                    "gamma_shape": gamma_shape if model == "gamma" else None,
                    "environment_index": environment_index,
                    **asdict(seeds),
                    "sweeps": sweeps,
                    "actual_separations": json.dumps(
                        [asdict(value) for value in separations], sort_keys=True
                    ),
                }
                case_row, rows = write_mixing_case(
                    comparison,
                    case_directory,
                    metadata=metadata,
                    max_lag_records=max_lag_records,
                )
                case_rows.append(case_row)
                observable_rows.extend(rows)
                pd.DataFrame(case_rows).to_csv(output_path / "case_summary.csv", index=False)
                pd.DataFrame(observable_rows).to_csv(
                    output_path / "observable_summary.csv", index=False
                )

    case_frame = pd.DataFrame(case_rows)
    observable_frame = pd.DataFrame(observable_rows)
    runtime_summary = (
        case_frame.groupby(["L", "model"], dropna=False, as_index=False)
        .agg(
            cases=("environment_id", "count"),
            total_chain_seconds=("total_chain_seconds", "sum"),
            mean_chain_pair_seconds=("total_chain_seconds", "mean"),
            mean_horizontal_flippable_rate=("horizontal_flippable_proposal_rate", "mean"),
            mean_vertical_flippable_rate=("vertical_flippable_proposal_rate", "mean"),
            mean_horizontal_change_rate=("horizontal_actual_change_rate", "mean"),
            mean_vertical_change_rate=("vertical_actual_change_rate", "mean"),
        )
    )
    runtime_summary.to_csv(output_path / "runtime_summary.csv", index=False)
    campaign = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sizes": list(sizes),
        "models": ["uniform", "gamma"],
        "gamma_shape": gamma_shape,
        "gamma_environments_per_size": gamma_environments,
        "uniform_replicates_per_size": uniform_replicates,
        "sweeps_per_chain": sweeps,
        "record_every_sweeps": record_every,
        "master_seed": master_seed,
        "cases": len(case_frame),
        "wall_seconds": perf_counter() - started,
        "warning": "Diagnostics only; numerical screens and trace overlap do not prove convergence.",
    }
    (output_path / "campaign_metadata.json").write_text(json.dumps(campaign, indent=2) + "\n")
    # Re-write final tables after every case was completed.
    case_frame.to_csv(output_path / "case_summary.csv", index=False)
    observable_frame.to_csv(output_path / "observable_summary.csv", index=False)
    return campaign
