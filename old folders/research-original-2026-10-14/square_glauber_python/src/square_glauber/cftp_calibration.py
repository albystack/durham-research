"""Bounded cost/censoring calibration for the certified CFTP sampler.

This module deliberately records computational certification horizons only.  It
does not estimate disorder covariance or fit any size-scaling model.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .environment import DimerEnvironment, sample_environment
from .experiment import _seed_uint64
from .height import center_height, face_heights
from .mixing_extension import environment_hardness, environment_hashes
from .perfect import CFTPNotCoalesced, MonotoneCFTPSampler, PerfectBackend


def _derived_seed(*parts: int) -> int:
    return _seed_uint64(np.random.SeedSequence([int(value) for value in parts]))


def _save_and_verify_environment(
    environment: DimerEnvironment,
    path: Path,
    *,
    environment_seed: int,
) -> dict[str, object]:
    """Persist an environment and require bitwise reload/regeneration equality."""
    np.savez_compressed(
        path,
        horizontal_weights=environment.horizontal_weights,
        vertical_weights=environment.vertical_weights,
    )
    with np.load(path, allow_pickle=False) as archive:
        stored_horizontal = np.asarray(archive["horizontal_weights"])
        stored_vertical = np.asarray(archive["vertical_weights"])
    regenerated = sample_environment(
        environment.L,
        np.random.Generator(np.random.PCG64(int(environment_seed))),
        environment.model,
        gamma_shape=(
            float(environment.gamma_shape)
            if environment.gamma_shape is not None
            else 1.0
        ),
    )
    reload_equal = bool(
        np.array_equal(stored_horizontal, environment.horizontal_weights)
        and np.array_equal(stored_vertical, environment.vertical_weights)
    )
    regeneration_equal = bool(
        np.array_equal(regenerated.horizontal_weights, environment.horizontal_weights)
        and np.array_equal(regenerated.vertical_weights, environment.vertical_weights)
    )
    if not reload_equal or not regeneration_equal:
        raise RuntimeError(f"environment verification failed for {path}")
    return {
        **environment_hashes(environment),
        "stored_npz": str(path),
        "reload_bitwise_equal": reload_equal,
        "seed_regeneration_bitwise_equal": regeneration_equal,
    }


def _group_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (phase, L, model), group in frame.groupby(
        ["phase", "L", "model"], sort=True
    ):
        horizons = group["coalescence_horizon_sweeps"].to_numpy(dtype=np.float64)
        completed = group.loc[group["certified"]]
        runtimes = group["runtime_seconds"].to_numpy(dtype=np.float64)
        attempts = group[
            ["lower_total_replayed_attempts", "upper_total_replayed_attempts"]
        ].sum(axis=1).to_numpy(dtype=np.float64)
        rows.append(
            {
                "phase": phase,
                "L": int(L),
                "model": model,
                "environment_count": int(group["environment_id"].nunique()),
                "sample_count": len(group),
                "certified_count": int(group["certified"].sum()),
                "censored_count": int(group["censored"].sum()),
                "censor_fraction": float(group["censored"].mean()),
                "horizon_min_sweeps": int(np.min(horizons)),
                "horizon_median_sweeps": float(np.median(horizons)),
                "horizon_q90_sweeps": float(np.quantile(horizons, 0.90)),
                "horizon_q95_sweeps": float(np.quantile(horizons, 0.95)),
                "horizon_max_sweeps": int(np.max(horizons)),
                "certified_horizon_median_sweeps": (
                    float(completed["coalescence_horizon_sweeps"].median())
                    if len(completed)
                    else np.nan
                ),
                "runtime_mean_seconds": float(np.mean(runtimes)),
                "runtime_median_seconds": float(np.median(runtimes)),
                "runtime_max_seconds": float(np.max(runtimes)),
                "runtime_total_seconds": float(np.sum(runtimes)),
                "replayed_attempts_per_second": float(
                    np.sum(attempts) / np.sum(runtimes)
                ),
            }
        )
    return pd.DataFrame(rows)


def _write_plots(samples: pd.DataFrame, output: Path) -> None:
    groups = list(
        samples.groupby(["phase", "L", "model"], sort=True)
    )
    figure, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=False)
    positions = np.arange(len(groups))
    labels: list[str] = []
    horizon_values: list[np.ndarray] = []
    runtime_values: list[np.ndarray] = []
    for (phase, L, model), group in groups:
        labels.append(f"{phase}\nL={L}\n{model}")
        horizon_values.append(
            group["coalescence_horizon_sweeps"].to_numpy(dtype=np.float64)
        )
        runtime_values.append(group["runtime_seconds"].to_numpy(dtype=np.float64))
    axes[0].boxplot(horizon_values, positions=positions, showfliers=True)
    axes[0].set_yscale("log", base=2)
    axes[0].set_ylabel("certifying/censoring horizon (sweeps)")
    axes[0].set_title("CFTP cost calibration (censored observations remain at cap)")
    axes[1].boxplot(runtime_values, positions=positions, showfliers=True)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("wall time per CFTP stream (seconds)")
    axes[1].set_xticks(positions, labels, rotation=30, ha="right")
    figure.tight_layout()
    figure.savefig(output / "cftp_horizon_and_runtime.png", dpi=170)
    plt.close(figure)


def run_cftp_cost_calibration(
    output: str | Path,
    *,
    standard_sizes: Iterable[int] = (8, 12, 16),
    standard_gamma_environments: int = 30,
    standard_uniform_environments: int = 1,
    standard_samples_per_environment: int = 2,
    probe_sizes: Iterable[int] = (24, 32),
    probe_gamma_environments: int = 5,
    probe_uniform_environments: int = 1,
    probe_samples_per_environment: int = 1,
    gamma_shape: float = 1.0,
    master_seed: int = 20260909,
    max_horizon_sweeps: int = 1_048_576,
    backend: PerfectBackend = "numba",
) -> dict[str, object]:
    """Run a bounded CFTP certification-cost study with fully frozen outputs."""
    standard_sizes = tuple(int(value) for value in standard_sizes)
    probe_sizes = tuple(int(value) for value in probe_sizes)
    if any(value <= 0 for value in standard_sizes + probe_sizes):
        raise ValueError("sizes must be positive")
    count_values = (
        standard_gamma_environments,
        standard_uniform_environments,
        standard_samples_per_environment,
        probe_gamma_environments,
        probe_uniform_environments,
        probe_samples_per_environment,
    )
    if any(value < 0 for value in count_values):
        raise ValueError("environment/sample counts must be nonnegative")
    if max_horizon_sweeps <= 0:
        raise ValueError("max_horizon_sweeps must be positive")

    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite nonempty output directory: {output_path}"
        )
    output_path.mkdir(parents=True, exist_ok=False)
    environment_directory = output_path / "environments"
    environment_directory.mkdir()

    design = (
        (
            "standard",
            standard_sizes,
            standard_gamma_environments,
            standard_uniform_environments,
            standard_samples_per_environment,
        ),
        (
            "probe",
            probe_sizes,
            probe_gamma_environments,
            probe_uniform_environments,
            probe_samples_per_environment,
        ),
    )
    environment_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    started = perf_counter()
    for phase, sizes, gamma_count, uniform_count, samples_per_environment in design:
        phase_code = 1 if phase == "probe" else 0
        for L in sizes:
            for model, environment_count in (
                ("uniform", uniform_count),
                ("gamma", gamma_count),
            ):
                model_code = 1 if model == "gamma" else 0
                for environment_index in range(environment_count):
                    environment_seed = _derived_seed(
                        master_seed,
                        phase_code,
                        L,
                        model_code,
                        environment_index,
                        0xE117,
                    )
                    environment = sample_environment(
                        L,
                        np.random.Generator(np.random.PCG64(environment_seed)),
                        model,
                        gamma_shape=gamma_shape,
                    )
                    environment_id = (
                        f"{phase}-{model}-L{L}-environment-{environment_index:03d}"
                    )
                    verification = _save_and_verify_environment(
                        environment,
                        environment_directory / f"{environment_id}.npz",
                        environment_seed=environment_seed,
                    )
                    hardness = environment_hardness(environment)
                    sampler = MonotoneCFTPSampler(environment, backend=backend)
                    minimum_heights = face_heights(sampler.minimum.state)
                    maximum_heights = face_heights(sampler.maximum.state)
                    environment_rows.append(
                        {
                            "phase": phase,
                            "environment_id": environment_id,
                            "environment_index": environment_index,
                            "L": L,
                            "model": model,
                            "gamma_shape": gamma_shape if model == "gamma" else np.nan,
                            "master_seed": master_seed,
                            "environment_seed": environment_seed,
                            **verification,
                            **hardness,
                            "minimum_directed_flips": sampler.minimum.directed_flips,
                            "maximum_directed_flips": sampler.maximum.directed_flips,
                            "minimum_height_sum": int(minimum_heights.sum()),
                            "maximum_height_sum": int(maximum_heights.sum()),
                        }
                    )
                    for sample_index in range(samples_per_environment):
                        cftp_seed = _derived_seed(
                            master_seed,
                            phase_code,
                            L,
                            model_code,
                            environment_index,
                            sample_index,
                            0xC0F7,
                        )
                        sample_started = perf_counter()
                        try:
                            result = sampler.sample(
                                cftp_seed,
                                max_horizon_sweeps=max_horizon_sweeps,
                                verify_order=True,
                            )
                            elapsed = perf_counter() - sample_started
                            sample_heights = face_heights(result.state)
                            lower_attempts = (
                                result.lower_update_counts.attempted_updates
                            )
                            upper_attempts = (
                                result.upper_update_counts.attempted_updates
                            )
                            row = {
                                "certified": True,
                                "censored": False,
                                "coalescence_horizon_sweeps": result.horizon_sweeps,
                                "doubling_iterations": len(result.iterations),
                                "maximum_height_gap_at_cap": 0,
                                "center_height": center_height(result.state),
                                "height_sum": int(sample_heights.sum()),
                                "lower_total_replayed_attempts": lower_attempts,
                                "upper_total_replayed_attempts": upper_attempts,
                                "runtime_seconds": elapsed,
                            }
                        except CFTPNotCoalesced as error:
                            elapsed = perf_counter() - sample_started
                            replay_sweeps = sum(
                                item.horizon_sweeps for item in error.iterations
                            )
                            attempts = replay_sweeps * (L - 1) ** 2
                            row = {
                                "certified": False,
                                "censored": True,
                                "coalescence_horizon_sweeps": error.horizon_sweeps,
                                "doubling_iterations": len(error.iterations),
                                "maximum_height_gap_at_cap": (
                                    error.iterations[-1]
                                    .maximum_pointwise_height_gap_at_time_zero
                                ),
                                "center_height": np.nan,
                                "height_sum": np.nan,
                                "lower_total_replayed_attempts": attempts,
                                "upper_total_replayed_attempts": attempts,
                                "runtime_seconds": elapsed,
                            }
                        sample_rows.append(
                            {
                                "phase": phase,
                                "environment_id": environment_id,
                                "environment_index": environment_index,
                                "sample_index": sample_index,
                                "L": L,
                                "model": model,
                                "gamma_shape": (
                                    gamma_shape if model == "gamma" else np.nan
                                ),
                                "master_seed": master_seed,
                                "environment_seed": environment_seed,
                                "cftp_seed": cftp_seed,
                                "environment_sha256": verification[
                                    "environment_sha256"
                                ],
                                "backend": backend,
                                "max_horizon_sweeps": max_horizon_sweeps,
                                **row,
                            }
                        )

    environments = pd.DataFrame(environment_rows)
    samples = pd.DataFrame(sample_rows)
    environments.to_csv(output_path / "environment_catalog.csv", index=False)
    samples.to_csv(output_path / "cftp_samples.csv", index=False)
    group_summary = _group_summary(samples)
    group_summary.to_csv(output_path / "group_summary.csv", index=False)
    censored = samples.loc[samples["censored"]].copy()
    censored.to_csv(output_path / "censored_samples_to_extend.csv", index=False)
    _write_plots(samples, output_path)

    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "bounded CFTP cost/censoring calibration; no scientific scaling fit",
        "backend": backend,
        "master_seed": master_seed,
        "gamma_shape": gamma_shape,
        "max_horizon_sweeps": max_horizon_sweeps,
        "standard_sizes": list(standard_sizes),
        "probe_sizes": list(probe_sizes),
        "environment_count": len(environments),
        "cftp_stream_count": len(samples),
        "certified_count": int(samples["certified"].sum()),
        "censored_count": int(samples["censored"].sum()),
        "all_environment_reloads_and_regenerations_match": bool(
            environments[
                ["reload_bitwise_equal", "seed_regeneration_bitwise_equal"]
            ].all(axis=None)
        ),
        "all_cftp_seeds_unique": bool(samples["cftp_seed"].is_unique),
        "all_environment_arrays_retained": bool(
            environments["stored_npz"].map(lambda value: Path(value).is_file()).all()
        ),
        "scientific_covariance_computed": False,
        "log_or_log_squared_fit_performed": False,
        "total_wall_seconds": perf_counter() - started,
        "group_summary": group_summary.to_dict(orient="records"),
        "next_action": (
            "extend every row in censored_samples_to_extend.csv with the identical "
            "environment and CFTP seed before any scientific use"
            if len(censored)
            else "no censoring at this cap; assess runtime tails before a paired pilot"
        ),
    }
    (output_path / "calibration_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def extend_censored_cftp_calibration(
    source: str | Path,
    output: str | Path,
    *,
    max_horizon_sweeps: int,
    backend: PerfectBackend = "numba",
) -> dict[str, object]:
    """Rerun only censored streams with identical arrays/seeds and a larger cap."""
    source_path = Path(source)
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite nonempty output directory: {output_path}"
        )
    output_path.mkdir(parents=True, exist_ok=False)
    old_rows = pd.read_csv(source_path / "censored_samples_to_extend.csv")
    if not len(old_rows):
        raise ValueError("source calibration has no censored streams to extend")
    old_cap = int(old_rows["max_horizon_sweeps"].max())
    if max_horizon_sweeps <= old_cap:
        raise ValueError("extension horizon must be larger than the source cap")
    environment_catalog = pd.read_csv(source_path / "environment_catalog.csv").set_index(
        "environment_id"
    )
    rows: list[dict[str, object]] = []
    started = perf_counter()
    for old in old_rows.itertuples(index=False):
        catalog = environment_catalog.loc[old.environment_id]
        archive_path = source_path / "environments" / f"{old.environment_id}.npz"
        with np.load(archive_path, allow_pickle=False) as archive:
            horizontal = np.asarray(archive["horizontal_weights"])
            vertical = np.asarray(archive["vertical_weights"])
        environment = DimerEnvironment(
            int(old.L),
            horizontal,
            vertical,
            model=str(old.model),
            gamma_shape=(float(old.gamma_shape) if str(old.model) == "gamma" else None),
        )
        hashes = environment_hashes(environment)
        if hashes["environment_sha256"] != str(old.environment_sha256):
            raise RuntimeError(f"saved environment hash changed for {old.environment_id}")
        regenerated = sample_environment(
            int(old.L),
            np.random.Generator(np.random.PCG64(int(old.environment_seed))),
            str(old.model),
            gamma_shape=(float(old.gamma_shape) if str(old.model) == "gamma" else 1.0),
        )
        if not (
            np.array_equal(environment.horizontal_weights, regenerated.horizontal_weights)
            and np.array_equal(environment.vertical_weights, regenerated.vertical_weights)
        ):
            raise RuntimeError(f"seed regeneration changed for {old.environment_id}")
        if str(catalog["environment_sha256"]) != hashes["environment_sha256"]:
            raise RuntimeError(f"catalog hash changed for {old.environment_id}")

        sampler = MonotoneCFTPSampler(environment, backend=backend)
        sample_started = perf_counter()
        result = None
        failure = None
        try:
            result = sampler.sample(
                int(old.cftp_seed),
                max_horizon_sweeps=max_horizon_sweeps,
                verify_order=True,
            )
            iterations = result.iterations
            certified = True
            horizon = result.horizon_sweeps
            lower_attempts = result.lower_update_counts.attempted_updates
            upper_attempts = result.upper_update_counts.attempted_updates
            heights = face_heights(result.state)
            center = center_height(result.state)
            height_sum = int(heights.sum())
            final_gap = 0
        except CFTPNotCoalesced as error:
            failure = error
            iterations = error.iterations
            certified = False
            horizon = error.horizon_sweeps
            replay_sweeps = sum(item.horizon_sweeps for item in iterations)
            lower_attempts = replay_sweeps * (int(old.L) - 1) ** 2
            upper_attempts = lower_attempts
            center = np.nan
            height_sum = np.nan
            final_gap = iterations[-1].maximum_pointwise_height_gap_at_time_zero
        elapsed = perf_counter() - sample_started
        old_iteration = next(
            (item for item in iterations if item.horizon_sweeps == int(old.max_horizon_sweeps)),
            None,
        )
        if old_iteration is None:
            raise RuntimeError(
                f"extended history lacks the old cap for {old.environment_id}"
            )
        if old_iteration.lower_upper_coalesced_at_time_zero:
            raise RuntimeError(
                f"stream unexpectedly coalesced at its former cap for {old.environment_id}"
            )
        if (
            old_iteration.maximum_pointwise_height_gap_at_time_zero
            != int(old.maximum_height_gap_at_cap)
        ):
            raise RuntimeError(
                f"replayed old-cap gap differs for {old.environment_id}"
            )
        rows.append(
            {
                "phase": old.phase,
                "environment_id": old.environment_id,
                "L": int(old.L),
                "model": old.model,
                "environment_seed": int(old.environment_seed),
                "cftp_seed": int(old.cftp_seed),
                "environment_sha256": hashes["environment_sha256"],
                "source_cap_sweeps": int(old.max_horizon_sweeps),
                "source_gap_replayed_exactly": True,
                "extension_cap_sweeps": max_horizon_sweeps,
                "certified": certified,
                "censored": not certified,
                "coalescence_horizon_sweeps": horizon,
                "doubling_iterations": len(iterations),
                "maximum_height_gap_at_extension_cap": final_gap,
                "center_height": center,
                "height_sum": height_sum,
                "lower_total_replayed_attempts": lower_attempts,
                "upper_total_replayed_attempts": upper_attempts,
                "runtime_seconds": elapsed,
                "backend": backend,
            }
        )
        del failure
    frame = pd.DataFrame(rows)
    frame.to_csv(output_path / "censored_stream_extensions.csv", index=False)
    remaining = frame.loc[frame["censored"]]
    remaining.to_csv(output_path / "still_censored_to_extend.csv", index=False)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(source_path),
        "backend": backend,
        "source_cap_sweeps": old_cap,
        "extension_cap_sweeps": max_horizon_sweeps,
        "streams_extended": len(frame),
        "certified_after_extension": int(frame["certified"].sum()),
        "still_censored": int(frame["censored"].sum()),
        "all_environment_hashes_match": True,
        "all_seed_regenerations_match": True,
        "all_old_cap_gaps_replayed_exactly": bool(
            frame["source_gap_replayed_exactly"].all()
        ),
        "total_wall_seconds": perf_counter() - started,
        "scientific_covariance_computed": False,
        "log_or_log_squared_fit_performed": False,
        "next_action": (
            "extend remaining saved streams again without changing their histories"
            if len(remaining)
            else "all initially censored streams now have CFTP certificates"
        ),
    }
    (output_path / "extension_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def consolidate_cftp_cost_calibration(
    source: str | Path,
    extension: str | Path,
    output: str | Path,
) -> dict[str, object]:
    """Combine capped and extended streams into one authoritative cost table."""
    source_path = Path(source)
    extension_path = Path(extension)
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite nonempty output directory: {output_path}"
        )
    output_path.mkdir(parents=True, exist_ok=False)
    base = pd.read_csv(source_path / "cftp_samples.csv")
    extended = pd.read_csv(extension_path / "censored_stream_extensions.csv")
    final = base.copy()
    for column in (
        "certified",
        "censored",
        "coalescence_horizon_sweeps",
        "runtime_seconds",
        "lower_total_replayed_attempts",
        "upper_total_replayed_attempts",
    ):
        final[f"initial_{column}"] = final[column]
    final["extension_runtime_seconds"] = 0.0
    final["calibration_compute_spent_seconds"] = final["runtime_seconds"]
    for item in extended.itertuples(index=False):
        mask = (final["environment_id"] == item.environment_id) & (
            final["cftp_seed"].astype(str) == str(item.cftp_seed)
        )
        if int(mask.sum()) != 1:
            raise RuntimeError(f"cannot uniquely match extension {item.environment_id}")
        if not bool(final.loc[mask, "initial_censored"].iloc[0]):
            raise RuntimeError(f"extension did not correspond to a censored row")
        for column in (
            "certified",
            "censored",
            "coalescence_horizon_sweeps",
            "runtime_seconds",
            "lower_total_replayed_attempts",
            "upper_total_replayed_attempts",
            "center_height",
            "height_sum",
        ):
            final.loc[mask, column] = getattr(item, column)
        final.loc[mask, "extension_runtime_seconds"] = item.runtime_seconds
        final.loc[mask, "calibration_compute_spent_seconds"] = (
            final.loc[mask, "initial_runtime_seconds"] + item.runtime_seconds
        )
    if not bool(final["certified"].all()):
        raise RuntimeError("consolidated table still contains uncertified streams")
    final.to_csv(output_path / "final_cftp_samples.csv", index=False)
    groups = _group_summary(final)
    groups.to_csv(output_path / "final_group_summary.csv", index=False)
    environment_costs = (
        final.groupby(["phase", "L", "model", "environment_id"], as_index=False)
        .agg(
            streams=("cftp_seed", "size"),
            environment_total_seconds=("runtime_seconds", "sum"),
            environment_max_horizon_sweeps=(
                "coalescence_horizon_sweeps",
                "max",
            ),
        )
    )
    environment_costs.to_csv(output_path / "environment_costs.csv", index=False)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(source_path),
        "extension": str(extension_path),
        "streams": len(final),
        "environments": int(final["environment_id"].nunique()),
        "all_streams_certified": bool(final["certified"].all()),
        "initially_censored": int(final["initial_censored"].sum()),
        "remaining_censored": int(final["censored"].sum()),
        "all_extension_histories_replayed_exactly": bool(
            extended["source_gap_replayed_exactly"].all()
        ),
        "total_initial_calibration_wall_seconds": float(
            json.loads((source_path / "calibration_summary.json").read_text())[
                "total_wall_seconds"
            ]
        ),
        "total_extension_wall_seconds": float(
            json.loads((extension_path / "extension_summary.json").read_text())[
                "total_wall_seconds"
            ]
        ),
        "scientific_covariance_computed": False,
        "log_or_log_squared_fit_performed": False,
        "group_summary": groups.to_dict(orient="records"),
    }
    (output_path / "final_calibration_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report
