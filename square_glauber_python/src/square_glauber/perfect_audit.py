"""Validation and bounded local calibration for monotone perfect sampling."""

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
from .exact import (
    EdgeMatching,
    canonical_edge,
    deterministic_test_weights,
    enumerate_perfect_matchings,
    gibbs_probabilities,
    rectangular_face_heights,
)
from .experiment import _seed_uint64
from .final_calibration import verify_authoritative_environment
from .height import center_height, face_heights
from .matching import DimerState
from .mixing_extension import environment_hashes
from .order_audit import matching_signature
from .perfect import CFTPNotCoalesced, MonotoneCFTPSampler


def _state_as_matching(state: DimerState) -> EdgeMatching:
    selected = []
    for i, j in np.argwhere(state.horizontal_occupied):
        selected.append(canonical_edge((int(i), int(j)), (int(i), int(j) + 1)))
    for i, j in np.argwhere(state.vertical_occupied):
        selected.append(canonical_edge((int(i), int(j)), (int(i) + 1, int(j))))
    return frozenset(selected)


def run_exact_perfect_sampling_diagnostic(
    output: str | Path,
    *,
    samples: int = 20_000,
    seed: int = 20260909,
    max_horizon_sweeps: int = 16_384,
) -> dict[str, object]:
    """Compare independent perfect 4x4 samples with the exact weighted law."""
    if samples <= 0:
        raise ValueError("samples must be positive")
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    horizontal, vertical = deterministic_test_weights(4, 4)
    environment = DimerEnvironment(4, horizontal, vertical, model="deterministic_test")
    matchings = enumerate_perfect_matchings(4, 4)
    probabilities = gibbs_probabilities(matchings, horizontal, vertical)
    matching_index = {matching: index for index, matching in enumerate(matchings)}
    exact_centers = np.asarray(
        [rectangular_face_heights(item, 4, 4)[1, 1] for item in matchings],
        dtype=np.int64,
    )
    exact_center_mean = float(probabilities @ exact_centers)
    exact_center_variance = float(
        probabilities @ (exact_centers - exact_center_mean) ** 2
    )
    sampler = MonotoneCFTPSampler(environment, backend="numba")
    rows: list[dict[str, object]] = []
    counts = np.zeros(len(matchings), dtype=np.int64)
    started = perf_counter()
    for sample_index in range(samples):
        cftp_seed = _seed_uint64(
            np.random.SeedSequence([seed, sample_index, 4, 0xECA7])
        )
        result = sampler.sample(
            cftp_seed,
            max_horizon_sweeps=max_horizon_sweeps,
            verify_order=False,
        )
        index = matching_index[_state_as_matching(result.state)]
        counts[index] += 1
        rows.append(
            {
                "sample_index": sample_index,
                "cftp_seed": cftp_seed,
                "matching_index": index,
                "center_height": center_height(result.state),
                "height_sum": int(face_heights(result.state).sum()),
                "coalescence_horizon_sweeps": result.horizon_sweeps,
                "doubling_iterations": len(result.iterations),
                "runtime_seconds": result.elapsed_seconds,
                "certified": result.certified,
            }
        )
    empirical = counts / counts.sum()
    state_table = pd.DataFrame(
        {
            "matching_index": np.arange(len(matchings)),
            "matching_signature": [matching_signature(item) for item in matchings],
            "center_height": exact_centers,
            "exact_probability": probabilities,
            "empirical_count": counts,
            "empirical_probability": empirical,
            "probability_error": empirical - probabilities,
        }
    )
    sample_table = pd.DataFrame(rows)
    sample_table.to_csv(output_path / "perfect_samples.csv.gz", index=False)
    state_table.to_csv(output_path / "exact_state_comparison.csv", index=False)
    center_values = np.unique(exact_centers)
    center_pmf = pd.DataFrame(
        {
            "center_height": center_values,
            "exact_probability": [
                float(probabilities[exact_centers == value].sum())
                for value in center_values
            ],
            "empirical_probability": [
                float(counts[exact_centers == value].sum() / counts.sum())
                for value in center_values
            ],
        }
    )
    center_pmf["probability_error"] = (
        center_pmf["empirical_probability"] - center_pmf["exact_probability"]
    )
    center_pmf.to_csv(output_path / "center_height_pmf.csv", index=False)
    empirical_centers = sample_table["center_height"].to_numpy(dtype=np.float64)
    empirical_center_mean = float(empirical_centers.mean())
    empirical_center_variance = float(empirical_centers.var(ddof=0))
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "grid": "4x4",
        "weights": "deterministic positive non-symmetric test weights",
        "matching_count": len(matchings),
        "samples": samples,
        "all_samples_cftp_certified": bool(sample_table["certified"].all()),
        "exact_center_mean": exact_center_mean,
        "empirical_center_mean": empirical_center_mean,
        "center_mean_error": empirical_center_mean - exact_center_mean,
        "center_mean_monte_carlo_se": float(
            np.sqrt(exact_center_variance / samples)
        ),
        "center_mean_standardized_error": float(
            (empirical_center_mean - exact_center_mean)
            / np.sqrt(exact_center_variance / samples)
        ),
        "exact_center_variance": exact_center_variance,
        "empirical_center_variance": empirical_center_variance,
        "center_variance_error": empirical_center_variance - exact_center_variance,
        "maximum_absolute_state_probability_error": float(
            np.max(np.abs(empirical - probabilities))
        ),
        "state_distribution_total_variation_distance": float(
            0.5 * np.sum(np.abs(empirical - probabilities))
        ),
        "empirical_states_with_zero_count": int(np.sum(counts == 0)),
        "exact_center_height_pmf": center_pmf[
            ["center_height", "exact_probability"]
        ].to_dict(orient="records"),
        "maximum_coalescence_horizon_sweeps": int(
            sample_table["coalescence_horizon_sweeps"].max()
        ),
        "median_coalescence_horizon_sweeps": float(
            sample_table["coalescence_horizon_sweeps"].median()
        ),
        "wall_seconds": perf_counter() - started,
        "interpretation": (
            "Empirical agreement is a diagnostic, not the proof of exactness. "
            "Exactness follows from the monotone CFTP construction and validated kernel."
        ),
    }
    (output_path / "exact_perfect_sampling_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def run_frozen_hard_cftp_calibration(
    source_case: str | Path,
    output: str | Path,
    *,
    samples: int = 50,
    seed: int = 20260909,
    max_horizon_sweeps: int = 1_048_576,
) -> dict[str, object]:
    """Measure true-bound CFTP horizons in the frozen hard L=16 environment."""
    if samples <= 0:
        raise ValueError("samples must be positive")
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    environment, source_summary, verification = verify_authoritative_environment(
        source_case
    )
    sampler = MonotoneCFTPSampler(environment, backend="numba")
    minimum_heights = face_heights(sampler.minimum.state)
    maximum_heights = face_heights(sampler.maximum.state)
    np.savez_compressed(
        output_path / "true_extremal_states.npz",
        minimum_horizontal_occupied=sampler.minimum.state.horizontal_occupied,
        minimum_vertical_occupied=sampler.minimum.state.vertical_occupied,
        maximum_horizontal_occupied=sampler.maximum.state.horizontal_occupied,
        maximum_vertical_occupied=sampler.maximum.state.vertical_occupied,
        minimum_heights=minimum_heights,
        maximum_heights=maximum_heights,
    )
    rows: list[dict[str, object]] = []
    censored = 0
    started = perf_counter()
    for sample_index in range(samples):
        cftp_seed = _seed_uint64(
            np.random.SeedSequence(
                [seed, int(source_summary["environment_seed"]), sample_index, 0xC0F7]
            )
        )
        try:
            result = sampler.sample(
                cftp_seed,
                max_horizon_sweeps=max_horizon_sweeps,
                verify_order=True,
            )
            rows.append(
                {
                    "sample_index": sample_index,
                    "cftp_seed": cftp_seed,
                    "certified": True,
                    "censored": False,
                    "coalescence_horizon_sweeps": result.horizon_sweeps,
                    "doubling_iterations": len(result.iterations),
                    "runtime_seconds": result.elapsed_seconds,
                    "center_height": center_height(result.state),
                    "height_sum": int(face_heights(result.state).sum()),
                    "lower_total_replayed_attempts": result.lower_update_counts.attempted_updates,
                    "upper_total_replayed_attempts": result.upper_update_counts.attempted_updates,
                }
            )
        except CFTPNotCoalesced as error:
            censored += 1
            rows.append(
                {
                    "sample_index": sample_index,
                    "cftp_seed": cftp_seed,
                    "certified": False,
                    "censored": True,
                    "coalescence_horizon_sweeps": error.horizon_sweeps,
                    "doubling_iterations": len(error.iterations),
                    "runtime_seconds": np.nan,
                    "center_height": np.nan,
                    "height_sum": np.nan,
                    "lower_total_replayed_attempts": np.nan,
                    "upper_total_replayed_attempts": np.nan,
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(output_path / "cftp_calibration_samples.csv", index=False)
    figure, axis = plt.subplots(figsize=(7.5, 4.5))
    completed = frame.loc[~frame["censored"]]
    if len(completed):
        values = completed["coalescence_horizon_sweeps"].to_numpy(dtype=np.int64)
        bins = np.unique(np.r_[values, 2 * values.max()])
        axis.hist(values, bins=bins, align="left", rwidth=0.85)
        axis.set_xscale("log", base=2)
    axis.set_xlabel("certifying CFTP past horizon (sweeps, powers of two)")
    axis.set_ylabel("samples")
    axis.set_title("Frozen hard L=16 true-bound CFTP calibration")
    figure.tight_layout()
    figure.savefig(output_path / "cftp_horizons.png", dpi=160)
    plt.close(figure)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "environment_id": source_summary["environment_id"],
        "environment_seed": int(source_summary["environment_seed"]),
        **environment_hashes(environment),
        "environment_verification": verification,
        "samples_requested": samples,
        "samples_certified": int((~frame["censored"]).sum()),
        "samples_censored": censored,
        "maximum_allowed_horizon_sweeps": max_horizon_sweeps,
        "minimum_construction": {
            "directed_flips": sampler.minimum.directed_flips,
            "processed_faces": sampler.minimum.processed_faces,
            "height_sum": int(minimum_heights.sum()),
        },
        "maximum_construction": {
            "directed_flips": sampler.maximum.directed_flips,
            "processed_faces": sampler.maximum.processed_faces,
            "height_sum": int(maximum_heights.sum()),
        },
        "horizon_min": (
            int(completed["coalescence_horizon_sweeps"].min())
            if len(completed)
            else None
        ),
        "horizon_median": (
            float(completed["coalescence_horizon_sweeps"].median())
            if len(completed)
            else None
        ),
        "horizon_max": (
            int(completed["coalescence_horizon_sweeps"].max())
            if len(completed)
            else None
        ),
        "mean_runtime_seconds": (
            float(completed["runtime_seconds"].mean()) if len(completed) else None
        ),
        "total_wall_seconds": perf_counter() - started,
        "warning": (
            "CFTP horizon is a certification horizon, not automatically a mixing-time estimate. "
            "No coupled output is treated as two independent samples."
        ),
    }
    (output_path / "cftp_calibration_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def run_perfect_pair_smoke(
    output: str | Path,
    *,
    sizes: Iterable[int] = (4, 8, 12, 16),
    seed: int = 20260909,
) -> dict[str, object]:
    """Exercise two independent perfect samples per fixed environment, without fits."""
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    rows: list[dict[str, object]] = []
    started = perf_counter()
    for L in sizes:
        for model in ("uniform", "gamma"):
            environment_seed = _seed_uint64(
                np.random.SeedSequence([seed, int(L), int(model == "gamma"), 0])
            )
            environment = sample_environment(
                int(L),
                np.random.Generator(np.random.PCG64(environment_seed)),
                model,
                gamma_shape=1.0,
            )
            first_seed = _seed_uint64(
                np.random.SeedSequence([seed, int(L), int(model == "gamma"), 1])
            )
            second_seed = _seed_uint64(
                np.random.SeedSequence([seed, int(L), int(model == "gamma"), 2])
            )
            sampler = MonotoneCFTPSampler(environment, backend="numba")
            pair = sampler.sample_pair(first_seed, second_seed, verify_order=True)
            rows.append(
                {
                    "L": int(L),
                    "model": model,
                    "gamma_shape": 1.0 if model == "gamma" else np.nan,
                    "environment_seed": environment_seed,
                    "environment_sha256": environment_hashes(environment)[
                        "environment_sha256"
                    ],
                    "first_cftp_seed": first_seed,
                    "second_cftp_seed": second_seed,
                    "seeds_distinct": first_seed != second_seed,
                    "first_certified": pair.first.certified,
                    "second_certified": pair.second.certified,
                    "first_horizon_sweeps": pair.first.horizon_sweeps,
                    "second_horizon_sweeps": pair.second.horizon_sweeps,
                    "first_center_height": center_height(pair.first.state),
                    "second_center_height": center_height(pair.second.state),
                    "first_runtime_seconds": pair.first.elapsed_seconds,
                    "second_runtime_seconds": pair.second.elapsed_seconds,
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(output_path / "perfect_pair_smoke.csv", index=False)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sizes": [int(value) for value in sizes],
        "models": ["uniform", "gamma"],
        "cases": len(frame),
        "all_pairs_use_distinct_seeds": bool(frame["seeds_distinct"].all()),
        "all_samples_certified": bool(
            frame[["first_certified", "second_certified"]].all(axis=None)
        ),
        "wall_seconds": perf_counter() - started,
        "scientific_scaling_analysis_performed": False,
    }
    (output_path / "perfect_pair_smoke_summary.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report
