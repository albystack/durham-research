"""Equivalence, benchmark and non-scientific accelerated-backend diagnostics."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .accelerated import (
    AcceleratedRandomStreams,
    NUMBA_AVAILABLE,
    compare_extremal_starts_accelerated,
    horizontal_probability_grid,
    require_numba,
    run_supplied_attempts_accelerated,
)
from .environment import DimerEnvironment, sample_environment
from .experiment import _seed_uint64
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, center_face_index, spatial_separations
from .glauber import UpdateCounts, run_supplied_attempts
from .height import face_heights, spatial_increments
from .matching import DimerState
from .mixing_campaign import write_mixing_case
from .mixing_extension import _load_frozen_environment, environment_hashes


def _states_equal(first: DimerState, second: DimerState) -> bool:
    return bool(
        np.array_equal(first.horizontal_occupied, second.horizontal_occupied)
        and np.array_equal(first.vertical_occupied, second.vertical_occupied)
    )


def run_backend_equivalence_audit(
    output: str | Path,
    *,
    sizes: Iterable[int] = (4, 6, 8, 12, 16),
    attempts: int = 20_000,
    checkpoint_attempts: int = 257,
    seed: int = 20260908,
) -> dict[str, object]:
    """Require exact state/count equality under identical supplied randomness."""
    require_numba()
    if attempts <= 0 or checkpoint_attempts <= 0:
        raise ValueError("attempts and checkpoint_attempts must be positive")
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    rows: list[dict[str, object]] = []
    started = perf_counter()
    for L in sizes:
        for model in ("uniform", "gamma"):
            for start_name in ("all_horizontal", "all_vertical"):
                identity = [seed, int(L), 0 if model == "uniform" else 1, 0 if start_name == "all_horizontal" else 1]
                environment_seed, stream_seed = (
                    _seed_uint64(item)
                    for item in np.random.SeedSequence(identity).spawn(2)
                )
                environment = sample_environment(
                    int(L),
                    np.random.Generator(np.random.PCG64(environment_seed)),
                    model,
                    gamma_shape=1.0,
                )
                stream_rng = np.random.Generator(np.random.PCG64(stream_seed))
                faces = stream_rng.integers(
                    0, (int(L) - 1) ** 2, size=attempts, dtype=np.int64
                )
                uniforms = stream_rng.random(attempts)
                constructor = getattr(DimerState, start_name)
                reference = constructor(int(L))
                accelerated = constructor(int(L))
                reference_counts = UpdateCounts()
                accelerated_counts = UpdateCounts()
                probabilities = horizontal_probability_grid(environment)
                checkpoints = 0
                for begin in range(0, attempts, checkpoint_attempts):
                    end = min(attempts, begin + checkpoint_attempts)
                    run_supplied_attempts(
                        reference,
                        environment,
                        faces[begin:end],
                        uniforms[begin:end],
                        reference_counts,
                    )
                    run_supplied_attempts_accelerated(
                        accelerated,
                        environment,
                        faces[begin:end],
                        uniforms[begin:end],
                        accelerated_counts,
                        horizontal_probabilities=probabilities,
                    )
                    checkpoints += 1
                    if not _states_equal(reference, accelerated):
                        raise AssertionError(
                            f"backend state mismatch for L={L}, {model}, {start_name}, attempt={end}"
                        )
                    if reference_counts != accelerated_counts:
                        raise AssertionError(
                            f"backend count mismatch for L={L}, {model}, {start_name}, attempt={end}"
                        )
                reference.validate()
                accelerated.validate()
                rows.append(
                    {
                        "L": int(L),
                        "model": model,
                        "start": start_name,
                        "environment_seed": environment_seed,
                        "stream_seed": stream_seed,
                        "attempts": attempts,
                        "checkpoints": checkpoints,
                        "attempted_updates": reference_counts.attempted_updates,
                        "flippable_proposals": reference_counts.flippable_proposals,
                        "actual_changes": reference_counts.actual_changes,
                        "states_identical": True,
                        "counts_identical": True,
                    }
                )
    pd.DataFrame(rows).to_csv(output_path / "equivalence_cases.csv", index=False)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "numba_available": NUMBA_AVAILABLE,
        "sizes": [int(value) for value in sizes],
        "models": ["uniform", "gamma"],
        "starts": ["all_horizontal", "all_vertical"],
        "attempts_per_case": attempts,
        "checkpoint_attempts": checkpoint_attempts,
        "case_count": len(rows),
        "all_states_and_counts_identical": True,
        "wall_seconds": perf_counter() - started,
    }
    (output_path / "equivalence_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def benchmark_backends(
    output: str | Path,
    *,
    sizes: Iterable[int] = (8, 12, 16, 24, 32),
    target_attempts: int = 2_000_000,
    repeats: int = 3,
    seed: int = 20260908,
) -> dict[str, object]:
    """Benchmark supplied-stream kernels, excluding random-array generation."""
    require_numba()
    if target_attempts <= 0 or repeats <= 0:
        raise ValueError("target_attempts and repeats must be positive")
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)

    # Compile outside all recorded timings.
    warm_environment = sample_environment(4, np.random.default_rng(seed), "gamma")
    warm_state = DimerState.all_horizontal(4)
    run_supplied_attempts_accelerated(
        warm_state,
        warm_environment,
        np.asarray([0], dtype=np.int64),
        np.asarray([0.5]),
    )

    rows: list[dict[str, object]] = []
    started = perf_counter()
    for L in sizes:
        face_count = (int(L) - 1) ** 2
        sweeps = max(1, int(np.ceil(target_attempts / face_count)))
        attempts = sweeps * face_count
        environment = sample_environment(
            int(L), np.random.default_rng(seed + int(L)), "gamma", gamma_shape=1.0
        )
        probabilities = horizontal_probability_grid(environment)
        for repeat in range(repeats):
            rng = np.random.default_rng(seed + 10_000 * int(L) + repeat)
            faces = rng.integers(0, face_count, size=attempts, dtype=np.int64)
            uniforms = rng.random(attempts)
            reference = DimerState.all_horizontal(int(L))
            accelerated = DimerState.all_horizontal(int(L))
            reference_started = perf_counter()
            reference_counts = run_supplied_attempts(
                reference, environment, faces, uniforms
            )
            reference_seconds = perf_counter() - reference_started
            accelerated_started = perf_counter()
            accelerated_counts = run_supplied_attempts_accelerated(
                accelerated,
                environment,
                faces,
                uniforms,
                horizontal_probabilities=probabilities,
            )
            accelerated_seconds = perf_counter() - accelerated_started
            if not _states_equal(reference, accelerated):
                raise AssertionError(f"benchmark state mismatch for L={L}, repeat={repeat}")
            if reference_counts != accelerated_counts:
                raise AssertionError(f"benchmark count mismatch for L={L}, repeat={repeat}")
            rows.append(
                {
                    "L": int(L),
                    "repeat": repeat,
                    "sweeps": sweeps,
                    "attempts": attempts,
                    "reference_seconds": reference_seconds,
                    "accelerated_seconds": accelerated_seconds,
                    "reference_attempts_per_second": attempts / reference_seconds,
                    "accelerated_attempts_per_second": attempts / accelerated_seconds,
                    "acceleration_factor": reference_seconds / accelerated_seconds,
                }
            )
    raw = pd.DataFrame(rows)
    raw.to_csv(output_path / "benchmark_repeats.csv", index=False)
    summary = (
        raw.groupby("L", as_index=False)
        .agg(
            sweeps=("sweeps", "first"),
            attempts=("attempts", "first"),
            reference_seconds_median=("reference_seconds", "median"),
            accelerated_seconds_median=("accelerated_seconds", "median"),
            reference_attempts_per_second_median=(
                "reference_attempts_per_second",
                "median",
            ),
            accelerated_attempts_per_second_median=(
                "accelerated_attempts_per_second",
                "median",
            ),
            acceleration_factor_median=("acceleration_factor", "median"),
            acceleration_factor_min=("acceleration_factor", "min"),
            acceleration_factor_max=("acceleration_factor", "max"),
        )
    )
    summary.to_csv(output_path / "benchmark_summary.csv", index=False)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sizes": [int(value) for value in sizes],
        "target_attempts_per_repeat": target_attempts,
        "repeats": repeats,
        "timing_scope": "supplied-stream update kernels; RNG-array generation excluded",
        "numba_compilation_excluded": True,
        "wall_seconds": perf_counter() - started,
    }
    (output_path / "benchmark_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def audit_common_randomness_ordering(
    output: str | Path,
    *,
    sizes: Iterable[int] = (4, 6, 8, 12, 16),
    attempts: int = 20_000,
    streams_per_model_size: int = 2,
    seed: int = 20260908,
) -> dict[str, object]:
    """Test pointwise face-height ordering after every common update."""
    if attempts <= 0 or streams_per_model_size <= 0:
        raise ValueError("attempts and streams_per_model_size must be positive")
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    rows: list[dict[str, object]] = []
    failure: dict[str, object] | None = None
    started = perf_counter()
    for L in sizes:
        for model in ("uniform", "gamma"):
            for stream_index in range(streams_per_model_size):
                identity = [seed, int(L), int(model == "gamma"), stream_index, 991]
                environment_sequence, stream_sequence = np.random.SeedSequence(identity).spawn(2)
                environment_seed = _seed_uint64(environment_sequence)
                stream_seed = _seed_uint64(stream_sequence)
                environment = sample_environment(
                    int(L),
                    np.random.Generator(np.random.PCG64(environment_seed)),
                    model,
                    gamma_shape=1.0,
                )
                rng = np.random.Generator(np.random.PCG64(stream_seed))
                faces = rng.integers(
                    0, (int(L) - 1) ** 2, size=attempts, dtype=np.int64
                )
                uniforms = rng.random(attempts)
                lower = DimerState.all_horizontal(int(L))
                upper = DimerState.all_vertical(int(L))
                coalescence_attempt: int | None = None
                checked = 0
                for index in range(attempts):
                    run_supplied_attempts(
                        lower, environment, faces[index : index + 1], uniforms[index : index + 1]
                    )
                    run_supplied_attempts(
                        upper, environment, faces[index : index + 1], uniforms[index : index + 1]
                    )
                    checked = index + 1
                    lower_heights = face_heights(lower)
                    upper_heights = face_heights(upper)
                    if not np.all(lower_heights <= upper_heights):
                        violating = np.argwhere(lower_heights > upper_heights)[0]
                        failure = {
                            "L": int(L),
                            "model": model,
                            "stream_index": stream_index,
                            "environment_seed": environment_seed,
                            "stream_seed": stream_seed,
                            "attempt": checked,
                            "face_index": int(faces[index]),
                            "uniform": float(uniforms[index]),
                            "violating_face": [int(value) for value in violating],
                            "lower_height": int(lower_heights[tuple(violating)]),
                            "upper_height": int(upper_heights[tuple(violating)]),
                        }
                        break
                    if coalescence_attempt is None and _states_equal(lower, upper):
                        coalescence_attempt = checked
                        break
                rows.append(
                    {
                        "L": int(L),
                        "model": model,
                        "stream_index": stream_index,
                        "environment_seed": environment_seed,
                        "stream_seed": stream_seed,
                        "attempts_checked": checked,
                        "ordering_preserved": failure is None,
                        "coalescence_attempt": coalescence_attempt,
                    }
                )
                if failure is not None:
                    break
            if failure is not None:
                break
        if failure is not None:
            break
    pd.DataFrame(rows).to_csv(output_path / "ordering_cases.csv", index=False)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "expected_order": "H_all_horizontal <= H_all_vertical pointwise",
        "checked_after_every_attempt": True,
        "requested_attempts_per_case": attempts,
        "streams_per_model_size": streams_per_model_size,
        "case_count": len(rows),
        "ordering_preserved_in_all_checked_cases": failure is None,
        "first_failure": failure,
        "wall_seconds": perf_counter() - started,
        "warning": "Exploratory finite testing is not a proof of attractiveness/monotonicity.",
    }
    (output_path / "ordering_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def run_accelerated_frozen_stress(
    source_case: str | Path,
    output: str | Path,
    *,
    sweeps: int = 1_000_000,
    record_every: int = 10,
    relative_separations: Iterable[float] = DEFAULT_RELATIVE_SEPARATIONS,
    max_lag_records: int = 20_000,
) -> dict[str, object]:
    """Run independent accelerated chains in one exactly saved environment."""
    require_numba()
    source_case_path = Path(source_case)
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    environment, source_summary, verification = _load_frozen_environment(source_case_path)
    metadata: dict[str, object] = {
        "environment_id": str(source_summary["environment_id"]),
        "L": environment.L,
        "model": "gamma",
        "gamma_shape": 1.0,
        "environment_index": int(source_summary["environment_index"]),
        "master_seed": int(source_summary["master_seed"]),
        "environment_seed": int(source_summary["environment_seed"]),
        "chain1_seed": int(source_summary["chain1_seed"]),
        "chain2_seed": int(source_summary["chain2_seed"]),
        "sweeps": sweeps,
        "record_every_sweeps": record_every,
        "source_case": str(source_case_path),
        "source_sweeps": int(source_summary["sweeps"]),
        "backend": "numba",
        "rng_stream_semantics": (
            "persistent independent PCG64 face-index and heat-bath-uniform "
            "substreams derived once from each recorded chain seed; one value "
            "from each stream per attempt and the uniform is ignored on null moves"
        ),
        **environment_hashes(environment),
    }
    started = perf_counter()
    comparison = compare_extremal_starts_accelerated(
        environment,
        np.random.Generator(np.random.PCG64(int(source_summary["chain1_seed"]))),
        np.random.Generator(np.random.PCG64(int(source_summary["chain2_seed"]))),
        sweeps=sweeps,
        record_every=record_every,
        separations=spatial_separations(environment.L, relative_separations),
    )
    case_row, observable_rows = write_mixing_case(
        comparison,
        output_path,
        metadata=metadata,
        max_lag_records=max_lag_records,
    )
    np.savez_compressed(
        output_path / "environment_weights.npz",
        horizontal_weights=environment.horizontal_weights,
        vertical_weights=environment.vertical_weights,
    )
    (output_path / "environment_verification.json").write_text(
        json.dumps(verification, indent=2) + "\n"
    )
    pd.DataFrame([case_row]).to_csv(output_path / "case_summary.csv", index=False)
    pd.DataFrame(observable_rows).to_csv(
        output_path / "observable_summary.csv", index=False
    )
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "environment_id": metadata["environment_id"],
        "environment_sha256": metadata["environment_sha256"],
        "sweeps_per_chain": sweeps,
        "record_every_sweeps": record_every,
        "wall_seconds": perf_counter() - started,
        "horizontal_chain_seconds": comparison.horizontal_start.elapsed_seconds,
        "vertical_chain_seconds": comparison.vertical_start.elapsed_seconds,
        "warning": "Stress/mixing diagnostic only; not a paired scaling experiment.",
    }
    (output_path / "stress_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def run_common_randomness_coupling(
    source_case: str | Path,
    output: str | Path,
    *,
    ordering_audit: str | Path,
    sweeps: int = 1_000_000,
    record_every: int = 10,
    seed: int = 20260908,
) -> dict[str, object]:
    """Exploratory ordered common-randomness coupling, separate from sampling."""
    require_numba()
    source_case_path = Path(source_case)
    output_path = Path(output)
    ordering_audit_path = Path(ordering_audit)
    ordering_report = json.loads(ordering_audit_path.read_text())
    if not ordering_report.get("ordering_preserved_in_all_checked_cases", False):
        raise RuntimeError("coupling refused because the ordering audit did not pass")
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    environment, source_summary, verification = _load_frozen_environment(source_case_path)
    coupling_seed = _seed_uint64(
        np.random.SeedSequence(
            [seed, environment.L, int(source_summary["environment_index"]), 0xC0A1ED]
        )
    )
    rng = np.random.Generator(np.random.PCG64(coupling_seed))
    streams = AcceleratedRandomStreams.from_generator(rng)
    lower = DimerState.all_horizontal(environment.L)
    upper = DimerState.all_vertical(environment.L)
    probabilities = horizontal_probability_grid(environment)
    lower_counts = UpdateCounts()
    upper_counts = UpdateCounts()
    face_count = (environment.L - 1) ** 2
    rows: list[dict[str, object]] = []
    coalescence_sweep: int | None = None
    started = perf_counter()

    def record(sweep: int) -> None:
        nonlocal coalescence_sweep
        lower_heights = face_heights(lower)
        upper_heights = face_heights(upper)
        ordered = bool(np.all(lower_heights <= upper_heights))
        if not ordered:
            raise AssertionError(f"common-randomness height ordering failed at sweep {sweep}")
        coalesced = _states_equal(lower, upper)
        if coalesced and coalescence_sweep is None:
            coalescence_sweep = sweep
        rows.append(
            {
                "sweep": sweep,
                "lower_center_height": int(
                    lower_heights[center_face_index(environment.L)]
                ),
                "upper_center_height": int(
                    upper_heights[center_face_index(environment.L)]
                ),
                "center_height_gap": int(
                    upper_heights[center_face_index(environment.L)]
                    - lower_heights[center_face_index(environment.L)]
                ),
                "lower_height_sum": int(lower_heights.sum()),
                "upper_height_sum": int(upper_heights.sum()),
                "height_sum_gap": int(upper_heights.sum() - lower_heights.sum()),
                "maximum_pointwise_height_gap": int(
                    np.max(upper_heights - lower_heights)
                ),
                "ordered": ordered,
                "coalesced": coalesced,
            }
        )

    record(0)
    completed = 0
    while completed < sweeps:
        block_sweeps = min(record_every, sweeps - completed)
        block_attempts = block_sweeps * face_count
        faces = streams.face_rng.integers(
            0, face_count, size=block_attempts, dtype=np.int64
        )
        uniforms = streams.uniform_rng.random(block_attempts)
        run_supplied_attempts_accelerated(
            lower,
            environment,
            faces,
            uniforms,
            lower_counts,
            horizontal_probabilities=probabilities,
        )
        run_supplied_attempts_accelerated(
            upper,
            environment,
            faces,
            uniforms,
            upper_counts,
            horizontal_probabilities=probabilities,
        )
        completed += block_sweeps
        record(completed)
        if coalescence_sweep is not None:
            break
    trace = pd.DataFrame(rows)
    trace.to_csv(output_path / "coupled_trace.csv.gz", index=False)
    figure, axes = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(trace["sweep"], trace["center_height_gap"], linewidth=0.9)
    axes[0].set_ylabel("centre gap")
    axes[1].plot(trace["sweep"], trace["height_sum_gap"], linewidth=0.9)
    axes[1].set_ylabel("height-sum gap")
    axes[2].plot(
        trace["sweep"], trace["maximum_pointwise_height_gap"], linewidth=0.9
    )
    axes[2].set_ylabel("maximum gap")
    axes[2].set_xlabel("sweep")
    figure.suptitle(
        "Common-randomness extremal coupling\n"
        "Exploratory only — coupled chains are not independent samples"
    )
    figure.tight_layout()
    figure.savefig(output_path / "coupled_gaps.png", dpi=160)
    plt.close(figure)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "environment_id": source_summary["environment_id"],
        "environment_sha256": environment_hashes(environment)["environment_sha256"],
        "coupling_seed": coupling_seed,
        "prerequisite_ordering_audit": str(ordering_audit_path),
        "prerequisite_ordering_audit_passed": True,
        "requested_sweeps": sweeps,
        "completed_sweeps": completed,
        "record_every_sweeps": record_every,
        "ordering_preserved_at_all_recorded_checkpoints": True,
        "coalescence_sweep_upper_bound": coalescence_sweep,
        "wall_seconds": perf_counter() - started,
        "environment_verification": verification,
        "warning": (
            "Exploratory common-randomness diagnostic only; finite ordering tests "
            "are not a proof and this coupled trace is not an independent-chain sample."
        ),
    }
    (output_path / "coupling_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report
