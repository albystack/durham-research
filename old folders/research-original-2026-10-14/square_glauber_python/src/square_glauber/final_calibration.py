"""Final frozen-environment mixing calibration and artifact comparisons.

This module orchestrates the validated accelerated backend.  It does not
implement or alter a transition kernel.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable, Mapping

import numpy as np
import pandas as pd

from .accelerated import compare_extremal_starts_accelerated, require_numba
from .diagnostics import integrated_autocorrelation_time
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, spatial_separations
from .mixing_campaign import write_mixing_case
from .mixing_extension import _load_frozen_environment, environment_hashes


EXPECTED_ENVIRONMENT_ID = "gamma-L16-diagnostic-02"
EXPECTED_ENVIRONMENT_SEED = 10195577682248251575
EXPECTED_ENVIRONMENT_SHA256 = (
    "e0defb728ab4ca8fae77ced3bb5aab00c464f8c412bf4339501a0254fe4d72c9"
)
EXPECTED_L = 16
EXPECTED_FACE_COUNT = 225


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_authoritative_environment(source_case: str | Path):
    """Load and strictly require the named, hashed frozen Gamma environment."""
    source_path = Path(source_case)
    environment, source_summary, verification = _load_frozen_environment(source_path)
    hashes = environment_hashes(environment)
    requirements = {
        "environment_id_matches": source_summary["environment_id"]
        == EXPECTED_ENVIRONMENT_ID,
        "L_matches": environment.L == EXPECTED_L,
        "horizontal_shape_matches": environment.horizontal_weights.shape
        == (EXPECTED_L, EXPECTED_L - 1),
        "vertical_shape_matches": environment.vertical_weights.shape
        == (EXPECTED_L - 1, EXPECTED_L),
        "environment_seed_matches": int(source_summary["environment_seed"])
        == EXPECTED_ENVIRONMENT_SEED,
        "environment_sha256_matches": hashes["environment_sha256"]
        == EXPECTED_ENVIRONMENT_SHA256,
        "stored_and_regenerated_arrays_match": bool(
            verification["horizontal_arrays_bitwise_equal_to_seed_regeneration"]
            and verification["vertical_arrays_bitwise_equal_to_seed_regeneration"]
            and verification["combined_hashes_match"]
        ),
    }
    if not all(requirements.values()):
        raise RuntimeError(
            "authoritative frozen-environment verification failed: "
            + json.dumps(requirements, sort_keys=True)
        )
    return environment, source_summary, verification | requirements


def _trace_prefix_verification(
    authoritative_one_million_case: Path, new_case: Path
) -> dict[str, object]:
    old = pd.read_csv(authoritative_one_million_case / "traces.csv.gz")
    new = pd.read_csv(new_case / "traces.csv.gz")
    old = old.sort_values(["start", "sweep"]).reset_index(drop=True)
    new_prefix = (
        new.loc[new["sweep"] <= 1_000_000, old.columns]
        .sort_values(["start", "sweep"])
        .reset_index(drop=True)
    )
    columns = {
        column: bool(
            np.array_equal(old[column].to_numpy(), new_prefix[column].to_numpy())
        )
        for column in old.columns
    }
    result = {
        "authoritative_one_million_directory": str(authoritative_one_million_case),
        "authoritative_rows": len(old),
        "new_prefix_rows": len(new_prefix),
        "columns_match_exactly": columns,
        "all_columns_and_rows_match_exactly": bool(
            len(old) == len(new_prefix) and all(columns.values())
        ),
    }
    if not result["all_columns_and_rows_match_exactly"]:
        raise RuntimeError("3m stream-safe trajectory does not match authoritative 1m prefix")
    return result


def _tail_summary(
    values: np.ndarray, fraction: float, record_every: int
) -> dict[str, object]:
    points = max(2, int(np.ceil(len(values) * fraction)))
    tail = np.asarray(values[-points:], dtype=np.float64)
    variance = float(tail.var(ddof=1))
    iact_records = float(integrated_autocorrelation_time(tail))
    ess = float(len(tail) / iact_records)
    if ess < 20:
        status = "UNSTABLE_OR_TOO_SHORT_ESS_BELOW_20"
    elif ess < 100:
        status = "LIMITED_ESS_BELOW_100"
    else:
        status = "ESTIMATE_REPORTED"
    return {
        "retained_fraction": fraction,
        "points": len(tail),
        "mean": float(tail.mean()),
        "variance": variance,
        "minimum": float(tail.min()),
        "maximum": float(tail.max()),
        "median": float(np.median(tail)),
        "iact_records": iact_records,
        "iact_sweeps": iact_records * record_every,
        "ess": ess,
        "mcse": float(np.sqrt(variance / ess)),
        "iact_status": status,
    }


def write_cumulative_tail_summary(
    traces: pd.DataFrame,
    output: str | Path,
    *,
    record_every: int,
    fractions: Iterable[float] = (0.1, 0.2, 0.3, 0.4, 0.5),
) -> pd.DataFrame:
    """Write retained-tail diagnostics for centre height and raw global sum."""
    rows: list[dict[str, object]] = []
    for start, group in traces.sort_values("sweep").groupby("start"):
        for observable in ("center_height", "height_sum"):
            values = group[observable].to_numpy(dtype=np.float64)
            for fraction in fractions:
                rows.append(
                    {"start": start, "observable": observable}
                    | _tail_summary(values, float(fraction), record_every)
                )
    frame = pd.DataFrame(rows)
    frame.to_csv(output, index=False)
    return frame


def _run_metrics(
    summary: Mapping[str, object], observable: str
) -> dict[str, float]:
    details = summary["observables"][observable]
    horizontal = details["horizontal_start"]
    vertical = details["vertical_start"]
    return {
        "horizontal_mean": float(horizontal["late_mean"]),
        "vertical_mean": float(vertical["late_mean"]),
        "horizontal_variance": float(horizontal["late_variance"]),
        "vertical_variance": float(vertical["late_variance"]),
        "horizontal_iact": float(
            horizontal["integrated_autocorrelation_time_sweeps"]
        ),
        "vertical_iact": float(
            vertical["integrated_autocorrelation_time_sweeps"]
        ),
        "horizontal_ess": float(horizontal["effective_sample_size_late"]),
        "vertical_ess": float(vertical["effective_sample_size_late"]),
        "split_rhat": float(details["split_rhat_late_half"]),
        "z_gap": float(details["z_gap_horizontal_minus_vertical"]),
        "absolute_start_mean_gap": abs(
            float(details["late_mean_gap_horizontal_minus_vertical"])
        ),
    }


def build_three_run_comparison(
    reference_200k_case: str | Path,
    stream_safe_1m_case: str | Path,
    stream_safe_3m_case: str | Path,
) -> pd.DataFrame:
    """Build the authoritative 200k/reference, 1m/Numba, 3m/Numba table."""
    cases = [
        ("200k reference", 200_000, Path(reference_200k_case)),
        ("1m stream-safe Numba", 1_000_000, Path(stream_safe_1m_case)),
        ("3m stream-safe Numba", 3_000_000, Path(stream_safe_3m_case)),
    ]
    rows: list[dict[str, object]] = []
    previous_by_observable: dict[str, dict[str, float]] = {}
    for run_name, expected_sweeps, case in cases:
        summary = json.loads((case / "summary.json").read_text())
        if int(summary["sweeps"]) != expected_sweeps:
            raise RuntimeError(f"unexpected sweep count in {case}")
        runtime = float(summary["horizontal_start"]["elapsed_seconds"]) + float(
            summary["vertical_start"]["elapsed_seconds"]
        )
        for observable in ("center_height", "height_sum"):
            metrics = _run_metrics(summary, observable)
            row: dict[str, object] = {
                "run": run_name,
                "sweeps": expected_sweeps,
                "observable": observable,
                **metrics,
                "runtime_seconds_two_chains": runtime,
            }
            previous = previous_by_observable.get(observable)
            if previous is not None:
                row.update(
                    {
                        "rhat_improvement_from_previous": previous["split_rhat"]
                        - metrics["split_rhat"],
                        "mean_gap_shrinkage_fraction_from_previous": (
                            1.0
                            - metrics["absolute_start_mean_gap"]
                            / previous["absolute_start_mean_gap"]
                            if previous["absolute_start_mean_gap"] > 0
                            else np.nan
                        ),
                        "horizontal_iact_ratio_to_previous": metrics["horizontal_iact"]
                        / previous["horizontal_iact"],
                        "vertical_iact_ratio_to_previous": metrics["vertical_iact"]
                        / previous["vertical_iact"],
                        "horizontal_ess_factor_from_previous": metrics["horizontal_ess"]
                        / previous["horizontal_ess"],
                        "vertical_ess_factor_from_previous": metrics["vertical_ess"]
                        / previous["vertical_ess"],
                        "horizontal_variance_relative_change": metrics[
                            "horizontal_variance"
                        ]
                        / previous["horizontal_variance"]
                        - 1.0,
                        "vertical_variance_relative_change": metrics[
                            "vertical_variance"
                        ]
                        / previous["vertical_variance"]
                        - 1.0,
                    }
                )
            rows.append(row)
            previous_by_observable[observable] = metrics
    frame = pd.DataFrame(rows)
    frame.to_csv(Path(stream_safe_3m_case) / "comparison_200k_1m_3m.csv", index=False)
    return frame


def backfill_exact_saved_counts(output: str | Path) -> dict[str, dict[str, int]]:
    """Add exact integer update counts to an artifact written by older summaries.

    The counts are exactly recoverable: actual changes equal attempts minus the
    recorded null moves, and flippable encounters are the unique integer whose
    ratio to attempts equals the recorded full-precision rate.
    """
    output_path = Path(output)
    summary_path = output_path / "summary.json"
    summary = json.loads(summary_path.read_text())
    recovered: dict[str, dict[str, int]] = {}
    for start_name in ("horizontal_start", "vertical_start"):
        details = summary[start_name]
        attempts = int(details["attempted_updates"])
        flippable = int(round(float(details["flippable_proposal_rate"]) * attempts))
        changes = attempts - int(details["null_moves"])
        if not np.isclose(
            flippable / attempts,
            float(details["flippable_proposal_rate"]),
            rtol=0.0,
            atol=np.finfo(float).eps,
        ):
            raise RuntimeError("recorded flippable rate does not recover an exact count")
        if not np.isclose(
            changes / attempts,
            float(details["actual_change_rate"]),
            rtol=0.0,
            atol=np.finfo(float).eps,
        ):
            raise RuntimeError("recorded change rate does not recover an exact count")
        details["flippable_proposals"] = flippable
        details["actual_changes"] = changes
        recovered[start_name] = {
            "attempted_updates": attempts,
            "flippable_proposals": flippable,
            "actual_changes": changes,
            "null_moves": int(details["null_moves"]),
        }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    artifact_path = output_path / "artifact_audit.json"
    if artifact_path.exists():
        artifact = json.loads(artifact_path.read_text())
        artifact["exact_update_counts"] = recovered
        artifact_path.write_text(json.dumps(artifact, indent=2) + "\n")
    case_path = output_path / "case_summary.csv"
    if case_path.exists():
        case = pd.read_csv(case_path)
        for label, start_name in (
            ("horizontal", "horizontal_start"),
            ("vertical", "vertical_start"),
        ):
            for field, value in recovered[start_name].items():
                case[f"{label}_{field}"] = value
        case.to_csv(case_path, index=False)
    return recovered


def audit_final_output_consistency(root: str | Path) -> dict[str, object]:
    """Cross-check the completed CSV/JSON final-calibration artifacts."""
    root_path = Path(root)
    hard = root_path / "hard_L16_3m"
    summary = json.loads((hard / "summary.json").read_text())
    artifact = json.loads((hard / "artifact_audit.json").read_text())
    observable = pd.read_csv(hard / "observable_summary.csv")
    comparison = pd.read_csv(hard / "comparison_200k_1m_3m.csv")
    windows = pd.read_csv(hard / "window_summary.csv")
    cumulative = pd.read_csv(hard / "cumulative_tail_summary.csv")
    traces = pd.read_csv(hard / "traces.csv.gz")
    acf = pd.read_csv(hard / "acf.csv.gz")

    observable_matches = True
    for row in observable.itertuples(index=False):
        details = summary["observables"][row.observable]
        first = details["horizontal_start"]
        second = details["vertical_start"]
        checks = (
            (row.horizontal_late_mean, first["late_mean"]),
            (row.vertical_late_mean, second["late_mean"]),
            (row.horizontal_late_variance, first["late_variance"]),
            (row.vertical_late_variance, second["late_variance"]),
            (
                row.horizontal_iact_sweeps,
                first["integrated_autocorrelation_time_sweeps"],
            ),
            (
                row.vertical_iact_sweeps,
                second["integrated_autocorrelation_time_sweeps"],
            ),
            (row.split_rhat_late_half, details["split_rhat_late_half"]),
            (row.z_gap_horizontal_minus_vertical, details["z_gap_horizontal_minus_vertical"]),
        )
        observable_matches &= all(
            np.isclose(float(left), float(right), rtol=1e-13, atol=1e-13)
            for left, right in checks
        )

    comparison_matches = True
    for row in comparison.loc[
        comparison["run"] == "3m stream-safe Numba"
    ].itertuples(index=False):
        metrics = _run_metrics(summary, row.observable)
        comparison_matches &= all(
            np.isclose(float(getattr(row, key)), value, rtol=1e-13, atol=1e-13)
            for key, value in metrics.items()
        )

    cumulative_matches = True
    for row in cumulative.loc[cumulative["retained_fraction"] == 0.5].itertuples(
        index=False
    ):
        start_key = (
            "horizontal_start" if row.start == "all_horizontal" else "vertical_start"
        )
        details = summary[start_key]["observables"][row.observable]
        cumulative_matches &= all(
            (
                np.isclose(row.mean, details["late_mean"]),
                np.isclose(row.variance, details["late_variance"]),
                np.isclose(
                    row.iact_sweeps,
                    details["integrated_autocorrelation_time_sweeps"],
                ),
                np.isclose(row.ess, details["effective_sample_size_late"]),
            )
        )

    expected_observables = set(summary["observables"])
    window_group_sizes = windows.groupby(["start", "observable"]).size()
    trace_group_sizes = traces.groupby("start").size()
    acf_group_sizes = acf.groupby(["start", "observable"]).size()
    checks: dict[str, object] = {
        "summary_observable_csv_matches": bool(observable_matches),
        "summary_comparison_csv_matches": bool(comparison_matches),
        "summary_cumulative_final_half_matches": bool(cumulative_matches),
        "window_groups_all_have_16_windows": bool((window_group_sizes == 16).all()),
        "window_group_count": int(len(window_group_sizes)),
        "trace_columns_include_all_spatial_increments": expected_observables
        == {"center_height", "height_sum"}
        | {column for column in traces if column.startswith("DeltaH_")},
        "trace_rows_per_start": {
            str(key): int(value) for key, value in trace_group_sizes.items()
        },
        "trace_sweeps_are_0_to_3m_by_10": bool(
            all(
                np.array_equal(
                    group.sort_values("sweep")["sweep"].to_numpy(),
                    np.arange(0, 3_000_001, 10),
                )
                for _, group in traces.groupby("start")
            )
        ),
        "acf_groups_all_have_100001_lags": bool((acf_group_sizes == 100_001).all()),
        "acf_group_count": int(len(acf_group_sizes)),
        "exact_count_payload_present": "exact_update_counts" in artifact,
        "cftp_optional_run_skipped_with_reason": bool(
            json.loads((root_path / "coupling" / "coupling_summary.json").read_text())[
                "run_performed"
            ]
            is False
        ),
        "cftp_prerequisite_note_present": (root_path / "CFTP_PREREQUISITES.md").is_file(),
    }
    boolean_checks = [value for value in checks.values() if isinstance(value, bool)]
    checks["all_consistency_checks_passed"] = all(boolean_checks)
    (root_path / "final_artifact_audit.json").write_text(
        json.dumps(checks, indent=2) + "\n"
    )
    if not checks["all_consistency_checks_passed"]:
        raise RuntimeError("cross-file final artifact audit failed")
    return checks


def run_final_hard_calibration(
    source_case: str | Path,
    authoritative_one_million_case: str | Path,
    output: str | Path,
    *,
    sweeps: int = 3_000_000,
    record_every: int = 10,
    relative_separations: Iterable[float] = DEFAULT_RELATIVE_SEPARATIONS,
    max_lag_records: int = 100_000,
    window_count: int = 16,
) -> dict[str, object]:
    """Run the two independent stream-safe 3m chains in the frozen hard case."""
    require_numba()
    if sweeps != 3_000_000:
        raise ValueError("this final calibration requires exactly 3,000,000 sweeps")
    if record_every <= 0 or sweeps % record_every:
        raise ValueError("record_every must be positive and divide 3,000,000")
    output_path = Path(output)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    output_path.mkdir(parents=True, exist_ok=False)
    source_path = Path(source_case)
    one_million_path = Path(authoritative_one_million_case)
    if not one_million_path.name.endswith("numba-1m-stream-safe"):
        raise ValueError("the authoritative 1m comparison must be the stream-safe directory")

    protected_files = [
        source_path / "summary.json",
        source_path / "traces.csv.gz",
        source_path / "environment_weights.npz",
        one_million_path / "summary.json",
        one_million_path / "traces.csv.gz",
    ]
    protected_hashes_before = {
        str(path): _file_sha256(path) for path in protected_files
    }
    environment, source_summary, verification = verify_authoritative_environment(
        source_path
    )
    verification.update(
        {
            "verified_before_chain_run": True,
            "expected_environment_sha256": EXPECTED_ENVIRONMENT_SHA256,
            "expected_environment_seed": EXPECTED_ENVIRONMENT_SEED,
            "expected_horizontal_shape": [EXPECTED_L, EXPECTED_L - 1],
            "expected_vertical_shape": [EXPECTED_L - 1, EXPECTED_L],
        }
    )
    (output_path / "environment_verification.json").write_text(
        json.dumps(verification, indent=2) + "\n"
    )
    metadata: dict[str, object] = {
        "environment_id": EXPECTED_ENVIRONMENT_ID,
        "L": EXPECTED_L,
        "model": "gamma",
        "gamma_shape": 1.0,
        "environment_index": int(source_summary["environment_index"]),
        "master_seed": int(source_summary["master_seed"]),
        "environment_seed": int(source_summary["environment_seed"]),
        "chain1_seed": int(source_summary["chain1_seed"]),
        "chain2_seed": int(source_summary["chain2_seed"]),
        "sweeps": sweeps,
        "record_every_sweeps": record_every,
        "source_case": str(source_path),
        "authoritative_one_million_case": str(one_million_path),
        "backend": "numba_stream_safe",
        "rng_stream_semantics": (
            "persistent independent PCG64 face-index and heat-bath-uniform "
            "substreams derived once from each recorded chain seed; one face and "
            "one uniform consumed per attempt; the uniform is ignored on null moves"
        ),
        **environment_hashes(environment),
    }
    total_started = perf_counter()
    comparison = compare_extremal_starts_accelerated(
        environment,
        np.random.Generator(np.random.PCG64(int(source_summary["chain1_seed"]))),
        np.random.Generator(np.random.PCG64(int(source_summary["chain2_seed"]))),
        sweeps=sweeps,
        record_every=record_every,
        separations=spatial_separations(environment.L, relative_separations),
    )
    chain_wall = perf_counter() - total_started
    case_row, observable_rows = write_mixing_case(
        comparison,
        output_path,
        metadata=metadata,
        max_lag_records=max_lag_records,
        late_fraction=0.5,
        window_count=window_count,
    )
    traces = pd.read_csv(output_path / "traces.csv.gz")
    cumulative = write_cumulative_tail_summary(
        traces,
        output_path / "cumulative_tail_summary.csv",
        record_every=record_every,
    )
    prefix = _trace_prefix_verification(one_million_path, output_path)
    (output_path / "one_million_prefix_verification.json").write_text(
        json.dumps(prefix, indent=2) + "\n"
    )
    protected_hashes_after = {
        str(path): _file_sha256(path) for path in protected_files
    }
    protected_unchanged = protected_hashes_before == protected_hashes_after
    expected_attempts = sweeps * EXPECTED_FACE_COUNT
    summary = json.loads((output_path / "summary.json").read_text())
    horizontal_attempts = int(summary["horizontal_start"]["attempted_updates"])
    vertical_attempts = int(summary["vertical_start"]["attempted_updates"])
    expected_rows = 2 * (sweeps // record_every + 1)
    artifact_checks = {
        "frozen_environment_hash_verified": verification[
            "environment_sha256_matches"
        ],
        "sweeps_per_chain": sweeps,
        "face_count": EXPECTED_FACE_COUNT,
        "expected_attempts_per_chain": expected_attempts,
        "horizontal_attempts": horizontal_attempts,
        "vertical_attempts": vertical_attempts,
        "attempt_counts_exact": horizontal_attempts == expected_attempts
        and vertical_attempts == expected_attempts,
        "expected_trace_rows": expected_rows,
        "actual_trace_rows": len(traces),
        "trace_row_count_exact": len(traces) == expected_rows,
        "one_million_prefix_exact": prefix["all_columns_and_rows_match_exactly"],
        "protected_prior_artifacts_unchanged": protected_unchanged,
        "authoritative_one_million_is_stream_safe": one_million_path.name.endswith(
            "numba-1m-stream-safe"
        ),
        "stream_safe_rng_metadata_present": "persistent independent PCG64"
        in str(summary["rng_stream_semantics"]),
        "scientific_sizes_run": [16],
        "paired_covariance_experiment_run": False,
        "hamilton_or_slurm_added": False,
    }
    if not all(
        artifact_checks[key]
        for key in (
            "frozen_environment_hash_verified",
            "attempt_counts_exact",
            "trace_row_count_exact",
            "one_million_prefix_exact",
            "protected_prior_artifacts_unchanged",
            "authoritative_one_million_is_stream_safe",
            "stream_safe_rng_metadata_present",
        )
    ):
        raise RuntimeError("final calibration artifact audit failed")
    (output_path / "artifact_audit.json").write_text(
        json.dumps(artifact_checks, indent=2) + "\n"
    )
    pd.DataFrame([case_row]).to_csv(output_path / "case_summary.csv", index=False)
    pd.DataFrame(observable_rows).to_csv(
        output_path / "observable_summary.csv", index=False
    )
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "chain_execution_wall_seconds": chain_wall,
        "total_wall_seconds_including_analysis_and_serialization": perf_counter()
        - total_started,
        "horizontal_chain_seconds": comparison.horizontal_start.elapsed_seconds,
        "vertical_chain_seconds": comparison.vertical_start.elapsed_seconds,
        "cumulative_tail_rows": len(cumulative),
        "acf_max_lag_records": max_lag_records,
        "acf_max_lag_sweeps": max_lag_records * record_every,
        "final_half_window_count": window_count,
        "warning": (
            "Independent-chain mixing diagnostic only. Classification requires "
            "numerical and visual review; coupled trajectories are not used."
        ),
    }
    (output_path / "run_metadata.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def write_final_results_summary(
    root: str | Path,
    *,
    classification: str,
    recommendation: str,
    late_window_assessment: str,
    visual_assessment: str,
    test_results: str,
) -> Path:
    """Finalize the human-reviewed scientific handoff after plots are inspected."""
    allowed_classifications = {
        "CREDIBLE FOR PILOT CALIBRATION",
        "BORDERLINE — EXTEND FURTHER",
        "NOT CREDIBLE — ORDINARY GLAUBER STILL TOO SLOW",
    }
    allowed_recommendations = {
        "A. proceed to multi-environment mixing calibration with accelerated ordinary Glauber",
        "B. prototype monotone CFTP before paired sampling",
        "C. continue investigating mixing / sampling before any scaling study",
    }
    if classification not in allowed_classifications:
        raise ValueError("invalid classification")
    if recommendation not in allowed_recommendations:
        raise ValueError("invalid recommendation")
    root_path = Path(root)
    hard = root_path / "hard_L16_3m"
    summary = json.loads((hard / "summary.json").read_text())
    runtime = json.loads((hard / "run_metadata.json").read_text())
    artifact = json.loads((hard / "artifact_audit.json").read_text())
    extremality = json.loads(
        (root_path / "extremality" / "extremality_summary.json").read_text()
    )
    monotonicity = json.loads(
        (root_path / "monotonicity" / "monotonicity_summary.json").read_text()
    )
    comparison = pd.read_csv(hard / "comparison_200k_1m_3m.csv")
    final_rows = comparison.loc[comparison["run"] == "3m stream-safe Numba"]
    lines = [
        "FINAL MIXING CALIBRATION + MONOTONE-COUPLING AUDIT",
        "Square-grid random-weight dimer Glauber project",
        "",
        f"Classification: {classification}",
        f"Final recommendation: {recommendation}",
        "",
        "3m independent-chain run",
        f"  environment: {summary['environment_id']}",
        f"  environment SHA-256: {summary['environment_sha256']}",
        f"  sweeps per chain: {summary['sweeps']:,}",
        f"  attempts per chain: {artifact['expected_attempts_per_chain']:,}",
        f"  horizontal chain runtime: {runtime['horizontal_chain_seconds']:.6f} s",
        f"  vertical chain runtime: {runtime['vertical_chain_seconds']:.6f} s",
        f"  chain execution wall: {runtime['chain_execution_wall_seconds']:.6f} s",
        f"  total wall incl. analysis/serialization: {runtime['total_wall_seconds_including_analysis_and_serialization']:.6f} s",
        "",
        "Final-half diagnostics",
    ]
    for row in final_rows.itertuples(index=False):
        lines.extend(
            [
                f"  {row.observable}:",
                f"    mean H/V = {row.horizontal_mean:.9g} / {row.vertical_mean:.9g}",
                f"    variance H/V = {row.horizontal_variance:.9g} / {row.vertical_variance:.9g}",
                f"    IACT sweeps H/V = {row.horizontal_iact:.9g} / {row.vertical_iact:.9g}",
                f"    ESS H/V = {row.horizontal_ess:.9g} / {row.vertical_ess:.9g}",
                f"    split R-hat = {row.split_rhat:.9g}",
                f"    z-gap = {row.z_gap:.9g}",
            ]
        )
    lines.extend(
        [
            "",
            "Late-window assessment",
            f"  {late_window_assessment}",
            "",
            "Visual assessment",
            f"  {visual_assessment}",
            "",
            "Height-order audit",
            "  all tested finite state spaces have unique legal pointwise extrema: "
            + str(
                extremality[
                    "all_tested_state_spaces_have_unique_legal_pointwise_extrema"
                ]
            ),
            "  implemented all-horizontal/all-vertical are extrema on every grid: "
            + str(
                extremality[
                    "implemented_horizontal_vertical_are_extrema_on_every_grid"
                ]
            ),
            "  exhaustive common-update small-state audit passed: "
            + str(
                monotonicity[
                    "all_exhaustive_small_state_checks_preserved_order"
                ]
            ),
            f"  exhaustive decision-regime checks: {monotonicity['total_regime_checks']:,}",
            "  mathematical status: finite evidence supports monotonicity, but the general height-lattice/local-admissibility proof is not in the repository.",
            "  CFTP status: not scientifically legitimate yet; the previously used horizontal/vertical states are not general global bounds.",
            "",
            "Tests",
            f"  {test_results}",
            "",
            "Integrity",
            "  Frozen hash, 675,000,000 attempts/chain, trace row count, exact 1m prefix, and preservation of prior key artifacts all passed.",
            "  No L=24/L=32 scientific run, paired covariance/scaling experiment, or Hamilton/Slurm work was performed.",
        ]
    )
    path = root_path / "RESULTS_SUMMARY.txt"
    path.write_text("\n".join(lines) + "\n")
    assessment = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "classification": classification,
        "recommendation": recommendation,
        "late_window_assessment": late_window_assessment,
        "visual_assessment": visual_assessment,
        "test_results": test_results,
    }
    (root_path / "assessment.json").write_text(
        json.dumps(assessment, indent=2) + "\n"
    )
    return path
