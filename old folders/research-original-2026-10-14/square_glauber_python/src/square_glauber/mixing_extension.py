"""Longer mixing calibration in previously saved frozen Gamma environments."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable, Mapping

import numpy as np
import pandas as pd

from .diagnostics import compare_extremal_starts
from .environment import DimerEnvironment, sample_environment
from .experiment import derive_seeds
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, spatial_separations, validate_L
from .mixing_campaign import trace_frame, write_mixing_case


EDGE_QUANTILES = (0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
ABS_LAMBDA_QUANTILES = (0.5, 0.75, 0.9, 0.95, 0.99)
ABS_LAMBDA_THRESHOLDS = (2.0, 4.0, 6.0, 8.0)


def _quantile_label(value: float) -> str:
    return f"q{int(round(100 * value)):02d}"


def array_sha256(array: np.ndarray) -> str:
    """Hash array identity including shape, dtype and C-order bytes."""
    contiguous = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(contiguous.dtype).encode("ascii"))
    digest.update(np.asarray(contiguous.shape, dtype=np.int64).tobytes())
    digest.update(contiguous.tobytes(order="C"))
    return digest.hexdigest()


def environment_hashes(environment: DimerEnvironment) -> dict[str, str]:
    """Return separate and combined SHA-256 identifiers for one environment."""
    horizontal_hash = array_sha256(environment.horizontal_weights)
    vertical_hash = array_sha256(environment.vertical_weights)
    combined = hashlib.sha256(
        f"horizontal:{horizontal_hash};vertical:{vertical_hash}".encode("ascii")
    ).hexdigest()
    return {
        "horizontal_weights_sha256": horizontal_hash,
        "vertical_weights_sha256": vertical_hash,
        "environment_sha256": combined,
    }


def environment_hardness(environment: DimerEnvironment) -> dict[str, float | int]:
    """Descriptive edge and face log-odds diagnostics; never a selection rule."""
    weights = np.concatenate(
        (environment.horizontal_weights.ravel(), environment.vertical_weights.ravel())
    )
    log_weights = np.log(weights)
    horizontal_log = np.log(environment.horizontal_weights)
    vertical_log = np.log(environment.vertical_weights)
    face_lambda = (
        horizontal_log[:-1, :]
        + horizontal_log[1:, :]
        - vertical_log[:, 1:]
        - vertical_log[:, :-1]
    )
    absolute_lambda = np.abs(face_lambda).ravel()
    result: dict[str, float | int] = {
        "edge_count": int(weights.size),
        "face_count": int(absolute_lambda.size),
        "edge_weight_min": float(weights.min()),
        "edge_weight_max": float(weights.max()),
        "edge_weight_mean": float(weights.mean()),
        "log_edge_weight_min": float(log_weights.min()),
        "log_edge_weight_max": float(log_weights.max()),
        "log_edge_weight_range": float(np.ptp(log_weights)),
        "abs_lambda_max": float(absolute_lambda.max()),
    }
    for quantile, value in zip(
        EDGE_QUANTILES, np.quantile(weights, EDGE_QUANTILES), strict=True
    ):
        result[f"edge_weight_{_quantile_label(quantile)}"] = float(value)
    for quantile, value in zip(
        ABS_LAMBDA_QUANTILES,
        np.quantile(absolute_lambda, ABS_LAMBDA_QUANTILES),
        strict=True,
    ):
        result[f"abs_lambda_{_quantile_label(quantile)}"] = float(value)
    for threshold in ABS_LAMBDA_THRESHOLDS:
        result[f"fraction_abs_lambda_gt_{int(threshold)}"] = float(
            np.mean(absolute_lambda > threshold)
        )
    return result


def _load_frozen_environment(
    case_directory: Path,
) -> tuple[DimerEnvironment, dict[str, object], dict[str, object]]:
    old_summary = json.loads((case_directory / "summary.json").read_text())
    if old_summary["model"] != "gamma" or float(old_summary["gamma_shape"]) != 1.0:
        raise ValueError(f"{case_directory} is not a Gamma(shape=1) case")
    L = validate_L(int(old_summary["L"]))
    with np.load(case_directory / "environment_weights.npz", allow_pickle=False) as archive:
        horizontal = np.array(archive["horizontal_weights"], copy=True)
        vertical = np.array(archive["vertical_weights"], copy=True)
    stored = DimerEnvironment(
        L,
        horizontal,
        vertical,
        model="gamma",
        gamma_shape=1.0,
    )

    expected = derive_seeds(
        int(old_summary["master_seed"]),
        L,
        "gamma",
        int(old_summary["environment_index"]),
    )
    seed_fields_match = all(
        int(old_summary[name]) == int(getattr(expected, name))
        for name in ("environment_seed", "chain1_seed", "chain2_seed")
    )
    regenerated = sample_environment(
        L,
        np.random.Generator(np.random.PCG64(int(old_summary["environment_seed"]))),
        "gamma",
        gamma_shape=1.0,
    )
    horizontal_equal = np.array_equal(
        stored.horizontal_weights, regenerated.horizontal_weights
    )
    vertical_equal = np.array_equal(stored.vertical_weights, regenerated.vertical_weights)
    stored_hashes = environment_hashes(stored)
    regenerated_hashes = environment_hashes(regenerated)
    verification: dict[str, object] = {
        "source_case_directory": str(case_directory),
        "environment_seed": int(old_summary["environment_seed"]),
        "recorded_seed_tuple_matches_current_derivation": seed_fields_match,
        "horizontal_arrays_bitwise_equal_to_seed_regeneration": horizontal_equal,
        "vertical_arrays_bitwise_equal_to_seed_regeneration": vertical_equal,
        **{f"stored_{key}": value for key, value in stored_hashes.items()},
        **{f"regenerated_{key}": value for key, value in regenerated_hashes.items()},
        "combined_hashes_match": (
            stored_hashes["environment_sha256"]
            == regenerated_hashes["environment_sha256"]
        ),
    }
    if not all((seed_fields_match, horizontal_equal, vertical_equal)):
        raise RuntimeError(f"frozen-environment verification failed for {case_directory}")
    return stored, old_summary, verification


def _verify_trace_prefix(old_case: Path, new_case: Path) -> dict[str, object]:
    old = pd.read_csv(old_case / "traces.csv.gz").sort_values(["start", "sweep"])
    new = pd.read_csv(new_case / "traces.csv.gz").sort_values(["start", "sweep"])
    old_max_sweep = int(old["sweep"].max())
    prefix = new.loc[new["sweep"] <= old_max_sweep, old.columns]
    old = old.reset_index(drop=True)
    prefix = prefix.reset_index(drop=True)
    column_matches = {
        column: bool(np.array_equal(old[column].to_numpy(), prefix[column].to_numpy()))
        for column in old.columns
    }
    result: dict[str, object] = {
        "old_max_sweep": old_max_sweep,
        "old_rows": len(old),
        "new_prefix_rows": len(prefix),
        "columns_match_exactly": column_matches,
        "all_trace_columns_match_exactly": (
            len(old) == len(prefix) and all(column_matches.values())
        ),
    }
    if not result["all_trace_columns_match_exactly"]:
        raise RuntimeError(f"long-run trace prefix differs from {old_case}")
    return result


def _old_or_new_observable_metrics(details: Mapping[str, object]) -> dict[str, float]:
    def numeric(value: object) -> float:
        return float("nan") if value is None else float(value)

    horizontal = details["horizontal_start"]
    vertical = details["vertical_start"]
    assert isinstance(horizontal, Mapping) and isinstance(vertical, Mapping)
    gap = float(details["late_mean_gap_horizontal_minus_vertical"])
    gap_se = float(details["late_mean_gap_standard_error"])
    if "z_gap_horizontal_minus_vertical" in details:
        z_gap = float(details["z_gap_horizontal_minus_vertical"])
    elif gap_se > 0:
        z_gap = gap / gap_se
    else:
        z_gap = 0.0 if gap == 0 else float(np.copysign(np.inf, gap))
    return {
        "iact_sweeps_max": max(
            float(horizontal["integrated_autocorrelation_time_sweeps"]),
            float(vertical["integrated_autocorrelation_time_sweeps"]),
        ),
        "ess_min": min(
            float(horizontal["effective_sample_size_late"]),
            float(vertical["effective_sample_size_late"]),
        ),
        "rhat": numeric(details["split_rhat_late_half"]),
        "z_gap": z_gap,
        "mean_gap": gap,
        "variance_gap": float(
            details["late_variance_difference_horizontal_minus_vertical"]
        ),
    }


def build_old_new_comparison(
    source: str | Path,
    output: str | Path,
) -> pd.DataFrame:
    """Build the requested old/new table, optionally merging manual assessments."""
    source_path = Path(source)
    output_path = Path(output)
    hardness = pd.read_csv(output_path / "environment_hardness.csv").set_index(
        "environment_id"
    )
    case_summary = pd.read_csv(output_path / "case_summary.csv").set_index(
        "environment_id"
    )
    rows: list[dict[str, object]] = []
    for environment_id in sorted(hardness.index, key=lambda item: (int(item.split("L")[1].split("-")[0]), item)):
        old = json.loads((source_path / environment_id / "summary.json").read_text())
        new = json.loads((output_path / environment_id / "summary.json").read_text())
        row: dict[str, object] = {
            "environment_id": environment_id,
            "L": int(new["L"]),
            "environment_seed": int(new["environment_seed"]),
            "chain1_seed": int(new["chain1_seed"]),
            "chain2_seed": int(new["chain2_seed"]),
            "environment_sha256": new["environment_sha256"],
            "old_sweeps": int(old["sweeps"]),
            "new_sweeps": int(new["sweeps"]),
        }
        for observable, short_name in (
            ("center_height", "center"),
            ("height_sum", "global"),
        ):
            for run_name, summary in (("old", old), ("new", new)):
                metrics = _old_or_new_observable_metrics(summary["observables"][observable])
                row.update(
                    {
                        f"{run_name}_{short_name}_{key}": value
                        for key, value in metrics.items()
                    }
                )
        row.update(
            {
                key: value
                for key, value in hardness.loc[environment_id].to_dict().items()
                if key != "L"
            }
        )
        row["new_pair_runtime_seconds"] = float(
            case_summary.loc[environment_id, "total_chain_seconds"]
        )
        row["new_case_wall_seconds"] = float(
            case_summary.loc[environment_id, "case_wall_seconds"]
        )
        row["start_state_separation_disappeared"] = "pending visual review"
        row["visual_assessment"] = "pending visual review"
        row["credible"] = "pending visual review"
        row["recommended_next_action"] = "pending visual review"
        rows.append(row)
    result = pd.DataFrame(rows)
    assessment_path = output_path / "visual_assessment.csv"
    if assessment_path.exists():
        assessments = pd.read_csv(assessment_path)
        assessment_columns = [
            "start_state_separation_disappeared",
            "visual_assessment",
            "credible",
            "recommended_next_action",
        ]
        result = result.drop(columns=assessment_columns).merge(
            assessments[["environment_id", *assessment_columns]],
            on="environment_id",
            how="left",
            validate="one_to_one",
        )
    result.to_csv(output_path / "old_new_comparison.csv", index=False)
    if assessment_path.exists():
        size_summary = (
            result.assign(credible_bool=result["credible"].eq("yes"))
            .groupby("L", as_index=False)
            .agg(
                environments=("environment_id", "count"),
                credible_environments=("credible_bool", "sum"),
                all_environments_credible=("credible_bool", "all"),
            )
        )
        size_summary.to_csv(output_path / "size_credibility_summary.csv", index=False)
    return result


def run_frozen_gamma_extension(
    *,
    source: str | Path,
    output: str | Path,
    sweeps_by_size: Mapping[int, int],
    record_every: int = 10,
    relative_separations: Iterable[float] = DEFAULT_RELATIVE_SEPARATIONS,
    max_lag_records: int = 5_000,
) -> dict[str, object]:
    """Run longer extremal-start chains in the exact saved Gamma environments."""
    source_path = Path(source)
    output_path = Path(output)
    if not source_path.is_dir():
        raise ValueError(f"source campaign directory does not exist: {source_path}")
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty output directory: {output_path}")
    if record_every <= 0 or max_lag_records < 0:
        raise ValueError("record_every must be positive and max_lag_records nonnegative")
    schedule = {validate_L(int(L)): int(sweeps) for L, sweeps in sweeps_by_size.items()}
    if any(sweeps <= 0 for sweeps in schedule.values()):
        raise ValueError("all sweep counts must be positive")
    relative_separations = tuple(float(value) for value in relative_separations)
    output_path.mkdir(parents=True, exist_ok=False)

    case_rows: list[dict[str, object]] = []
    observable_rows: list[dict[str, object]] = []
    hardness_rows: list[dict[str, object]] = []
    verification_rows: list[dict[str, object]] = []
    campaign_started = perf_counter()
    source_metadata = json.loads((source_path / "campaign_metadata.json").read_text())

    source_cases: list[Path] = []
    for L in sorted(schedule):
        cases = sorted(source_path.glob(f"gamma-L{L}-diagnostic-*"))
        if len(cases) != 3:
            raise ValueError(f"expected exactly 3 saved Gamma environments for L={L}, found {len(cases)}")
        source_cases.extend(cases)

    for case_number, old_case in enumerate(source_cases, start=1):
        environment, old_summary, verification = _load_frozen_environment(old_case)
        L = environment.L
        environment_id = str(old_summary["environment_id"])
        print(
            f"[{case_number}/{len(source_cases)}] {environment_id}: "
            f"verified frozen environment; running {schedule[L]} sweeps per chain",
            flush=True,
        )
        case_started = perf_counter()
        comparison = compare_extremal_starts(
            environment,
            np.random.Generator(np.random.PCG64(int(old_summary["chain1_seed"]))),
            np.random.Generator(np.random.PCG64(int(old_summary["chain2_seed"]))),
            sweeps=schedule[L],
            record_every=record_every,
            separations=spatial_separations(L, relative_separations),
        )
        case_directory = output_path / environment_id
        case_directory.mkdir(parents=True, exist_ok=False)
        np.savez_compressed(
            case_directory / "environment_weights.npz",
            horizontal_weights=environment.horizontal_weights,
            vertical_weights=environment.vertical_weights,
        )
        hashes = environment_hashes(environment)
        hardness = environment_hardness(environment)
        metadata: dict[str, object] = {
            "environment_id": environment_id,
            "L": L,
            "model": "gamma",
            "gamma_shape": 1.0,
            "environment_index": int(old_summary["environment_index"]),
            "master_seed": int(old_summary["master_seed"]),
            "environment_seed": int(old_summary["environment_seed"]),
            "chain1_seed": int(old_summary["chain1_seed"]),
            "chain2_seed": int(old_summary["chain2_seed"]),
            "sweeps": schedule[L],
            "source_sweeps": int(old_summary["sweeps"]),
            "source_campaign": str(source_path),
            "actual_separations": old_summary["actual_separations"],
            **hashes,
        }
        case_row, rows = write_mixing_case(
            comparison,
            case_directory,
            metadata=metadata,
            max_lag_records=max_lag_records,
        )
        prefix_verification = _verify_trace_prefix(old_case, case_directory)
        (case_directory / "environment_verification.json").write_text(
            json.dumps(verification, indent=2) + "\n"
        )
        (case_directory / "environment_hardness.json").write_text(
            json.dumps(hardness, indent=2) + "\n"
        )
        (case_directory / "trace_prefix_verification.json").write_text(
            json.dumps(prefix_verification, indent=2) + "\n"
        )
        case_row.update(hardness)
        case_row["case_wall_seconds"] = perf_counter() - case_started
        case_row["source_trace_prefix_matches_exactly"] = prefix_verification[
            "all_trace_columns_match_exactly"
        ]
        case_rows.append(case_row)
        observable_rows.extend(rows)
        hardness_rows.append({"environment_id": environment_id, "L": L, **hardness})
        verification_rows.append({"environment_id": environment_id, "L": L, **verification})

        pd.DataFrame(case_rows).to_csv(output_path / "case_summary.csv", index=False)
        pd.DataFrame(observable_rows).to_csv(
            output_path / "observable_summary.csv", index=False
        )
        pd.DataFrame(hardness_rows).to_csv(
            output_path / "environment_hardness.csv", index=False
        )
        pd.DataFrame(verification_rows).to_csv(
            output_path / "environment_verification.csv", index=False
        )
        print(
            f"[{case_number}/{len(source_cases)}] {environment_id}: "
            f"complete in {case_row['case_wall_seconds']:.3f} s",
            flush=True,
        )

    case_frame = pd.DataFrame(case_rows)
    runtime_summary = (
        case_frame.groupby("L", as_index=False)
        .agg(
            cases=("environment_id", "count"),
            total_chain_seconds=("total_chain_seconds", "sum"),
            mean_chain_pair_seconds=("total_chain_seconds", "mean"),
            total_case_wall_seconds=("case_wall_seconds", "sum"),
            mean_case_wall_seconds=("case_wall_seconds", "mean"),
            mean_horizontal_flippable_rate=(
                "horizontal_flippable_proposal_rate",
                "mean",
            ),
            mean_vertical_flippable_rate=("vertical_flippable_proposal_rate", "mean"),
            mean_horizontal_change_rate=("horizontal_actual_change_rate", "mean"),
            mean_vertical_change_rate=("vertical_actual_change_rate", "mean"),
        )
    )
    runtime_summary.to_csv(output_path / "runtime_summary.csv", index=False)
    campaign = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_campaign": str(source_path),
        "source_campaign_master_seed": int(source_metadata["master_seed"]),
        "source_sweeps_per_chain": int(source_metadata["sweeps_per_chain"]),
        "model": "gamma",
        "gamma_shape": 1.0,
        "sweeps_by_size": {str(L): sweeps for L, sweeps in schedule.items()},
        "record_every_sweeps": record_every,
        "acf_max_lag_records": max_lag_records,
        "post_transient_fraction": 0.5,
        "post_transient_windows": 4,
        "environment_count": len(case_frame),
        "wall_seconds": perf_counter() - campaign_started,
        "all_seed_and_array_verifications_passed": bool(
            pd.DataFrame(verification_rows)[
                [
                    "recorded_seed_tuple_matches_current_derivation",
                    "horizontal_arrays_bitwise_equal_to_seed_regeneration",
                    "vertical_arrays_bitwise_equal_to_seed_regeneration",
                    "combined_hashes_match",
                ]
            ].all(axis=None)
        ),
        "all_20k_trace_prefixes_match_exactly": bool(
            case_frame["source_trace_prefix_matches_exactly"].all()
        ),
        "warning": (
            "Exploratory diagnostics only. Numerical thresholds are not by "
            "themselves a convergence declaration, and no environment was removed."
        ),
    }
    (output_path / "campaign_metadata.json").write_text(
        json.dumps(campaign, indent=2) + "\n"
    )
    build_old_new_comparison(source_path, output_path)
    return campaign
