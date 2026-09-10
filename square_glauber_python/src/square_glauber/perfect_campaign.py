"""Restart-safe paired perfect-sampling campaign for Hamilton production.

Every scientific row is one frozen edge environment with two conditionally
independent, CFTP-certified dimer samples.  A pair is never emitted until both
samples certify.  Environments and individual sample states are persisted
before aggregation so censored streams can be extended without redrawing or
selection on computational hardness.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable

import numpy as np
import pandas as pd

from .environment import DimerEnvironment, sample_environment
from .experiment import MODEL_CODES, SeedBundle, derive_seeds
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, spatial_separations, validate_L
from .height import face_heights, spatial_increments
from .matching import DimerState
from .mixing_extension import array_sha256, environment_hardness, environment_hashes
from .perfect import CFTPNotCoalesced, MonotoneCFTPSampler, PerfectBackend, PerfectSample


MANIFEST_COLUMNS = (
    "task_id",
    "campaign_id",
    "master_seed",
    "model",
    "gamma_shape",
    "L",
    "environment_start",
    "environment_stop",
    "environment_count",
    "max_horizon_sweeps",
    "backend",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    temporary.replace(path)


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def create_campaign_manifest(
    schedule: str | Path,
    output: str | Path,
    *,
    campaign_id: str,
    master_seed: int,
    backend: PerfectBackend = "numba",
) -> pd.DataFrame:
    """Expand a size/model schedule into deterministic disjoint task rows."""
    if not campaign_id or any(character.isspace() for character in campaign_id):
        raise ValueError("campaign_id must be a nonempty identifier without whitespace")
    if master_seed < 0:
        raise ValueError("master_seed must be nonnegative")
    if backend not in ("reference", "numba"):
        raise ValueError("backend must be 'reference' or 'numba'")
    source = pd.read_csv(schedule)
    required = {
        "L",
        "model",
        "gamma_shape",
        "environments",
        "environments_per_task",
        "max_horizon_sweeps",
    }
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"schedule is missing columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    coverage: set[tuple[str, int, int]] = set()
    for item in source.itertuples(index=False):
        L = validate_L(int(item.L))
        model = str(item.model)
        if model not in MODEL_CODES:
            raise ValueError(f"unknown model {model!r}")
        gamma_shape = float(item.gamma_shape)
        if model == "gamma" and (not np.isfinite(gamma_shape) or gamma_shape <= 0):
            raise ValueError("Gamma shape must be finite and positive")
        total = int(item.environments)
        batch = int(item.environments_per_task)
        maximum = int(item.max_horizon_sweeps)
        if total <= 0 or batch <= 0 or maximum <= 0:
            raise ValueError("environment counts, batch sizes, and horizons must be positive")
        for start in range(0, total, batch):
            stop = min(start + batch, total)
            for environment_index in range(start, stop):
                key = (model, L, environment_index)
                if key in coverage:
                    raise ValueError(f"schedule contains duplicate environment {key}")
                coverage.add(key)
            rows.append(
                {
                    "task_id": len(rows),
                    "campaign_id": campaign_id,
                    "master_seed": int(master_seed),
                    "model": model,
                    "gamma_shape": gamma_shape if model == "gamma" else np.nan,
                    "L": L,
                    "environment_start": start,
                    "environment_stop": stop,
                    "environment_count": stop - start,
                    "max_horizon_sweeps": maximum,
                    "backend": backend,
                }
            )
    if not rows:
        raise ValueError("schedule produced no tasks")
    environment_seeds: set[int] = set()
    cftp_seeds: set[int] = set()
    for model, L, environment_index in coverage:
        seeds = derive_seeds(master_seed, L, model, environment_index)
        if seeds.environment_seed in environment_seeds:
            raise RuntimeError("derived environment-seed collision")
        if (
            seeds.chain1_seed == seeds.chain2_seed
            or seeds.chain1_seed in cftp_seeds
            or seeds.chain2_seed in cftp_seeds
        ):
            raise RuntimeError("derived CFTP-seed collision")
        environment_seeds.add(seeds.environment_seed)
        cftp_seeds.update((seeds.chain1_seed, seeds.chain2_seed))
    frame = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    target = Path(output)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing manifest: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(target, index=False)
    metadata = {
        "created_utc": _utc_now(),
        "campaign_id": campaign_id,
        "master_seed": int(master_seed),
        "backend": backend,
        "schedule": str(Path(schedule)),
        "schedule_sha256": file_sha256(schedule),
        "manifest": str(target),
        "manifest_sha256": file_sha256(target),
        "task_count": len(frame),
        "environment_count": int(frame["environment_count"].sum()),
        "perfect_sample_count": int(2 * frame["environment_count"].sum()),
        "all_environment_seeds_unique": len(environment_seeds) == len(coverage),
        "all_cftp_seeds_unique": len(cftp_seeds) == 2 * len(coverage),
        "statistical_unit": "one frozen environment with two independent certified CFTP samples",
    }
    _atomic_json(target.with_suffix(".json"), metadata)
    return frame


def _case_path(root: Path, model: str, L: int, environment_index: int) -> Path:
    return root / "cases" / model / f"L{L:03d}" / f"environment-{environment_index:06d}"


def _load_environment(path: Path, *, L: int, model: str, gamma_shape: float) -> DimerEnvironment:
    with np.load(path, allow_pickle=False) as archive:
        horizontal = np.array(archive["horizontal_weights"], copy=True)
        vertical = np.array(archive["vertical_weights"], copy=True)
    return DimerEnvironment(
        L,
        horizontal,
        vertical,
        model=model,
        gamma_shape=gamma_shape if model == "gamma" else None,
    )


def _prepare_environment(
    case: Path,
    *,
    L: int,
    model: str,
    gamma_shape: float,
    environment_index: int,
    seeds: SeedBundle,
    campaign_id: str,
) -> tuple[DimerEnvironment, dict[str, object]]:
    weights_path = case / "environment_weights.npz"
    metadata_path = case / "environment.json"
    regenerated = sample_environment(
        L,
        np.random.Generator(np.random.PCG64(seeds.environment_seed)),
        model,
        gamma_shape=gamma_shape,
    )
    if weights_path.exists():
        environment = _load_environment(
            weights_path, L=L, model=model, gamma_shape=gamma_shape
        )
        if not (
            np.array_equal(environment.horizontal_weights, regenerated.horizontal_weights)
            and np.array_equal(environment.vertical_weights, regenerated.vertical_weights)
        ):
            raise RuntimeError(f"stored environment does not regenerate bitwise: {case}")
    else:
        environment = regenerated
        _atomic_npz(
            weights_path,
            horizontal_weights=environment.horizontal_weights,
            vertical_weights=environment.vertical_weights,
        )
    verification = {
        "campaign_id": campaign_id,
        "environment_id": f"{model}-L{L}-env{environment_index:06d}",
        "environment_index": environment_index,
        "L": L,
        "model": model,
        "gamma_shape": gamma_shape if model == "gamma" else None,
        **asdict(seeds),
        **environment_hashes(environment),
        **environment_hardness(environment),
        "weights_npz_sha256": file_sha256(weights_path),
        "reload_bitwise_equal": True,
        "seed_regeneration_bitwise_equal": True,
    }
    if metadata_path.exists():
        stored = json.loads(metadata_path.read_text())
        invariant_keys = (
            "campaign_id",
            "environment_id",
            "environment_seed",
            "environment_sha256",
            "L",
            "model",
        )
        if any(stored.get(key) != verification.get(key) for key in invariant_keys):
            raise RuntimeError(f"environment metadata conflict in {metadata_path}")
    else:
        _atomic_json(metadata_path, verification)
    return environment, verification


def _state_from_npz(path: Path, L: int) -> DimerState:
    with np.load(path, allow_pickle=False) as archive:
        state = DimerState(
            L,
            np.asarray(archive["horizontal_occupied"], dtype=np.bool_),
            np.asarray(archive["vertical_occupied"], dtype=np.bool_),
        )
    state.validate()
    return state


def _sample_metadata(
    result: PerfectSample,
    *,
    environment_sha256: str,
    sample_index: int,
    separations: tuple,
) -> dict[str, object]:
    heights = face_heights(result.state)
    increments = spatial_increments(heights, separations)
    lower = result.lower_update_counts
    upper = result.upper_update_counts
    lower_attempts = lower.attempted_updates
    upper_attempts = upper.attempted_updates
    return {
        **result.metadata(),
        "sample_index": sample_index,
        "environment_sha256": environment_sha256,
        "center_height": int(heights[(result.state.L - 2) // 2, (result.state.L - 2) // 2]),
        "height_sum": int(heights.sum()),
        "height_mean": float(heights.mean()),
        "spatial_increments": {str(key): int(value) for key, value in increments.items()},
        "horizontal_occupied_sha256": array_sha256(result.state.horizontal_occupied),
        "vertical_occupied_sha256": array_sha256(result.state.vertical_occupied),
        "lower_flippable_proposals": lower.flippable_proposals,
        "lower_actual_changes": lower.actual_changes,
        "lower_flippable_rate": lower.flippable_proposals / lower_attempts,
        "lower_change_rate": lower.actual_changes / lower_attempts,
        "upper_flippable_proposals": upper.flippable_proposals,
        "upper_actual_changes": upper.actual_changes,
        "upper_flippable_rate": upper.flippable_proposals / upper_attempts,
        "upper_change_rate": upper.actual_changes / upper_attempts,
    }


def _load_certified_sample(
    case: Path,
    *,
    sample_index: int,
    seed: int,
    L: int,
    environment_sha256: str,
    separations: tuple,
) -> tuple[DimerState, dict[str, object]] | None:
    metadata_path = case / f"sample{sample_index}.json"
    state_path = case / f"sample{sample_index}_state.npz"
    if not metadata_path.exists():
        return None
    if not state_path.exists():
        raise RuntimeError(f"certified metadata exists without state: {metadata_path}")
    metadata = json.loads(metadata_path.read_text())
    if (
        not metadata.get("certified_perfect_sample")
        or int(metadata["cftp_seed"]) != int(seed)
        or metadata["environment_sha256"] != environment_sha256
    ):
        raise RuntimeError(f"certified sample metadata conflict: {metadata_path}")
    state = _state_from_npz(state_path, L)
    if (
        array_sha256(state.horizontal_occupied) != metadata["horizontal_occupied_sha256"]
        or array_sha256(state.vertical_occupied) != metadata["vertical_occupied_sha256"]
    ):
        raise RuntimeError(f"certified sample state hash conflict: {state_path}")
    heights = face_heights(state)
    increments = {
        str(key): value for key, value in spatial_increments(heights, separations).items()
    }
    center = (L - 2) // 2
    if (
        int(metadata["center_height"]) != int(heights[center, center])
        or int(metadata["height_sum"]) != int(heights.sum())
        or metadata["spatial_increments"] != increments
    ):
        raise RuntimeError(f"saved observables do not match certified state: {state_path}")
    return state, metadata


def _obtain_sample(
    sampler: MonotoneCFTPSampler,
    case: Path,
    *,
    sample_index: int,
    seed: int,
    max_horizon_sweeps: int,
    environment_sha256: str,
    separations: tuple,
) -> tuple[DimerState | None, dict[str, object]]:
    existing = _load_certified_sample(
        case,
        sample_index=sample_index,
        seed=seed,
        L=sampler.environment.L,
        environment_sha256=environment_sha256,
        separations=separations,
    )
    if existing is not None:
        return existing
    censor_directory = case / "censoring"
    censor_path = censor_directory / f"sample{sample_index}_cap-{max_horizon_sweeps}.json"
    if censor_path.exists():
        return None, json.loads(censor_path.read_text())
    try:
        result = sampler.sample(
            seed,
            max_horizon_sweeps=max_horizon_sweeps,
            verify_order=True,
        )
    except CFTPNotCoalesced as error:
        metadata = {
            "sample_index": sample_index,
            "cftp_seed": int(seed),
            "certified_perfect_sample": False,
            "censored": True,
            "environment_sha256": environment_sha256,
            "max_horizon_sweeps": max_horizon_sweeps,
            "horizon_sweeps": error.horizon_sweeps,
            "iterations": [asdict(value) for value in error.iterations],
            "instruction": "extend this identical seed and environment; never discard the case",
        }
        _atomic_json(censor_path, metadata)
        return None, metadata
    metadata = _sample_metadata(
        result,
        environment_sha256=environment_sha256,
        sample_index=sample_index,
        separations=separations,
    )
    _atomic_npz(
        case / f"sample{sample_index}_state.npz",
        horizontal_occupied=result.state.horizontal_occupied,
        vertical_occupied=result.state.vertical_occupied,
    )
    _atomic_json(case / f"sample{sample_index}.json", metadata)
    return result.state, metadata


def _pair_metadata(
    environment: dict[str, object],
    first: dict[str, object],
    second: dict[str, object],
    *,
    separations: tuple,
) -> dict[str, object]:
    if int(first["cftp_seed"]) == int(second["cftp_seed"]):
        raise RuntimeError("paired perfect samples unexpectedly share a seed")
    return {
        **environment,
        "created_utc": _utc_now(),
        "sampling_method": "monotone CFTP with true pointwise height extrema",
        "kernel": "random scan of all bounded faces with null moves and fixed-environment heat bath",
        "pair_conditionally_independent_given_environment": True,
        "both_samples_cftp_certified": True,
        "sample1": first,
        "sample2": second,
        "separations": [asdict(item) for item in separations],
    }


def run_campaign_task(
    manifest: str | Path,
    task_id: int,
    output: str | Path,
    *,
    max_horizon_override: int | None = None,
) -> dict[str, object]:
    """Run or resume one manifest task without ever redrawing a saved case."""
    manifest_path = Path(manifest)
    frame = pd.read_csv(manifest_path)
    selected = frame.loc[frame["task_id"] == int(task_id)]
    if len(selected) != 1:
        raise ValueError(f"task_id {task_id} appears {len(selected)} times in manifest")
    row = selected.iloc[0]
    model = str(row["model"])
    L = validate_L(int(row["L"]))
    gamma_shape = float(row["gamma_shape"]) if model == "gamma" else 1.0
    maximum = (
        int(max_horizon_override)
        if max_horizon_override is not None
        else int(row["max_horizon_sweeps"])
    )
    if maximum <= 0:
        raise ValueError("maximum CFTP horizon must be positive")
    backend = str(row["backend"])
    if backend not in ("reference", "numba"):
        raise ValueError(f"invalid backend {backend!r}")
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    separations = spatial_separations(L, DEFAULT_RELATIVE_SEPARATIONS)
    started = perf_counter()
    cases: list[dict[str, object]] = []
    for environment_index in range(
        int(row["environment_start"]), int(row["environment_stop"])
    ):
        case_started = perf_counter()
        case = _case_path(root, model, L, environment_index)
        case.mkdir(parents=True, exist_ok=True)
        seeds = derive_seeds(int(row["master_seed"]), L, model, environment_index)
        environment, environment_metadata = _prepare_environment(
            case,
            L=L,
            model=model,
            gamma_shape=gamma_shape,
            environment_index=environment_index,
            seeds=seeds,
            campaign_id=str(row["campaign_id"]),
        )
        sampler = MonotoneCFTPSampler(environment, backend=backend)  # type: ignore[arg-type]
        first_state, first = _obtain_sample(
            sampler,
            case,
            sample_index=1,
            seed=seeds.chain1_seed,
            max_horizon_sweeps=maximum,
            environment_sha256=str(environment_metadata["environment_sha256"]),
            separations=separations,
        )
        second_state, second = _obtain_sample(
            sampler,
            case,
            sample_index=2,
            seed=seeds.chain2_seed,
            max_horizon_sweeps=maximum,
            environment_sha256=str(environment_metadata["environment_sha256"]),
            separations=separations,
        )
        complete = first_state is not None and second_state is not None
        if complete:
            first_state.validate()
            second_state.validate()
            pair_path = case / "pair.json"
            pair = _pair_metadata(
                environment_metadata, first, second, separations=separations
            )
            if pair_path.exists():
                stored = json.loads(pair_path.read_text())
                if (
                    stored["environment_sha256"] != pair["environment_sha256"]
                    or stored["sample1"]["cftp_seed"] != pair["sample1"]["cftp_seed"]
                    or stored["sample2"]["cftp_seed"] != pair["sample2"]["cftp_seed"]
                ):
                    raise RuntimeError(f"pair metadata conflict: {pair_path}")
            else:
                _atomic_json(pair_path, pair)
        cases.append(
            {
                "environment_id": environment_metadata["environment_id"],
                "complete": complete,
                "sample1_certified": first_state is not None,
                "sample2_certified": second_state is not None,
                "case_elapsed_seconds": perf_counter() - case_started,
            }
        )
    report = {
        "created_utc": _utc_now(),
        "task_id": int(task_id),
        "campaign_id": str(row["campaign_id"]),
        "manifest_sha256": file_sha256(manifest_path),
        "L": L,
        "model": model,
        "environment_start": int(row["environment_start"]),
        "environment_stop": int(row["environment_stop"]),
        "max_horizon_sweeps": maximum,
        "backend": backend,
        "environment_count": len(cases),
        "complete_pair_count": sum(bool(item["complete"]) for item in cases),
        "incomplete_pair_count": sum(not bool(item["complete"]) for item in cases),
        "elapsed_seconds": perf_counter() - started,
        "cases": cases,
    }
    _atomic_json(root / "tasks" / f"task-{int(task_id):06d}.json", report)
    return report


def _expected_cases(manifest: pd.DataFrame) -> Iterable[dict[str, object]]:
    for row in manifest.itertuples(index=False):
        for environment_index in range(int(row.environment_start), int(row.environment_stop)):
            yield {
                "campaign_id": str(row.campaign_id),
                "master_seed": int(row.master_seed),
                "model": str(row.model),
                "gamma_shape": float(row.gamma_shape) if str(row.model) == "gamma" else np.nan,
                "L": int(row.L),
                "environment_index": environment_index,
                "environment_id": f"{row.model}-L{int(row.L)}-env{environment_index:06d}",
                "task_id": int(row.task_id),
            }


def _slurm_array_spec(values: Iterable[int]) -> str:
    ordered = sorted(set(int(value) for value in values))
    if not ordered:
        return ""
    groups: list[str] = []
    start = previous = ordered[0]
    for value in ordered[1:]:
        if value == previous + 1:
            previous = value
            continue
        groups.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    groups.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(groups)


def aggregate_campaign(
    manifest: str | Path,
    output: str | Path,
    *,
    require_complete: bool = True,
) -> dict[str, object]:
    """Audit all expected cases and emit environment-blocked analysis tables."""
    manifest_path = Path(manifest)
    root = Path(output)
    manifest_frame = pd.read_csv(manifest_path)
    center_rows: list[dict[str, object]] = []
    spatial_rows: list[dict[str, object]] = []
    diagnostic_rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    for expected in _expected_cases(manifest_frame):
        model = str(expected["model"])
        L = int(expected["L"])
        environment_index = int(expected["environment_index"])
        case = _case_path(root, model, L, environment_index)
        pair_path = case / "pair.json"
        status = "missing"
        if pair_path.exists():
            pair = json.loads(pair_path.read_text())
            expected_id = str(expected["environment_id"])
            if (
                pair["environment_id"] != expected_id
                or pair["campaign_id"] != expected["campaign_id"]
                or int(pair["master_seed"]) != int(expected["master_seed"])
                or int(pair["L"]) != L
                or pair["model"] != model
                or not pair["both_samples_cftp_certified"]
            ):
                raise RuntimeError(f"pair invariant failed: {pair_path}")
            stored_environment = json.loads((case / "environment.json").read_text())
            stored_first = json.loads((case / "sample1.json").read_text())
            stored_second = json.loads((case / "sample2.json").read_text())
            if (
                pair["environment_sha256"] != stored_environment["environment_sha256"]
                or pair["sample1"] != stored_first
                or pair["sample2"] != stored_second
            ):
                raise RuntimeError(f"pair does not match its component artifacts: {pair_path}")
            status = "complete"
            first = pair["sample1"]
            second = pair["sample2"]
            base = {
                "master_seed": int(pair["master_seed"]),
                "environment_seed": int(pair["environment_seed"]),
                "chain1_seed": int(pair["chain1_seed"]),
                "chain2_seed": int(pair["chain2_seed"]),
                "L": L,
                "model": model,
                "gamma_shape": pair["gamma_shape"],
                "environment_index": environment_index,
                "environment_id": pair["environment_id"],
                "environment_sha256": pair["environment_sha256"],
                "chain1_start": "certified_cftp",
                "chain2_start": "certified_cftp",
            }
            center_rows.append(
                base
                | {
                    "measurement_index": 0,
                    "H1_center": int(first["center_height"]),
                    "H2_center": int(second["center_height"]),
                    "H1_global_sum": int(first["height_sum"]),
                    "H2_global_sum": int(second["height_sum"]),
                }
            )
            for separation in pair["separations"]:
                key = str(int(separation["r"]))
                spatial_rows.append(
                    base
                    | {
                        "measurement_index": 0,
                        "requested_rho": float(separation["requested_rho"]),
                        "r": int(separation["r"]),
                        "DeltaH1": int(first["spatial_increments"][key]),
                        "DeltaH2": int(second["spatial_increments"][key]),
                    }
                )
            diagnostic = base | {
                key: value
                for key, value in pair.items()
                if key.startswith("edge_")
                or key.startswith("log_edge_")
                or key.startswith("abs_lambda_")
                or key.startswith("fraction_abs_lambda_")
            }
            for prefix, sample in (("chain1", first), ("chain2", second)):
                for source_key, target_key in (
                    ("coalescence_horizon_sweeps", "cftp_horizon_sweeps"),
                    ("elapsed_seconds", "runtime_seconds"),
                    ("lower_total_replayed_attempts", "lower_replayed_attempts"),
                    ("upper_total_replayed_attempts", "upper_replayed_attempts"),
                    ("lower_flippable_rate", "lower_flippable_rate"),
                    ("lower_change_rate", "lower_change_rate"),
                    ("upper_flippable_rate", "upper_flippable_rate"),
                    ("upper_change_rate", "upper_change_rate"),
                ):
                    diagnostic[f"{prefix}_{target_key}"] = sample[source_key]
            diagnostic_rows.append(diagnostic)
        elif case.exists() and any((case / "censoring").glob("*.json")):
            status = "censored"
        status_rows.append(expected | {"status": status, "case_path": str(case)})

    center = pd.DataFrame(center_rows)
    spatial = pd.DataFrame(spatial_rows)
    diagnostics = pd.DataFrame(diagnostic_rows)
    statuses = pd.DataFrame(status_rows)
    root.mkdir(parents=True, exist_ok=True)
    center.to_csv(root / "center_measurements.csv", index=False)
    spatial.to_csv(root / "spatial_measurements.csv", index=False)
    diagnostics.to_csv(root / "chain_diagnostics.csv", index=False)
    statuses.to_csv(root / "campaign_status.csv", index=False)
    status_counts = statuses["status"].value_counts().to_dict()
    incomplete = statuses.loc[statuses["status"] != "complete"]
    incomplete_task_ids = sorted(set(incomplete["task_id"].astype(int)))
    (root / "incomplete_task_ids.txt").write_text(
        "".join(f"{value}\n" for value in incomplete_task_ids)
    )
    (root / "incomplete_slurm_array_spec.txt").write_text(
        _slurm_array_spec(incomplete_task_ids) + ("\n" if incomplete_task_ids else "")
    )
    (root / "censored_environment_ids.txt").write_text(
        "".join(
            f"{value}\n"
            for value in statuses.loc[
                statuses["status"] == "censored", "environment_id"
            ]
        )
    )
    report = {
        "created_utc": _utc_now(),
        "campaign_id": str(manifest_frame.iloc[0]["campaign_id"]),
        "manifest": str(manifest_path),
        "manifest_sha256": file_sha256(manifest_path),
        "expected_environments": len(statuses),
        "complete_environments": int(status_counts.get("complete", 0)),
        "censored_environments": int(status_counts.get("censored", 0)),
        "missing_environments": int(status_counts.get("missing", 0)),
        "incomplete_task_count": len(incomplete_task_ids),
        "incomplete_slurm_array_spec": _slurm_array_spec(incomplete_task_ids),
        "center_rows": len(center),
        "spatial_rows": len(spatial),
        "all_pairs_cftp_certified": len(center) == len(statuses),
        "statistical_unit": "whole frozen environment; one conditionally independent certified pair",
        "selection_policy": "no environment is discarded; censored cases must be extended with identical arrays and seeds",
    }
    _atomic_json(root / "campaign_summary.json", report)
    if require_complete and len(center) != len(statuses):
        raise RuntimeError(
            f"campaign incomplete: {report['complete_environments']} complete, "
            f"{report['censored_environments']} censored, {report['missing_environments']} missing"
        )
    return report
