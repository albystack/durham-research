"""Shared-environment paired-chain pilot and reproducible output contract."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Iterable

import numpy as np
import pandas as pd

from .environment import DimerEnvironment, sample_environment
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, SpatialSeparation, spatial_separations, validate_L
from .glauber import ChainResult, run_chain
from .matching import DimerState

MODEL_CODES = {"uniform": 0, "gamma": 1}


@dataclass(frozen=True)
class SeedBundle:
    master_seed: int
    environment_seed: int
    chain1_seed: int
    chain2_seed: int


def _seed_uint64(sequence: np.random.SeedSequence) -> int:
    return int(sequence.generate_state(1, dtype=np.uint64)[0])


def derive_seeds(master_seed: int, L: int, model: str, environment_index: int) -> SeedBundle:
    """Derive stable, independent PCG64 seeds without sharing RNG objects."""
    validate_L(L)
    if model not in MODEL_CODES:
        raise ValueError(f"unknown model {model!r}")
    if master_seed < 0 or environment_index < 0:
        raise ValueError("master_seed and environment_index must be nonnegative")
    root = np.random.SeedSequence([master_seed, L, MODEL_CODES[model], environment_index])
    environment_sequence, chain1_sequence, chain2_sequence = root.spawn(3)
    return SeedBundle(
        int(master_seed),
        _seed_uint64(environment_sequence),
        _seed_uint64(chain1_sequence),
        _seed_uint64(chain2_sequence),
    )


@dataclass
class PairedRun:
    environment: DimerEnvironment
    seeds: SeedBundle
    chain1: ChainResult
    chain2: ChainResult
    separations: tuple[SpatialSeparation, ...]


def run_paired_environment(
    *,
    L: int,
    model: str,
    gamma_shape: float,
    environment_index: int,
    master_seed: int,
    burnin_sweeps: int,
    measurement_gap_sweeps: int,
    num_measurements: int,
    relative_separations: Iterable[float] = DEFAULT_RELATIVE_SEPARATIONS,
) -> PairedRun:
    """Run two independent chains against the exact same immutable arrays."""
    seeds = derive_seeds(master_seed, L, model, environment_index)
    environment_rng = np.random.Generator(np.random.PCG64(seeds.environment_seed))
    chain1_rng = np.random.Generator(np.random.PCG64(seeds.chain1_seed))
    chain2_rng = np.random.Generator(np.random.PCG64(seeds.chain2_seed))
    environment = sample_environment(
        L, environment_rng, model, gamma_shape=gamma_shape
    )
    separations = spatial_separations(L, relative_separations)
    common_arguments = dict(
        environment=environment,
        burnin_sweeps=burnin_sweeps,
        measurement_gap_sweeps=measurement_gap_sweeps,
        num_measurements=num_measurements,
        separations=separations,
    )
    chain1 = run_chain(
        initial_matching=DimerState.all_horizontal(L),
        rng=chain1_rng,
        **common_arguments,
    )
    chain2 = run_chain(
        initial_matching=DimerState.all_vertical(L),
        rng=chain2_rng,
        **common_arguments,
    )
    if chain1_rng is chain2_rng:
        raise AssertionError("independent chains unexpectedly share one RNG object")
    if not (
        chain1.final_state.L == environment.L == chain2.final_state.L
        and environment.horizontal_weights.flags.writeable is False
        and environment.vertical_weights.flags.writeable is False
    ):
        raise AssertionError("paired frozen-environment invariant failed")
    return PairedRun(environment, seeds, chain1, chain2, separations)


def _base_row(
    paired: PairedRun,
    model: str,
    gamma_shape: float,
    L: int,
    environment_index: int,
) -> dict[str, object]:
    return {
        "master_seed": paired.seeds.master_seed,
        "environment_seed": paired.seeds.environment_seed,
        "chain1_seed": paired.seeds.chain1_seed,
        "chain2_seed": paired.seeds.chain2_seed,
        "L": L,
        "model": model,
        "gamma_shape": gamma_shape if model == "gamma" else np.nan,
        "environment_index": environment_index,
        "environment_id": f"{model}-L{L}-env{environment_index:06d}",
        "chain1_start": "all_horizontal",
        "chain2_start": "all_vertical",
    }


def run_pilot(
    *,
    sizes: Iterable[int],
    models: Iterable[str],
    gamma_shape: float,
    environments: int,
    burnin_sweeps: int,
    measurement_gap_sweeps: int,
    num_measurements: int,
    master_seed: int,
    output: str | Path,
    relative_separations: Iterable[float] = DEFAULT_RELATIVE_SEPARATIONS,
) -> dict[str, object]:
    """Run a deliberately sequential local pilot and write auditable tables."""
    size_values = tuple(validate_L(int(value)) for value in sizes)
    model_values = tuple(str(value) for value in models)
    if not size_values:
        raise ValueError("at least one size is required")
    if not model_values or any(model not in MODEL_CODES for model in model_values):
        raise ValueError("models must be a nonempty subset of: uniform, gamma")
    if environments <= 0:
        raise ValueError("environments must be positive")
    if num_measurements <= 0:
        raise ValueError("num_measurements must be positive")

    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    measurement_rows: list[dict[str, object]] = []
    spatial_rows: list[dict[str, object]] = []
    diagnostic_rows: list[dict[str, object]] = []
    started = perf_counter()
    for model in model_values:
        for L in size_values:
            for environment_index in range(environments):
                paired = run_paired_environment(
                    L=L,
                    model=model,
                    gamma_shape=gamma_shape,
                    environment_index=environment_index,
                    master_seed=master_seed,
                    burnin_sweeps=burnin_sweeps,
                    measurement_gap_sweeps=measurement_gap_sweeps,
                    num_measurements=num_measurements,
                    relative_separations=relative_separations,
                )
                base = _base_row(paired, model, gamma_shape, L, environment_index)
                for measurement_index, (height1, height2) in enumerate(
                    zip(paired.chain1.center_heights, paired.chain2.center_heights, strict=True)
                ):
                    measurement_rows.append(
                        base
                        | {
                            "measurement_index": measurement_index,
                            "H1_center": height1,
                            "H2_center": height2,
                        }
                    )
                    for separation in paired.separations:
                        spatial_rows.append(
                            base
                            | {
                                "measurement_index": measurement_index,
                                "requested_rho": separation.requested_rho,
                                "r": separation.r,
                                "DeltaH1": paired.chain1.spatial_increment_measurements[
                                    separation.r
                                ][measurement_index],
                                "DeltaH2": paired.chain2.spatial_increment_measurements[
                                    separation.r
                                ][measurement_index],
                            }
                        )
                diagnostic_rows.append(
                    base
                    | {
                        **{f"chain1_{key}": value for key, value in paired.chain1.diagnostics().items()},
                        **{f"chain2_{key}": value for key, value in paired.chain2.diagnostics().items()},
                    }
                )
    elapsed = perf_counter() - started
    pd.DataFrame(measurement_rows).to_csv(output_path / "center_measurements.csv", index=False)
    pd.DataFrame(spatial_rows).to_csv(output_path / "spatial_measurements.csv", index=False)
    pd.DataFrame(diagnostic_rows).to_csv(output_path / "chain_diagnostics.csv", index=False)
    metadata: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "exploratory": True,
        "mixing_warning": "Burn-in and spacing are user inputs; output does not certify convergence.",
        "statistical_unit": "environment_id; repeated measurements within an environment are not independent environments",
        "seed_policy": "SeedSequence([master_seed,L,model_code,environment_index]).spawn(3), then stored uint64 PCG64 seeds",
        "master_seed": master_seed,
        "sizes": list(size_values),
        "models": list(model_values),
        "gamma_shape": gamma_shape,
        "environments_per_model_size": environments,
        "burnin_sweeps": burnin_sweeps,
        "measurement_gap_sweeps": measurement_gap_sweeps,
        "num_measurements_per_chain_environment": num_measurements,
        "relative_separations": [float(value) for value in relative_separations],
        "elapsed_seconds": elapsed,
    }
    (output_path / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata
