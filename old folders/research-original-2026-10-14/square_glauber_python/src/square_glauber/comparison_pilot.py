"""Small, resumable feasibility pilot; never a production-campaign launcher.

Usage: python -m square_glauber.comparison_pilot --config CONFIG --output DIR
The cap can be raised on restart, but identities, weights and random maps cannot
change. Every requested stream has a status; only complete pairs have observables.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import platform
import subprocess
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from .comparison import (
    ComparisonModel,
    comparison_environment,
    comparison_seeds,
    measure_probes,
    symmetric_probes,
)
from .height import face_heights
from .matching import DimerState
from .perfect import CFTPNotCoalesced, MonotoneCFTPSampler
from .perfect_campaign import _atomic_json, _atomic_npz
from .restart import sample_with_checkpoint


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def law_id(model: ComparisonModel) -> str:
    # Round-trip precision prevents distinct nearby parameters sharing a label.
    alpha, beta, scale = (
        repr(float(x)).removesuffix(".0")
        for x in (model.alpha, model.beta, model.scale)
    )
    return f"{model.name}:alpha={alpha}:beta={beta}:scale={scale}:a-north,b-east"


def validate_config(config: dict) -> tuple[list[ComparisonModel], list[int]]:
    """Hard limits make this entry point unsuitable for accidental large runs."""
    sizes = config["sizes"]
    count = config["environments_per_size"]
    cap = config["max_horizon_sweeps"]
    if (
        not sizes
        or len(sizes) != len(set(sizes))
        or any(type(L) is not int or L < 4 or L > 32 or L % 2 for L in sizes)
    ):
        raise ValueError("pilot sizes must be distinct even integers in [4,32]")
    if type(count) is not int or not 1 <= count <= 64:
        raise ValueError("pilot requires 1..64 environments per size/model")
    if type(cap) is not int or cap <= 0 or cap & (cap - 1) or cap > 2**20:
        raise ValueError("pilot cap must be a power of two, at most 2**20 sweeps")
    models = [ComparisonModel(**m) for m in config["models"]]
    if not models or len({law_id(m) for m in models}) != len(models):
        raise ValueError("models must be nonempty and distinct")
    if len(models) * len(sizes) * count > 1024:
        raise ValueError("pilot exceeds 1024 environments")
    return models, sizes


def _save_csv(path: Path, rows: list[dict]) -> None:
    temporary = path.with_suffix(".tmp")
    pd.DataFrame(rows).to_csv(temporary, index=False)
    temporary.replace(path)


def run_pilot(config: dict, output: Path, *, backend="numba") -> dict:
    """Single-writer, restartable pilot. Censored environments are retained."""
    models, sizes = validate_config(config)
    output.mkdir(parents=True, exist_ok=True)
    # A local/Slurm task must own this output root. Advisory locking also works
    # for crash recovery: the OS releases the lock when the process exits.
    with (output / ".writer.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        identity = {k: v for k, v in config.items() if k != "max_horizon_sweeps"}
        identity["schema"] = 1
        manifest_path = output / "manifest.json"
        if manifest_path.exists() and json.loads(manifest_path.read_text()) != identity:
            raise ValueError("pilot identity changed; use a fresh output directory")
        _atomic_json(manifest_path, identity)
        try:
            revision = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip()
        except (OSError, subprocess.CalledProcessError):
            revision = "unavailable"
        _atomic_json(
            output / "run_metadata.json",
            {
                "config": config,
                "backend": backend,
                "python": platform.python_version(),
                "numpy": np.__version__,
                "git_revision": revision,
                "source_sha256": {
                    p.name: file_hash(p) for p in Path(__file__).parent.glob("*.py")
                },
                "purpose": "feasibility only; censored cases prohibit complete-campaign scaling claims",
            },
        )
        statuses, measurements = [], []
        used_seeds = set()
        for model_index, model in enumerate(models):
            for L in sizes:
                probes = symmetric_probes(L)
                _atomic_json(
                    output / f"probes-L{L}.json",
                    [
                        {**asdict(p), "distance": p.distance, "label": p.label}
                        for p in probes
                    ],
                )
                for index in range(config["environments_per_size"]):
                    seeds = comparison_seeds(config["master_seed"], L, model, index)
                    if used_seeds.intersection(seeds):
                        raise RuntimeError(
                            "seed collision across requested pilot streams"
                        )
                    used_seeds.update(seeds)
                    case = output / "cases" / f"model{model_index}-L{L}-env{index:06d}"
                    case.mkdir(parents=True, exist_ok=True)
                    environment = comparison_environment(
                        L, np.random.default_rng(seeds[0]), model
                    )
                    weights_path = case / "weights.npz"
                    if weights_path.exists():
                        with np.load(weights_path, allow_pickle=False) as saved:
                            for key, actual in (
                                ("horizontal", environment.horizontal_weights),
                                ("vertical", environment.vertical_weights),
                            ):
                                if not np.array_equal(saved[key], actual):
                                    raise ValueError(
                                        "saved weights differ from seeded regeneration"
                                    )
                    else:
                        _atomic_npz(
                            weights_path,
                            horizontal=environment.horizontal_weights,
                            vertical=environment.vertical_weights,
                        )
                    env_hash = file_hash(weights_path)
                    sampler = MonotoneCFTPSampler(environment, backend=backend)
                    paired = []
                    for replica, seed in enumerate(seeds[1:], 1):
                        state_path = case / f"sample{replica}_state.npz"
                        metadata_path = case / f"sample{replica}.json"
                        base = dict(
                            domain="square",
                            size=L,
                            model=model.name,
                            law_id=law_id(model),
                            environment_id=case.name,
                            environment_seed=str(seeds[0]),
                            replica=replica,
                            cftp_seed=str(seed),
                            environment_sha256=env_hash,
                        )
                        if metadata_path.exists():
                            meta = json.loads(metadata_path.read_text())
                            if (
                                meta["cftp_seed"] != str(seed)
                                or meta["environment_sha256"] != env_hash
                                or not meta["certified"]
                                or meta["state_sha256"] != file_hash(state_path)
                            ):
                                raise ValueError(
                                    "saved certificate identity/hash mismatch"
                                )
                            with np.load(state_path, allow_pickle=False) as saved:
                                state = DimerState(
                                    L,
                                    saved["horizontal_occupied"],
                                    saved["vertical_occupied"],
                                )
                            state.validate()
                        else:
                            try:
                                sample = sample_with_checkpoint(
                                    sampler,
                                    seed,
                                    checkpoint_path=case
                                    / f"sample{replica}_checkpoint.json",
                                    max_horizon_sweeps=config["max_horizon_sweeps"],
                                )
                            except CFTPNotCoalesced as error:
                                checkpoint = json.loads(
                                    (
                                        case / f"sample{replica}_checkpoint.json"
                                    ).read_text()
                                )
                                status = {
                                    **base,
                                    "certified": False,
                                    "censored": True,
                                    "horizon_sweeps": error.horizon_sweeps,
                                    "elapsed_seconds": checkpoint["payload"][
                                        "elapsed_seconds"
                                    ],
                                }
                                _atomic_json(
                                    case / f"sample{replica}_status.json", status
                                )
                                statuses.append(status)
                                continue
                            state = sample.state
                            _atomic_npz(
                                state_path,
                                horizontal_occupied=state.horizontal_occupied,
                                vertical_occupied=state.vertical_occupied,
                            )
                            meta = {
                                **base,
                                **sample.metadata(),
                                "cftp_seed": str(seed),
                                "certified": True,
                                "censored": False,
                                "horizon_sweeps": sample.horizon_sweeps,
                                "elapsed_seconds": sample.elapsed_seconds,
                                "state_sha256": file_hash(state_path),
                            }
                            _atomic_json(metadata_path, meta)
                        statuses.append(meta)
                        _atomic_json(case / f"sample{replica}_status.json", meta)
                        heights = face_heights(state)
                        paired.append(
                            (replica, heights, measure_probes(heights, probes))
                        )
                    if len(paired) == 2:
                        for probe in probes:
                            measurements.append(
                                {
                                    "domain": "square",
                                    "size": L,
                                    "model": model.name,
                                    "law_id": law_id(model),
                                    "environment_id": case.name,
                                    "probe": probe.label,
                                    "axis": probe.axis,
                                    "r": probe.r,
                                    "distance": probe.distance,
                                    "height_scale": 1,
                                    "start_row": probe.start[0],
                                    "start_column": probe.start[1],
                                    "end_row": probe.end[0],
                                    "end_column": probe.end[1],
                                    "value1": paired[0][2][probe.label],
                                    "value2": paired[1][2][probe.label],
                                }
                            )
                    # Publish progress even if later environments get interrupted.
                    _save_csv(output / "stream_status.csv", statuses)
                    _save_csv(output / "measurements.csv", measurements)
        frame = pd.DataFrame(statuses)
        summary_rows = []
        for (model, L), group in frame.groupby(["law_id", "size"]):
            done = group[group["certified"]]
            row = dict(
                law_id=model,
                size=int(L),
                requested_streams=len(group),
                certified_streams=len(done),
                censored_streams=int(group["censored"].sum()),
                complete_pairs=int(
                    group.groupby("environment_id")["certified"].sum().eq(2).sum()
                ),
                total_recorded_seconds=float(group["elapsed_seconds"].sum()),
            )
            for quantity in ("elapsed_seconds", "horizon_sweeps"):
                for q in (0.5, 0.9, 0.99):
                    row[f"completed_{quantity}_q{q:g}"] = (
                        float(done[quantity].quantile(q)) if len(done) else None
                    )
            summary_rows.append(row)
        report = {
            "all_pairs_certified": bool(frame["certified"].all()),
            "requested_environments": len(frame) // 2,
            "measurements_sha256": file_hash(output / "measurements.csv"),
            "note": "Quantiles condition on certification; censored horizons are lower bounds. No scaling fit.",
            "groups": summary_rows,
        }
        _atomic_json(output / "summary.json", report)
        _save_csv(output / "summary.csv", summary_rows)
        return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--backend", choices=("reference", "numba"), default="numba")
    args = parser.parse_args()
    report = run_pilot(
        json.loads(args.config.read_text()), args.output, backend=args.backend
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
