"""Optional one-sample workers for already-frozen campaigns.

Run only after all legacy campaign tasks have exited: their writer does not
participate in these advisory locks. Workers never create/redraw environments,
write pair.json, or update task reports. A separate coordinator owns pair.json.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict
import fcntl
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .experiment import derive_seeds
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, spatial_separations
from .mixing_extension import environment_hashes
from .perfect import CFTPNotCoalesced, MonotoneCFTPSampler
from .perfect_campaign import (
    _atomic_json, _atomic_npz, _case_path, _load_certified_sample,
    _load_environment, _pair_metadata, _sample_metadata, file_sha256,
)


@contextmanager
def exclusive_lock(path: Path):
    """Nonblocking POSIX advisory lock; never unlink a live lock inode."""
    with path.open("a+") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"another restart worker owns {path}") from error
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def frozen_case(manifest, task_id, output, environment_index):
    """Load immutable existing inputs, checking manifest, seeds and both hashes."""
    manifest = Path(manifest)
    manifest_metadata = json.loads(manifest.with_suffix(".json").read_text())
    if file_sha256(manifest) != manifest_metadata["manifest_sha256"]:
        raise RuntimeError("frozen manifest hash mismatch")
    frame = pd.read_csv(manifest)
    selected = frame.loc[frame["task_id"] == int(task_id)]
    if len(selected) != 1:
        raise ValueError("task_id must identify exactly one frozen manifest row")
    row = selected.iloc[0]
    if not int(row.environment_start) <= environment_index < int(row.environment_stop):
        raise ValueError("environment index is outside the selected task")
    L, model = int(row.L), str(row.model)
    case = _case_path(Path(output), model, L, environment_index)
    metadata = json.loads((case / "environment.json").read_text())
    seeds = derive_seeds(int(row.master_seed), L, model, environment_index)
    expected = {
        **asdict(seeds), "campaign_id": str(row.campaign_id), "L": L,
        "model": model, "environment_index": environment_index,
        "environment_id": f"{model}-L{L}-env{environment_index:06d}",
        "gamma_shape": float(row.gamma_shape) if model == "gamma" else None,
    }
    if any(metadata.get(key) != value for key, value in expected.items()):
        raise RuntimeError("frozen environment seed/identity metadata mismatch")
    weights = case / "environment_weights.npz"
    if file_sha256(weights) != metadata["weights_npz_sha256"]:
        raise RuntimeError("frozen environment NPZ hash mismatch")
    environment = _load_environment(
        weights, L=L, model=model,
        gamma_shape=float(row.gamma_shape) if model == "gamma" else 1.0,
    )
    if any(metadata.get(key) != value for key, value in environment_hashes(environment).items()):
        raise RuntimeError("frozen environment array hash mismatch")
    return case, environment, metadata, seeds, str(row.backend)


def run_sample_worker(manifest, task_id, output, environment_index, sample_index,
                      *, max_horizon_sweeps, chunk_attempts=4_000_000,
                      legacy_censor_path=None):
    if sample_index not in (1, 2):
        raise ValueError("sample_index must be 1 or 2")
    if max_horizon_sweeps <= 0 or chunk_attempts <= 0:
        raise ValueError("horizon and chunk_attempts must be positive")
    case, environment, metadata, seeds, backend = frozen_case(
        manifest, task_id, output, environment_index)
    separations = spatial_separations(environment.L, DEFAULT_RELATIVE_SEPARATIONS)
    seed = seeds.chain1_seed if sample_index == 1 else seeds.chain2_seed
    with exclusive_lock(case / f".sample{sample_index}.restart.lock"):
        existing = _load_certified_sample(
            case, sample_index=sample_index, seed=seed, L=environment.L,
            environment_sha256=metadata["environment_sha256"], separations=separations)
        if existing is not None:
            return existing[1]
        if (case / "pair.json").exists():
            raise RuntimeError("pair exists but selected certified sample is absent")
        # Lazy import keeps coordinator/inspection usable independently of restart.
        from .restart import sample_with_checkpoint
        sampler = MonotoneCFTPSampler(environment, backend=backend,
                                      chunk_attempts=chunk_attempts)
        try:
            result = sample_with_checkpoint(
                sampler, seed, checkpoint_path=case / f"sample{sample_index}_restart.json",
                max_horizon_sweeps=max_horizon_sweeps,
                legacy_censor_path=legacy_censor_path)
        except CFTPNotCoalesced as error:
            censor = {
                "sample_index": sample_index, "cftp_seed": int(seed),
                "certified_perfect_sample": False, "censored": True,
                "environment_sha256": metadata["environment_sha256"],
                "max_horizon_sweeps": max_horizon_sweeps,
                "horizon_sweeps": error.horizon_sweeps,
                "iterations": [asdict(item) for item in error.iterations],
                "instruction": "extend identical frozen environment and seed; never discard",
                "restart_info": getattr(error, "restart_info", {}),
            }
            target = case / "censoring" / f"sample{sample_index}_restart_cap-{max_horizon_sweeps}.json"
            if not target.exists():
                _atomic_json(target, censor)
            return censor
        value = _sample_metadata(result, environment_sha256=metadata["environment_sha256"],
                                 sample_index=sample_index, separations=separations)
        value["restart_info"] = getattr(result, "restart_info", {})
        value["update_counts_and_rates_scope"] = value["restart_info"].get(
            "update_count_scope", "all_replayed_updates")
        state_path = case / f"sample{sample_index}_state.npz"
        arrays = dict(horizontal_occupied=result.state.horizontal_occupied,
                      vertical_occupied=result.state.vertical_occupied)
        if state_path.exists():
            with np.load(state_path, allow_pickle=False) as saved:
                if any(not np.array_equal(saved[key], array) for key, array in arrays.items()):
                    raise RuntimeError("orphan saved sample state differs; refusing overwrite")
        else:
            _atomic_npz(state_path, **arrays)
        _atomic_json(case / f"sample{sample_index}.json", value)
        return value


def consolidate_case(manifest, task_id, output, environment_index):
    case, environment, metadata, seeds, _ = frozen_case(
        manifest, task_id, output, environment_index)
    separations = spatial_separations(environment.L, DEFAULT_RELATIVE_SEPARATIONS)
    with exclusive_lock(case / ".pair.restart.lock"), \
         exclusive_lock(case / ".sample1.restart.lock"), \
         exclusive_lock(case / ".sample2.restart.lock"):
        samples = [_load_certified_sample(
            case, sample_index=index, seed=seed, L=environment.L,
            environment_sha256=metadata["environment_sha256"], separations=separations)
            for index, seed in ((1, seeds.chain1_seed), (2, seeds.chain2_seed))]
        if any(sample is None for sample in samples):
            if (case / "pair.json").exists():
                raise RuntimeError("pair exists without two verified samples")
            return {"complete": False, "environment_id": metadata["environment_id"]}
        pair = _pair_metadata(metadata, samples[0][1], samples[1][1], separations=separations)
        target = case / "pair.json"
        if target.exists():
            saved = json.loads(target.read_text())
            if any(saved.get(key) != value for key, value in pair.items() if key != "created_utc"):
                raise RuntimeError("existing pair conflicts with verified samples")
        else:
            _atomic_json(target, pair)
        return {"complete": True, "environment_id": metadata["environment_id"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("sample", "consolidate"))
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--task-id", required=True, type=int)
    parser.add_argument("--environment-index", required=True, type=int)
    parser.add_argument("--sample-index", choices=(1, 2), type=int)
    parser.add_argument("--max-horizon-sweeps", type=int)
    parser.add_argument("--chunk-attempts", type=int, default=4_000_000)
    parser.add_argument("--legacy-censor-path")
    args = parser.parse_args()
    common = (args.manifest, args.task_id, args.output, args.environment_index)
    if args.action == "sample":
        if args.sample_index is None or args.max_horizon_sweeps is None:
            parser.error("sample requires --sample-index and --max-horizon-sweeps")
        result = run_sample_worker(*common, args.sample_index,
            max_horizon_sweeps=args.max_horizon_sweeps,
            chunk_attempts=args.chunk_attempts, legacy_censor_path=args.legacy_censor_path)
    else:
        result = consolidate_case(*common)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
