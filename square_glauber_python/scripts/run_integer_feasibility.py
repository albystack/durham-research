"""Fixed-cohort direct sampling; original arrays required for authenticated replay.

The environment seeds are retained; sampler seeds have a distinct versioned
namespace. These samples do not complete or overwrite old CFTP certificates.
"""

import argparse
import shutil
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

from square_glauber.comparison import (
    ComparisonModel,
    comparison_environment,
    comparison_seeds,
    measure_probes,
    symmetric_probes,
)
from square_glauber.height import face_heights
from square_glauber.integer_sampling import (
    IntegerFrontierSampler,
    VERSION,
    load_verified_environment,
)
from square_glauber.comparison_pilot import law_id

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "square_glauber_python/outputs/integer_feasibility_20261001"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--original-cases", type=Path)
    source.add_argument(
        "--allow-seed-regeneration",
        action="store_true",
        help="Diagnostic only; byte identity must still be audited",
    )
    parser.add_argument(
        "--pilot-status",
        type=Path,
        default=ROOT
        / "square_glauber_python/outputs/hamilton_matched_20260930/feasibility/stream_status.csv",
    )
    args = parser.parse_args()
    output = args.output
    old = pd.read_csv(
        args.pilot_status, dtype={"environment_seed": str, "cftp_seed": str}
    )
    config = json.loads(
        (ROOT / "docs/validation/numerical_20261001/plan.json").read_text()
    )["square"]
    if output.exists():
        raise ValueError("output exists; preserve it and choose a new directory")
    output.mkdir(parents=True)
    statuses = []
    measurements = []
    for name in ("aztec_gamma", "aztec_mean", "uniform"):
        model = ComparisonModel(name=name)
        cohort = (
            config["random_cohort"]
            if name == "aztec_gamma"
            else config["control_cohort"]
        )
        for size, count in cohort.items():
            L = int(size)
            for index in range(count):
                original_seeds = comparison_seeds(
                    config["original_master_seed"], L, model, index
                )
                record = old[
                    (old.model == name)
                    & (old["size"] == L)
                    & old.environment_seed.eq(str(original_seeds[0]))
                ]
                expected_hash = (
                    record.environment_sha256.iloc[0] if len(record) else None
                )
                original_path = None
                if args.original_cases is not None and L >= 8:
                    model_index = ("aztec_gamma", "aztec_mean", "uniform").index(name)
                    original_path = (
                        args.original_cases
                        / f"model{model_index}-L{L}-env{index:06d}"
                        / "weights.npz"
                    )
                    if expected_hash is None:
                        raise ValueError(
                            "original environment missing from frozen pilot status"
                        )
                    env = load_verified_environment(
                        original_path, L, name, expected_hash
                    )
                else:
                    env = comparison_environment(
                        L, np.random.default_rng(original_seeds[0]), model
                    )
                case = output / "cases" / f"{name}-L{L}-env{index:06d}"
                case.mkdir(parents=True)
                np.savez_compressed(
                    case / "weights.npz",
                    horizontal=env.horizontal_weights,
                    vertical=env.vertical_weights,
                )
                if original_path is not None:
                    shutil.copyfile(original_path, case / "weights.npz")
                archive_hash = hashlib.sha256(
                    (case / "weights.npz").read_bytes()
                ).hexdigest()
                identity_verified = (
                    expected_hash is not None and archive_hash == expected_hash
                )
                start = perf_counter()
                sampler = IntegerFrontierSampler(
                    env, max_states=config["max_states_per_environment"]
                )
                setup = perf_counter() - start
                seeds = [
                    int.from_bytes(
                        hashlib.sha256(f"{VERSION}:{s}".encode()).digest()[:16], "big"
                    )
                    for s in original_seeds[1:]
                ]
                assert len(set(seeds)) == 2
                values = []
                sampling_start = perf_counter()
                for replica, seed in enumerate(seeds, 1):
                    state = sampler.sample(seed)
                    np.savez_compressed(
                        case / f"sample{replica}.npz",
                        horizontal=state.horizontal_occupied,
                        vertical=state.vertical_occupied,
                    )
                    values.append(
                        measure_probes(face_heights(state), symmetric_probes(L))
                    )
                meta = dict(
                    model=name,
                    L=L,
                    index=index,
                    environment_seed=str(original_seeds[0]),
                    sampler=VERSION,
                    original_archive_verified=identity_verified,
                    original_archive_expected_sha256=expected_hash,
                    input_provenance="loaded_verified_archive"
                    if original_path
                    else "seed_regenerated_local",
                    original_chain_seeds=[str(s) for s in original_seeds[1:]],
                    sampler_seeds=[str(s) for s in seeds],
                    partition_bits=sampler.total.bit_length(),
                    cached_states=sampler.partition.cache_info().currsize,
                    partition_sha256=hashlib.sha256(
                        sampler.total.to_bytes(
                            (sampler.total.bit_length() + 7) // 8, "big"
                        )
                    ).hexdigest(),
                    setup_seconds=setup,
                    pair_seconds=perf_counter() - sampling_start,
                    complete_pair=True,
                    exactness="integer Gibbs law conditional on stored dyadic weights and ideal random bits; not CFTP",
                    files_sha256={
                        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in case.glob("*.npz")
                    },
                )
                (case / "certificate.json").write_text(
                    json.dumps(meta, indent=2) + "\n"
                )
                statuses.append(meta)
                for probe in symmetric_probes(L):
                    measurements.append(
                        dict(
                            domain="square",
                            size=L,
                            model=name,
                            law_id=law_id(model),
                            environment_id=case.name,
                            probe=probe.label,
                            axis=probe.axis,
                            r=probe.r,
                            distance=probe.distance,
                            height_scale=1,
                            value1=values[0][probe.label],
                            value2=values[1][probe.label],
                        )
                    )
                print(
                    name,
                    L,
                    index,
                    round(setup, 4),
                    sampler.partition.cache_info().currsize,
                    flush=True,
                )
    pd.DataFrame(statuses).to_csv(output / "status.csv", index=False)
    pd.DataFrame(measurements).to_csv(output / "measurements.csv", index=False)
    summary = dict(
        version=VERSION,
        complete_pairs=len(statuses),
        all_pairs_certified=True,
        measurements_sha256=hashlib.sha256(
            (output / "measurements.csv").read_bytes()
        ).hexdigest(),
        note="Bounded fixed cohort; no scaling fit; new streams do not resume old CFTP",
        groups=[],
    )
    for (model, L), group in pd.DataFrame(statuses).groupby(["model", "L"]):
        summary["groups"].append(
            dict(
                model=model,
                L=int(L),
                pairs=len(group),
                setup_seconds=float(group.setup_seconds.sum()),
                sampling_seconds=float(group.pair_seconds.sum()),
                max_cached_states=int(group.cached_states.max()),
            )
        )
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
