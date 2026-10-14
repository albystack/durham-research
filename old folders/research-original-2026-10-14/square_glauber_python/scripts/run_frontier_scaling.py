"""Run one frozen feasibility case in a fresh, resource-limited HPC worker.

This writes all outcomes, including explicit budget exhaustion. It never
replaces a difficult environment or turns an incomplete partition into a draw.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
from time import perf_counter

import numpy as np

from square_glauber.budgeted_frontier import (
    BudgetedFrontierSampler,
    FrontierBudgetExceeded,
    VERSION,
    peak_rss_bytes,
)
from square_glauber.comparison import (
    ComparisonModel,
    comparison_environment,
    measure_probes,
    symmetric_probes,
)
from square_glauber.height import face_heights
from square_glauber.integer_sampling import load_verified_environment


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def run_case(plan_path, index, output, original_cases=None):
    plan = json.loads(plan_path.read_text())
    if plan["sampler_version"] != VERSION:
        raise ValueError("wrong sampler plan")
    case = plan["cases"][index]
    if case["case_index"] != index:
        raise ValueError("case index mismatch")
    directory = output / f"case-{index:03d}"
    directory.mkdir(parents=True, exist_ok=False)
    start = perf_counter()
    status = dict(
        case=case,
        plan_sha256=digest(plan_path),
        sampler_version=VERSION,
        state="STARTED",
        complete_pair=False,
        python=platform.python_version(),
        numpy=np.__version__,
        budgets=plan["budgets"],
    )
    write_json(directory / "status.json", status)
    try:
        if case["input_source"] == "original_archive":
            if original_cases is None:
                raise ValueError("original archive root required")
            path = original_cases / case["original_relative_path"]
            environment = load_verified_environment(
                path, case["L"], case["model"], case["original_sha256"]
            )
            shutil.copyfile(path, directory / "weights.npz")
        else:
            environment = comparison_environment(
                case["L"],
                np.random.default_rng(int(case["environment_seed"])),
                ComparisonModel(name=case["model"]),
            )
            np.savez_compressed(
                directory / "weights.npz",
                horizontal=environment.horizontal_weights,
                vertical=environment.vertical_weights,
            )
        status["weights_sha256"] = digest(directory / "weights.npz")
        write_json(directory / "status.json", status)
        sampler = BudgetedFrontierSampler(environment, **plan["budgets"])
        status.update(sampler.stats())
        status["setup_seconds"] = sampler.setup_seconds
        states = [sampler.sample(int(seed)) for seed in case["replica_seeds"]]
        probes = symmetric_probes(
            case["L"], fractions=(0.125, 0.25), axes=("diagonal",)
        )
        values = [measure_probes(face_heights(s), probes) for s in states]
        for replica, state in enumerate(states, 1):
            state.validate()
            np.savez_compressed(
                directory / f"sample{replica}.npz",
                horizontal=state.horizontal_occupied,
                vertical=state.vertical_occupied,
            )
        measurements = [
            dict(
                probe=p.label,
                r=p.r,
                distance=p.distance,
                start=p.start,
                end=p.end,
                height_scale=1,
                value1=values[0][p.label],
                value2=values[1][p.label],
            )
            for p in probes
        ]
        write_json(directory / "measurements.json", measurements)
        status.update(
            state="COMPLETE",
            complete_pair=True,
            partition_bits=sampler.total.bit_length(),
            partition_sha256=hashlib.sha256(
                sampler.total.to_bytes((sampler.total.bit_length() + 7) // 8, "big")
            ).hexdigest(),
            files_sha256={
                p.name: digest(p)
                for p in directory.iterdir()
                if p.name != "status.json"
            },
            exactness="Gibbs for stored dyadic weights, conditional on ideal random bits; not CFTP",
        )
    except FrontierBudgetExceeded as error:
        status.update(error.stats)
        status.update(state="BUDGET_EXHAUSTED", reason=error.reason)
    except Exception as error:
        status.update(state="ERROR", reason=f"{type(error).__name__}: {error}")
        raise
    finally:
        status["total_seconds"] = perf_counter() - start
        status["peak_rss_bytes"] = peak_rss_bytes()
        write_json(directory / "status.json", status)
    print(json.dumps(status), flush=True)
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--case", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--original-cases", type=Path)
    args = parser.parse_args()
    run_case(args.plan, args.case, args.output, args.original_cases)
