"""Run one predeclared model/size/domain cell; never discard difficult inputs.

A fresh output directory is mandatory. Every square input and both matchings
are saved before a completion certificate. Failures leave explicit status and
prevent analysis of a selected subset. Julia supplies the independent Aztec
implementation, with the same physical law and probe conventions.
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
from square_glauber.budgeted_frontier import peak_rss_bytes
from square_glauber.comparison import (
    ComparisonModel,
    comparison_environment,
    measure_probes,
    symmetric_probes,
)
from square_glauber.embedded_shuffle import VERSION, EmbeddedSquareSampler
from square_glauber.height import face_heights

NAMESPACE = "matched-shuffle-20261004"
FIELDS = [
    "domain",
    "size",
    "model",
    "law_id",
    "environment_id",
    "probe",
    "axis",
    "r",
    "distance",
    "height_scale",
    "value1",
    "value2",
    "environment_seed",
    "replica1_seed",
    "replica2_seed",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seed(master, domain, model, size, index, role):
    text = f"{NAMESPACE}:{master}:{domain}:{model}:{size}:{index}:{role}"
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "big")


def write_json(path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)


def run_square(plan_path, cell, directory, validation=False):
    plan = json.loads(plan_path.read_text())
    count = plan["validation_pairs"] if validation else cell["pairs"]
    master = plan["validation_seed"] if validation else plan["master_seed"]
    directory.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    status = {
        "state": "STARTED",
        "completed": 0,
        "cell": cell,
        "pairs": count,
        "validation_only": validation,
        "plan_sha256": digest(plan_path),
        "sampler_version": VERSION,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "master_seed": master,
    }
    write_json(directory / "manifest.json", status)
    size, model = cell["size"], cell["model"]
    probes = symmetric_probes(size, fractions=(0.125, 0.25), axes=("diagonal",))
    law = f"{model}:alpha=0.2:beta=0.25:scale=1:a-north,b-east"
    temporary = directory / "measurements.csv.tmp"
    try:
        with (
            temporary.open("w", newline="") as out,
            (directory / "certificates.jsonl").open("w") as certificates,
        ):
            writer = csv.DictWriter(out, fieldnames=FIELDS)
            writer.writeheader()
            for index in range(1, count + 1):
                seeds = [
                    seed(master, "square", model, size, index, role)
                    for role in ("environment", "replica1", "replica2")
                ]
                env = comparison_environment(
                    size, np.random.default_rng(seeds[0]), ComparisonModel(name=model)
                )
                case = directory / "cases" / f"env-{index:05d}"
                case.mkdir(parents=True)
                np.savez_compressed(
                    case / "weights.npz",
                    horizontal=env.horizontal_weights,
                    vertical=env.vertical_weights,
                )
                # Reduction errors and zero-cell failures abort this cell. No
                # replacement environment, epsilon, burn-in or censoring rule.
                sampler = EmbeddedSquareSampler(env, cast=float)
                states = [sampler.sample(s) for s in seeds[1:]]
                values = [measure_probes(face_heights(s), probes) for s in states]
                for replica, state in enumerate(states, 1):
                    state.validate()
                    np.savez_compressed(
                        case / f"sample{replica}.npz",
                        horizontal=state.horizontal_occupied,
                        vertical=state.vertical_occupied,
                    )
                rows = []
                for probe in probes:
                    rows.append(
                        {
                            "domain": "square",
                            "size": size,
                            "model": model,
                            "law_id": law,
                            "environment_id": f"square-{model}-n{size}-env{index}",
                            "probe": probe.label,
                            "axis": "diagonal",
                            "r": probe.r,
                            "distance": probe.distance,
                            "height_scale": 1,
                            "value1": values[0][probe.label],
                            "value2": values[1][probe.label],
                            "environment_seed": str(seeds[0]),
                            "replica1_seed": str(seeds[1]),
                            "replica2_seed": str(seeds[2]),
                        }
                    )
                write_json(case / "rows.json", rows)
                certificate = {
                    "index": index,
                    "files": {
                        name: digest(case / name)
                        for name in (
                            "weights.npz",
                            "sample1.npz",
                            "sample2.npz",
                            "rows.json",
                        )
                    },
                }
                certificates.write(json.dumps(certificate) + "\n")
                writer.writerows(rows)
                status["completed"] = index
                if index % 64 == 0:
                    out.flush()
                    certificates.flush()
                    write_json(directory / "status.json", status)
                    print(f"square {model} L={size}: {index}/{count}", flush=True)
                del sampler, states
        temporary.replace(directory / "measurements.csv")
        status.update(
            state="COMPLETE",
            measurements_sha256=digest(directory / "measurements.csv"),
            certificates_sha256=digest(directory / "certificates.jsonl"),
            arithmetic="Float64 reductions; exact rational creation coins for represented weights; finite numerical audit, not universal exactness",
        )
    except Exception as error:
        status.update(state="ERROR", reason=f"{type(error).__name__}: {error}")
        raise
    finally:
        status.update(
            seconds=time.perf_counter() - started, peak_rss_bytes=peak_rss_bytes()
        )
        write_json(directory / "status.json", status)
    print(json.dumps(status), flush=True)
    return status


def run_cell(plan_path, cell_index, output, validation=False):
    plan = json.loads(plan_path.read_text())
    if plan["namespace"] != NAMESPACE:
        raise ValueError("unexpected campaign namespace")
    cell = plan["cells"][cell_index]
    if cell["index"] != cell_index:
        raise ValueError("cell index mismatch")
    destination = output / f"cell-{cell_index:02d}"
    if cell["domain"] == "square":
        return run_square(plan_path, cell, destination, validation)
    root = Path(__file__).resolve().parents[2]
    command = [
        os.environ.get("RESEARCH_JULIA", "julia"),
        "--startup-file=no",
        f"--project={os.environ.get('RESEARCH_JULIA_PROJECT', root / 'aztec')}",
        str(root / "aztec/scripts/run_matched_shuffle.jl"),
        "--output",
        str(destination),
        "--size",
        str(cell["size"]),
        "--model",
        cell["model"],
        "--environments",
        str(plan["validation_pairs"] if validation else cell["pairs"]),
        "--seed",
        str(plan["validation_seed"] if validation else plan["master_seed"]),
        "--plan-sha256",
        digest(plan_path),
        "--validation",
        str(validation).lower(),
    ]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--cell", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--validation-only", action="store_true")
    args = parser.parse_args()
    run_cell(args.plan, args.cell, args.output, args.validation_only)
