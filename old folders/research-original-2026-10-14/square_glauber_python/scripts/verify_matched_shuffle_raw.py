"""Revalidate archived weights, matchings, observables and file certificates.

Run on Hamilton where all raw files are retained. The compact-table analysis
is a separate gate; it never claims to have read raw files it did not receive.
"""

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
from analyse_matched_shuffle import load_complete
from run_matched_shuffle import digest, write_json
from square_glauber.comparison import (
    ComparisonModel,
    comparison_environment,
    family_masks,
    measure_probes,
    symmetric_probes,
)
from square_glauber.environment import DimerEnvironment
from square_glauber.height import face_heights
from square_glauber.matching import DimerState


def verify(root, plan_path, output, validation=False):
    frame, hashes = load_complete(root, plan_path, validation)
    plan = json.loads(plan_path.read_text())
    files, pairs, observations = 0, 0, 0
    for cell in plan["cells"]:
        directory = root / f"cell-{cell['index']:02d}"
        n = cell["size"]
        count = plan["validation_pairs"] if validation else cell["pairs"]
        if cell["domain"] == "square":
            certificate = directory / "certificates.jsonl"
            status = json.loads((directory / "status.json").read_text())
            if digest(certificate) != status["certificates_sha256"]:
                raise ValueError("certificate checksum mismatch")
            records = [
                json.loads(line) for line in certificate.read_text().splitlines()
            ]
            if [r["index"] for r in records] != list(range(1, count + 1)):
                raise ValueError("certificate coverage mismatch")
            for record in records:
                case = directory / "cases" / f"env-{record['index']:05d}"
                if set(record["files"]) != {
                    "weights.npz",
                    "sample1.npz",
                    "sample2.npz",
                    "rows.json",
                }:
                    raise ValueError("wrong certified files")
                for name, expected in record["files"].items():
                    if digest(case / name) != expected:
                        raise ValueError("raw checksum mismatch")
                    files += 1
                with np.load(case / "weights.npz", allow_pickle=False) as w:
                    env = DimerEnvironment(n, w["horizontal"], w["vertical"])
                bmask, amask = family_masks(n)
                if (
                    not (env.horizontal_weights[~bmask] == 1).all()
                    or not (env.vertical_weights[~amask] == 1).all()
                ):
                    raise ValueError("wrong fixed edge families")
                if cell["model"] != "aztec_gamma":
                    control = comparison_environment(
                        n, np.random.default_rng(0), ComparisonModel(name=cell["model"])
                    )
                    if not np.array_equal(
                        env.horizontal_weights, control.horizontal_weights
                    ) or not np.array_equal(
                        env.vertical_weights, control.vertical_weights
                    ):
                        raise ValueError("wrong deterministic weights")
                rows = json.loads((case / "rows.json").read_text())
                probes = symmetric_probes(
                    n, fractions=(0.125, 0.25), axes=("diagonal",)
                )
                for replica in (1, 2):
                    with np.load(
                        case / f"sample{replica}.npz", allow_pickle=False
                    ) as s:
                        state = DimerState(n, s["horizontal"], s["vertical"])
                    state.validate()
                    values = measure_probes(face_heights(state), probes)
                    for row in rows:
                        if values[row["probe"]] != row[f"value{replica}"]:
                            raise ValueError("raw observable mismatch")
                        observations += 1
                for row in rows:
                    selected = frame[
                        (frame.environment_id == row["environment_id"])
                        & (frame.probe == row["probe"])
                    ]
                    if len(selected) != 1 or any(
                        str(selected.iloc[0][k]) != str(row[k])
                        for k in (
                            "value1",
                            "value2",
                            "environment_seed",
                            "replica1_seed",
                            "replica2_seed",
                        )
                    ):
                        raise ValueError("case/table mismatch")
                pairs += 1
        else:
            reconstructed = []
            for index in range(1, count + 1):
                case = directory / "cases" / f"env-{index:05d}"
                names = set()
                for line in (case / "sha256.txt").read_text().splitlines():
                    expected, name = line.split()
                    if (
                        name
                        not in {
                            "a.txt",
                            "b.txt",
                            "sample1.txt",
                            "sample2.txt",
                            "rows.csv",
                        }
                        or name in names
                    ):
                        raise ValueError("unexpected certificate name")
                    if digest(case / name) != expected:
                        raise ValueError("Aztec raw checksum mismatch")
                    names.add(name)
                    files += 1
                if len(names) != 5:
                    raise ValueError("incomplete Aztec certificate")
                for name, mean in (("a.txt", 0.2), ("b.txt", 0.25)):
                    weights = np.loadtxt(case / name)
                    if (
                        weights.shape != (n, n)
                        or not np.isfinite(weights).all()
                        or not (weights > 0).all()
                    ):
                        raise ValueError("invalid family weights")
                    if (
                        cell["model"] != "aztec_gamma"
                        and not (
                            weights == (mean if cell["model"] == "aztec_mean" else 1)
                        ).all()
                    ):
                        raise ValueError("wrong Aztec control law")
                reconstructed.append((case / "rows.csv").read_text())
            table = (directory / "measurements.csv").read_text().split("\n", 1)[1]
            if table != "".join(reconstructed):
                raise ValueError("Aztec case/table mismatch")
            repository = Path(__file__).resolve().parents[2]
            subprocess.run(
                [
                    os.environ.get("RESEARCH_JULIA", "julia"),
                    "--startup-file=no",
                    f"--project={os.environ.get('RESEARCH_JULIA_PROJECT', repository / 'aztec')}",
                    str(repository / "aztec/scripts/verify_matched_shuffle_raw.jl"),
                    str(directory),
                    str(n),
                    str(count),
                ],
                check=True,
            )
            pairs += count
            observations += 4 * count
    report = {
        "state": "RAW_VERIFIED",
        "validation_only": validation,
        "plan_sha256": digest(plan_path),
        "tables_sha256": hashes,
        "pairs": pairs,
        "matchings": 2 * pairs,
        "observations": observations,
        "files_verified": files,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, report)
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--validation-only", action="store_true")
    args = parser.parse_args()
    verify(args.input, args.plan, args.output, args.validation_only)
