"""Verify and analyse the six-size corrected Aztec refinement after completion.

Refuses missing sizes, wrong counts, repeated seeds, incomplete probe vectors,
version mismatches and altered tables. Uses the preselected probes/windows and
10,000 whole-environment bootstrap draws; no automatic sample-size escalation.
"""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from square_glauber.followup_analysis import bootstrap_table, curvature


def load_refinement(root, expected_pairs=2048):
    frames = []
    source_versions = set()
    for n in (16, 32, 48, 64, 96, 128):
        directory = root / f"n{n}"
        summary = json.loads((directory / "summary.json").read_text())
        manifest = json.loads((directory / "manifest.json").read_text())
        source_versions.add((manifest["source_hashes"], manifest["julia"]))
        path = directory / "measurements.csv"
        if (
            not summary["all_pairs_certified"]
            or summary["environments"] != expected_pairs
        ):
            raise ValueError("incomplete size or wrong environment count")
        if (
            summary["sampler_version"] != "gamma-interval-coins-v1"
            or manifest["master_seed"] != 2026100102
        ):
            raise ValueError("wrong sampler/seed namespace")
        if manifest["order"] != n or manifest["environments"] != expected_pairs:
            raise ValueError("manifest mismatch")
        if (
            hashlib.sha256(path.read_bytes()).hexdigest()
            != summary["measurements_sha256"]
        ):
            raise ValueError("measurement checksum mismatch")
        frame = pd.read_csv(path, dtype={"environment_seed": str, "coin_seed": str})
        if (
            not frame["size"].eq(n).all()
            or frame.environment_id.nunique() != expected_pairs
            or len(frame) != 2 * expected_pairs
        ):
            raise ValueError("wrong size or incomplete pairs")
        if (
            not frame.domain.eq("aztec").all()
            or not frame.height_scale.eq(1).all()
            or not frame.axis.eq("diagonal").all()
        ):
            raise ValueError("wrong geometry/height convention")
        if set(frame.law_id) != {
            "aztec_gamma:alpha=0.2:beta=0.25:scale=1:a-north,b-east"
        }:
            raise ValueError("wrong edge law")
        if frame.duplicated(["environment_id", "probe"]).any() or set(frame.probe) != {
            "diagonal:rho=0.125",
            "diagonal:rho=0.25",
        }:
            raise ValueError("duplicate/missing probes")
        for probe, rho in [("diagonal:rho=0.125", 0.125), ("diagonal:rho=0.25", 0.25)]:
            group = frame[frame.probe == probe]
            if (
                not group.r.eq(n * rho).all()
                or not (abs(group.distance - 2**0.5 * n * rho) < 1e-12).all()
            ):
                raise ValueError("incorrect physical distance")
        frames.append(frame)
    if len(source_versions) != 1:
        raise ValueError("mixed source snapshots or Julia versions")
    all_data = pd.concat(frames, ignore_index=True)
    for _, group in all_data.groupby("environment_id"):
        if group.environment_seed.nunique() != 1 or group.coin_seed.nunique() != 1:
            raise ValueError("inconsistent replica environment identity")
    seeds = (
        all_data.drop_duplicates("environment_id")[["environment_seed", "coin_seed"]]
        .to_numpy()
        .ravel()
    )
    if len(set(seeds)) != len(seeds):
        raise ValueError("global seed collision")
    return all_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--validation-only",
        action="store_true",
        help="Four-pair integration test, never the scientific refinement",
    )
    args = parser.parse_args()
    data = load_refinement(args.input, 4 if args.validation_only else 2048)
    count = 100 if args.validation_only else 10000
    table, draws = bootstrap_table(data, count, 2026100103)
    fits = []
    for window in [(16, 128), (32, 128), (48, 128)]:
        result, _ = curvature(
            table,
            draws,
            "aztec_gamma",
            ["diagonal:rho=0.125", "diagonal:rho=0.25"],
            window,
        )
        fits.append(result)
    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "intervals.csv", index=False)
    report = dict(
        validation_only=args.validation_only,
        bootstrap=count,
        seed=2026100103,
        source_tables_sha256={
            f"n{n}": hashlib.sha256(
                (args.input / f"n{n}/measurements.csv").read_bytes()
            ).hexdigest()
            for n in (16, 32, 48, 64, 96, 128)
        },
        fits=fits,
        interpretation="Fixed-window exploratory reassessment; Float64 reductions with robust coins. Not a matched square contrast or asymptotic proof.",
    )
    (args.output / "analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                validation_only=args.validation_only,
                primary_c=fits[0]["c"],
                interval=fits[0]["interval"],
            )
        )
    )


if __name__ == "__main__":
    main()
