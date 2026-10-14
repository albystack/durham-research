"""Restartable batches for the fixed larger Aztec/square comparison.

The scientific samplers are unchanged. Each batch has a distinct recorded
master seed and at most 512 pairs. Completed batches are checked and reused;
an interrupted attempt is preserved and requires explicit --retry. A retry
uses the same plan and seeds, never replacement environments. Raw verification
must pass before a batch receives its completion receipt. Analysis requires
every declared batch and resamples whole environment vectors across batches.
"""

import argparse
import fcntl
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from analyse_matched_shuffle import PROBES, load_complete
from run_matched_shuffle import digest, run_square, write_json
from square_glauber.followup_analysis import bootstrap_table, curvature
from verify_matched_shuffle_raw import verify

ROOT = Path(__file__).resolve().parents[2]


def read_plan(path):
    plan = json.loads(path.read_text())
    batches = plan["batches"]
    if [b["index"] for b in batches] != list(range(len(batches))):
        raise ValueError("batch indices must be consecutive")
    if len({b["master_seed"] for b in batches}) != len(batches):
        raise ValueError("duplicate batch master seed")
    for b in batches:
        if (
            b["domain"] not in ("square", "aztec")
            or b["model"] not in plan["models"]
            or b["size"] not in plan["sizes"]
            or not 2 <= b["pairs"] <= 512
        ):
            raise ValueError("invalid batch")
    for domain in ("square", "aztec"):
        for model, count in plan["models"].items():
            for size in plan["sizes"]:
                actual = sum(
                    b["pairs"]
                    for b in batches
                    if (b["domain"], b["model"], b["size"]) == (domain, model, size)
                )
                if actual != count:
                    raise ValueError("incomplete domain/model/size allocation")
    return plan


def batch_plan(plan, batch):
    """Adapt one batch to the previously validated table/seed/raw verifier."""
    return {
        "namespace": "matched-shuffle-20261004",
        "master_seed": batch["master_seed"],
        "validation_seed": batch["master_seed"],
        "validation_pairs": batch["pairs"],
        "sizes": [batch["size"]],
        "probes": [0.125, 0.25],
        "models": {batch["model"]: batch["pairs"]},
        "purpose": plan["purpose"],
        "cells": [dict(batch, index=0)],
    }


def checked_batch(plan_path, batch, root):
    """Validate the receipt and compact evidence; raw checking occurred in worker."""
    plan = read_plan(plan_path)
    directory = root / f"batch-{batch['index']:04d}"
    receipt = json.loads((directory / "complete.json").read_text())
    if receipt["plan_sha256"] != digest(plan_path) or receipt["batch"] != batch:
        raise ValueError("receipt belongs to another plan or batch")
    if not re.fullmatch(r"attempt-[0-9]{4}", receipt["attempt"]):
        raise ValueError("invalid attempt path")
    attempt = directory / receipt["attempt"]
    for name, expected in receipt["files"].items():
        if name not in (
            "plan.json",
            "raw_verification.json",
            "cells/cell-00/measurements.csv",
        ):
            raise ValueError("unexpected receipt file")
        if digest(attempt / name) != expected:
            raise ValueError("completed batch evidence changed")
    if set(receipt["files"]) != {
        "plan.json",
        "raw_verification.json",
        "cells/cell-00/measurements.csv",
    }:
        raise ValueError("incomplete receipt")
    if json.loads((attempt / "plan.json").read_text()) != batch_plan(plan, batch):
        raise ValueError("batch plan changed")
    frame, hashes = load_complete(attempt / "cells", attempt / "plan.json")
    raw = json.loads((attempt / "raw_verification.json").read_text())
    if (
        raw["state"] != "RAW_VERIFIED"
        or raw["pairs"] != batch["pairs"]
        or raw["matchings"] != 2 * batch["pairs"]
        or raw["observations"] != 4 * batch["pairs"]
        or raw["tables_sha256"] != hashes
        or raw["plan_sha256"] != digest(attempt / "plan.json")
    ):
        raise ValueError("inconsistent raw-verification certificate")
    return frame, receipt


def run_batch(plan_path, index, root, retry=False):
    plan = read_plan(plan_path)
    batch = plan["batches"][index]
    directory = root / f"batch-{index:04d}"
    directory.mkdir(parents=True, exist_ok=True)
    # The operating system releases this advisory lock even after a scheduler kill.
    with (directory / "writer.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        identity = {"plan_sha256": digest(plan_path), "batch": batch}
        identity_path = directory / "identity.json"
        if identity_path.exists() and json.loads(identity_path.read_text()) != identity:
            raise ValueError("refusing to change batch identity")
        write_json(identity_path, identity)
        if (directory / "complete.json").exists():
            checked_batch(plan_path, batch, root)
            print(f"batch {index}: verified existing completion", flush=True)
            return
        previous = sorted(directory.glob("attempt-*"))
        if previous and not retry:
            raise ValueError("incomplete attempt preserved; explicit --retry required")
        attempt = directory / f"attempt-{len(previous) + 1:04d}"
        attempt.mkdir()
        bp = attempt / "plan.json"
        write_json(bp, batch_plan(plan, batch))
        destination = attempt / "cells/cell-00"
        if batch["domain"] == "square":
            run_square(bp, batch_plan(plan, batch)["cells"][0], destination)
        else:
            command = [
                os.environ.get("RESEARCH_JULIA", "julia"),
                "--startup-file=no",
                f"--project={os.environ.get('RESEARCH_JULIA_PROJECT', ROOT / 'aztec')}",
                str(ROOT / "aztec/scripts/run_larger_matched_batch.jl"),
                "--output",
                str(destination),
                "--size",
                str(batch["size"]),
                "--model",
                batch["model"],
                "--environments",
                str(batch["pairs"]),
                "--seed",
                str(batch["master_seed"]),
                "--plan-sha256",
                digest(bp),
                "--validation",
                "false",
            ]
            subprocess.run(command, check=True)
        # A completed record in an interrupted attempt must replay identically.
        for old in previous:
            for prior in (old / "cells/cell-00/cases").glob("env-*"):
                cert = prior / (
                    "sha256.txt" if batch["domain"] == "aztec" else "rows.json"
                )
                if cert.exists():
                    fresh = destination / "cases" / prior.name
                    for name in (
                        ("a.txt", "b.txt", "sample1.txt", "sample2.txt", "rows.csv")
                        if batch["domain"] == "aztec"
                        else ("weights.npz", "sample1.npz", "sample2.npz", "rows.json")
                    ):
                        if not (fresh / name).exists() or digest(
                            prior / name
                        ) != digest(fresh / name):
                            raise ValueError(
                                "retry changed a previously completed case"
                            )
        verify(attempt / "cells", bp, attempt / "raw_verification.json")
        receipt = dict(
            identity,
            attempt=attempt.name,
            files={
                name: digest(attempt / name)
                for name in (
                    "plan.json",
                    "raw_verification.json",
                    "cells/cell-00/measurements.csv",
                )
            },
        )
        write_json(directory / "complete.json", receipt)
        checked_batch(plan_path, batch, root)
        print(f"batch {index}: all {batch['pairs']} pairs raw-verified", flush=True)


def collect(plan_path, root):
    plan = read_plan(plan_path)
    frames, receipts, seeds = [], {}, set()
    for batch in plan["batches"]:
        frame, receipt = checked_batch(plan_path, batch, root)
        # Environment labels inside each legacy batch start at one. Prefixing
        # with the batch identity makes them unique without changing raw data.
        prefix = f"batch-{batch['index']:04d}:"
        frame["environment_id"] = prefix + frame.environment_id
        roles = ["environment_seed"] + (
            ["replica1_seed", "replica2_seed"]
            if batch["domain"] == "square"
            else ["coin_seed"]
        )
        for _, group in frame.groupby("environment_id", sort=False):
            for role in roles:
                value = str(group[role].iloc[0])
                if value in seeds:
                    raise ValueError("cross-batch role seed collision")
                seeds.add(value)
        frames.append(frame)
        receipts[str(batch["index"])] = receipt
    return pd.concat(frames, ignore_index=True), receipts


def control_curvature_test(data, table, windows, permutations, seed):
    """Predeclared maximum-standardised test of all six control contrasts."""
    labels = [
        (model, window) for model in ("aztec_mean", "uniform") for window in windows
    ]
    groups = []
    null_covariance = np.zeros((len(labels), len(labels)))
    observed = np.zeros(len(labels))
    for (domain, model, size), group in data[data.model != "aztec_gamma"].groupby(
        ["domain", "model", "size"], sort=True
    ):
        a = group.pivot(index="environment_id", columns="probe", values="value1")[
            PROBES
        ].to_numpy(float)
        b = group.pivot(index="environment_id", columns="probe", values="value2")[
            PROBES
        ].to_numpy(float)
        a, b = a - a.mean(0), b - b.mean(0)
        w = np.zeros((len(labels), 2))
        for j, (target, window) in enumerate(labels):
            if target != model or not window[0] <= size <= window[1]:
                continue
            selected = table[
                (table.domain == domain)
                & (table.model == model)
                & table["size"].between(*window)
            ]
            design = []
            for row in selected.itertuples():
                z = np.log(row.distance)
                v = np.zeros(5)
                k = PROBES.index(row.probe)
                v[2 * k : 2 * k + 2] = [1, z]
                v[-1] = z * z
                design.append(v)
            influence = np.linalg.pinv(design)[-1]
            for row, value in zip(selected.itertuples(), influence):
                if row.size == size:
                    w[j, PROBES.index(row.probe)] = value * (
                        1 if domain == "aztec" else -1
                    )
        observed += w @ ((a * b).sum(0) / (len(a) - 1))
        null_covariance += w @ ((a.T @ a) * (b.T @ b) / (len(a) - 1) ** 3) @ w.T
        groups.append((a, b, w))
    sd = np.sqrt(np.diag(null_covariance))
    if np.any(sd <= 0):
        raise ValueError("degenerate control permutation variance")
    rng = np.random.default_rng(seed)
    maxima = np.empty(permutations)
    for k in range(permutations):
        value = sum(
            w @ ((a * b[rng.permutation(len(b))]).sum(0) / (len(b) - 1))
            for a, b, w in groups
        )
        maxima[k] = np.max(abs(value / sd))
    p = float((1 + (maxima >= np.max(abs(observed / sd))).sum()) / (permutations + 1))
    return {
        "predeclared": True,
        "permutations": permutations,
        "seed": seed,
        "joint_p": p,
        "contrasts": [
            {
                "model": m,
                "window": window,
                "delta_c": float(observed[j]),
                "null_sd": float(sd[j]),
            }
            for j, (m, window) in enumerate(labels)
        ],
    }


def analyse(plan_path, root, output):
    plan = read_plan(plan_path)
    data, receipts = collect(plan_path, root)
    settings = plan["analysis"]
    tables, fits, estimates = [], [], {}
    for domain in ("square", "aztec"):
        table, draws = bootstrap_table(
            data[data.domain == domain],
            settings["bootstrap"],
            settings["bootstrap_seed"] + (domain == "aztec"),
        )
        for model in plan["models"]:
            for window in settings["windows"]:
                fit, boot = curvature(table, draws, model, PROBES, window)
                fit["domain"] = domain
                fits.append(fit)
                estimates[domain, model, tuple(window)] = (fit["c"], boot)
        table["domain"] = domain
        tables.append(table)
    contrasts = []
    for model in plan["models"]:
        for window in settings["windows"]:
            a, ab = estimates["aztec", model, tuple(window)]
            s, sb = estimates["square", model, tuple(window)]
            contrasts.append(
                {
                    "model": model,
                    "window": window,
                    "delta_c": a - s,
                    "interval": np.quantile(ab - sb, [0.025, 0.975]).tolist(),
                    "primary": model == "aztec_gamma"
                    and window == settings["primary_window"],
                }
            )
    table = pd.concat(tables, ignore_index=True)
    controls = control_curvature_test(
        data,
        table,
        settings["windows"],
        settings["permutations"],
        settings["control_seed"],
    )
    for c in controls["contrasts"]:
        other = next(
            x
            for x in contrasts
            if x["model"] == c["model"] and x["window"] == c["window"]
        )
        if abs(c["delta_c"] - other["delta_c"]) > 1e-10:
            raise ValueError("independent control coefficient disagreement")
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "intervals.csv", index=False)
    write_json(output / "batch_receipts.json", receipts)
    write_json(
        output / "analysis.json",
        {
            "state": "VALIDATION_ONLY" if plan["validation_only"] else "ANALYSED",
            "plan_sha256": digest(plan_path),
            "pairs": len(data) // 2,
            "batches": len(receipts),
            "batch_receipts_sha256": digest(output / "batch_receipts.json"),
            "fits": fits,
            "contrasts": contrasts,
            "control_curvature_test": controls,
            "interpretation": "Validation data: no scientific inference"
            if plan["validation_only"]
            else "Finite-window matched comparison; an interval containing zero does not prove absence or equivalence. No adaptive extension.",
        },
    )
    print(
        json.dumps(
            {
                "pairs": len(data) // 2,
                "batches": len(receipts),
                "validation_only": plan["validation_only"],
            }
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["run", "analyse"])
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--batch", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--retry", action="store_true")
    args = parser.parse_args()
    if args.action == "run":
        if args.batch is None:
            parser.error("--batch required")
        run_batch(args.plan.resolve(), args.batch, args.root.resolve(), args.retry)
    else:
        if args.output is None:
            parser.error("--output required")
        analyse(args.plan.resolve(), args.root.resolve(), args.output.resolve())
