"""Verify every planned case and report feasibility, never a scaling fit.

Missing, interrupted and budget-exhausted cases remain in the denominator.
Completed pairs require input hashes, matching constraints and independently
recomputed height observables to agree with their saved records.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

from square_glauber.budgeted_frontier import VERSION
from square_glauber.comparison import measure_probes, symmetric_probes
from square_glauber.height import face_heights
from square_glauber.matching import DimerState


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(plan_path, root):
    plan = json.loads(plan_path.read_text())
    if plan["sampler_version"] != VERSION:
        raise ValueError("wrong sampler plan")
    rows = []
    for case in plan["cases"]:
        directory = root / f"case-{case['case_index']:03d}"
        path = directory / "status.json"
        row = dict(
            case_index=case["case_index"],
            L=case["L"],
            model=case["model"],
            state="MISSING",
            complete_pair=False,
        )
        if not path.exists():
            rows.append(row)
            continue
        status = json.loads(path.read_text())
        if (
            status["case"] != case
            or status["plan_sha256"] != digest(plan_path)
            or status["sampler_version"] != VERSION
            or status["budgets"] != plan["budgets"]
        ):
            raise ValueError("case/plan/version/budget mismatch")
        state = status["state"]
        if state not in {"COMPLETE", "BUDGET_EXHAUSTED", "STARTED", "ERROR"}:
            raise ValueError("unknown outcome")
        if bool(status["complete_pair"]) != (state == "COMPLETE"):
            raise ValueError("inconsistent completion flag")
        if state in {"COMPLETE", "BUDGET_EXHAUSTED"}:
            actual = digest(directory / "weights.npz")
            if actual != status["weights_sha256"]:
                raise ValueError("input checksum mismatch")
            if (
                case["input_source"] == "original_archive"
                and actual != case["original_sha256"]
            ):
                raise ValueError("original input identity mismatch")
        if state == "BUDGET_EXHAUSTED":
            if status["reason"] not in {"time_limit", "memory_limit", "state_limit"}:
                raise ValueError("unknown budget reason")
            if (
                list(directory.glob("sample*.npz"))
                or (directory / "measurements.json").exists()
            ):
                raise ValueError("censored case contains partial sample outputs")
        if state == "COMPLETE":
            required = {
                "weights.npz",
                "sample1.npz",
                "sample2.npz",
                "measurements.json",
            }
            if set(status["files_sha256"]) != required:
                raise ValueError("incomplete file certificate")
            for name, expected in status["files_sha256"].items():
                if digest(directory / name) != expected:
                    raise ValueError("sample/table checksum mismatch")
            probes = symmetric_probes(
                case["L"], fractions=(0.125, 0.25), axes=("diagonal",)
            )
            values = []
            for replica in (1, 2):
                with np.load(
                    directory / f"sample{replica}.npz", allow_pickle=False
                ) as a:
                    sample = DimerState(case["L"], a["horizontal"], a["vertical"])
                    sample.validate()
                    values.append(measure_probes(face_heights(sample), probes))
            expected = [
                dict(
                    probe=p.label,
                    r=p.r,
                    distance=p.distance,
                    start=list(p.start),
                    end=list(p.end),
                    height_scale=1,
                    value1=values[0][p.label],
                    value2=values[1][p.label],
                )
                for p in probes
            ]
            if json.loads((directory / "measurements.json").read_text()) != expected:
                raise ValueError("saved observables differ from remeasurement")
        row.update(
            {
                k: status[k]
                for k in [
                    "state",
                    "complete_pair",
                    "reason",
                    "cached_states",
                    "setup_seconds",
                    "total_seconds",
                    "peak_rss_bytes",
                ]
                if k in status
            }
        )
        row["status_sha256"] = digest(path)
        rows.append(row)
    outcomes = Counter(r["state"] for r in rows)
    return dict(
        expected_cases=len(rows),
        outcomes=dict(outcomes),
        complete_pairs=outcomes["COMPLETE"],
        all_pairs_complete=outcomes["COMPLETE"] == len(rows),
        all_outcomes_recorded=not any(outcomes[s] for s in ("MISSING", "STARTED")),
        plan_sha256=digest(plan_path),
        cases=rows,
        interpretation="Feasibility only; no curvature fit or omission of failed cases",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = summarize(args.plan, args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "cases"}))
