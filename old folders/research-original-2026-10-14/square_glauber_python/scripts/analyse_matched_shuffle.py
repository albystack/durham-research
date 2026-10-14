"""Verify the complete frozen comparison, then estimate paired-domain curvature.

No fit is produced for an incomplete cohort. Confidence intervals resample
whole environments independently between domains. The primary estimand is
c_Aztec-c_square, not a comparison of two separate significance decisions.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from run_matched_shuffle import digest, seed, write_json
from square_glauber.control_audit import joint_covariance_test
from square_glauber.followup_analysis import bootstrap_table, curvature

PROBES = ["diagonal:rho=0.125", "diagonal:rho=0.25"]


def load_complete(root, plan_path, validation=False):
    plan = json.loads(plan_path.read_text())
    master = plan["validation_seed"] if validation else plan["master_seed"]
    tables, hashes, allseeds = [], {}, []
    for cell in plan["cells"]:
        directory = root / f"cell-{cell['index']:02d}"
        status = json.loads((directory / "status.json").read_text())
        manifest = json.loads((directory / "manifest.json").read_text())
        count = plan["validation_pairs"] if validation else cell["pairs"]
        if status["state"] != "COMPLETE" or status["completed"] != count:
            raise ValueError(f"incomplete cell {cell['index']}")
        if (
            manifest["plan_sha256"] != digest(plan_path)
            or manifest["master_seed"] != master
            or manifest["validation_only"] != validation
            or manifest["pairs"] != count
        ):
            raise ValueError("manifest plan/count/seed mismatch")
        version = (
            "embedded-square-shuffle-v1"
            if cell["domain"] == "square"
            else "gamma-interval-coins-v1"
        )
        if manifest["sampler_version"] != version:
            raise ValueError("sampler mismatch")
        path = directory / "measurements.csv"
        hashes[directory.name] = digest(path)
        if hashes[directory.name] != status["measurements_sha256"]:
            raise ValueError("table checksum mismatch")
        frame = pd.read_csv(
            path,
            dtype={
                key: str
                for key in [
                    "environment_seed",
                    "replica1_seed",
                    "replica2_seed",
                    "coin_seed",
                ]
            },
        )
        law = f"{cell['model']}:alpha=0.2:beta=0.25:scale=1:a-north,b-east"
        if (
            len(frame) != 2 * count
            or frame.environment_id.nunique() != count
            or frame.duplicated(["environment_id", "probe"]).any()
        ):
            raise ValueError("missing or duplicate environments/probes")
        for column, value in [
            ("domain", cell["domain"]),
            ("size", cell["size"]),
            ("model", cell["model"]),
            ("law_id", law),
            ("axis", "diagonal"),
            ("height_scale", 1),
        ]:
            if not frame[column].eq(value).all():
                raise ValueError(f"wrong {column}")
        for index in range(1, count + 1):
            name = f"{cell['domain']}-{cell['model']}-n{cell['size']}-env{index}"
            group = frame[frame.environment_id == name]
            if len(group) != 2 or set(group.probe) != set(PROBES):
                raise ValueError("incomplete environment vector")
            roles = [("environment_seed", "environment")]
            roles += (
                [("replica1_seed", "replica1"), ("replica2_seed", "replica2")]
                if cell["domain"] == "square"
                else [("coin_seed", "paired-coins")]
            )
            for column, role in roles:
                expected = str(
                    seed(
                        master, cell["domain"], cell["model"], cell["size"], index, role
                    )
                )
                if not group[column].eq(expected).all():
                    raise ValueError("seed provenance mismatch")
                allseeds.append(expected)
        for probe, rho in zip(PROBES, [0.125, 0.25]):
            g = frame[frame.probe == probe]
            if not g.r.eq(cell["size"] * rho).all() or not np.allclose(
                g.distance, np.sqrt(2) * cell["size"] * rho, rtol=0, atol=1e-12
            ):
                raise ValueError("physical probe mismatch")
        if not np.isfinite(frame[["value1", "value2"]]).all().all():
            raise ValueError("non-finite observables")
        tables.append(frame)
    if len(set(allseeds)) != len(allseeds):
        raise ValueError("seed collision")
    return pd.concat(tables, ignore_index=True), hashes


def analyse(root, plan_path, output, validation=False):
    plan = json.loads(plan_path.read_text())
    data, hashes = load_complete(root, plan_path, validation)
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "validation_only": validation,
        "plan_sha256": digest(plan_path),
        "tables_sha256": hashes,
        "complete_cells": len(plan["cells"]),
        "complete_pairs": len(data) // 2,
    }
    if validation:
        report["state"] = "VALIDATED_NO_SCIENTIFIC_FIT"
    else:
        tables, fits, coefficients = [], [], {}
        for domain in ["square", "aztec"]:
            # Distinct seeds give independent domain resamples; each row of
            # their difference is one joint bootstrap replicate.
            bs_seed = plan["analysis"]["bootstrap_seed"] + (domain == "aztec")
            table, draws = bootstrap_table(
                data[data.domain == domain], plan["analysis"]["bootstrap"], bs_seed
            )
            for model in plan["models"]:
                for window in plan["analysis"]["windows"]:
                    fit, boot = curvature(table, draws, model, PROBES, window)
                    fit["domain"] = domain
                    fits.append(fit)
                    coefficients[(domain, model, tuple(window))] = (fit["c"], boot)
            table["domain"] = domain
            tables.append(table)
        controls, control_labels = [], []
        for (domain, model, size), group in data[data.model != "aztec_gamma"].groupby(
            ["domain", "model", "size"], sort=True
        ):
            first = group.pivot(
                index="environment_id", columns="probe", values="value1"
            )[PROBES]
            second = group.pivot(
                index="environment_id", columns="probe", values="value2"
            ).reindex(index=first.index, columns=PROBES)
            controls.append((first.to_numpy(), second.to_numpy()))
            control_labels.extend(
                [
                    {"domain": domain, "model": model, "size": int(size), "probe": p}
                    for p in PROBES
                ]
            )
        control_test = joint_covariance_test(
            controls, permutations=9999, seed=2026100411
        )
        control_test["labels"] = control_labels
        contrasts = []
        for model in plan["models"]:
            for window in plan["analysis"]["windows"]:
                a, ab = coefficients[("aztec", model, tuple(window))]
                s, sb = coefficients[("square", model, tuple(window))]
                contrasts.append(
                    {
                        "model": model,
                        "window": window,
                        "delta_c": a - s,
                        "interval": np.quantile(ab - sb, [0.025, 0.975]).tolist(),
                        "primary": model == "aztec_gamma" and window == [8, 64],
                    }
                )
        pd.concat(tables, ignore_index=True).to_csv(
            output / "intervals.csv", index=False
        )
        report.update(
            state="ANALYSED",
            fits=fits,
            contrasts=contrasts,
            joint_control_test=control_test,
            interpretation="Finite-window estimates; zero-compatible square curvature does not prove absence of a squared-log term. Cutoff windows and controls are sensitivity analyses, not alternative primary tests.",
        )
    write_json(output / "analysis.json", report)
    print(
        json.dumps(
            {k: v for k, v in report.items() if k not in ("tables_sha256", "fits")}
        )
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--validation-only", action="store_true")
    args = parser.parse_args()
    analyse(args.input, args.plan, args.output, args.validation_only)
