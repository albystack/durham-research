"""Prospectively selected, whole-environment bootstrap of curvature differences.

Input is the common long-form schema written by comparison_pilot and the Julia
bridge. No old rightward data are silently relabelled. Production windows and
probe lists must be specified explicitly. Censored campaigns are refused.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def prepare_blocks(frame: pd.DataFrame, probes: list[str], minimum: int, maximum: int):
    """Reject partial vectors and rank-deficient fits instead of dropping rows."""
    selected = frame[
        frame["size"].between(minimum, maximum) & frame.probe.isin(probes)
    ].copy()
    if selected.empty or selected.duplicated(["size", "environment_id", "probe"]).any():
        raise ValueError("empty selection or duplicate environment/probe")
    blocks = []
    design = []
    descriptions = []
    for size, group in selected.groupby("size", sort=True):
        first = group.pivot(
            index="environment_id", columns="probe", values="value1"
        ).reindex(columns=probes)
        second = group.pivot(
            index="environment_id", columns="probe", values="value2"
        ).reindex(index=first.index, columns=probes)
        if len(first) < 3 or first.isna().any().any() or second.isna().any().any():
            raise ValueError(
                "each size needs >=3 complete environment vectors for every declared probe"
            )
        if (
            not np.isfinite(first.to_numpy()).all()
            or not np.isfinite(second.to_numpy()).all()
        ):
            raise ValueError("nonfinite measurements")
        blocks.append((first.to_numpy(float), second.to_numpy(float)))
        for index, probe in enumerate(probes):
            rows = group[group.probe.eq(probe)]
            if rows.distance.nunique() != 1 or float(rows.distance.iloc[0]) <= 0:
                raise ValueError("inconsistent/nonpositive physical distance")
            z = np.log(float(rows.distance.iloc[0]))
            row = np.zeros(2 * len(probes) + 1)
            row[2 * index : 2 * index + 2] = [1, z]
            row[-1] = z * z
            design.append(row)
            descriptions.append(
                dict(size=int(size), probe=probe, distance=float(rows.distance.iloc[0]))
            )
    matrix = np.asarray(design)
    if len(blocks) < 4 or np.linalg.matrix_rank(matrix) != matrix.shape[1]:
        raise ValueError("need >=4 sizes and a full-rank pooled log-quadratic design")
    # Four distinct physical distances per probe avoid nominal-size replication
    # masquerading as evidence of curvature at a rounded small separation.
    for probe in probes:
        if len({r["distance"] for r in descriptions if r["probe"] == probe}) < 4:
            raise ValueError("each declared probe needs >=4 distinct actual distances")
    return blocks, matrix, descriptions


def covariance_vector(blocks, rng=None):
    result = []
    for a, b in blocks:
        if rng is not None:
            indices = rng.integers(0, len(a), len(a))
            a, b = a[indices], b[indices]
        result.extend(((a - a.mean(0)) * (b - b.mean(0))).sum(0) / (len(a) - 1))
    return np.asarray(result)


def compare_curvature(
    aztec: pd.DataFrame,
    square: pd.DataFrame,
    *,
    probes: list[str],
    aztec_window: tuple[int, int],
    square_window: tuple[int, int],
    bootstrap=10000,
    seed=2026092904,
) -> dict:
    """Independent campaign bootstraps, preserving joint probes within each."""
    if bootstrap < 2 or not probes or len(set(probes)) != len(probes):
        raise ValueError("need distinct declared probes and >=2 bootstrap repetitions")
    laws = set(aztec.law_id.unique()) | set(square.law_id.unique())
    if len(laws) != 1:
        raise ValueError("comparison requires the same single weight law")
    if not set(aztec.domain) == {"aztec"} or not set(square.domain) == {"square"}:
        raise ValueError("domain labels are missing or swapped")
    for data in (aztec, square):
        if not data.height_scale.eq(1).all() or not data.axis.eq("diagonal").all():
            raise ValueError(
                "comparison requires diagonal probes in common unit height scale"
            )
    prepared = [
        prepare_blocks(data, probes, *window)
        for data, window in ((aztec, aztec_window), (square, square_window))
    ]
    inverses = [np.linalg.pinv(item[1]) for item in prepared]
    point = []
    fits = []
    for (blocks, matrix, description), inverse in zip(prepared, inverses):
        y = covariance_vector(blocks)
        coefficients = inverse @ y
        log_matrix = matrix[:, :-1]
        log_coefficients = np.linalg.lstsq(log_matrix, y, rcond=None)[0]
        point.append(float(coefficients[-1]))
        fits.append(
            dict(
                rows=description,
                covariances=y.tolist(),
                coefficients=coefficients.tolist(),
                log_sse=float(np.sum((y - log_matrix @ log_coefficients) ** 2)),
                quadratic_sse=float(np.sum((y - matrix @ coefficients) ** 2)),
                residuals=(y - matrix @ coefficients).tolist(),
                condition_number=float(np.linalg.cond(matrix)),
            )
        )
    sequences = np.random.SeedSequence(seed).spawn(2)
    draws = []
    for (blocks, _, _), inverse, sequence in zip(prepared, inverses, sequences):
        rng = np.random.default_rng(sequence)
        draws.append(
            np.array(
                [
                    float((inverse @ covariance_vector(blocks, rng))[-1])
                    for _ in range(bootstrap)
                ]
            )
        )
    difference = draws[0] - draws[1]

    def interval(values):
        return np.quantile(values, [0.025, 0.975]).tolist()

    return dict(
        law_id=next(iter(laws)),
        probes=probes,
        aztec_window=aztec_window,
        square_window=square_window,
        bootstrap_replicates=bootstrap,
        seed=seed,
        aztec_c=point[0],
        square_c=point[1],
        delta_c=point[0] - point[1],
        aztec_interval=interval(draws[0]),
        square_interval=interval(draws[1]),
        delta_interval=interval(difference),
        bootstrap_fraction_delta_positive=float(np.mean(difference > 0)),
        fit_details=dict(aztec=fits[0], square=fits[1]),
        interpretation="Finite-window curvature contrast, not an asymptotic or boundary-causality proof. Size labels are domain-specific; fit uses recorded physical distance.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for domain in ("aztec", "square"):
        parser.add_argument(f"--{domain}", type=Path, required=True)
        parser.add_argument(f"--{domain}-summary", type=Path, required=True)
        parser.add_argument(f"--{domain}-window", type=int, nargs=2, required=True)
    parser.add_argument("--law-id", required=True)
    parser.add_argument("--probes", nargs="+", required=True)
    parser.add_argument("--bootstrap", type=int, default=10000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = []
    for domain in ("aztec", "square"):
        summary = json.loads(getattr(args, domain + "_summary").read_text())
        if not summary.get("all_pairs_certified", False):
            raise ValueError(
                f"{domain} campaign incomplete; no automatic complete-campaign contrast"
            )
        data_path = getattr(args, domain)
        digest = hashlib.sha256(data_path.read_bytes()).hexdigest()
        if summary.get("measurements_sha256") != digest:
            raise ValueError(
                f"{domain} summary does not authenticate the measurement table"
            )
        frame = pd.read_csv(data_path)
        data.append(frame[frame.law_id.eq(args.law_id) & frame.axis.eq("diagonal")])
    report = compare_curvature(
        *data,
        probes=args.probes,
        aztec_window=tuple(args.aztec_window),
        square_window=tuple(args.square_window),
        bootstrap=args.bootstrap,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("delta_c", "delta_interval")}, indent=2))


if __name__ == "__main__":
    main()
