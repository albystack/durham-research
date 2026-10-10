"""Exploratory environment-block bootstrap for recovered and corrected pilots.

Within a size/model, one multinomial resample moves the complete environment
vector: all probes and both replicas together. Different sizes/models are
independently resampled. Percentile intervals are pointwise, not simultaneous.
Fits describe finite windows and cannot repair missing/censored environments.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def paired_moments(a, b, counts=None):
    """Unbiased sample T,D,V; optional multinomial counts of total N per row."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or a.ndim != 2 or len(a) < 2:
        raise ValueError("matching N-by-probe matrices, N>=2 required")
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError("finite data required")
    n = len(a)
    if counts is None:
        counts = np.ones((1, n))
    counts = np.asarray(counts)
    if (
        counts.ndim != 2
        or counts.shape[1] != n
        or np.any(counts < 0)
        or not np.all(counts.sum(1) == n)
    ):
        raise ValueError("counts must be nonnegative resamples of N environments")
    # Fixed centring reduces cancellation and leaves every resampled moment invariant.
    a, b = a - a.mean(0), b - b.mean(0)
    sa, sb = counts @ a, counts @ b
    va = (counts @ (a * a) - sa * sa / n) / (n - 1)
    vb = (counts @ (b * b) - sb * sb / n) / (n - 1)
    d = (counts @ (a * b) - sa * sb / n) / (n - 1)
    v = (va + vb) / 2
    return np.stack((v - d, d, v), axis=-1)


def bootstrap_table(frame, repetitions, seed):
    """Return moments table and aligned bootstrap disorder covariance columns."""
    if (
        repetitions < 2
        or frame.duplicated(["model", "size", "environment_id", "probe"]).any()
    ):
        raise ValueError("need >=2 repetitions and unique environment/probe rows")
    rng = np.random.default_rng(seed)
    rows, draws = [], []
    for (model, size), group in frame.groupby(["model", "size"], sort=True):
        a = group.pivot(index="environment_id", columns="probe", values="value1")
        b = group.pivot(
            index="environment_id", columns="probe", values="value2"
        ).reindex(index=a.index, columns=a.columns)
        if a.isna().any().any() or b.isna().any().any():
            raise ValueError("partial probe vectors may not be dropped")
        x, y = a.to_numpy(), b.to_numpy()
        point = paired_moments(x, y)[0]
        boot = np.empty((repetitions, len(a.columns), 3))
        for start in range(0, repetitions, 100):
            weights = rng.multinomial(
                len(a), np.full(len(a), 1 / len(a)), size=min(100, repetitions - start)
            )
            boot[start : start + len(weights)] = paired_moments(x, y, weights)
        for j, probe in enumerate(a.columns):
            g = group[group.probe == probe]
            if g.distance.nunique() != 1:
                raise ValueError("inconsistent physical distances")
            row = dict(
                model=model,
                size=int(size),
                probe=probe,
                pairs=len(a),
                distance=float(g.distance.iloc[0]),
            )
            for k, name in enumerate(("thermal", "disorder", "total")):
                row[name] = float(point[j, k])
                row[name + "_low"], row[name + "_high"] = np.quantile(
                    boot[:, j, k], [0.025, 0.975]
                ).tolist()
            rows.append(row)
        draws.append(boot[:, :, 1])
    return pd.DataFrame(rows), np.concatenate(draws, axis=1)


def curvature(table, draws, model, probes, window):
    """Probe-specific intercept/log slope and one shared squared-log coefficient."""
    selected = table[
        (table.model == model)
        & table.probe.isin(probes)
        & table["size"].between(*window)
    ]
    if selected["size"].nunique() < 4:
        raise ValueError("fewer than four sizes")
    if any(set(group.probe) != set(probes) for _, group in selected.groupby("size")):
        raise ValueError("incomplete probe coverage across sizes")
    design = []
    for _, row in selected.iterrows():
        z = np.log(row.distance)
        x = np.zeros(2 * len(probes) + 1)
        j = probes.index(row.probe)
        x[2 * j : 2 * j + 2] = [1, z]
        x[-1] = z * z
        design.append(x)
    for probe in probes:
        if selected[selected.probe == probe].distance.nunique() < 4:
            raise ValueError("fewer than four distinct distances for a probe")
    design = np.asarray(design)
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("rank-deficient design")
    inverse = np.linalg.pinv(design)
    y = selected.disorder.to_numpy()
    coeff = inverse @ y
    boot = draws[:, selected.index] @ inverse[-1]
    residual = y - design @ coeff
    linear = design[:, :-1]
    linear_coeff = np.linalg.lstsq(linear, y, rcond=None)[0]
    result = dict(
        model=model,
        probes=probes,
        window=list(window),
        c=float(coeff[-1]),
        interval=np.quantile(boot, [0.025, 0.975]).tolist(),
        condition_number=float(np.linalg.cond(design)),
        log_sse=float(np.sum((y - linear @ linear_coeff) ** 2)),
        quadratic_sse=float(residual @ residual),
        rows=selected[["size", "probe", "distance", "disorder"]].to_dict("records"),
        residuals=residual.tolist(),
    )
    return result, boot


def recovered_frame(root):
    inputs = (
        root
        / "archive/history/square_glauber_python/reports/provisional_20260915/input"
    )
    recovered = (
        root
        / "square_glauber_python/outputs/hamilton_matched_20260930/september_remeasurement/exploratory_measurements.csv"
    )
    new = pd.read_csv(recovered).rename(columns={"L": "size"})
    old = pd.read_csv(inputs / "spatial_measurements.csv").rename(
        columns={"L": "size", "DeltaH1": "value1", "DeltaH2": "value2"}
    )
    old["probe"] = old.requested_rho.map(lambda x: f"rightward:rho={x:.12g}")
    old["distance"] = old.r
    centre = pd.read_csv(inputs / "center_measurements.csv").rename(
        columns={"L": "size", "H1_center": "value1", "H2_center": "value2"}
    )
    centre["probe"] = "center"
    centre["distance"] = centre["size"]
    assert (
        set(new.environment_id) == set(centre.environment_id) == set(old.environment_id)
    )
    columns = [
        "model",
        "size",
        "environment_id",
        "probe",
        "distance",
        "value1",
        "value2",
    ]
    return pd.concat(
        [new[columns], old[columns], centre[columns]], ignore_index=True
    ), [
        recovered,
        inputs / "spatial_measurements.csv",
        inputs / "center_measurements.csv",
    ]


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["recovered", "aztec"])
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[3]
    )
    args = parser.parse_args()
    root = args.root
    out = root / "docs/validation/followup_20260930"
    plan = json.loads((out / "analysis_plan.json").read_text())
    if args.mode == "recovered":
        frame, inputs = recovered_frame(root)
        config = plan["recovered_analysis"]
        seed = config["seed"]
        models = sorted(frame.model.unique())
        probe_sets = [
            [f"{axis}:rho={rho:.12g}" for rho in config["probes_rho"]]
            for axis in config["axes"]
        ]
        probe_sets += [[f"{axis}:rho=0.25"] for axis in config["axes"]] + [["center"]]
    else:
        directory = root / "aztec/output/corrected_reassessment_20260930"
        inputs = [directory / "measurements.csv"]
        frame = pd.read_csv(inputs[0])
        config = plan["aztec_pilot"]
        seed = config["bootstrap_seed"]
        summary = json.loads((directory / "summary.json").read_text())
        assert (
            summary["measurements_sha256"]
            == hashlib.sha256(inputs[0].read_bytes()).hexdigest()
        )
        models = ["aztec_gamma"]
        probe_sets = [[f"diagonal:rho={rho:.12g}" for rho in config["probes_rho"]]]
        # Explicit exploratory amendment: these two probes cover all six sizes.
        probe_sets.append(["diagonal:rho=0.125", "diagonal:rho=0.25"])
    table, draws = bootstrap_table(frame, config["bootstrap"], seed)
    table.to_csv(out / f"{args.mode}_intervals.csv", index=False)
    fits, refused, paired = {}, [], []
    for model in models:
        for probes in probe_sets:
            for window in config["windows"]:
                key = f"{model}|{','.join(probes)}|{window[0]}-{window[1]}"
                try:
                    result, boot = curvature(table, draws, model, probes, window)
                    fits[key] = (result, boot)
                except ValueError as error:
                    refused.append(
                        dict(
                            model=model, probes=probes, window=window, reason=str(error)
                        )
                    )
    if args.mode == "recovered":
        for key, (fit, boot) in fits.items():
            if "horizontal:" not in key:
                continue
            other = key.replace("horizontal:", "rightward:")
            if other in fits:
                reference, reference_boot = fits[other]
                paired.append(
                    dict(
                        model=fit["model"],
                        probes=fit["probes"],
                        window=fit["window"],
                        delta_c=fit["c"] - reference["c"],
                        interval=np.quantile(
                            boot - reference_boot, [0.025, 0.975]
                        ).tolist(),
                        meaning="symmetric horizontal minus original rightward, paired environment bootstrap",
                    )
                )
    report = dict(
        mode=args.mode,
        bootstrap=config["bootstrap"],
        seed=seed,
        input_sha256={
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in inputs
        },
        fits=[f for f, _ in fits.values()],
        refused_fits=refused,
        paired_geometry_contrasts=paired,
        interpretation="Exploratory finite-window OLS curvature; pointwise percentile intervals; no asymptotic claim. September intervals do not account for its unresolved environment.",
    )
    (out / f"{args.mode}_analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"{args.mode}: {len(table)} moment rows, {len(fits)} fits, {len(refused)} refused designs"
    )


if __name__ == "__main__":
    main()
