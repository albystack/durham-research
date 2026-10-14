"""Reproduce the completed spatial-height analysis from the public measurements.

The default checks the data, recomputes every point estimate, and uses the
recorded Hamilton bootstrap intervals. --bootstrap 10000 independently reruns
the fixed whole-environment bootstrap and the joint control permutation test.
Run from any directory; output defaults to the repository's output/analysis.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [
    str(ROOT / "square_glauber_python/src"),
    str(ROOT / "square_glauber_python/scripts"),
]
from square_glauber.followup_analysis import bootstrap_table, curvature, paired_moments

DOMAINS = ("square", "aztec")
PROBES = ["diagonal:rho=0.125", "diagonal:rho=0.25"]


def load_data(directory):
    """Reject changed files, missing pairs and invalid probe geometry."""
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, checksum in manifest["sha256"].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != checksum:
            raise ValueError(f"Checksum mismatch: {name}")
    data = pd.read_csv(
        directory / "measurements.csv.gz",
        dtype={
            key: str
            for key in (
                "environment_seed",
                "replica1_seed",
                "replica2_seed",
                "coin_seed",
            )
        },
    )
    plan = json.loads((directory / "plan.json").read_text())
    if (
        len(data) != 2 * manifest["pairs"]
        or data.environment_id.nunique() != manifest["pairs"]
    ):
        raise ValueError("Incomplete paired-environment dataset")
    if data.duplicated(["environment_id", "probe"]).any():
        raise ValueError("Duplicate environment/probe record")
    if not np.isfinite(data[["value1", "value2", "distance"]]).all().all():
        raise ValueError("Non-finite measurements")
    if set(data.domain) != set(DOMAINS) or set(data.model) != set(plan["models"]):
        raise ValueError("Unexpected domains or models")
    if set(data.probe) != set(PROBES) or set(data["size"]) != set(plan["sizes"]):
        raise ValueError("Unexpected probes or sizes")
    if not data.axis.eq("diagonal").all() or not data.height_scale.eq(1).all():
        raise ValueError("Wrong height convention")
    for (domain, model, size), group in data.groupby(["domain", "model", "size"]):
        if len(group) != 2 * plan["models"][model]:
            raise ValueError(f"Incomplete cell: {domain}, {model}, {size}")
        for probe, fraction in zip(PROBES, (0.125, 0.25)):
            selected = group[group.probe == probe]
            if len(selected) != plan["models"][model] or not np.allclose(
                selected.distance, np.sqrt(2) * size * fraction, rtol=0, atol=1e-12
            ):
                raise ValueError("Missing probes or inconsistent physical distance")
    return data, plan


def point_table(data):
    """One row per domain/model/size/probe: T, D and their sum V."""
    rows = []
    for (domain, model, size), group in data.groupby(
        ["domain", "model", "size"], sort=True
    ):
        first = group.pivot(index="environment_id", columns="probe", values="value1")[
            PROBES
        ]
        second = group.pivot(
            index="environment_id", columns="probe", values="value2"
        ).reindex(index=first.index, columns=PROBES)
        if first.isna().any().any() or second.isna().any().any():
            raise ValueError("Partial environment vectors cannot be dropped")
        moments = paired_moments(first.to_numpy(), second.to_numpy())[0]
        for j, probe in enumerate(PROBES):
            rows.append(
                {
                    "domain": domain,
                    "model": model,
                    "size": int(size),
                    "probe": probe,
                    "pairs": len(first),
                    "distance": float(group[group.probe == probe].distance.iloc[0]),
                    "thermal": float(moments[j, 0]),
                    "disorder": float(moments[j, 1]),
                    "total": float(moments[j, 2]),
                }
            )
    return pd.DataFrame(rows)


def fit_coefficients(table, plan, recorded):
    """Check the coefficients against the released analysis; retain its intervals."""
    fits = []
    for domain in DOMAINS:
        selected = table[table.domain == domain].reset_index(drop=True)
        # Two identical draws let the shared fitting routine compute point estimates.
        draws = np.tile(selected.disorder.to_numpy(), (2, 1))
        for model in plan["models"]:
            for window in plan["analysis"]["windows"]:
                fit, _ = curvature(selected, draws, model, PROBES, window)
                original = next(
                    f
                    for f in recorded["fits"]
                    if (f["domain"], f["model"], f["window"]) == (domain, model, window)
                )
                if abs(fit["c"] - original["c"]) > 1e-10:
                    raise ValueError(
                        "Recomputed coefficient disagrees with released result"
                    )
                fit.update(domain=domain, interval=original["interval"])
                fits.append(fit)
    return fits


def rerun_uncertainty(data, plan, repetitions):
    """Resample whole paired environments, keeping both probes together."""
    from larger_campaign import control_curvature_test

    settings = plan["analysis"]
    fits, tables, estimates = [], [], {}
    for domain in DOMAINS:
        table, draws = bootstrap_table(
            data[data.domain == domain],
            repetitions,
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
    table = pd.concat(tables, ignore_index=True)
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
    controls = control_curvature_test(
        data,
        table,
        settings["windows"],
        settings["permutations"],
        settings["control_seed"],
    )
    return table, fits, contrasts, controls


def plot_results(table, fits, destination):
    """Show both Gamma covariance curves in the same physical units."""
    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(
        1, 2, figsize=(10, 3.8), sharey=True, constrained_layout=True
    )
    for ax, domain, title in zip(
        axes, ("aztec", "square"), ("Aztec diamond", "Ordinary square grid")
    ):
        fit = next(
            f
            for f in fits
            if f["domain"] == domain
            and f["model"] == "aztec_gamma"
            and f["window"] == [16, 128]
        )
        selected = table[(table.domain == domain) & (table.model == "aztec_gamma")]
        design, values = [], []
        for row in selected.itertuples():
            x = np.zeros(5)
            j = PROBES.index(row.probe)
            z = np.log(row.distance)
            x[2 * j : 2 * j + 2] = [1, z]
            x[-1] = z * z
            design.append(x)
            values.append(row.disorder)
        beta = np.linalg.lstsq(design, values, rcond=None)[0]
        for j, (probe, color) in enumerate(zip(PROBES, ("#156d83", "#c7682b"))):
            g = selected[selected.probe == probe].sort_values("distance")
            yerr = np.vstack(
                (g.disorder - g.disorder_low, g.disorder_high - g.disorder)
            )
            ax.errorbar(
                np.log(g.distance),
                g.disorder,
                yerr=yerr,
                fmt="o",
                ms=4,
                color=color,
                capsize=2,
                label=f"probe fraction = 1/{8 if j == 0 else 4}",
            )
            z = np.linspace(np.log(g.distance.min()), np.log(g.distance.max()), 150)
            ax.plot(
                z,
                beta[2 * j] + beta[2 * j + 1] * z + beta[-1] * z * z,
                color=color,
                lw=1.5,
            )
        lo, hi = fit["interval"]
        ax.set_title(
            f"{title}\nc = {fit['c']:.3f}  (95% CI {lo:.3f}–{hi:.3f})", fontsize=11
        )
        ax.set_xlabel("log of physical separation")
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    axes[0].set_ylabel("Disorder covariance D")
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "output/analysis")
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help="0 uses recorded intervals; 10000 reproduces the fixed resampling",
    )
    args = parser.parse_args()
    if args.bootstrap == 1 or args.bootstrap < 0:
        parser.error("--bootstrap must be zero or at least two")
    data, plan = load_data(args.data)
    recorded = json.loads((args.data / "analysis.json").read_text())
    points = point_table(data)
    saved = pd.read_csv(args.data / "intervals.csv")
    merged = points.merge(
        saved,
        on=["domain", "model", "size", "probe"],
        suffixes=("", "_saved"),
        validate="one_to_one",
    )
    for name in ("thermal", "disorder", "total"):
        if not np.allclose(merged[name], merged[name + "_saved"], rtol=0, atol=1e-10):
            raise ValueError(
                f"Recomputed {name} estimates disagree with released table"
            )
    fits = fit_coefficients(points, plan, recorded)
    table, contrasts, controls = (
        saved,
        recorded["contrasts"],
        recorded["control_curvature_test"],
    )
    if args.bootstrap:
        print(
            f"Recomputing {args.bootstrap:,} environment bootstrap draws and the fixed control test…",
            flush=True,
        )
        table, fits, contrasts, controls = rerun_uncertainty(data, plan, args.bootstrap)
    args.output.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(
        [
            {
                "domain": f["domain"],
                "model": f["model"],
                "minimum_size": f["window"][0],
                "maximum_size": f["window"][1],
                "c": f["c"],
                "ci_low": f["interval"][0],
                "ci_high": f["interval"][1],
            }
            for f in fits
        ]
    )
    summary.to_csv(args.output / "coefficients.csv", index=False)
    table.to_csv(args.output / "moments.csv", index=False)
    report = {
        "pairs": len(data) // 2,
        "fits": fits,
        "contrasts": contrasts,
        "control_curvature_test": controls,
        "uncertainty": "Recomputed from measurements"
        if args.bootstrap
        else "Recorded Hamilton intervals; point estimates recomputed",
        "bootstrap_repetitions": args.bootstrap,
    }
    (args.output / "analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_results(table, fits, args.output / "gamma_covariance.png")
    print(
        summary[(summary.model == "aztec_gamma") & (summary.minimum_size == 16)][
            ["domain", "c", "ci_low", "ci_high"]
        ].to_string(index=False)
    )
    print(f"Verified {len(data) // 2:,} paired environments. Output: {args.output}")


if __name__ == "__main__":
    main()
