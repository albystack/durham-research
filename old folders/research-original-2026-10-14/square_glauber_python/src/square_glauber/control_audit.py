"""Read-only pairing audit and joint all-one covariance permutation test.

Permutation moves a complete second-replica vector (centre and every probe)
within its size group. The test statistic is the largest absolute covariance
standardised by its exact permutation-null standard deviation. This preserves
within-replica spatial dependence, unlike permuting each observable separately.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .experiment import derive_seeds

SEEDS = ["master_seed", "environment_seed", "chain1_seed", "chain2_seed"]


def joint_covariance_test(
    groups: list[tuple[np.ndarray, np.ndarray]], *, permutations=9999, seed=2026092903
) -> dict:
    """Monte Carlo randomisation p-value with the standard plus-one correction."""
    if permutations < 1 or not groups:
        raise ValueError("need groups and at least one permutation")
    centered, observed, covariance = [], [], []
    for first, second in groups:
        first, second = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
        if (
            first.shape != second.shape
            or first.ndim != 2
            or first.shape[0] < 3
            or not np.isfinite(first).all()
            or not np.isfinite(second).all()
        ):
            raise ValueError(
                "finite paired matrices of equal shape with >=3 rows required"
            )
        n = len(first)
        a, b = first - first.mean(axis=0), second - second.mean(axis=0)
        # Var_pi[sum_i a_i b_pi(i)/(n-1)] = ||a||² ||b||²/(n-1)^3.
        sd = np.sqrt((a * a).sum(axis=0) * (b * b).sum(axis=0) / (n - 1) ** 3)
        cov = (a * b).sum(axis=0) / (n - 1)
        z = np.divide(cov, sd, out=np.zeros_like(cov), where=sd > 0)
        covariance.extend(cov.tolist())
        observed.extend(z.tolist())
        centered.append((a, b, sd))
    observed = np.abs(np.asarray(observed))
    rng = np.random.default_rng(seed)
    maxima = np.empty(permutations)
    for k in range(permutations):
        maximum = 0.0
        for a, b, sd in centered:
            cov = (a * b[rng.permutation(len(b))]).sum(axis=0) / (len(b) - 1)
            z = np.divide(cov, sd, out=np.zeros_like(cov), where=sd > 0)
            maximum = max(maximum, float(np.abs(z).max()))
        maxima[k] = maximum
    return dict(
        permutations=permutations,
        seed=seed,
        statistic=float(observed.max()),
        joint_p_value=float(
            (1 + (maxima >= observed.max()).sum()) / (permutations + 1)
        ),
        covariance=covariance,
        absolute_standardized_covariance=observed.tolist(),
        familywise_adjusted_p_values=[
            float((1 + (maxima >= v).sum()) / (permutations + 1)) for v in observed
        ],
        null_maximum_95_percentile=float(np.quantile(maxima, 0.95)),
    )


def audit_snapshot(input_dir: Path, *, permutations=9999, seed=2026092903) -> dict:
    center = pd.read_csv(
        input_dir / "center_measurements.csv", dtype={k: str for k in SEEDS}
    )
    spatial = pd.read_csv(
        input_dir / "spatial_measurements.csv", dtype={k: str for k in SEEDS}
    )
    if center["environment_id"].duplicated().any():
        raise ValueError("duplicate centre environment rows")
    if spatial.duplicated(["environment_id", "requested_rho"]).any():
        raise ValueError("duplicate spatial environment/probe rows")
    seen = set()
    for row in center.itertuples():
        expected = derive_seeds(
            int(row.master_seed), int(row.L), row.model, int(row.environment_index)
        )
        for key in SEEDS:
            if int(getattr(row, key)) != getattr(expected, key):
                raise ValueError(
                    f"seed derivation mismatch: {row.environment_id}/{key}"
                )
        for key in SEEDS[1:]:
            value = int(getattr(row, key))
            if value in seen:
                raise ValueError("seed reused across environments or replica streams")
            seen.add(value)
    keys = ["environment_id", *SEEDS, "environment_sha256", "L", "model"]
    joined = spatial.merge(
        center[keys],
        on="environment_id",
        suffixes=("", "_center"),
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    if not joined["_merge"].eq("both").all():
        raise ValueError("spatial row has no centre environment")
    for key in keys[1:]:
        if not joined[key].eq(joined[key + "_center"]).all():
            raise ValueError(f"spatial/centre pairing mismatch: {key}")
    groups, labels = [], []
    for L, group in center[center.model.eq("uniform")].groupby("L", sort=True):
        group = group.sort_values("environment_id").set_index("environment_id")
        rows = spatial[spatial.model.eq("uniform") & spatial.L.eq(L)]
        first = [group.H1_center.to_numpy()]
        second = [group.H2_center.to_numpy()]
        labels.append(dict(L=int(L), observable="center"))
        for rho, probe in rows.groupby("requested_rho", sort=True):
            probe = probe.set_index("environment_id").reindex(group.index)
            if (
                probe[["DeltaH1", "DeltaH2"]].isna().any().any()
                or probe.r.nunique() != 1
            ):
                raise ValueError("incomplete or inconsistent probe coverage")
            first.append(probe.DeltaH1.to_numpy())
            second.append(probe.DeltaH2.to_numpy())
            labels.append(
                dict(
                    L=int(L),
                    observable="spatial",
                    requested_rho=float(rho),
                    r=int(probe.r.iloc[0]),
                )
            )
        groups.append((np.column_stack(first), np.column_stack(second)))
    test = joint_covariance_test(groups, permutations=permutations, seed=seed)
    test["observables"] = labels
    summary_path = input_dir / "campaign_summary.json"
    campaign = json.loads(summary_path.read_text()) if summary_path.exists() else None
    return dict(
        source="September frozen snapshot; no new samples generated",
        input_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (
                input_dir / "center_measurements.csv",
                input_dir / "spatial_measurements.csv",
            )
        },
        environment_pairs=len(center),
        spatial_rows=len(spatial),
        seeds_verified=True,
        replica_seeds_globally_unique=True,
        joins_verified=True,
        note="Permutation tests pairing independence under the all-one model; it cannot validate single-replica Gibbs sampling. Historical state files were not required or certified by this table audit.",
        campaign_summary=campaign,
        joint_control_test=test,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=9999)
    args = parser.parse_args()
    report = audit_snapshot(args.input, permutations=args.permutations)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: report["joint_control_test"][k]
                for k in ("statistic", "joint_p_value", "permutations")
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
