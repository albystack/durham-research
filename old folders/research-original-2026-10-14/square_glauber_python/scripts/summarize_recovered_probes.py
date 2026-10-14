"""Descriptive paired moments for recovered September probes; no scaling fit.

Run with INPUT.csv OUTPUT.csv. Values retain the input height convention.
These are exploratory measurements of the historical all-edge Gamma(1) and
uniform laws, not the new two-family matched law. The original unresolved
environment remains missing, so these moments describe the recovered sample.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def summarize(source: Path, destination: Path) -> None:
    frame = pd.read_csv(source)
    if frame.duplicated(["environment_id", "probe"]).any():
        raise ValueError("duplicate environment/probe rows")
    rows = []
    for key, group in frame.groupby(["model", "L", "probe"], sort=True):
        if len(group) < 2:
            raise ValueError("at least two paired environments required")
        if group[["r", "distance"]].drop_duplicates().shape[0] != 1:
            raise ValueError("probe geometry varies within a group")
        x = group.value1.to_numpy(dtype=float)
        y = group.value2.to_numpy(dtype=float)
        if not (np.isfinite(x).all() and np.isfinite(y).all()):
            raise ValueError("non-finite measurements")
        thermal = float(np.var(x - y, ddof=1) / 2)
        disorder = float(np.cov(x, y, ddof=1)[0, 1])
        total = float((np.var(x, ddof=1) + np.var(y, ddof=1)) / 2)
        if not np.isclose(thermal + disorder, total):
            raise AssertionError("paired variance identity failed")
        rows.append(
            dict(
                model=key[0],
                L=int(key[1]),
                probe=key[2],
                pairs=len(group),
                r=int(group.r.iloc[0]),
                distance=float(group.distance.iloc[0]),
                thermal=thermal,
                disorder=disorder,
                total=total,
                note="Exploratory recovered sample; no uncertainty interval or scaling claim",
            )
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(destination, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    summarize(args.source, args.destination)
