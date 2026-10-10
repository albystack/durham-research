"""Read certified historical states and add explicitly exploratory probes.

Only states whose environment hash and both seeds match the frozen centre
snapshot are accepted. Missing states are listed; they are never regenerated or
silently replaced. The September primary measurements are not modified.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .comparison import measure_probes, symmetric_probes
from .geometry import spatial_separations
from .height import face_heights
from .perfect_campaign import _load_certified_sample


def remeasure(snapshot: Path, state_root: Path, output: Path) -> dict:
    frame = pd.read_csv(
        snapshot,
        dtype={k: str for k in ("environment_seed", "chain1_seed", "chain2_seed")},
    )
    if frame.environment_id.duplicated().any():
        raise ValueError("duplicate snapshot environments")
    expected = frame.set_index("environment_id")
    found = {}
    for path in sorted(state_root.rglob("environment.json")):
        meta = json.loads(path.read_text())
        key = meta.get("environment_id")
        if key not in expected.index:
            continue
        row = expected.loc[key]
        if meta.get("environment_sha256") != row.environment_sha256 or any(
            int(meta.get(k, -1)) != int(row[k])
            for k in ("environment_seed", "chain1_seed", "chain2_seed")
        ):
            continue
        states = []
        for replica in (1, 2):
            loaded = _load_certified_sample(
                path.parent,
                sample_index=replica,
                seed=int(row[f"chain{replica}_seed"]),
                L=int(row.L),
                environment_sha256=row.environment_sha256,
                separations=spatial_separations(int(row.L)),
            )
            if loaded is None:
                break
            states.append(loaded[0])
        if len(states) != 2:
            continue
        # The snapshot itself provides another independent join check.
        heights = [face_heights(state) for state in states]
        c = (int(row.L) - 2) // 2
        if any(
            int(h[c, c]) != int(row[f"H{k}_center"]) for k, h in enumerate(heights, 1)
        ):
            raise ValueError(f"state does not reproduce snapshot central height: {key}")
        probes = symmetric_probes(int(row.L))
        values = [measure_probes(h, probes) for h in heights]
        rows = [
            dict(
                environment_id=key,
                L=int(row.L),
                model=row.model,
                provenance="exploratory_remeasurement_of_certified_September_states",
                probe=p.label,
                r=p.r,
                distance=p.distance,
                start_row=p.start[0],
                start_column=p.start[1],
                end_row=p.end[0],
                end_column=p.end[1],
                value1=values[0][p.label],
                value2=values[1][p.label],
            )
            for p in probes
        ]
        if key in found and found[key] != rows:
            raise ValueError(f"conflicting certified state copies: {key}")
        found[key] = rows
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([r for rows in found.values() for r in rows]).to_csv(
        output / "exploratory_measurements.csv", index=False
    )
    missing = frame[~frame.environment_id.isin(found)]
    missing[
        [
            "environment_id",
            "L",
            "model",
            "environment_sha256",
            "chain1_seed",
            "chain2_seed",
        ]
    ].to_csv(output / "missing_states.csv", index=False)
    report = dict(
        requested_snapshot_pairs=len(frame),
        available_certified_pairs=len(found),
        missing_pairs=len(missing),
        complete=len(missing) == 0,
        note="Availability is not a sampling design. Do not fit the available subset as a complete campaign.",
    )
    (output / "availability.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snapshot", type=Path, required=True)
    p.add_argument("--states", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(remeasure(a.snapshot, a.states, a.output), indent=2))


if __name__ == "__main__":
    main()
