"""Bounded diagnostics using unchanged pilot environments and random maps.

Regenerates weights from recorded seeds (no new disorder draw), checks pilot
seed identities, and replays the first four L=8 random cases from true extrema.
The cap is 2^18 sweeps; this never edits the original pilot or its checkpoints.
Endpoint snapshots are diagnostics, NOT certified samples unless bounds agree.
"""

import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

from square_glauber.comparison import (
    ComparisonModel,
    comparison_environment,
    comparison_seeds,
)
from square_glauber.glauber import UpdateCounts
from square_glauber.height import face_heights
from square_glauber.perfect import MonotoneCFTPSampler, _random_block

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/validation/followup_20260930"


def main():
    plan = json.loads((OUT / "analysis_plan.json").read_text())["square_diagnostic"]
    old = pd.read_csv(
        ROOT
        / "square_glauber_python/outputs/hamilton_matched_20260930/feasibility/stream_status.csv",
        dtype={"environment_seed": str, "chain_seed": str, "cftp_seed": str},
    )
    rows, histories, maps = [], [], []
    for name in ("aztec_gamma", "aztec_mean", "uniform"):
        model = ComparisonModel(name=name)
        for L in (8, 12, 16, 24, 32):
            for index in range(16):
                seeds = comparison_seeds(plan["master_seed"], L, model, index)
                env = comparison_environment(L, np.random.default_rng(seeds[0]), model)
                group = old[
                    (old.model == name)
                    & (old["size"] == L)
                    & old.environment_id.str.endswith(f"env{index:06d}")
                ]
                assert len(group) == 2 and all(
                    int(v) == seeds[0] for v in group.environment_seed
                )
                for replica in (1, 2):
                    assert (
                        int(group[group.replica == replica].cftp_seed.iloc[0])
                        == seeds[replica]
                    )
                # Stable log odds and minority probability, even if Float64 pH rounds to one.
                logh = np.log(env.horizontal_weights[:-1]) + np.log(
                    env.horizontal_weights[1:]
                )
                logv = np.log(env.vertical_weights[:, :-1]) + np.log(
                    env.vertical_weights[:, 1:]
                )
                odds = logh - logv
                logminor = -np.logaddexp(0, np.abs(odds))
                sampler = MonotoneCFTPSampler(env)
                p = sampler.horizontal_probabilities
                rows.append(
                    dict(
                        model=name,
                        L=L,
                        index=index,
                        environment_seed=str(seeds[0]),
                        min_weight=float(
                            min(
                                env.horizontal_weights.min(), env.vertical_weights.min()
                            )
                        ),
                        max_abs_log_odds=float(abs(odds).max()),
                        min_minority_probability=float(np.exp(logminor.min())),
                        rounded_zero=int((p == 0).sum()),
                        rounded_one=int((p == 1).sum()),
                        faces_minority_below_inverse_cap=int(
                            (logminor < -np.log(65536)).sum()
                        ),
                        faces=p.size,
                    )
                )
                if (
                    name != "aztec_gamma"
                    or L != plan["local_replay_L"]
                    or index not in plan["environment_indices"]
                ):
                    continue
                for replica in plan["replicas"]:
                    seed = seeds[replica]
                    for horizon in plan["horizons"]:
                        blocks = [
                            _random_block(seed, L, i, 1 if i == 0 else 2 ** (i - 1))
                            for i in reversed(range(horizon.bit_length()))
                        ]
                        lo, hi = (
                            sampler.minimum.state.copy(),
                            sampler.maximum.state.copy(),
                        )
                        lc, hc = UpdateCounts(), UpdateCounts()
                        start = perf_counter()
                        sampler._apply_random_blocks(
                            lo, hi, blocks, lc, hc, verify_order=True
                        )
                        gap = face_heights(hi) - face_heights(lo)
                        assert np.all(gap >= 0)
                        if horizon == 1:
                            reference = MonotoneCFTPSampler(env, backend="reference")
                            rlo, rhi = (
                                reference.minimum.state.copy(),
                                reference.maximum.state.copy(),
                            )
                            reference._apply_random_blocks(
                                rlo,
                                rhi,
                                blocks,
                                UpdateCounts(),
                                UpdateCounts(),
                                verify_order=True,
                            )
                            assert np.array_equal(
                                face_heights(lo), face_heights(rlo)
                            ) and np.array_equal(face_heights(hi), face_heights(rhi))
                        histories.append(
                            dict(
                                index=index,
                                replica=replica,
                                chain_seed=str(seed),
                                horizon=horizon,
                                coalesced=bool(not gap.any()),
                                max_gap=int(gap.max()),
                                disagreeing_faces=int((gap > 0).sum()),
                                sum_gap=int(gap.sum()),
                                seconds=perf_counter() - start,
                                lower_flippable=lc.flippable_proposals,
                                lower_changes=lc.actual_changes,
                                upper_flippable=hc.flippable_proposals,
                                upper_changes=hc.actual_changes,
                                attempts_per_bound=lc.attempted_updates,
                            )
                        )
                        if horizon == plan["horizons"][-1]:
                            for i, j in np.ndindex(gap.shape):
                                maps.append(
                                    dict(
                                        index=index,
                                        replica=replica,
                                        row=i,
                                        column=j,
                                        gap=int(gap[i, j]),
                                        horizontal_probability=float(p[i, j]),
                                        minority_probability=float(
                                            np.exp(logminor[i, j])
                                        ),
                                    )
                                )
    pd.DataFrame(rows).to_csv(OUT / "weight_diagnostics.csv", index=False)
    pd.DataFrame(histories).to_csv(OUT / "cftp_replay.csv", index=False)
    pd.DataFrame(maps).to_csv(OUT / "cftp_endpoint_maps.csv", index=False)
    print("Saved bounded diagnostics; no original checkpoints changed.", flush=True)


if __name__ == "__main__":
    main()
