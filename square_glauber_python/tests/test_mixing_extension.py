import json

import numpy as np
import pandas as pd
import pytest

from square_glauber.environment import DimerEnvironment
from square_glauber.mixing_campaign import run_mixing_campaign
from square_glauber.mixing_extension import (
    array_sha256,
    environment_hardness,
    environment_hashes,
    run_frozen_gamma_extension,
)


def test_environment_hash_and_hardness_are_deterministic() -> None:
    horizontal = np.asarray(
        [[1.0, 2.0, 4.0], [0.5, 1.5, 3.0], [2.5, 0.25, 1.25], [0.75, 5.0, 2.25]]
    )
    vertical = np.asarray(
        [[1.25, 0.5, 2.0, 1.0], [4.0, 1.0, 0.75, 2.5], [0.25, 3.0, 1.5, 0.5]]
    )
    environment = DimerEnvironment(4, horizontal, vertical)
    duplicate = DimerEnvironment(4, horizontal.copy(), vertical.copy())
    assert environment_hashes(environment) == environment_hashes(duplicate)
    assert array_sha256(horizontal) != array_sha256(vertical)
    metrics = environment_hardness(environment)
    all_weights = np.concatenate((horizontal.ravel(), vertical.ravel()))
    assert metrics["edge_count"] == 24
    assert metrics["face_count"] == 9
    assert metrics["edge_weight_min"] == float(all_weights.min())
    assert metrics["edge_weight_max"] == float(all_weights.max())
    assert np.isclose(metrics["edge_weight_mean"], all_weights.mean())
    assert 0 <= metrics["fraction_abs_lambda_gt_2"] <= 1


def test_frozen_extension_reuses_weights_seeds_and_trace_prefix(tmp_path) -> None:
    source = tmp_path / "source"
    output = tmp_path / "extended"
    run_mixing_campaign(
        sizes=(4,),
        gamma_environments=3,
        uniform_replicates=1,
        gamma_shape=1.0,
        sweeps=8,
        record_every=2,
        master_seed=20260908,
        output=source,
        max_lag_records=3,
    )
    report = run_frozen_gamma_extension(
        source=source,
        output=output,
        sweeps_by_size={4: 16},
        record_every=2,
        max_lag_records=4,
    )
    assert report["environment_count"] == 3
    assert report["all_seed_and_array_verifications_passed"]
    assert report["all_20k_trace_prefixes_match_exactly"]
    cases = pd.read_csv(output / "case_summary.csv")
    assert cases["source_trace_prefix_matches_exactly"].all()
    assert set(cases["sweeps"]) == {16}
    assert (output / "environment_hardness.csv").exists()
    assert (output / "old_new_comparison.csv").exists()
    for environment_id in cases["environment_id"]:
        old_summary = json.loads((source / environment_id / "summary.json").read_text())
        new_summary = json.loads((output / environment_id / "summary.json").read_text())
        assert old_summary["environment_seed"] == new_summary["environment_seed"]
        verification = json.loads(
            (output / environment_id / "environment_verification.json").read_text()
        )
        assert verification["combined_hashes_match"]
        assert (output / environment_id / "window_summary.csv").exists()

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run_frozen_gamma_extension(
            source=source,
            output=output,
            sweeps_by_size={4: 20},
            record_every=2,
        )
