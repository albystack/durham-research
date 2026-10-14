import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from square_glauber.analysis import analyze_pilot
from square_glauber.perfect_campaign import (
    aggregate_campaign,
    create_campaign_manifest,
    file_sha256,
    run_campaign_task,
)


def _write_schedule(path) -> None:
    pd.DataFrame(
        [
            {
                "L": 4,
                "model": "gamma",
                "gamma_shape": 1.0,
                "environments": 2,
                "environments_per_task": 1,
                "max_horizon_sweeps": 16384,
            },
            {
                "L": 4,
                "model": "uniform",
                "gamma_shape": np.nan,
                "environments": 1,
                "environments_per_task": 1,
                "max_horizon_sweeps": 16384,
            },
        ]
    ).to_csv(path, index=False)


def test_manifest_task_resume_and_complete_aggregation(tmp_path) -> None:
    schedule = tmp_path / "schedule.csv"
    manifest = tmp_path / "manifest.csv"
    output = tmp_path / "campaign"
    _write_schedule(schedule)
    frame = create_campaign_manifest(
        schedule,
        manifest,
        campaign_id="unit-perfect-pairs",
        master_seed=99017,
        backend="reference",
    )
    assert len(frame) == 3
    assert frame["environment_count"].sum() == 3
    assert frame["task_id"].tolist() == [0, 1, 2]
    assert json.loads(manifest.with_suffix(".json").read_text())["perfect_sample_count"] == 6

    for task_id in frame["task_id"]:
        report = run_campaign_task(manifest, int(task_id), output)
        assert report["complete_pair_count"] == 1
        assert report["incomplete_pair_count"] == 0

    environment = output / "cases/gamma/L004/environment-000000/environment_weights.npz"
    pair = output / "cases/gamma/L004/environment-000000/pair.json"
    environment_hash_before = file_sha256(environment)
    pair_before = json.loads(pair.read_text())
    resumed = run_campaign_task(manifest, 0, output)
    pair_after = json.loads(pair.read_text())
    assert resumed["complete_pair_count"] == 1
    assert file_sha256(environment) == environment_hash_before
    assert pair_after == pair_before

    aggregate = aggregate_campaign(manifest, output)
    assert aggregate["expected_environments"] == 3
    assert aggregate["complete_environments"] == 3
    assert aggregate["all_pairs_cftp_certified"]
    center = pd.read_csv(output / "center_measurements.csv")
    spatial = pd.read_csv(output / "spatial_measurements.csv")
    diagnostics = pd.read_csv(output / "chain_diagnostics.csv")
    assert len(center) == 3
    assert len(spatial) == 3
    assert len(diagnostics) == 3
    assert center["environment_id"].nunique() == 3
    assert (center["chain1_seed"] != center["chain2_seed"]).all()
    assert center["environment_sha256"].str.len().eq(64).all()
    assert diagnostics["chain1_cftp_horizon_sweeps"].gt(0).all()


def test_incomplete_aggregation_retains_status_without_scientific_pair(tmp_path) -> None:
    schedule = tmp_path / "schedule.csv"
    pd.DataFrame(
        [
            {
                "L": 4,
                "model": "gamma",
                "gamma_shape": 1.0,
                "environments": 1,
                "environments_per_task": 1,
                "max_horizon_sweeps": 1,
            }
        ]
    ).to_csv(schedule, index=False)
    manifest = tmp_path / "manifest.csv"
    create_campaign_manifest(
        schedule,
        manifest,
        campaign_id="unit-censoring",
        master_seed=772244,
        backend="reference",
    )
    report = run_campaign_task(manifest, 0, tmp_path / "campaign")
    assert report["complete_pair_count"] == 0
    assert report["incomplete_pair_count"] == 1
    aggregate = aggregate_campaign(
        manifest, tmp_path / "campaign", require_complete=False
    )
    assert aggregate["complete_environments"] == 0
    assert aggregate["censored_environments"] == 1
    assert aggregate["center_rows"] == 0
    assert aggregate["incomplete_task_count"] == 1
    assert (tmp_path / "campaign/incomplete_slurm_array_spec.txt").read_text() == "0\n"
    with pytest.raises(ValueError, match="incomplete perfect-sampling campaign"):
        analyze_pilot(tmp_path / "campaign", tmp_path / "analysis")


def test_frozen_production_manifest_has_unique_complete_coverage() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest_path = (
        root / "configs/perfect_covariance_campaign_20260909_launch_manifest.csv"
    )
    metadata = json.loads(manifest_path.with_suffix(".json").read_text())
    frame = pd.read_csv(manifest_path)
    assert len(frame) == 934
    assert frame["environment_count"].sum() == 23_800
    assert metadata["perfect_sample_count"] == 47_600
    assert metadata["master_seed"] == 2_026_090_901
    assert metadata["manifest_sha256"] == file_sha256(manifest_path)
    assert frame["task_id"].tolist() == list(range(934))
    expected: set[tuple[str, int, int]] = set()
    for row in frame.itertuples(index=False):
        assert row.environment_stop - row.environment_start == row.environment_count
        for index in range(row.environment_start, row.environment_stop):
            key = (row.model, row.L, index)
            assert key not in expected
            expected.add(key)
    assert len(expected) == 23_800
