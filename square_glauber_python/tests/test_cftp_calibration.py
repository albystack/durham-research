import json

import pandas as pd
import pytest

from square_glauber.cftp_calibration import (
    consolidate_cftp_cost_calibration,
    extend_censored_cftp_calibration,
    run_cftp_cost_calibration,
)


def test_tiny_reference_cost_calibration_is_frozen_and_auditable(tmp_path) -> None:
    output = tmp_path / "calibration"
    report = run_cftp_cost_calibration(
        output,
        standard_sizes=(4,),
        standard_gamma_environments=2,
        standard_uniform_environments=1,
        standard_samples_per_environment=1,
        probe_sizes=(),
        probe_gamma_environments=0,
        probe_uniform_environments=0,
        probe_samples_per_environment=0,
        master_seed=7701,
        max_horizon_sweeps=16_384,
        backend="reference",
    )
    assert report["environment_count"] == 3
    assert report["cftp_stream_count"] == 3
    assert report["certified_count"] == 3
    assert report["censored_count"] == 0
    assert report["all_environment_reloads_and_regenerations_match"]
    assert report["all_cftp_seeds_unique"]
    assert report["all_environment_arrays_retained"]
    assert not report["scientific_covariance_computed"]
    assert not report["log_or_log_squared_fit_performed"]

    environments = pd.read_csv(output / "environment_catalog.csv")
    samples = pd.read_csv(output / "cftp_samples.csv")
    assert len(environments) == 3
    assert len(samples) == 3
    assert environments["environment_sha256"].str.len().eq(64).all()
    assert samples["certified"].all()
    assert (output / "censored_samples_to_extend.csv").is_file()
    saved = json.loads((output / "calibration_summary.json").read_text())
    assert saved["cftp_stream_count"] == 3

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run_cftp_cost_calibration(
            output,
            standard_sizes=(4,),
            standard_gamma_environments=1,
            standard_uniform_environments=0,
            standard_samples_per_environment=1,
            probe_sizes=(),
            backend="reference",
        )


def test_censored_stream_extension_replays_the_old_history(tmp_path) -> None:
    source = tmp_path / "source"
    initial = run_cftp_cost_calibration(
        source,
        standard_sizes=(4,),
        standard_gamma_environments=1,
        standard_uniform_environments=0,
        standard_samples_per_environment=1,
        probe_sizes=(),
        probe_gamma_environments=0,
        probe_uniform_environments=0,
        probe_samples_per_environment=0,
        master_seed=8801,
        max_horizon_sweeps=1,
        backend="reference",
    )
    assert initial["censored_count"] == 1
    output = tmp_path / "extended"
    extended = extend_censored_cftp_calibration(
        source,
        output,
        max_horizon_sweeps=16_384,
        backend="reference",
    )
    assert extended["streams_extended"] == 1
    assert extended["certified_after_extension"] == 1
    assert extended["still_censored"] == 0
    assert extended["all_old_cap_gaps_replayed_exactly"]
    consolidated = consolidate_cftp_cost_calibration(
        source, output, tmp_path / "consolidated"
    )
    assert consolidated["all_streams_certified"]
    assert consolidated["initially_censored"] == 1
    assert consolidated["remaining_censored"] == 0
