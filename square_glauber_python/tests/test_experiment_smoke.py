import json

import pandas as pd

from square_glauber.analysis import analyze_pilot
from square_glauber.experiment import run_pilot


def test_uniform_and_gamma_pilot_and_analysis_smoke(tmp_path) -> None:
    pilot = tmp_path / "pilot"
    metadata = run_pilot(
        sizes=[4, 6, 8],
        models=["uniform", "gamma"],
        gamma_shape=1.0,
        environments=3,
        burnin_sweeps=4,
        measurement_gap_sweeps=2,
        num_measurements=2,
        master_seed=20260907,
        output=pilot,
    )
    assert metadata["environments_per_model_size"] == 3
    center = pd.read_csv(pilot / "center_measurements.csv")
    spatial = pd.read_csv(pilot / "spatial_measurements.csv")
    diagnostics = pd.read_csv(pilot / "chain_diagnostics.csv")
    assert len(center) == 2 * 3 * 3 * 2
    assert center["environment_id"].nunique() == 2 * 3 * 3
    required = {
        "master_seed",
        "environment_seed",
        "chain1_seed",
        "chain2_seed",
        "L",
        "model",
        "environment_id",
        "H1_center",
        "H2_center",
    }
    assert required <= set(center.columns)
    assert (center["chain1_seed"] != center["chain2_seed"]).all()
    assert {"requested_rho", "r", "DeltaH1", "DeltaH2"} <= set(spatial.columns)
    assert {
        "chain1_flippable_proposal_rate",
        "chain1_actual_change_rate",
        "chain1_null_moves",
    } <= set(diagnostics.columns)

    analysis = pilot / "analysis"
    report = analyze_pilot(
        pilot, analysis, bootstrap_replicates=10, bootstrap_seed=20260908
    )
    assert report["statistical_unit"].startswith("whole environment")
    summary = pd.read_csv(analysis / "center_summary.csv")
    row = summary[(summary["model"] == "gamma") & (summary["L"] == 4)].iloc[0]
    environment_means = pd.read_csv(analysis / "environment_center_means.csv")
    group = environment_means[
        (environment_means["model"] == "gamma") & (environment_means["L"] == 4)
    ]
    expected_connected = 0.5 * (group.H1_environment_mean - group.H2_environment_mean).var(ddof=1)
    expected_covariance = group.H1_environment_mean.cov(group.H2_environment_mean)
    assert row.connected_component == expected_connected
    assert row.disorder_covariance == expected_covariance
    assert (analysis / "center_covariance_fits.png").exists()
    assert json.loads((analysis / "analysis_metadata.json").read_text())["exploratory"]
