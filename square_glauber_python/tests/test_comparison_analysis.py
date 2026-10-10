import numpy as np
import pandas as pd
import pytest
from square_glauber.comparison_analysis import compare_curvature, prepare_blocks
from square_glauber.control_audit import joint_covariance_test


def test_joint_control_test_detects_deliberately_shared_replica_stream():
    rng = np.random.default_rng(810)
    values = rng.normal(size=(80, 4))
    test = joint_covariance_test([(values, values)], permutations=199, seed=20)
    assert test["joint_p_value"] == 1 / 200
    assert max(test["familywise_adjusted_p_values"]) == 1 / 200


def synthetic(domain):
    rng = np.random.default_rng(24)
    rows = []
    for size in (16, 32, 64, 128, 256):
        for env in range(40):
            shared = rng.normal() * np.log(size)
            for probe, rho in (("p", 0.125), ("q", 0.25)):
                rows.append(
                    dict(
                        domain=domain,
                        size=size,
                        environment_id=str(env),
                        probe=probe,
                        distance=size * rho,
                        law_id="same",
                        axis="diagonal",
                        height_scale=1,
                        value1=shared + rng.normal(),
                        value2=shared + rng.normal(),
                    )
                )
    return pd.DataFrame(rows)


def test_direct_contrast_and_block_bootstrap():
    aztec, square = synthetic("aztec"), synthetic("square")
    result = compare_curvature(
        aztec,
        square,
        probes=["p", "q"],
        aztec_window=(16, 256),
        square_window=(16, 256),
        bootstrap=100,
        seed=290,
    )
    assert result["delta_c"] == pytest.approx(0)
    assert result["delta_interval"][0] < 0 < result["delta_interval"][1]
    assert result["aztec_c"] > 0
    # Independent resampling, rather than sharing campaign bootstrap indices.
    assert result["delta_interval"][1] - result["delta_interval"][0] > 0


def test_incomplete_vectors_and_collapsed_distances_refused():
    data = synthetic("square")
    with pytest.raises(ValueError, match="complete"):
        prepare_blocks(data.iloc[1:], ["p", "q"], 16, 256)
    data["distance"] = 1
    with pytest.raises(ValueError, match="rank"):
        prepare_blocks(data, ["p", "q"], 16, 256)


def test_unmatched_laws_refused():
    a, b = synthetic("aztec"), synthetic("square")
    b["law_id"] = "other"
    with pytest.raises(ValueError, match="same single weight law"):
        compare_curvature(
            a,
            b,
            probes=["p", "q"],
            aztec_window=(16, 256),
            square_window=(16, 256),
            bootstrap=10,
        )
