import numpy as np
import pandas as pd
import pytest

from square_glauber.followup_analysis import paired_moments, bootstrap_table, curvature


def test_weighted_resampling_matches_explicit_environment_duplication():
    a = np.array([[1, 8], [3, 4], [-2, 7], [9, 2.0]])
    b = np.array([[4, 2], [-1, 8], [6, 3], [2, 5.0]])
    counts = np.array([[2, 0, 1, 1], [0, 3, 0, 1]])
    actual = paired_moments(a, b, counts)
    for k, row in enumerate(counts):
        idx = np.repeat(np.arange(4), row)
        expected_d = np.array(
            [np.cov(a[idx, j], b[idx, j], ddof=1)[0, 1] for j in range(2)]
        )
        expected_t = np.var(a[idx] - b[idx], axis=0, ddof=1) / 2
        assert np.allclose(actual[k, :, 0], expected_t)
        assert np.allclose(actual[k, :, 1], expected_d)
        assert np.allclose(actual[k, :, 2], expected_t + expected_d)


def test_bootstrap_moves_all_probes_and_replicas_together():
    rows = [
        dict(
            model="x",
            size=8,
            environment_id=str(i),
            probe=p,
            distance=2,
            value1=i * i,
            value2=3 * i + 1,
        )
        for i in range(6)
        for p in ("a", "b")
    ]
    table, draws = bootstrap_table(pd.DataFrame(rows), 100, 12)
    assert np.array_equal(draws[:, 0], draws[:, 1])
    assert table.disorder.iloc[0] == table.disorder.iloc[1]
    with pytest.raises(ValueError, match="partial"):
        bootstrap_table(pd.DataFrame(rows[:-1]), 100, 12)


def test_curvature_refuses_rounded_distance_pseudoreplication():
    table = pd.DataFrame(
        dict(
            model=["x"] * 5,
            size=[8, 12, 16, 20, 24],
            probe=["p"] * 5,
            distance=[1, 1, 2, 2, 3],
            disorder=[0] * 5,
        )
    )
    with pytest.raises(ValueError, match="distinct distances"):
        curvature(table, np.zeros((10, 5)), "x", ["p"], [8, 24])


def test_curvature_refuses_partial_probe_coverage():
    table = (
        pd.DataFrame(
            [
                dict(model="x", size=size, probe=probe, distance=size / 4, disorder=0)
                for size in (8, 12, 16, 20, 24)
                for probe in ("p", "q")
            ]
        )
        .iloc[1:]
        .reset_index(drop=True)
    )
    with pytest.raises(ValueError, match="incomplete probe coverage"):
        curvature(table, np.zeros((10, len(table))), "x", ["p", "q"], [8, 24])
