"""Scientific-contract tests for the new law and physical observables."""

import numpy as np
import pytest
from square_glauber.accelerated import NUMBA_AVAILABLE
from square_glauber.comparison import (
    ComparisonModel,
    aztec_face_xy,
    comparison_environment,
    comparison_seeds,
    family_masks,
    measure_probes,
    symmetric_probes,
)
from square_glauber.comparison_pilot import run_pilot
from square_glauber.exact import validate_exact_kernel
from square_glauber.height import face_heights
from square_glauber.perfect import MonotoneCFTPSampler


def test_family_placement_has_paper_order_one_probability():
    env = comparison_environment(
        2, np.random.default_rng(1), ComparisonModel("aztec_mean")
    )
    # The occupied top/bottom pair has yellow b; left/right pair has red a.
    np.testing.assert_array_equal(env.horizontal_weights, [[1.0], [0.25]])
    np.testing.assert_array_equal(env.vertical_weights, [[0.2, 1.0]])
    assert env.horizontal_probability(0, 0) == pytest.approx(0.25 / 0.45)


@pytest.mark.parametrize("L", [2, 4])
@pytest.mark.parametrize("name", ["aztec_gamma", "aztec_mean", "uniform"])
def test_matched_law_exact_detailed_balance(L, name):
    env = comparison_environment(L, np.random.default_rng(297), ComparisonModel(name))
    report = validate_exact_kernel(L, L, env.horizontal_weights, env.vertical_weights)
    assert report.stationarity_residual < 1e-12
    assert report.detailed_balance_residual < 1e-12
    assert report.illegal_transition_count == 0
    h_mask, v_mask = family_masks(L)
    assert np.all(env.horizontal_weights[~h_mask] == 1)
    assert np.all(env.vertical_weights[~v_mask] == 1)


def test_family_normalisation_is_not_silently_gauge_equivalent():
    original = comparison_environment(
        2, np.random.default_rng(1), ComparisonModel("aztec_mean")
    )
    assert original.horizontal_probability(0, 0) != pytest.approx(0.5)


def test_physical_probe_coordinates_and_height_offset_invariance():
    for L in (4, 8, 12, 32):
        probes = symmetric_probes(L)
        assert len({(p.axis, p.r) for p in probes}) == len(probes)
        for p in probes:
            assert p.end[1] - p.start[1] == p.r
            assert p.end[0] - p.start[0] == (p.r if p.axis == "diagonal" else 0)
            assert 0 <= p.start[0] <= p.end[0] < L - 1
        heights = np.arange((L - 1) ** 2).reshape((L - 1, L - 1))
        assert measure_probes(heights, probes) == measure_probes(heights + 100, probes)
    for n in (4, 8, 16):
        for r in range(1, n):
            left = (n - r) // 2 + 1
            x1, y1 = aztec_face_xy(n, n + 1, left)
            x2, y2 = aztec_face_xy(n, n + 1, left + r)
            assert (x2 - x1, y2 - y1) == (r, -r)


def test_model_parameter_seed_namespaces_do_not_overlap():
    settings = [
        ComparisonModel(),
        ComparisonModel(alpha=0.3),
        ComparisonModel("uniform"),
    ]
    seeds = [
        s for m in settings for i in range(8) for s in comparison_seeds(29, 4, m, i)
    ]
    assert len(seeds) == len(set(seeds))


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="optional Numba backend unavailable")
def test_matched_cftp_reference_matches_accelerated():
    env = comparison_environment(4, np.random.default_rng(91), ComparisonModel())
    first = MonotoneCFTPSampler(env, backend="reference").sample(
        92, max_horizon_sweeps=4096
    )
    second = MonotoneCFTPSampler(env, backend="numba").sample(
        92, max_horizon_sweeps=4096
    )
    np.testing.assert_array_equal(face_heights(first.state), face_heights(second.state))
    assert first.horizon_sweeps == second.horizon_sweeps


def test_pilot_censor_restart_preserves_identity_and_both_streams(tmp_path):
    config = dict(
        sizes=[4],
        environments_per_size=2,
        max_horizon_sweeps=1,
        master_seed=19,
        models=[dict(name="uniform")],
    )
    initial = run_pilot(config, tmp_path, backend="reference")
    assert initial["requested_environments"] == 2
    assert not initial["all_pairs_certified"]
    config["max_horizon_sweeps"] = 1024
    final = run_pilot(config, tmp_path, backend="reference")
    assert final["all_pairs_certified"]
    old = (tmp_path / "measurements.csv").read_bytes()
    assert run_pilot(config, tmp_path, backend="reference")["all_pairs_certified"]
    assert (tmp_path / "measurements.csv").read_bytes() == old
    config["master_seed"] += 1
    with pytest.raises(ValueError, match="identity"):
        run_pilot(config, tmp_path, backend="reference")
