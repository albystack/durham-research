import numpy as np

from square_glauber.diagnostics import (
    autocorrelation,
    compare_extremal_starts,
    consecutive_window_summaries,
    integrated_autocorrelation_time,
    split_rhat,
)
from square_glauber.environment import sample_environment
from square_glauber.geometry import spatial_separations


def test_autocorrelation_and_iact_are_finite() -> None:
    values = [0, 1, 0, -1] * 10
    result = autocorrelation(values, max_lag=5)
    assert result.shape == (6,)
    assert result[0] == 1
    assert integrated_autocorrelation_time(values) >= 1
    assert integrated_autocorrelation_time([3, 3, 3]) == 1
    assert split_rhat([values, values]) < 1.05


def test_fft_autocorrelation_matches_direct_definition() -> None:
    values = np.asarray([2.0, -1.0, 4.0, 0.5, 3.0])
    centered = values - values.mean()
    expected = np.correlate(centered, centered, mode="full")[len(values) - 1 :]
    expected /= centered @ centered
    np.testing.assert_allclose(autocorrelation(values), expected, atol=1e-14)


def test_post_transient_consecutive_windows_cover_tail() -> None:
    rows = consecutive_window_summaries(
        list(range(9)), list(range(0, 90, 10)), post_transient_fraction=0.5, window_count=4
    )
    assert [row["sweep_start"] for row in rows] == [40, 60, 70, 80]
    assert sum(int(row["points"]) for row in rows) == 5
    assert rows[-1]["mean"] == 8.0
    assert rows[-1]["minimum"] == 8.0
    assert rows[-1]["maximum"] == 8.0
    assert rows[-1]["median"] == 8.0


def test_mixing_comparison_uses_distinct_starts() -> None:
    environment = sample_environment(4, np.random.default_rng(401), "uniform")
    result = compare_extremal_starts(
        environment,
        np.random.default_rng(402),
        np.random.default_rng(403),
        sweeps=20,
        record_every=4,
        separations=spatial_separations(4),
    )
    assert result.horizontal_start.start == "all_horizontal"
    assert result.vertical_start.start == "all_vertical"
    assert result.horizontal_start.sweeps == [0, 4, 8, 12, 16, 20]
    assert len(result.horizontal_start.height_sums) == 6
    summary = result.summary()
    assert "warning" in summary
    assert {"center_height", "height_sum"} <= set(summary["observables"])
    assert "split_rhat_height_sum_late_half" in summary
    center = summary["observables"]["center_height"]
    assert "z_gap_horizontal_minus_vertical" in center
    assert "late_mean_mcse" in center["horizontal_start"]
    assert len(center["horizontal_start"]["post_transient_windows"]) == 3
