"""Local, deliberately non-conclusive mixing diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from numpy.typing import NDArray

from .environment import DimerEnvironment
from .geometry import SpatialSeparation, center_face_index
from .glauber import UpdateCounts, run_attempts
from .height import face_heights, spatial_increments
from .matching import DimerState


def autocorrelation(values: list[int] | NDArray[np.float64], max_lag: int | None = None) -> NDArray[np.float64]:
    """Biased-normalization autocorrelation computed by an exact FFT convolution."""
    data = np.asarray(values, dtype=np.float64)
    if data.ndim != 1 or data.size == 0:
        raise ValueError("autocorrelation requires a nonempty one-dimensional trace")
    maximum = data.size - 1 if max_lag is None else min(int(max_lag), data.size - 1)
    if maximum < 0:
        raise ValueError("max_lag must be nonnegative")
    centered = data - data.mean()
    variance_sum = float(centered @ centered)
    if variance_sum == 0:
        result = np.zeros(maximum + 1, dtype=np.float64)
        result[0] = 1.0
        return result
    # Zero padding to at least 2*n-1 makes the circular FFT convolution equal
    # to the desired linear autocovariance sums.  This has the same biased
    # normalization as ``np.correlate`` but remains practical for long traces.
    transform_length = 1 << (2 * data.size - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=transform_length)
    correlations = np.fft.irfft(
        spectrum * np.conjugate(spectrum), n=transform_length
    )[: maximum + 1]
    return np.asarray(correlations / variance_sum, dtype=np.float64)


def integrated_autocorrelation_time(values: list[int] | NDArray[np.float64]) -> float:
    """Geyer-style initial-positive-pair IACT estimate, lower bounded by one."""
    data = np.asarray(values, dtype=np.float64)
    if data.size < 2 or np.var(data) == 0:
        return 1.0
    acf = autocorrelation(data, max_lag=data.size // 2)
    positive_sum = 0.0
    for first_lag in range(1, len(acf), 2):
        pair = float(acf[first_lag])
        if first_lag + 1 < len(acf):
            pair += float(acf[first_lag + 1])
        if not np.isfinite(pair) or pair <= 0:
            break
        positive_sum += pair
    return max(1.0, 1.0 + 2.0 * positive_sum)


def split_rhat(
    chains: list[list[int] | NDArray[np.float64]], *, tail_fraction: float = 0.5
) -> float:
    """Basic split-R-hat across chain tails, without rank normalization."""
    if len(chains) < 2:
        raise ValueError("split-R-hat requires at least two chains")
    if not 0 < tail_fraction <= 1:
        raise ValueError("tail_fraction must lie in (0,1]")
    arrays = [np.asarray(chain, dtype=np.float64) for chain in chains]
    tail_length = min(max(1, int(np.ceil(len(chain) * tail_fraction))) for chain in arrays)
    half_length = tail_length // 2
    if half_length < 2:
        return float("nan")
    split_chains = []
    for chain in arrays:
        tail = chain[-tail_length:]
        split_chains.extend((tail[:half_length], tail[-half_length:]))
    matrix = np.stack(split_chains)
    within = float(np.mean(np.var(matrix, axis=1, ddof=1)))
    between = float(half_length * np.var(np.mean(matrix, axis=1), ddof=1))
    if within == 0:
        return 1.0 if between == 0 else float("inf")
    variance_estimate = ((half_length - 1) / half_length) * within + between / half_length
    return float(np.sqrt(max(variance_estimate / within, 0.0)))


def consecutive_window_summaries(
    values: list[int] | NDArray[np.float64],
    sweeps: list[int] | NDArray[np.int64],
    *,
    post_transient_fraction: float = 0.5,
    window_count: int = 4,
) -> list[dict[str, float | int]]:
    """Summarize consecutive windows covering the post-transient trace tail."""
    if not 0 < post_transient_fraction <= 1:
        raise ValueError("post_transient_fraction must lie in (0,1]")
    if window_count <= 0:
        raise ValueError("window_count must be positive")
    data = np.asarray(values, dtype=np.float64)
    sweep_array = np.asarray(sweeps, dtype=np.int64)
    if data.ndim != 1 or sweep_array.ndim != 1 or len(data) != len(sweep_array):
        raise ValueError("values and sweeps must be one-dimensional and equally sized")
    if len(data) == 0:
        raise ValueError("window summaries require a nonempty trace")
    tail_points = max(1, int(np.ceil(len(data) * post_transient_fraction)))
    start_index = len(data) - tail_points
    windows = np.array_split(
        np.arange(start_index, len(data), dtype=np.int64), min(window_count, tail_points)
    )
    rows: list[dict[str, float | int]] = []
    for number, indices in enumerate(windows, start=1):
        window = data[indices]
        rows.append(
            {
                "window": number,
                "points": len(window),
                "sweep_start": int(sweep_array[indices[0]]),
                "sweep_end": int(sweep_array[indices[-1]]),
                "mean": float(window.mean()),
                "variance": float(window.var(ddof=1)) if len(window) > 1 else 0.0,
                "minimum": float(window.min()),
                "maximum": float(window.max()),
                "median": float(np.median(window)),
            }
        )
    return rows


def _observable_summary(
    values: list[int],
    sweeps: list[int],
    late_fraction: float,
    record_every_sweeps: int,
    window_count: int,
) -> dict[str, object]:
    start_index = max(0, len(values) - max(1, int(np.ceil(len(values) * late_fraction))))
    late = np.asarray(values[start_index:], dtype=np.float64)
    late_acf = autocorrelation(late, max_lag=min(1, len(late) - 1))
    iact = integrated_autocorrelation_time(late)
    variance = float(late.var(ddof=1)) if len(late) > 1 else 0.0
    effective_sample_size = float(len(late) / iact)
    return {
        "late_window_points": len(late),
        "late_mean": float(late.mean()),
        "late_variance": variance,
        "late_lag1_autocorrelation": float(late_acf[1]) if len(late_acf) > 1 else np.nan,
        "integrated_autocorrelation_time_records": iact,
        "integrated_autocorrelation_time_sweeps": iact * record_every_sweeps,
        "effective_sample_size_late": effective_sample_size,
        "late_mean_mcse": float(np.sqrt(variance / effective_sample_size)),
        "post_transient_iact_multiples": effective_sample_size,
        "post_transient_windows": consecutive_window_summaries(
            values,
            sweeps,
            post_transient_fraction=late_fraction,
            window_count=window_count,
        ),
    }


@dataclass
class MixingChainTrace:
    start: str
    sweeps: list[int]
    center_heights: list[int]
    height_sums: list[int]
    spatial_increments: dict[int, list[int]]
    counts: UpdateCounts
    elapsed_seconds: float

    def observable_series(self) -> dict[str, list[int]]:
        result = {
            "center_height": self.center_heights,
            "height_sum": self.height_sums,
        }
        result.update({f"DeltaH_r{r}": values for r, values in self.spatial_increments.items()})
        return result

    def summary(
        self,
        late_fraction: float = 0.5,
        record_every_sweeps: int = 1,
        window_count: int = 4,
    ) -> dict[str, object]:
        if not 0 < late_fraction <= 1:
            raise ValueError("late_fraction must lie in (0,1]")
        observable_summaries = {
            name: _observable_summary(
                values,
                self.sweeps,
                late_fraction,
                record_every_sweeps,
                window_count,
            )
            for name, values in self.observable_series().items()
        }
        center = observable_summaries["center_height"]
        return {
            "start": self.start,
            "trace_points": len(self.center_heights),
            # Preserve the original centre-height summary fields for callers.
            **center,
            "observables": observable_summaries,
            "flippable_proposal_rate": self.counts.flippable_proposal_rate,
            "actual_change_rate": self.counts.actual_change_rate,
            "attempted_updates": self.counts.attempted_updates,
            "flippable_proposals": self.counts.flippable_proposals,
            "actual_changes": self.counts.actual_changes,
            "null_moves": self.counts.null_moves,
            "elapsed_seconds": self.elapsed_seconds,
        }


@dataclass
class MixingComparison:
    horizontal_start: MixingChainTrace
    vertical_start: MixingChainTrace
    record_every_sweeps: int

    def summary(
        self, *, late_fraction: float = 0.5, window_count: int = 4
    ) -> dict[str, object]:
        horizontal = self.horizontal_start.summary(
            late_fraction=late_fraction,
            record_every_sweeps=self.record_every_sweeps,
            window_count=window_count,
        )
        vertical = self.vertical_start.summary(
            late_fraction=late_fraction,
            record_every_sweeps=self.record_every_sweeps,
            window_count=window_count,
        )
        observable_comparisons = {}
        horizontal_series = self.horizontal_start.observable_series()
        vertical_series = self.vertical_start.observable_series()
        for name in horizontal_series:
            first = horizontal["observables"][name]
            second = vertical["observables"][name]
            gap = float(first["late_mean"] - second["late_mean"])
            standard_error = float(
                np.sqrt(
                    float(first["late_mean_mcse"]) ** 2
                    + float(second["late_mean_mcse"]) ** 2
                )
            )
            z_gap = (
                gap / standard_error
                if standard_error > 0
                else (0.0 if gap == 0 else float(np.copysign(np.inf, gap)))
            )
            observable_comparisons[name] = {
                "horizontal_start": first,
                "vertical_start": second,
                "late_mean_gap_horizontal_minus_vertical": gap,
                "late_variance_difference_horizontal_minus_vertical": float(
                    first["late_variance"] - second["late_variance"]
                ),
                "late_variance_ratio_horizontal_over_vertical": (
                    float(first["late_variance"] / second["late_variance"])
                    if second["late_variance"]
                    else None
                ),
                "late_mean_gap_standard_error": standard_error,
                "z_gap_horizontal_minus_vertical": z_gap,
                "absolute_standardized_late_mean_gap": (
                    abs(z_gap)
                ),
                "split_rhat_late_half": split_rhat(
                    [horizontal_series[name], vertical_series[name]],
                    tail_fraction=late_fraction,
                ),
            }
        center_comparison = observable_comparisons["center_height"]
        return {
            "warning": "Exploratory mixing diagnostic only; agreement does not prove convergence.",
            "record_every_sweeps": self.record_every_sweeps,
            "horizontal_start": horizontal,
            "vertical_start": vertical,
            "observables": observable_comparisons,
            "late_mean_gap_horizontal_minus_vertical": center_comparison[
                "late_mean_gap_horizontal_minus_vertical"
            ],
            "late_variance_ratio_horizontal_over_vertical": center_comparison[
                "late_variance_ratio_horizontal_over_vertical"
            ],
            "split_rhat_center_late_half": center_comparison["split_rhat_late_half"],
            "split_rhat_height_sum_late_half": observable_comparisons["height_sum"][
                "split_rhat_late_half"
            ],
        }


def _run_trace(
    environment: DimerEnvironment,
    state: DimerState,
    start_name: str,
    rng: np.random.Generator,
    sweeps: int,
    record_every: int,
    separations: tuple[SpatialSeparation, ...],
) -> MixingChainTrace:
    attempts_per_sweep = (environment.L - 1) ** 2
    counts = UpdateCounts()
    sweep_values = [0]
    heights = face_heights(state)
    center_values = [int(heights[center_face_index(environment.L)])]
    height_sum_values = [int(heights.sum())]
    increment_values = {item.r: [value] for item, value in zip(
        separations, spatial_increments(heights, separations).values(), strict=True
    )}
    started = perf_counter()
    completed = 0
    while completed < sweeps:
        block = min(record_every, sweeps - completed)
        run_attempts(state, environment, rng, block * attempts_per_sweep, counts)
        completed += block
        heights = face_heights(state)
        sweep_values.append(completed)
        center_values.append(int(heights[center_face_index(environment.L)]))
        height_sum_values.append(int(heights.sum()))
        current = spatial_increments(heights, separations)
        for r, value in current.items():
            increment_values[r].append(value)
    state.validate()
    return MixingChainTrace(
        start_name,
        sweep_values,
        center_values,
        height_sum_values,
        increment_values,
        counts,
        perf_counter() - started,
    )


def compare_extremal_starts(
    environment: DimerEnvironment,
    rng_horizontal: np.random.Generator,
    rng_vertical: np.random.Generator,
    *,
    sweeps: int,
    record_every: int,
    separations: tuple[SpatialSeparation, ...] = (),
) -> MixingComparison:
    """Compare independent all-horizontal and all-vertical starts."""
    if sweeps <= 0:
        raise ValueError("sweeps must be positive")
    if record_every <= 0:
        raise ValueError("record_every must be positive")
    horizontal = _run_trace(
        environment,
        DimerState.all_horizontal(environment.L),
        "all_horizontal",
        rng_horizontal,
        sweeps,
        record_every,
        separations,
    )
    vertical = _run_trace(
        environment,
        DimerState.all_vertical(environment.L),
        "all_vertical",
        rng_vertical,
        sweeps,
        record_every,
        separations,
    )
    return MixingComparison(horizontal, vertical, record_every)
