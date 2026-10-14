"""Optional Numba backend for the validated random-scan heat-bath kernel."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Final

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .diagnostics import MixingChainTrace, MixingComparison
from .environment import DimerEnvironment
from .geometry import SpatialSeparation, center_face_index
from .glauber import UpdateCounts, _validated_supplied_stream
from .height import face_heights, spatial_increments
from .matching import DimerState

try:
    import numba as _numba
except ImportError:  # pragma: no cover - exercised in environments without the extra
    _numba = None

NUMBA_AVAILABLE: Final = _numba is not None


@dataclass
class AcceleratedRandomStreams:
    """Persistent independent streams for faces and heat-bath uniforms."""

    face_rng: np.random.Generator
    uniform_rng: np.random.Generator

    @classmethod
    def from_generator(cls, rng: np.random.Generator) -> AcceleratedRandomStreams:
        seeds = np.asarray(rng.bit_generator.random_raw(2), dtype=np.uint64)
        return cls(
            np.random.Generator(np.random.PCG64(int(seeds[0]))),
            np.random.Generator(np.random.PCG64(int(seeds[1]))),
        )


def require_numba() -> None:
    if not NUMBA_AVAILABLE:
        raise RuntimeError(
            "the accelerated backend requires Numba; install square-glauber[accel]"
        )


def horizontal_probability_grid(environment: DimerEnvironment) -> NDArray[np.float64]:
    """Precompute the authoritative stable heat-bath probability for each face."""
    side = environment.L - 1
    probabilities = np.empty((side, side), dtype=np.float64)
    for i in range(side):
        for j in range(side):
            probabilities[i, j] = environment.horizontal_probability(i, j)
    return probabilities


if NUMBA_AVAILABLE:

    @_numba.njit(cache=True)
    def _numba_supplied_kernel(
        horizontal_occupied: NDArray[np.bool_],
        vertical_occupied: NDArray[np.bool_],
        horizontal_probabilities: NDArray[np.float64],
        face_indices: NDArray[np.int64],
        uniforms: NDArray[np.float64],
    ) -> tuple[int, int]:
        side = horizontal_probabilities.shape[0]
        flippable = 0
        changed = 0
        for index in range(face_indices.size):
            face = face_indices[index]
            i = face // side
            j = face - i * side
            top = horizontal_occupied[i, j]
            bottom = horizontal_occupied[i + 1, j]
            left = vertical_occupied[i, j]
            right = vertical_occupied[i, j + 1]
            is_horizontal = top and bottom and not left and not right
            is_vertical = left and right and not top and not bottom
            if not is_horizontal and not is_vertical:
                continue
            flippable += 1
            choose_horizontal = uniforms[index] < horizontal_probabilities[i, j]
            if choose_horizontal == is_horizontal:
                continue
            horizontal_occupied[i, j] = choose_horizontal
            horizontal_occupied[i + 1, j] = choose_horizontal
            vertical_occupied[i, j] = not choose_horizontal
            vertical_occupied[i, j + 1] = not choose_horizontal
            changed += 1
        return flippable, changed

else:
    _numba_supplied_kernel = None


def run_supplied_attempts_accelerated(
    state: DimerState,
    environment: DimerEnvironment,
    face_indices: ArrayLike,
    uniforms: ArrayLike,
    counts: UpdateCounts | None = None,
    *,
    horizontal_probabilities: NDArray[np.float64] | None = None,
) -> UpdateCounts:
    """Run Numba updates in-place using an externally supplied random stream."""
    require_numba()
    if state.L != environment.L:
        raise ValueError("matching and environment have different L")
    faces, values = _validated_supplied_stream(state.L, face_indices, uniforms)
    probabilities = (
        horizontal_probability_grid(environment)
        if horizontal_probabilities is None
        else np.asarray(horizontal_probabilities, dtype=np.float64)
    )
    expected_shape = (state.L - 1, state.L - 1)
    if probabilities.shape != expected_shape:
        raise ValueError(
            f"horizontal_probabilities must have shape {expected_shape}, got {probabilities.shape}"
        )
    probabilities = np.ascontiguousarray(probabilities)
    flippable, changed = _numba_supplied_kernel(
        state.horizontal_occupied,
        state.vertical_occupied,
        probabilities,
        faces,
        values,
    )
    result = UpdateCounts() if counts is None else counts
    result.attempted_updates += len(faces)
    result.flippable_proposals += int(flippable)
    result.actual_changes += int(changed)
    return result


def run_attempts_accelerated(
    state: DimerState,
    environment: DimerEnvironment,
    rng: np.random.Generator | AcceleratedRandomStreams,
    attempts: int,
    counts: UpdateCounts | None = None,
    *,
    chunk_attempts: int = 1_000_000,
    horizontal_probabilities: NDArray[np.float64] | None = None,
) -> UpdateCounts:
    """Ordinary RNG-driven accelerated API with the same transition kernel.

    Faces are sampled uniformly from all bounded faces.  One uniform is drawn
    per attempted update and ignored on null moves.  This consumes RNG values
    differently from the reference convenience API, but implements the same
    Markov transition law.  Use the supplied-stream APIs for bitwise backend
    equivalence tests.
    """
    require_numba()
    if attempts < 0:
        raise ValueError("attempts must be nonnegative")
    if chunk_attempts <= 0:
        raise ValueError("chunk_attempts must be positive")
    result = UpdateCounts() if counts is None else counts
    streams = (
        rng
        if isinstance(rng, AcceleratedRandomStreams)
        else AcceleratedRandomStreams.from_generator(rng)
    )
    probabilities = (
        horizontal_probability_grid(environment)
        if horizontal_probabilities is None
        else horizontal_probabilities
    )
    remaining = attempts
    face_count = (state.L - 1) ** 2
    while remaining:
        block = min(remaining, chunk_attempts)
        faces = streams.face_rng.integers(0, face_count, size=block, dtype=np.int64)
        uniforms = streams.uniform_rng.random(block)
        run_supplied_attempts_accelerated(
            state,
            environment,
            faces,
            uniforms,
            result,
            horizontal_probabilities=probabilities,
        )
        remaining -= block
    return result


def _run_trace_accelerated(
    environment: DimerEnvironment,
    state: DimerState,
    start_name: str,
    rng: np.random.Generator,
    sweeps: int,
    record_every: int,
    separations: tuple[SpatialSeparation, ...],
) -> MixingChainTrace:
    """Accelerated chain with heights derived only at recording boundaries."""
    attempts_per_sweep = (environment.L - 1) ** 2
    probabilities = horizontal_probability_grid(environment)
    streams = AcceleratedRandomStreams.from_generator(rng)
    counts = UpdateCounts()
    sweep_values = [0]
    heights = face_heights(state)
    center_values = [int(heights[center_face_index(environment.L)])]
    height_sum_values = [int(heights.sum())]
    initial_increments = spatial_increments(heights, separations)
    increment_values = {r: [value] for r, value in initial_increments.items()}
    started = perf_counter()
    completed = 0
    while completed < sweeps:
        block = min(record_every, sweeps - completed)
        run_attempts_accelerated(
            state,
            environment,
            streams,
            block * attempts_per_sweep,
            counts,
            horizontal_probabilities=probabilities,
        )
        completed += block
        heights = face_heights(state)
        sweep_values.append(completed)
        center_values.append(int(heights[center_face_index(environment.L)]))
        height_sum_values.append(int(heights.sum()))
        for r, value in spatial_increments(heights, separations).items():
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


def compare_extremal_starts_accelerated(
    environment: DimerEnvironment,
    rng_horizontal: np.random.Generator,
    rng_vertical: np.random.Generator,
    *,
    sweeps: int,
    record_every: int,
    separations: tuple[SpatialSeparation, ...] = (),
) -> MixingComparison:
    """Independent extremal-start diagnostic using the optional Numba backend."""
    if sweeps <= 0 or record_every <= 0:
        raise ValueError("sweeps and record_every must be positive")
    horizontal = _run_trace_accelerated(
        environment,
        DimerState.all_horizontal(environment.L),
        "all_horizontal",
        rng_horizontal,
        sweeps,
        record_every,
        separations,
    )
    vertical = _run_trace_accelerated(
        environment,
        DimerState.all_vertical(environment.L),
        "all_vertical",
        rng_vertical,
        sweeps,
        record_every,
        separations,
    )
    return MixingComparison(horizontal, vertical, record_every)
