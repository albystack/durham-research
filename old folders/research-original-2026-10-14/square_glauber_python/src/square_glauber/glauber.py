"""Literal random-all-faces single-face heat-bath dynamics."""

from __future__ import annotations

from dataclasses import dataclass, field
from time import perf_counter

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .environment import DimerEnvironment
from .geometry import SpatialSeparation, center_face_index
from .height import face_heights, spatial_increments
from .matching import DimerState, FaceOrientation


@dataclass(frozen=True)
class UpdateResult:
    face: tuple[int, int]
    flippable: bool
    changed: bool
    old_orientation: FaceOrientation | None
    new_orientation: FaceOrientation | None
    probability_horizontal: float | None


@dataclass
class UpdateCounts:
    attempted_updates: int = 0
    flippable_proposals: int = 0
    actual_changes: int = 0

    @property
    def null_moves(self) -> int:
        """All unchanged attempts, including same-orientation heat-bath redraws."""
        return self.attempted_updates - self.actual_changes

    @property
    def flippable_proposal_rate(self) -> float:
        return self.flippable_proposals / self.attempted_updates if self.attempted_updates else 0.0

    @property
    def actual_change_rate(self) -> float:
        return self.actual_changes / self.attempted_updates if self.attempted_updates else 0.0

    def add(self, result: UpdateResult) -> None:
        self.attempted_updates += 1
        self.flippable_proposals += int(result.flippable)
        self.actual_changes += int(result.changed)


@dataclass(frozen=True)
class TraceRecord:
    sweep: int
    phase: str
    center_height: int
    spatial_increments: dict[int, int]


@dataclass
class ChainResult:
    center_heights: list[int]
    spatial_increment_measurements: dict[int, list[int]]
    final_state: DimerState
    counts: UpdateCounts
    elapsed_seconds: float
    trace: list[TraceRecord] = field(default_factory=list)

    def diagnostics(self) -> dict[str, float | int]:
        return {
            "attempted_updates": self.counts.attempted_updates,
            "flippable_proposals": self.counts.flippable_proposals,
            "actual_changes": self.counts.actual_changes,
            "null_moves": self.counts.null_moves,
            "flippable_proposal_rate": self.counts.flippable_proposal_rate,
            "actual_change_rate": self.counts.actual_change_rate,
            "elapsed_seconds": self.elapsed_seconds,
        }


def heatbath_update(
    state: DimerState,
    environment: DimerEnvironment,
    rng: np.random.Generator,
) -> UpdateResult:
    """Attempt one face chosen uniformly from all ``(L-1)^2`` faces."""
    if state.L != environment.L:
        raise ValueError("matching and environment have different L")
    i = int(rng.integers(0, state.L - 1))
    j = int(rng.integers(0, state.L - 1))
    old = state.face_orientation(i, j)
    if old is None:
        return UpdateResult((i, j), False, False, None, None, None)
    probability_horizontal = environment.horizontal_probability(i, j)
    new: FaceOrientation = "horizontal" if rng.random() < probability_horizontal else "vertical"
    changed = state.set_face_orientation(i, j, new)
    return UpdateResult((i, j), True, changed, old, new, probability_horizontal)


def run_attempts(
    state: DimerState,
    environment: DimerEnvironment,
    rng: np.random.Generator,
    attempts: int,
    counts: UpdateCounts | None = None,
) -> UpdateCounts:
    if attempts < 0:
        raise ValueError("attempts must be nonnegative")
    result_counts = UpdateCounts() if counts is None else counts
    for _ in range(attempts):
        result_counts.add(heatbath_update(state, environment, rng))
    return result_counts


def _validated_supplied_stream(
    L: int,
    face_indices: ArrayLike,
    uniforms: ArrayLike,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Validate a deterministic row-major face/uniform update stream."""
    faces = np.asarray(face_indices)
    random_values = np.asarray(uniforms, dtype=np.float64)
    if faces.ndim != 1 or random_values.ndim != 1:
        raise ValueError("face_indices and uniforms must be one-dimensional")
    if len(faces) != len(random_values):
        raise ValueError("face_indices and uniforms must have equal length")
    if not np.issubdtype(faces.dtype, np.integer):
        raise TypeError("face_indices must have an integer dtype")
    faces = np.asarray(faces, dtype=np.int64)
    face_count = (L - 1) ** 2
    if np.any(faces < 0) or np.any(faces >= face_count):
        raise ValueError(f"face_indices must lie in [0, {face_count})")
    if not np.all(np.isfinite(random_values)) or np.any(random_values < 0) or np.any(
        random_values >= 1
    ):
        raise ValueError("uniforms must be finite values in [0,1)")
    return np.ascontiguousarray(faces), np.ascontiguousarray(random_values)


def heatbath_update_supplied(
    state: DimerState,
    environment: DimerEnvironment,
    face_index: int,
    uniform: float,
) -> UpdateResult:
    """Reference update using externally supplied randomness.

    ``face_index`` is a row-major index in ``[0, (L-1)^2)``.  The supplied
    uniform is ignored on a nonflippable face, exactly as the heat-bath choice
    is irrelevant on a null move.  The ordinary RNG-driven API is unchanged.
    """
    if state.L != environment.L:
        raise ValueError("matching and environment have different L")
    faces, values = _validated_supplied_stream(
        state.L,
        np.asarray([face_index], dtype=np.int64),
        np.asarray([uniform], dtype=np.float64),
    )
    side = state.L - 1
    i, j = divmod(int(faces[0]), side)
    old = state.face_orientation(i, j)
    if old is None:
        return UpdateResult((i, j), False, False, None, None, None)
    probability_horizontal = environment.horizontal_probability(i, j)
    new: FaceOrientation = (
        "horizontal" if values[0] < probability_horizontal else "vertical"
    )
    changed = state.set_face_orientation(i, j, new)
    return UpdateResult((i, j), True, changed, old, new, probability_horizontal)


def run_supplied_attempts(
    state: DimerState,
    environment: DimerEnvironment,
    face_indices: ArrayLike,
    uniforms: ArrayLike,
    counts: UpdateCounts | None = None,
) -> UpdateCounts:
    """Run the authoritative Python backend on a deterministic update stream."""
    if state.L != environment.L:
        raise ValueError("matching and environment have different L")
    faces, values = _validated_supplied_stream(state.L, face_indices, uniforms)
    result_counts = UpdateCounts() if counts is None else counts
    side = state.L - 1
    for face_index, uniform in zip(faces, values, strict=True):
        i, j = divmod(int(face_index), side)
        old = state.face_orientation(i, j)
        if old is None:
            result_counts.attempted_updates += 1
            continue
        probability_horizontal = environment.horizontal_probability(i, j)
        new: FaceOrientation = (
            "horizontal" if uniform < probability_horizontal else "vertical"
        )
        changed = state.set_face_orientation(i, j, new)
        result_counts.attempted_updates += 1
        result_counts.flippable_proposals += 1
        result_counts.actual_changes += int(changed)
    return result_counts


def run_chain(
    environment: DimerEnvironment,
    initial_matching: DimerState,
    rng: np.random.Generator,
    burnin_sweeps: int,
    measurement_gap_sweeps: int,
    num_measurements: int,
    record_trace: bool = False,
    *,
    separations: tuple[SpatialSeparation, ...] = (),
) -> ChainResult:
    """Run a reusable chain; one sweep is exactly ``(L-1)^2`` attempts.

    The first measurement is made immediately after burn-in. Later
    measurements are separated by ``measurement_gap_sweeps`` full sweeps.
    If ``record_trace`` is true, height observables are additionally recorded
    at the initial state and after every executed sweep.
    """
    if initial_matching.L != environment.L:
        raise ValueError("initial matching and environment have different L")
    for name, value in (
        ("burnin_sweeps", burnin_sweeps),
        ("measurement_gap_sweeps", measurement_gap_sweeps),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if isinstance(num_measurements, bool) or not isinstance(num_measurements, int):
        raise ValueError("num_measurements must be a positive integer")
    if num_measurements <= 0:
        raise ValueError("num_measurements must be positive")

    state = initial_matching.copy()
    state.validate()
    attempts_per_sweep = (state.L - 1) ** 2
    counts = UpdateCounts()
    trace: list[TraceRecord] = []
    completed_sweeps = 0
    started = perf_counter()

    def record(phase: str) -> None:
        heights = face_heights(state)
        trace.append(
            TraceRecord(
                completed_sweeps,
                phase,
                int(heights[center_face_index(state.L)]),
                spatial_increments(heights, separations),
            )
        )

    def execute_sweeps(number: int, phase: str) -> None:
        nonlocal completed_sweeps
        for _ in range(number):
            run_attempts(state, environment, rng, attempts_per_sweep, counts)
            completed_sweeps += 1
            if record_trace:
                record(phase)

    if record_trace:
        record("initial")
    execute_sweeps(burnin_sweeps, "burnin")

    centers: list[int] = []
    increments = {item.r: [] for item in separations}
    for measurement in range(num_measurements):
        if measurement:
            execute_sweeps(measurement_gap_sweeps, "measurement_gap")
        heights = face_heights(state)
        centers.append(int(heights[center_face_index(state.L)]))
        values = spatial_increments(heights, separations)
        for r, value in values.items():
            increments[r].append(value)
    elapsed = perf_counter() - started
    state.validate()
    return ChainResult(centers, increments, state, counts, elapsed, trace)
