"""Monotone coupling-from-the-past for exact weighted dimer samples.

The random map is exactly the validated all-face single-face heat-bath update.
Every attempted update consumes one uniformly sampled bounded-face index and one
uniform variate; nonflippable choices remain null moves.  No active-face clock
or alternative transition kernel is used.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter
from typing import Literal

import numpy as np

from .accelerated import (
    horizontal_probability_grid,
    require_numba,
    run_supplied_attempts_accelerated,
)
from .environment import DimerEnvironment
from .experiment import _seed_uint64
from .extrema import true_extremal_states
from .glauber import UpdateCounts, run_supplied_attempts
from .height import face_heights
from .matching import DimerState

PerfectBackend = Literal["reference", "numba"]


@dataclass(frozen=True)
class CFTPRandomBlock:
    """Replayable independent randomness assigned to one past-time block."""

    block_index: int
    sweeps: int
    face_seed: int
    uniform_seed: int


@dataclass(frozen=True)
class CFTPIteration:
    horizon_sweeps: int
    block_count: int
    lower_upper_coalesced_at_time_zero: bool
    maximum_pointwise_height_gap_at_time_zero: int


@dataclass
class PerfectSample:
    state: DimerState
    cftp_seed: int
    backend: PerfectBackend
    certified: bool
    horizon_sweeps: int
    iterations: tuple[CFTPIteration, ...]
    blocks: tuple[CFTPRandomBlock, ...]
    lower_update_counts: UpdateCounts
    upper_update_counts: UpdateCounts
    elapsed_seconds: float
    minimum_directed_flips: int
    maximum_directed_flips: int

    def metadata(self) -> dict[str, object]:
        return {
            "cftp_seed": self.cftp_seed,
            "backend": self.backend,
            "certified_perfect_sample": self.certified,
            "coalescence_horizon_sweeps": self.horizon_sweeps,
            "doubling_iterations": len(self.iterations),
            "iterations": [asdict(value) for value in self.iterations],
            "random_blocks": [asdict(value) for value in self.blocks],
            "lower_total_replayed_attempts": self.lower_update_counts.attempted_updates,
            "upper_total_replayed_attempts": self.upper_update_counts.attempted_updates,
            "elapsed_seconds": self.elapsed_seconds,
            "minimum_directed_flips": self.minimum_directed_flips,
            "maximum_directed_flips": self.maximum_directed_flips,
            "kernel": (
                "uniform random scan over all bounded faces; null move if not "
                "flippable; fixed-environment ac/(ac+bd) heat bath otherwise"
            ),
        }


@dataclass(frozen=True)
class PerfectSamplePair:
    """Two conditionally independent perfect samples in one environment."""

    first: PerfectSample
    second: PerfectSample

    def __post_init__(self) -> None:
        if self.first.cftp_seed == self.second.cftp_seed:
            raise ValueError("paired perfect samples require distinct CFTP seeds")
        if not self.first.certified or not self.second.certified:
            raise ValueError("both paired samples must be CFTP-certified")


class CFTPNotCoalesced(RuntimeError):
    """Raised instead of returning a biased/uncertified state at the horizon cap."""

    def __init__(self, horizon_sweeps: int, iterations: tuple[CFTPIteration, ...]):
        super().__init__(
            f"CFTP bounds did not coalesce by {horizon_sweeps} past sweeps; "
            "extend the identical seed/history rather than using a state"
        )
        self.horizon_sweeps = horizon_sweeps
        self.iterations = iterations


def _random_block(cftp_seed: int, L: int, block_index: int, sweeps: int) -> CFTPRandomBlock:
    face_sequence, uniform_sequence = np.random.SeedSequence(
        [int(cftp_seed), int(L), int(block_index), 0xC0F7]
    ).spawn(2)
    return CFTPRandomBlock(
        block_index=block_index,
        sweeps=sweeps,
        face_seed=_seed_uint64(face_sequence),
        uniform_seed=_seed_uint64(uniform_sequence),
    )


def _states_equal(first: DimerState, second: DimerState) -> bool:
    return bool(
        np.array_equal(first.horizontal_occupied, second.horizontal_occupied)
        and np.array_equal(first.vertical_occupied, second.vertical_occupied)
    )


class MonotoneCFTPSampler:
    """Reusable exact sampler for one immutable frozen edge environment."""

    def __init__(
        self,
        environment: DimerEnvironment,
        *,
        backend: PerfectBackend = "numba",
        chunk_attempts: int = 1_000_000,
    ) -> None:
        if backend not in ("reference", "numba"):
            raise ValueError("backend must be 'reference' or 'numba'")
        if backend == "numba":
            require_numba()
        if chunk_attempts <= 0:
            raise ValueError("chunk_attempts must be positive")
        self.environment = environment
        self.backend = backend
        self.chunk_attempts = chunk_attempts
        self.minimum, self.maximum = true_extremal_states(environment.L)
        self.horizontal_probabilities = horizontal_probability_grid(environment)

    def _apply_random_blocks(
        self,
        lower: DimerState,
        upper: DimerState,
        blocks: list[CFTPRandomBlock],
        lower_counts: UpdateCounts,
        upper_counts: UpdateCounts,
        *,
        verify_order: bool,
    ) -> None:
        face_count = (self.environment.L - 1) ** 2
        for block in blocks:
            faces_rng = np.random.Generator(np.random.PCG64(block.face_seed))
            uniforms_rng = np.random.Generator(np.random.PCG64(block.uniform_seed))
            remaining = block.sweeps * face_count
            while remaining:
                count = min(remaining, self.chunk_attempts)
                faces = faces_rng.integers(0, face_count, size=count, dtype=np.int64)
                uniforms = uniforms_rng.random(count)
                if self.backend == "numba":
                    run_supplied_attempts_accelerated(
                        lower,
                        self.environment,
                        faces,
                        uniforms,
                        lower_counts,
                        horizontal_probabilities=self.horizontal_probabilities,
                    )
                    run_supplied_attempts_accelerated(
                        upper,
                        self.environment,
                        faces,
                        uniforms,
                        upper_counts,
                        horizontal_probabilities=self.horizontal_probabilities,
                    )
                else:
                    run_supplied_attempts(
                        lower, self.environment, faces, uniforms, lower_counts
                    )
                    run_supplied_attempts(
                        upper, self.environment, faces, uniforms, upper_counts
                    )
                remaining -= count
                if verify_order and not np.all(
                    face_heights(lower) <= face_heights(upper)
                ):
                    raise AssertionError("common-randomness map violated height order")

    def sample(
        self,
        cftp_seed: int,
        *,
        initial_horizon_sweeps: int = 1,
        max_horizon_sweeps: int = 1_048_576,
        verify_order: bool = True,
    ) -> PerfectSample:
        """Return a certified time-zero sample or raise without returning a state.

        On each unsuccessful iteration, an independent block is prepended in the
        past and every previously assigned block is replayed unchanged.  This is
        the consistent deterministic-history contract needed by CFTP.
        """
        if initial_horizon_sweeps <= 0:
            raise ValueError("initial_horizon_sweeps must be positive")
        if max_horizon_sweeps < initial_horizon_sweeps:
            raise ValueError("max_horizon_sweeps is smaller than the initial horizon")
        started = perf_counter()
        blocks = [
            _random_block(
                int(cftp_seed), self.environment.L, 0, initial_horizon_sweeps
            )
        ]
        iteration_rows: list[CFTPIteration] = []
        lower_total = UpdateCounts()
        upper_total = UpdateCounts()
        while True:
            lower = self.minimum.state.copy()
            upper = self.maximum.state.copy()
            iteration_lower = UpdateCounts()
            iteration_upper = UpdateCounts()
            self._apply_random_blocks(
                lower,
                upper,
                blocks,
                iteration_lower,
                iteration_upper,
                verify_order=verify_order,
            )
            lower_total.attempted_updates += iteration_lower.attempted_updates
            lower_total.flippable_proposals += iteration_lower.flippable_proposals
            lower_total.actual_changes += iteration_lower.actual_changes
            upper_total.attempted_updates += iteration_upper.attempted_updates
            upper_total.flippable_proposals += iteration_upper.flippable_proposals
            upper_total.actual_changes += iteration_upper.actual_changes
            lower_heights = face_heights(lower)
            upper_heights = face_heights(upper)
            ordered = bool(np.all(lower_heights <= upper_heights))
            if not ordered:
                raise AssertionError("CFTP bounds lost pointwise height order")
            coalesced = _states_equal(lower, upper)
            horizon = sum(block.sweeps for block in blocks)
            iteration_rows.append(
                CFTPIteration(
                    horizon_sweeps=horizon,
                    block_count=len(blocks),
                    lower_upper_coalesced_at_time_zero=coalesced,
                    maximum_pointwise_height_gap_at_time_zero=int(
                        np.max(upper_heights - lower_heights)
                    ),
                )
            )
            if coalesced:
                lower.validate()
                return PerfectSample(
                    state=lower,
                    cftp_seed=int(cftp_seed),
                    backend=self.backend,
                    certified=True,
                    horizon_sweeps=horizon,
                    iterations=tuple(iteration_rows),
                    blocks=tuple(blocks),
                    lower_update_counts=lower_total,
                    upper_update_counts=upper_total,
                    elapsed_seconds=perf_counter() - started,
                    minimum_directed_flips=self.minimum.directed_flips,
                    maximum_directed_flips=self.maximum.directed_flips,
                )
            if horizon >= max_horizon_sweeps:
                raise CFTPNotCoalesced(horizon, tuple(iteration_rows))
            extension = min(horizon, max_horizon_sweeps - horizon)
            new_block = _random_block(
                int(cftp_seed),
                self.environment.L,
                len(blocks),
                extension,
            )
            blocks.insert(0, new_block)

    def sample_pair(
        self,
        first_cftp_seed: int,
        second_cftp_seed: int,
        *,
        initial_horizon_sweeps: int = 1,
        max_horizon_sweeps: int = 1_048_576,
        verify_order: bool = True,
    ) -> PerfectSamplePair:
        """Generate independent CFTP samples conditional on this environment."""
        if int(first_cftp_seed) == int(second_cftp_seed):
            raise ValueError("paired perfect samples require distinct CFTP seeds")
        first = self.sample(
            int(first_cftp_seed),
            initial_horizon_sweeps=initial_horizon_sweeps,
            max_horizon_sweeps=max_horizon_sweeps,
            verify_order=verify_order,
        )
        second = self.sample(
            int(second_cftp_seed),
            initial_horizon_sweeps=initial_horizon_sweeps,
            max_horizon_sweeps=max_horizon_sweeps,
            verify_order=verify_order,
        )
        return PerfectSamplePair(first, second)
