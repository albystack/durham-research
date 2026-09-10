"""Frozen positive edge environments for the random-bond dimer model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .geometry import validate_L

FloatArray = NDArray[np.float64]
EnvironmentSampler = Callable[[np.random.Generator, tuple[int, ...], float], FloatArray]


def _uniform_weights(
    rng: np.random.Generator, shape: tuple[int, ...], parameter: float
) -> FloatArray:
    del rng, parameter
    return np.ones(shape, dtype=np.float64)


def _gamma_weights(
    rng: np.random.Generator, shape: tuple[int, ...], parameter: float
) -> FloatArray:
    if not np.isfinite(parameter) or parameter <= 0:
        raise ValueError("gamma_shape must be finite and strictly positive")
    return np.asarray(rng.gamma(parameter, 1.0 / parameter, size=shape), dtype=np.float64)


ENVIRONMENT_SAMPLERS: dict[str, EnvironmentSampler] = {
    "uniform": _uniform_weights,
    "gamma": _gamma_weights,
}


@dataclass(frozen=True)
class DimerEnvironment:
    """One immutable realization of undirected grid-edge weights."""

    L: int
    horizontal_weights: FloatArray
    vertical_weights: FloatArray
    model: str = "custom"
    gamma_shape: float | None = None

    def __post_init__(self) -> None:
        L = validate_L(self.L)
        horizontal = np.array(self.horizontal_weights, dtype=np.float64, copy=True)
        vertical = np.array(self.vertical_weights, dtype=np.float64, copy=True)
        if horizontal.shape != (L, L - 1):
            raise ValueError(
                f"horizontal_weights must have shape {(L, L - 1)}, got {horizontal.shape}"
            )
        if vertical.shape != (L - 1, L):
            raise ValueError(
                f"vertical_weights must have shape {(L - 1, L)}, got {vertical.shape}"
            )
        if not np.all(np.isfinite(horizontal)) or not np.all(horizontal > 0):
            raise ValueError("horizontal edge weights must be finite and positive")
        if not np.all(np.isfinite(vertical)) or not np.all(vertical > 0):
            raise ValueError("vertical edge weights must be finite and positive")
        horizontal.setflags(write=False)
        vertical.setflags(write=False)
        object.__setattr__(self, "horizontal_weights", horizontal)
        object.__setattr__(self, "vertical_weights", vertical)

    def face_weights(self, i: int, j: int) -> tuple[float, float, float, float]:
        """Return ``(top, right, bottom, left)`` around bounded face ``(i,j)``."""
        if not (0 <= i < self.L - 1 and 0 <= j < self.L - 1):
            raise IndexError(f"face {(i, j)} is outside the bounded-face array")
        return (
            float(self.horizontal_weights[i, j]),
            float(self.vertical_weights[i, j + 1]),
            float(self.horizontal_weights[i + 1, j]),
            float(self.vertical_weights[i, j]),
        )

    def horizontal_probability(self, i: int, j: int) -> float:
        """Stable evaluation of ``ac / (ac + bd)`` for a face."""
        a, b, c, d = self.face_weights(i, j)
        log_ratio = (np.log(b) + np.log(d)) - (np.log(a) + np.log(c))
        if log_ratio >= 0:
            factor = np.exp(-log_ratio)
            return float(factor / (1.0 + factor))
        factor = np.exp(log_ratio)
        return float(1.0 / (1.0 + factor))


def sample_environment(
    L: int,
    rng: np.random.Generator,
    model: str = "gamma",
    *,
    gamma_shape: float = 1.0,
) -> DimerEnvironment:
    """Sample one environment; its arrays remain fixed for every chain update."""
    validate_L(L)
    try:
        sampler = ENVIRONMENT_SAMPLERS[model]
    except KeyError as error:
        choices = ", ".join(sorted(ENVIRONMENT_SAMPLERS))
        raise ValueError(f"unknown environment model {model!r}; choose {choices}") from error
    parameter = float(gamma_shape)
    horizontal = sampler(rng, (L, L - 1), parameter)
    vertical = sampler(rng, (L - 1, L), parameter)
    return DimerEnvironment(
        L,
        horizontal,
        vertical,
        model=model,
        gamma_shape=parameter if model == "gamma" else None,
    )
