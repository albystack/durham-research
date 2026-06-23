"""Central configuration for the random dimer/LERW experiments.

This is the one file to edit when the research question changes.  It contains
the lattice directions, the random edge-weight models, and the experiment
presets used by the command-line runner.

Important convention:
    All random models below randomise every nearest-neighbour edge.  The old
    two-random-weight template (N=E=1, S=u, W=v) is deliberately not present.
"""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import dataclass


# Directions are ordered clockwise around a vertex.  The simulation uses this
# order both for transition weights and for computing turns in the LERW path.
DIRECTIONS = ("N", "E", "S", "W")

STEP_VECTOR: dict[str, tuple[int, int]] = {
    "N": (0, 1),
    "E": (1, 0),
    "S": (0, -1),
    "W": (-1, 0),
}


@dataclass(frozen=True)
class ModelSpec:
    """A named positive edge-weight distribution.

    A model describes the law of one undirected lattice-edge weight.  In a
    trial, every edge encountered by the walk receives one fixed sample from
    this law and keeps that value for the rest of the trial.  That is the
    quenched random environment.
    """

    name: str
    parameter_name: str | None
    distribution: str
    description: str
    color: str


@dataclass(frozen=True)
class ExperimentPreset:
    """A complete simulation recipe used by `python3 main.py run NAME`."""

    description: str
    sizes: tuple[int, ...]
    samples_per_size: int
    models: tuple[str, ...]
    base_seed: int = 20260618
    max_steps: int | None = None
    workers: int = 1
    bootstrap_reps: int = 250


MODEL_LIBRARY: dict[str, ModelSpec] = {
    "symmetric": ModelSpec(
        name="symmetric",
        parameter_name=None,
        distribution="constant",
        description="Simple random walk baseline: every edge has weight 1.",
        color="#111827",
    ),
    "gamma_edges": ModelSpec(
        name="gamma_edges",
        parameter_name="shape",
        distribution="gamma",
        description="IID Gamma(shape, scale=1/shape) edge weights with mean 1.",
        color="#0f766e",
    ),
    "exp_edges": ModelSpec(
        name="exp_edges",
        parameter_name=None,
        distribution="exponential",
        description="IID Exponential(1) edge weights.",
        color="#2563eb",
    ),
    "lognormal_edges": ModelSpec(
        name="lognormal_edges",
        parameter_name="sigma",
        distribution="lognormal",
        description="IID lognormal edge weights normalised to mean 1.",
        color="#7c3aed",
    ),
    "pareto_edges": ModelSpec(
        name="pareto_edges",
        parameter_name="alpha",
        distribution="pareto",
        description="IID scaled Pareto(alpha) edge weights, alpha > 1.",
        color="#c2410c",
    ),
    "uniform_edges": ModelSpec(
        name="uniform_edges",
        parameter_name="a",
        distribution="uniform",
        description="IID Uniform(1-a, 1+a) edge weights, 0 < a < 1.",
        color="#047857",
    ),
    "beta_edges": ModelSpec(
        name="beta_edges",
        parameter_name="a",
        distribution="beta",
        description="IID Beta(a,a) weights rescaled to mean 1.",
        color="#a21caf",
    ),
    "weibull_edges": ModelSpec(
        name="weibull_edges",
        parameter_name="shape",
        distribution="weibull",
        description="IID Weibull edge weights rescaled to mean 1.",
        color="#0891b2",
    ),
    "inverse_gamma_edges": ModelSpec(
        name="inverse_gamma_edges",
        parameter_name="alpha",
        distribution="inverse_gamma",
        description="IID inverse-Gamma edge weights with mean 1, alpha > 1.",
        color="#be123c",
    ),
    "bernoulli_edges": ModelSpec(
        name="bernoulli_edges",
        parameter_name="a",
        distribution="bernoulli",
        description="IID two-point weights 1-a or 1+a with equal probability.",
        color="#4d7c0f",
    ),
    "triangular_edges": ModelSpec(
        name="triangular_edges",
        parameter_name="a",
        distribution="triangular",
        description="IID triangular weights on [1-a, 1+a] with mode 1.",
        color="#9333ea",
    ),
}


# The main preset deliberately contains roughly ten different disorder laws.
# Smoke is tiny and is meant only to check that the pipeline still runs.
EXPERIMENTS: dict[str, ExperimentPreset] = {
    "smoke": ExperimentPreset(
        description="Tiny run for checking the command line and output files.",
        sizes=(8, 16),
        samples_per_size=5,
        models=(
            "symmetric",
            "gamma_edges:1.0",
            "exp_edges",
            "lognormal_edges:0.6",
        ),
        workers=1,
        bootstrap_reps=50,
    ),
    "pilot": ExperimentPreset(
        description="Medium run for checking finite-size trends before a long job.",
        sizes=(16, 32, 64, 128),
        samples_per_size=200,
        models=(
            "symmetric",
            "gamma_edges:0.7",
            "gamma_edges:1.0",
            "exp_edges",
            "lognormal_edges:0.7",
            "pareto_edges:2.2",
            "uniform_edges:0.6",
            "beta_edges:0.7",
            "weibull_edges:0.8",
            "inverse_gamma_edges:2.5",
            "bernoulli_edges:0.7",
            "triangular_edges:0.8",
        ),
        workers=0,
        bootstrap_reps=250,
    ),
    "super_rough_probe": ExperimentPreset(
        description="Main probe for Var(h(0)) proxy scaling: log L vs (log L)^2.",
        sizes=(16, 32, 64, 128, 256, 512),
        samples_per_size=1500,
        models=(
            "symmetric",
            "gamma_edges:0.5",
            "gamma_edges:1.0",
            "gamma_edges:2.0",
            "exp_edges",
            "lognormal_edges:0.5",
            "lognormal_edges:1.0",
            "pareto_edges:2.0",
            "pareto_edges:3.0",
            "uniform_edges:0.7",
            "beta_edges:0.5",
            "weibull_edges:0.7",
            "inverse_gamma_edges:2.2",
            "bernoulli_edges:0.8",
            "triangular_edges:0.8",
        ),
        workers=0,
        bootstrap_reps=500,
    ),
    "super_rough_large": ExperimentPreset(
        description="Long continuation of the main probe restricted to L=512 and L=1024.",
        sizes=(512, 1024),
        samples_per_size=1500,
        models=(
            "symmetric",
            "gamma_edges:0.5",
            "gamma_edges:1.0",
            "gamma_edges:2.0",
            "exp_edges",
            "lognormal_edges:0.5",
            "lognormal_edges:1.0",
            "pareto_edges:2.0",
            "pareto_edges:3.0",
            "uniform_edges:0.7",
            "beta_edges:0.5",
            "weibull_edges:0.7",
            "inverse_gamma_edges:2.2",
            "bernoulli_edges:0.8",
            "triangular_edges:0.8",
        ),
        workers=0,
        bootstrap_reps=500,
    ),
}


def parse_model(text: str) -> tuple[str, float | None]:
    """Parse strings such as `gamma_edges:1.0` into a model and parameter."""

    if ":" in text:
        model_name, raw_parameter = text.split(":", 1)
        parameter = float(raw_parameter)
    else:
        model_name = text
        parameter = None
    validate_model(model_name, parameter)
    return model_name, parameter


def validate_model(model_name: str, parameter: float | None) -> None:
    """Fail early when a model string cannot define a positive distribution."""

    if model_name not in MODEL_LIBRARY:
        available = ", ".join(sorted(MODEL_LIBRARY))
        raise ValueError(f"unknown model {model_name!r}. Available models: {available}")

    spec = MODEL_LIBRARY[model_name]
    if spec.parameter_name is None and parameter is not None:
        raise ValueError(f"{model_name} does not take a parameter")
    if spec.parameter_name is not None and parameter is None:
        raise ValueError(f"{model_name} needs parameter `{spec.parameter_name}`")

    if parameter is not None and parameter <= 0:
        raise ValueError(f"{model_name} needs a positive parameter")
    if spec.distribution == "pareto" and parameter is not None and parameter <= 1:
        raise ValueError("pareto_edges needs alpha > 1")
    if spec.distribution == "inverse_gamma" and parameter is not None and parameter <= 1:
        raise ValueError("inverse_gamma_edges needs alpha > 1")
    if spec.distribution in {"uniform", "bernoulli", "triangular"} and parameter is not None:
        if parameter >= 1:
            raise ValueError(f"{model_name} needs 0 < a < 1")


def model_label(model_name: str, parameter: float | str | None) -> str:
    """Human-readable label used in tables and plots."""

    if parameter in (None, "", "None"):
        return model_name
    spec = MODEL_LIBRARY[model_name]
    return f"{model_name} ({spec.parameter_name}={float(parameter):g})"


def safe_model_name(model_name: str, parameter: float | str | None) -> str:
    """Filename-safe model label."""

    if parameter in (None, "", "None"):
        return model_name
    return f"{model_name}_{float(parameter):g}".replace(".", "p").replace("-", "m")


def model_color(model_name: str) -> str:
    """Preferred plot colour for a model family."""

    return MODEL_LIBRARY[model_name].color


def stable_seed(*parts: object) -> int:
    """Create a reproducible integer seed from arbitrary identifying parts."""

    text = ":".join(str(part) for part in parts)
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")


def _require_parameter(distribution: str, parameter: float | None) -> float:
    if parameter is None:
        raise ValueError(f"{distribution} needs a parameter")
    return parameter


def sample_edge_weight(rng: random.Random, model_name: str, parameter: float | None) -> float:
    """Sample one positive edge weight from the chosen model.

    The returned weights are usually normalised to have mean one.  That keeps
    the overall time scale comparable across models; the disorder shape, not
    the average edge strength, is what changes.
    """

    validate_model(model_name, parameter)
    distribution = MODEL_LIBRARY[model_name].distribution

    if distribution == "constant":
        return 1.0
    if distribution == "gamma":
        shape = _require_parameter(distribution, parameter)
        return rng.gammavariate(shape, 1.0 / shape)
    if distribution == "exponential":
        return rng.expovariate(1.0)
    if distribution == "lognormal":
        sigma = _require_parameter(distribution, parameter)
        return rng.lognormvariate(-0.5 * sigma * sigma, sigma)
    if distribution == "pareto":
        alpha = _require_parameter(distribution, parameter)
        return ((alpha - 1.0) / alpha) * rng.paretovariate(alpha)
    if distribution == "uniform":
        a = _require_parameter(distribution, parameter)
        return rng.uniform(1.0 - a, 1.0 + a)
    if distribution == "beta":
        a = _require_parameter(distribution, parameter)
        return 2.0 * rng.betavariate(a, a)
    if distribution == "weibull":
        shape = _require_parameter(distribution, parameter)
        scale = 1.0 / math.gamma(1.0 + 1.0 / shape)
        return rng.weibullvariate(scale, shape)
    if distribution == "inverse_gamma":
        alpha = _require_parameter(distribution, parameter)
        return (alpha - 1.0) / rng.gammavariate(alpha, 1.0)
    if distribution == "bernoulli":
        a = _require_parameter(distribution, parameter)
        return 1.0 + a if rng.random() < 0.5 else 1.0 - a
    if distribution == "triangular":
        a = _require_parameter(distribution, parameter)
        return rng.triangular(1.0 - a, 1.0 + a, 1.0)

    raise ValueError(f"unknown distribution {distribution!r}")


def available_models_text() -> str:
    """Compact model catalogue for the CLI."""

    lines = []
    for name, spec in MODEL_LIBRARY.items():
        if spec.parameter_name:
            signature = f"{name}:<{spec.parameter_name}>"
        else:
            signature = name
        lines.append(f"{signature:30s} {spec.description}")
    return "\n".join(lines)
