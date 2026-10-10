"""Exhaustive distributions versus independent square matching weights."""

from collections import defaultdict
from fractions import Fraction

import numpy as np
import pytest
from square_glauber.comparison import ComparisonModel, comparison_environment
from square_glauber.embedded_shuffle import (
    EmbeddedSquareSampler,
    coin_probability,
    diamond,
    move_and_holes,
    reduce_level,
    validate_host,
)
from square_glauber.exact import enumerate_perfect_matchings, matching_occupancy


def key(h, v):
    return h.tobytes(), v.tobytes()


def full_distribution(sampler):
    distribution = {frozenset(): Fraction(1)}
    for n in range(1, sampler.order + 1):
        next_distribution = defaultdict(Fraction)
        for old, mass in distribution.items():
            moved, holes = move_and_holes(n, old)
            branches = [(moved, mass)]
            for ids in holes:
                p = coin_probability(sampler.levels[n], ids)
                branches = [
                    (edges | set(pair), weight * q)
                    for edges, weight in branches
                    for pair, q in [((ids[0], ids[2]), p), ((ids[1], ids[3]), 1 - p)]
                    if q
                ]
            for edges, weight in branches:
                validate_host(n, edges, sampler.levels[n])
                next_distribution[frozenset(edges)] += weight
        distribution = next_distribution
    result = defaultdict(Fraction)
    for edges, mass in distribution.items():
        s = sampler.extract(edges)
        result[key(s.horizontal_occupied, s.vertical_occupied)] += mass
    return result


@pytest.mark.parametrize("L", [2, 4])
@pytest.mark.parametrize("model", ["aztec_gamma", "aztec_mean", "uniform"])
def test_every_shuffle_outcome_equals_square_gibbs(L, model):
    env = comparison_environment(
        L, np.random.default_rng(813), ComparisonModel(name=model)
    )
    sampler = EmbeddedSquareSampler(env)
    actual = full_distribution(sampler)
    raw = {}
    for matching in enumerate_perfect_matchings(L):
        h, v = matching_occupancy(matching, L, L)
        weight = Fraction(1)
        for a, occupied in [(env.horizontal_weights, h), (env.vertical_weights, v)]:
            for w in a[occupied]:
                weight *= Fraction(float(w))
        raw[key(h, v)] = weight
    total = sum(raw.values())
    assert actual == {k: weight / total for k, weight in raw.items()}
    assert sum(actual.values()) == 1


@pytest.mark.parametrize("L", [6, 8, 12])
def test_float_embedding_forbids_crossings_and_returns_valid_squares(L):
    env = comparison_environment(L, np.random.default_rng(917 + L), ComparisonModel())
    sampler = EmbeddedSquareSampler(env, cast=float)
    for seed in range(5):
        sampler.sample(seed).validate()
    assert len(sampler.exterior_witness) * 2 == len(diamond(L - 1).vertices) - L * L


def test_uniform_six_square_entire_distribution():
    env = comparison_environment(
        6, np.random.default_rng(1), ComparisonModel(name="uniform")
    )
    actual = full_distribution(EmbeddedSquareSampler(env))
    assert len(actual) == 6728
    assert set(actual.values()) == {Fraction(1, 6728)}


def test_singular_cells_and_positive_underflow_are_never_silently_repaired():
    with pytest.raises(ValueError, match="singular"):
        reduce_level(2, [0.0] * len(diamond(2).edges))
    with pytest.raises(ArithmeticError, match="underflow"):
        reduce_level(2, [1e-250] * len(diamond(2).edges))


def test_arbitrary_all_edge_weights_distribution():
    from square_glauber.environment import DimerEnvironment
    from square_glauber.integer_sampling import IntegerFrontierSampler

    rng = np.random.default_rng(192)
    env = DimerEnvironment(4, rng.uniform(0.1, 2, (4, 3)), rng.uniform(0.1, 2, (3, 4)))
    actual = full_distribution(EmbeddedSquareSampler(env))
    oracle = IntegerFrontierSampler(env)
    from square_glauber.matching import DimerState

    for matching in enumerate_perfect_matchings(4):
        h, v = matching_occupancy(matching, 4, 4)
        assert actual[key(h, v)] == oracle.matching_probability(DimerState(4, h, v))
