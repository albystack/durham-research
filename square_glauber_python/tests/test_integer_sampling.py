from fractions import Fraction
import random
import hashlib

import numpy as np
import pytest

from square_glauber.environment import DimerEnvironment
from square_glauber.comparison import ComparisonModel, comparison_environment
from square_glauber.exact import enumerate_perfect_matchings, matching_occupancy
from square_glauber.extrema import true_extremal_states
from square_glauber.height import face_heights
from square_glauber.integer_sampling import (
    IntegerFrontierSampler,
    bernoulli_ratio,
    face_ratio,
    coupled_integer_update,
    randbelow_bits,
    load_verified_environment,
)
from square_glauber.matching import DimerState


class ScriptedBits:
    def __init__(self, values):
        self.values = iter(values)

    def getrandbits(self, k):
        value = next(self.values)
        assert 0 <= value < 2**k
        return value


def test_original_environment_loading_requires_archive_identity(tmp_path):
    env = comparison_environment(4, np.random.default_rng(31), ComparisonModel())
    path = tmp_path / "weights.npz"
    np.savez_compressed(
        path, horizontal=env.horizontal_weights, vertical=env.vertical_weights
    )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    loaded = load_verified_environment(path, 4, env.model, digest)
    assert np.array_equal(loaded.horizontal_weights, env.horizontal_weights)
    assert np.array_equal(loaded.vertical_weights, env.vertical_weights)
    with pytest.raises(ValueError, match="checksum"):
        load_verified_environment(path, 4, env.model, "wrong")


def test_rejection_and_extremely_rare_outcomes_remain_possible():
    assert randbelow_bits(ScriptedBits([7, 6, 4]), 5) == 4
    n = 2**200
    assert bernoulli_ratio(ScriptedBits([0]), 1, n)
    assert not bernoulli_ratio(ScriptedBits([n - 1]), 1, n)
    assert not bernoulli_ratio(ScriptedBits([n - 1]), n - 1, n)
    with pytest.raises(ValueError):
        bernoulli_ratio(random.Random(1), 2, 1)


def test_rounded_one_is_not_one_in_integer_kernel():
    env = DimerEnvironment(2, np.ones((2, 1)), np.array([[2.0**-100, 2.0**-100]]))
    assert env.horizontal_probability(0, 0) == 1.0
    n, d = face_ratio(env, 0, 0)
    assert Fraction(n, d) == Fraction(2**200, 2**200 + 1)
    assert not bernoulli_ratio(ScriptedBits([2**200]), n, d)


@pytest.mark.parametrize("L", [2, 4])
def test_all_tiny_matching_masses_equal_independent_rational_gibbs(L):
    env = comparison_environment(L, np.random.default_rng(782 + L), ComparisonModel())
    sampler = IntegerFrontierSampler(env)
    states = []
    products = []
    for matching in enumerate_perfect_matchings(L):
        h, v = matching_occupancy(matching, L, L)
        states.append(DimerState(L, h, v))
        product = Fraction(1)
        for w in env.horizontal_weights[h]:
            product *= Fraction(float(w))
        for w in env.vertical_weights[v]:
            product *= Fraction(float(w))
        products.append(product)
    z = sum(products)
    assert Fraction(sampler.total) / (sampler.edge_scale ** (L * L // 2)) == z
    assert all(
        sampler.matching_probability(s) == p / z for s, p in zip(states, products)
    )
    assert sum(sampler.matching_probability(s) for s in states) == 1
    for seed in range(20):
        sampler.sample(seed).validate()


def test_all_conditional_branch_masses_are_normalized():
    env = comparison_environment(4, np.random.default_rng(44), ComparisonModel())
    sampler = IntegerFrontierSampler(env)

    def visit(pos, mask, prob):
        if pos == 16:
            return prob
        if mask & 1:
            return visit(pos + 1, mask >> 1, prob)
        total = Fraction(0)
        for _, weight, nextmask in sampler.branches(pos, mask):
            mass = weight * sampler.partition(pos + 1, nextmask)
            if mass:
                total += visit(
                    pos + 1,
                    nextmask,
                    prob * Fraction(mass, sampler.partition(pos, mask)),
                )
        return total

    assert visit(0, 0, Fraction(1)) == 1


def test_new_common_maps_preserve_height_order_and_seed_replay():
    env = comparison_environment(4, np.random.default_rng(17), ComparisonModel())
    minimum, maximum = true_extremal_states(4)
    lo, hi = minimum.state.copy(), maximum.state.copy()
    rng = random.Random(19)
    for _ in range(200):
        coupled_integer_update(lo, hi, env, rng)
        assert np.all(face_heights(lo) <= face_heights(hi))
    sampler = IntegerFrontierSampler(env)
    a, b = sampler.sample(21), sampler.sample(21)
    assert np.array_equal(a.horizontal_occupied, b.horizontal_occupied)
    with pytest.raises(ValueError, match="L<=12"):
        IntegerFrontierSampler(
            comparison_environment(14, np.random.default_rng(1), ComparisonModel())
        )


def test_integer_heatbath_detailed_balance_on_every_tiny_flip():
    env = comparison_environment(4, np.random.default_rng(319), ComparisonModel())
    sampler = IntegerFrontierSampler(env)
    for matching in enumerate_perfect_matchings(4):
        h, v = matching_occupancy(matching, 4, 4)
        state = DimerState(4, h, v)
        for i in range(3):
            for j in range(3):
                orientation = state.face_orientation(i, j)
                if orientation is None:
                    continue
                n, d = face_ratio(env, i, j)
                p = Fraction(n, d)
                flipped = state.copy()
                flipped.flip_face(i, j)
                forward = 1 - p if orientation == "horizontal" else p
                reverse = p if orientation == "horizontal" else 1 - p
                assert (
                    sampler.matching_probability(state) * forward
                    == sampler.matching_probability(flipped) * reverse
                )
