"""Independent laws, exact gauge invariance and honest resource censoring."""

from fractions import Fraction

import numpy as np
import pytest

from square_glauber.budgeted_frontier import (
    BudgetedFrontierSampler,
    FrontierBudgetExceeded,
)
from square_glauber.comparison import ComparisonModel, comparison_environment
from square_glauber.environment import DimerEnvironment
from square_glauber.exact import enumerate_perfect_matchings, matching_occupancy
from square_glauber.integer_sampling import IntegerFrontierSampler
from square_glauber.matching import DimerState


@pytest.mark.parametrize("model", ["aztec_gamma", "aztec_mean", "uniform"])
def test_gauge_and_every_small_matching_mass(model):
    env = comparison_environment(
        4, np.random.default_rng(47), ComparisonModel(name=model)
    )
    new = BudgetedFrontierSampler(env)
    old = IntegerFrontierSampler(env)
    products, states = [], []
    for matching in enumerate_perfect_matchings(4):
        h, v = matching_occupancy(matching, 4, 4)
        state = DimerState(4, h, v)
        states.append(state)
        product = Fraction(1)
        for a, occupied in [(env.horizontal_weights, h), (env.vertical_weights, v)]:
            for w in a[occupied]:
                product *= Fraction(float(w))
        products.append(product)
    assert Fraction(new.total) / new.matching_scale == sum(products)
    for state, weight in zip(states, products):
        assert new.matching_probability(state) == weight / sum(products)
        assert new.matching_probability(state) == old.matching_probability(state)
    for seed in range(12):
        a, b = new.sample(seed), new.sample(seed)
        assert np.array_equal(a.horizontal_occupied, b.horizontal_occupied)
        assert np.array_equal(a.vertical_occupied, b.vertical_occupied)

    # Gauge factors change integer bit lengths and RNG consumption; compare
    # conditional probabilities, not bitwise trajectories across versions.
    pending, seen = [(0, 0)], set()
    while pending:
        pos, mask = pending.pop()
        if (pos, mask) in seen or pos == 16:
            continue
        seen.add((pos, mask))
        if mask & 1:
            pending.append((pos + 1, mask >> 1))
            continue
        if not old.partition(pos, mask):
            continue
        for (_, ow, child), (_, nw, same_child) in zip(
            old.branches(pos, mask), new.branches(pos, mask)
        ):
            assert child == same_child
            assert Fraction(
                ow * old.partition(pos + 1, child), old.partition(pos, mask)
            ) == Fraction(nw * new.partition(pos + 1, child), new.partition(pos, mask))
            pending.append((pos + 1, child))


def test_deeper_partition_and_subnormal_positive_weights():
    env = comparison_environment(8, np.random.default_rng(59), ComparisonModel())
    old, new = IntegerFrontierSampler(env), BudgetedFrontierSampler(env)
    assert (
        Fraction(new.total) / new.matching_scale
        == Fraction(old.total) / old.edge_scale**32
    )
    assert new.stats()["cached_states"] < old.partition.cache_info().currsize
    tiny = np.nextafter(0.0, 1.0)
    env = DimerEnvironment(2, np.array([[tiny], [1.0]]), np.array([[1.0, 1.0]]))
    sampler = BudgetedFrontierSampler(env)
    assert sampler.matching_probability(DimerState.all_horizontal(2)) == Fraction(
        tiny
    ) / (1 + Fraction(tiny))


@pytest.mark.parametrize(
    "budget,reason",
    [
        ({"max_states": 1}, "state_limit"),
        ({"max_seconds": 1e-15}, "time_limit"),
        ({"max_rss_bytes": 1}, "memory_limit"),
    ],
)
def test_resource_exhaustion_never_returns_a_sampler(budget, reason):
    env = comparison_environment(4, np.random.default_rng(1), ComparisonModel())
    with pytest.raises(FrontierBudgetExceeded) as caught:
        BudgetedFrontierSampler(env, **budget)
    assert caught.value.reason == reason
    assert caught.value.stats["cached_states"] >= 0


def test_larger_widths_remain_guarded():
    env = comparison_environment(20, np.random.default_rng(1), ComparisonModel())
    with pytest.raises(ValueError, match="L<=18"):
        BudgetedFrontierSampler(env)
