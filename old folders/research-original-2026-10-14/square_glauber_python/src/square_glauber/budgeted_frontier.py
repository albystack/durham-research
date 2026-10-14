"""Bounded wider-square experiment, exact for the stored dyadic edge weights.

Two memory reductions leave the Gibbs law unchanged: scale edges separately at
each black vertex (every matching uses one such edge), and cache only free
vertices after skipping already occupied ones. A postorder stack avoids Python
recursion limits. Width remains exponential; this is a feasibility tool.

The historical reference sampler and CFTP stream identities are unchanged.
"""

from fractions import Fraction
from math import gcd, isfinite
import resource
import sys
from time import perf_counter

import numpy as np

from .integer_sampling import IntegerFrontierSampler

VERSION = "integer-frontier-budgeted-v1"


def peak_rss_bytes():
    """Process high-water RSS: bytes on macOS, KiB on Linux/Hamilton."""
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


class FrontierBudgetExceeded(RuntimeError):
    """A censored computation, never a certified draw or a zero partition."""

    def __init__(self, reason, stats):
        super().__init__(reason)
        self.reason = reason
        self.stats = stats


class BudgetedFrontierSampler(IntegerFrontierSampler):
    """Same exact conditional law as the reference, with an explicit budget.

    Sampling and local branch enumeration are inherited from the tested
    reference. Partition queries after construction only read the completed
    cache. Every instance belongs in a fresh worker for interpretable peak RSS.
    """

    def __init__(
        self,
        environment,
        *,
        max_states=3_000_000,
        max_seconds=120,
        max_rss_bytes=6 * 1024**3,
    ):
        if environment.L > 18:
            raise ValueError("bounded feasibility restricted to L<=18")
        if (
            max_states < 1
            or not isfinite(max_seconds)
            or max_seconds <= 0
            or max_rss_bytes < 1
        ):
            raise ValueError("positive finite budgets required")
        self.environment = environment
        self.L = environment.L
        self.max_states = int(max_states)
        self.max_seconds = float(max_seconds)
        self.max_rss_bytes = int(max_rss_bytes)
        self.started = perf_counter()
        self.memo = {}
        self._check_budget()
        self._integer_weights()
        self._build_partition()
        self.total = self.partition(0, 0)
        if self.total <= 0:
            raise ValueError("zero partition function")
        self.setup_seconds = perf_counter() - self.started

    def _integer_weights(self):
        """One independent exact gauge factor per odd-parity black vertex.

        Multiplying all incident edges at vertex b by g_b multiplies every
        matching by product_b g_b. Opposite endpoints are white, so no edge is
        assigned twice. Dividing by a local gcd reduces integer bit lengths.
        """
        env, L = self.environment, self.L
        self.horizontal = np.empty(env.horizontal_weights.shape, dtype=object)
        self.vertical = np.empty(env.vertical_weights.shape, dtype=object)
        self.matching_scale = Fraction(1)
        for i in range(L):
            for j in range(L):
                if (i + j) % 2 != 1:
                    continue
                edges = []
                if j:
                    edges.append((self.horizontal, env.horizontal_weights, i, j - 1))
                if j + 1 < L:
                    edges.append((self.horizontal, env.horizontal_weights, i, j))
                if i:
                    edges.append((self.vertical, env.vertical_weights, i - 1, j))
                if i + 1 < L:
                    edges.append((self.vertical, env.vertical_weights, i, j))
                ratios = [
                    float(source[r, c]).as_integer_ratio() for _, source, r, c in edges
                ]
                denominator = max(d for _, d in ratios)
                integers = [n * (denominator // d) for n, d in ratios]
                common = 0
                for w in integers:
                    common = gcd(common, w)
                for (target, _, r, c), w in zip(edges, integers):
                    target[r, c] = w // common
                self.matching_scale *= Fraction(denominator, common)

    def _key(self, position, mask):
        # Skipped vertices have no weight factor: their edge was counted when
        # its earlier endpoint was processed.
        while mask & 1:
            position += 1
            mask >>= 1
        return (position << self.L) | mask

    def partition(self, position, mask):
        return self.memo[self._key(position, mask)]

    def stats(self):
        return dict(
            cached_states=len(self.memo),
            elapsed_seconds=perf_counter() - self.started,
            peak_rss_bytes=peak_rss_bytes(),
        )

    def _check_budget(self):
        stats = self.stats()
        if stats["elapsed_seconds"] >= self.max_seconds:
            raise FrontierBudgetExceeded("time_limit", stats)
        if stats["peak_rss_bytes"] >= self.max_rss_bytes:
            raise FrontierBudgetExceeded("memory_limit", stats)

    def _build_partition(self):
        stack = [(self._key(0, 0), False)]
        mask_bits = (1 << self.L) - 1
        iterations = 0
        while stack:
            iterations += 1
            if iterations % 1024 == 0:
                self._check_budget()
            key, expanded = stack.pop()
            if key in self.memo:
                continue
            position, mask = key >> self.L, key & mask_bits
            if position == self.L * self.L:
                value = int(mask == 0)
            else:
                children = [
                    (weight, self._key(position + 1, nextmask))
                    for _, weight, nextmask in self.branches(position, mask)
                ]
                if not expanded:
                    stack.append((key, True))
                    stack.extend(
                        (child, False)
                        for _, child in children
                        if child not in self.memo
                    )
                    continue
                value = sum(weight * self.memo[child] for weight, child in children)
            if len(self.memo) >= self.max_states:
                raise FrontierBudgetExceeded("state_limit", self.stats())
            self.memo[key] = value
        self._check_budget()
