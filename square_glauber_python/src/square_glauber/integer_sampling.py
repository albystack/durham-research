"""Versioned exact-integer sampling for positive stored Float64 edge weights.

A Float64 is a dyadic rational. Multiplying every edge by a common power of two
makes all weights integers and multiplies every matching weight by the same
factor. Integer partition functions and random-bit rejection therefore avoid
rounding probabilities to zero/one. Exactness is conditional on the stored
weights and the ideal independent-bit model; it does not make Gamma generation
or a pseudorandom generator mathematically exact. Historical samplers unchanged.
"""

from functools import cache
from fractions import Fraction
from math import gcd
import random
import hashlib
from pathlib import Path

import numpy as np

from .environment import DimerEnvironment
from .matching import DimerState

VERSION = "integer-frontier-v1"


def load_verified_environment(path, L, model, expected_sha256):
    """Load the actual frozen arrays; seeded regeneration is not proof of identity."""
    path = Path(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError("original environment archive checksum mismatch")
    with np.load(path, allow_pickle=False) as data:
        return DimerEnvironment(L, data["horizontal"], data["vertical"], model=model)


def randbelow_bits(rng, bound: int) -> int:
    """Unbiased integer in [0,bound), assuming rng.getrandbits gives fair bits."""
    if not isinstance(bound, int) or bound <= 0:
        raise ValueError("positive integer bound required")
    if bound == 1:
        return 0
    bits = (bound - 1).bit_length()
    while True:
        value = rng.getrandbits(bits)
        if value < bound:
            return value


def bernoulli_ratio(rng, numerator: int, denominator: int) -> bool:
    """No Float64 division; even a probability smaller than 2^-1074 survives."""
    if (
        not isinstance(numerator, int)
        or not isinstance(denominator, int)
        or not 0 <= numerator <= denominator
        or denominator <= 0
    ):
        raise ValueError("require integer 0 <= numerator <= denominator > 0")
    return randbelow_bits(rng, denominator) < numerator


def face_ratio(environment, i, j):
    """Exact ac/(ac+bd) for the represented input weights, without products in float."""
    a, b, c, d = map(Fraction, environment.face_weights(i, j))
    probability = a * c / (a * c + b * d)
    return probability.numerator, probability.denominator


def coupled_integer_update(lower, upper, environment, rng):
    """One new-version common heat-bath map; draws a coin even on null faces.

    Uniform face selection and the state-independent rational coin preserve the
    existing heat-order monotone coupling. This is a correctness reference, not
    an accelerated or checkpoint-compatible replacement for the old CFTP code.
    """
    if lower.L != environment.L or upper.L != environment.L:
        raise ValueError("inconsistent grid sizes")
    side = environment.L - 1
    i, j = divmod(randbelow_bits(rng, side * side), side)
    n, d = face_ratio(environment, i, j)
    orientation = "horizontal" if bernoulli_ratio(rng, n, d) else "vertical"
    for state in (lower, upper):
        if state.face_orientation(i, j) is not None:
            state.set_face_orientation(i, j, orientation)
    return i, j, orientation


class IntegerFrontierSampler:
    """Exact weighted rectangle recursion, restricted to even squares <=12.

    Scan vertices in row-major order. Bit zero of mask says whether the current
    vertex is already matched from an earlier vertex. If it is free, its only
    possible new partners are right/down. Cache the total integer weight of
    completions Z(position,mask), then choose branches proportional to w*Z.
    Width is exponential: this class deliberately rejects larger squares.
    """

    def __init__(self, environment: DimerEnvironment, *, max_states=2_000_000):
        if environment.L > 12:
            raise ValueError("integer frontier reference restricted to L<=12")
        if max_states <= 0:
            raise ValueError("positive state budget required")
        self.environment = environment
        self.L = environment.L
        arrays = [environment.horizontal_weights, environment.vertical_weights]
        ratios = [float(w).as_integer_ratio() for a in arrays for w in a.flat]
        denominator = max(d for _, d in ratios)  # All denominators are powers of two.
        weights = [n * (denominator // d) for n, d in ratios]
        common = 0
        for w in weights:
            common = gcd(common, w)
        weights = [w // common for w in weights]
        split = environment.horizontal_weights.size
        self.horizontal = np.array(weights[:split], dtype=object).reshape(
            arrays[0].shape
        )
        self.vertical = np.array(weights[split:], dtype=object).reshape(arrays[1].shape)
        self.edge_scale = Fraction(denominator, common)
        self._visited = 0

        @cache
        def partition(position, mask):
            self._visited += 1
            if self._visited > max_states:
                raise RuntimeError(
                    "frontier state budget exhausted; no sample returned"
                )
            if position == self.L * self.L:
                return int(mask == 0)
            if mask & 1:
                return partition(position + 1, mask >> 1)
            return sum(
                weight * partition(position + 1, newmask)
                for _, weight, newmask in self.branches(position, mask)
            )

        self.partition = partition
        self.total = partition(0, 0)
        if self.total <= 0:
            raise ValueError("zero partition function")

    def branches(self, position, mask):
        i, j = divmod(position, self.L)
        result = []
        if j + 1 < self.L and not mask & 2:
            result.append(("horizontal", int(self.horizontal[i, j]), (mask | 2) >> 1))
        if i + 1 < self.L:
            result.append(
                (
                    "vertical",
                    int(self.vertical[i, j]),
                    (mask >> 1) | (1 << (self.L - 1)),
                )
            )
        return result

    def sample(self, seed: int):
        """Fresh independent versioned seed per replica, never an old CFTP continuation."""
        rng = random.Random(int(seed))
        h = np.zeros((self.L, self.L - 1), bool)
        v = np.zeros((self.L - 1, self.L), bool)
        mask = 0
        for position in range(self.L * self.L):
            if mask & 1:
                mask >>= 1
                continue
            choices = [
                (orientation, weight * self.partition(position + 1, nextmask), nextmask)
                for orientation, weight, nextmask in self.branches(position, mask)
            ]
            target = randbelow_bits(rng, sum(weight for _, weight, _ in choices))
            i, j = divmod(position, self.L)
            for orientation, weight, nextmask in choices:
                if target < weight:
                    (h if orientation == "horizontal" else v)[i, j] = True
                    mask = nextmask
                    break
                target -= weight
            else:
                raise AssertionError("invalid categorical draw")
        state = DimerState(self.L, h, v)
        state.validate()
        return state

    def matching_probability(self, state):
        """Exact rational mass; used to compare ALL tiny tilings with enumeration."""
        state.validate()
        if state.L != self.L:
            raise ValueError("grid mismatch")
        weight = 1
        for w in self.horizontal[state.horizontal_occupied]:
            weight *= int(w)
        for w in self.vertical[state.vertical_occupied]:
            weight *= int(w)
        return Fraction(weight, self.total)
