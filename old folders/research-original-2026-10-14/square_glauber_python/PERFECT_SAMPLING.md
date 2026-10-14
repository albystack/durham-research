# Certified perfect sampling for the frozen random-bond dimer model

## Scope

This implementation samples the same finite square-grid perfect matchings, in
the same frozen positive edge environment, as the validated Glauber code.  It
does not introduce a different Gibbs law or transition kernel.  Its random map
is one ordinary update attempt:

1. choose uniformly from all `(L-1)^2` bounded faces;
2. do nothing when the face is not flippable;
3. otherwise choose horizontal with `ac/(ac+bd)` and vertical with
   `bd/(ac+bd)`.

The method is monotone coupling from the past (CFTP).  “Perfect” means that a
returned state has no finite burn-in bias: the code returns it only after all
possible past starting states have been certified to give the same time-zero
state.  It is still a numerical implementation using PCG64 pseudo-randomness,
IEEE-754 weights and probabilities.

## Height order and true bounds

All height fields use the existing exterior-zero/fixed-cut convention.  At a
fixed face, the heights of any two legal tilings differ by a multiple of four.
Define

`h <= g` if and only if `h(f) <= g(f)` at every bounded face `f`.

For domino tilings of a simply connected region, normalized height functions
form a finite distributive lattice.  Local flips change one face by `+4` or
`-4`; repeatedly applying all available downward (respectively upward) flips
therefore reaches the unique minimum (respectively maximum).  The rectangle in
this project is simply connected.  `extrema.py` constructs these two states and
checks that no directed flip remains.

This distinction matters: the convenient all-horizontal and all-vertical
tilings are not generally the lattice bounds.  Exact enumeration already
refutes that assumption on 4x4.  They remain useful mixing starts, but CFTP uses
the newly constructed true minimum and maximum.

References:

- Propp and Wilson, *Exact sampling with coupled Markov chains and applications
  to statistical mechanics* (1996), especially the tiling-height discussion:
  <https://www.math.cmu.edu/~af1p/Teaching/MCC17/Papers/propp_wilson.pdf>
- Desreux and Rémila, *An optimal algorithm to generate tilings* (2005), on
  flip accessibility and the distributive lattice of domino tilings:
  <https://www.sciencedirect.com/science/article/pii/S1570866705000092>

## Why the common update is monotone

Only the selected face height can change.  A perfect matching can show exactly
seven occupancy masks around one face: no occupied edge, one of four single
occupied edges, the horizontal opposite pair, or the vertical opposite pair.
The last two are the flippable cases.

For two ordered states, compare the selected face and its four neighbours.  The
two selected-face heights differ by `0, 4, 8, ...`.  A single update moves either
height by at most four, so gaps larger than eight are automatically safe.  The
code exhausts both face parities, all 49 pairs of local masks, gaps `0, 4, 8`,
and both possible common heat-bath outcomes, retaining precisely the locally
ordered cases.  All 496 admissible cases preserve order.  This is a finite,
domain-independent local proof, rather than an extrapolation from a few global
simulations.

The frozen environment is essential.  At a given face,

`p_f = ac/(ac+bd)`

depends only on its four fixed weights, not on the matching.  Thus both
trajectories use the same threshold whenever the face is flippable in both.
The local truth table also covers when only one or neither state is flippable.

## CFTP algorithm

For a requested deterministic CFTP seed:

1. assign independent PCG64 face-index and heat-bath-uniform streams to a
   replayable past block;
2. start at the true minimum and maximum at time `-T`;
3. apply exactly the same stored random maps to both bounds through time zero;
4. if their occupancies agree at time zero, return that common state;
5. otherwise prepend an independently seeded block, double the past horizon,
   restart at the true bounds and replay every later block unchanged.

Monotonicity sandwiches every legal starting state between the two bounding
trajectories.  Therefore coalescence of the true bounds implies coalescence of
the entire state space.  Propp--Wilson CFTP then gives the unique stationary
law.  The validated detailed-balance calculation identifies that law as

`pi_omega(M) proportional to product(w_e for e in M)`.

The state space is finite.  Rectangular domino tilings are connected by local
flips, so the positive-weight heat-bath chain is irreducible.  It is aperiodic
because an orientation can be redrawn unchanged (and nonflippable face choices
are null moves).  These facts supply the remaining CFTP prerequisites.

If the maximum horizon is reached without coalescence, the sampler raises
`CFTPNotCoalesced` and returns no state.  Such an environment or seed must be
extended using the same replayable history; it must not be discarded or
silently replaced.

Two calls with distinct CFTP seeds give two independent samples conditional on
the same immutable environment.  That is the correct future input to the
shared-environment covariance analysis.  Coalescence horizons are diagnostics,
not themselves independent equilibrium observations or automatic mixing-time
estimates.

## Validation commands

Install the optional backend in a supported Python environment:

```bash
python -m pip install '.[test,accel]'
```

Then run:

```bash
python -m square_glauber.cli local-monotonicity-audit \
  --output outputs/perfect_sampling_validation/local_monotonicity

python -m square_glauber.cli exact-perfect-sampling \
  --samples 20000 --seed 20260909 \
  --output outputs/perfect_sampling_validation/exact_weighted_4x4

python -m square_glauber.cli cftp-hard-calibration \
  --source-case outputs/mixing_validation_20260908_extended/gamma-L16-diagnostic-02 \
  --samples 50 --seed 20260909 \
  --output outputs/perfect_sampling_validation/hard_L16

python -m square_glauber.cli perfect-pair-smoke \
  --sizes 4,8,12,16 --seed 20260909 \
  --output outputs/perfect_sampling_validation/pair_smoke
```

These are sampler validation and computational calibration commands.  They do
not fit `log L` or `(log L)^2`, and they do not start a scientific scaling run.

The subsequent paired production implementation lives in
`perfect_campaign.py`. It stores the frozen environment before sampling,
persists each certified state separately, and publishes a paired scientific
row only after two distinct CFTP streams certify. See the Hamilton campaign
section of [`PROJECT_GUIDE.md`](PROJECT_GUIDE.md) for the frozen manifest and
launch procedure.
