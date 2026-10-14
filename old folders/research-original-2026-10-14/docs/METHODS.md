# Model and analysis conventions

## Weight law and geometries

The completed experiment uses independent Gamma edge-weight families with
shape parameters `alpha=0.2`, `beta=0.25`, **scale one**, and the remaining
families fixed at one. Thus the random family means are 0.2 and 0.25.
Family placement is fixed by the corrected coordinate adapters, with `a`
on north and `b` on east edges in the documented square-cell convention.

For an environment ω, a perfect matching M has conditional probability
proportional to the product of its occupied-edge weights. Every vertex is
covered once. The environment stays fixed while two independent replicas
are generated, and is replaced for the next pair.

Aztec `size` is the diamond order; square `size` is the even vertex side length.
The graphs and their induced boundaries differ. At a common size label, the
selected diagonal probes have the same physical endpoints relative to the
centre and the same Euclidean separation. The two fractions are 1/8 and 1/4:
`r = size × fraction`, `distance = sqrt(2) × r` in primal-edge units.

Heights live on faces. Crossing an occupied edge changes height by magnitude
three; an empty edge gives magnitude one, with signs determined by the vertex
colour and crossing direction. The fixed boundary reference is common to both
replicas. Each observable is a symmetric spatial increment, in height units one.

## Sampling used in this campaign

Aztec tilings use weighted domino shuffling. A 29 September audit found an
orientation mismatch between the historical Gamma reduction and tiling arrays.
The corrected adapter reflects the family arrays into reduction coordinates,
then reflects creation tables back. Small exact-distribution tests verify this
mapping. Details are in [the orientation audit](history/MATCHED_COMPARISON.md).

The ordinary square region is embedded in a larger Aztec host. Crossing edges
are assigned zero weight, the outside region has a fixed compatible matching,
and unreachable zero-support states are pruned. The restriction to the square
therefore has the desired conditional product-weight law. The Python
implementation is `square_glauber/embedded_shuffle.py`.

Creation coins use robust arithmetic for represented weights. Reduction tables
use Float64 arithmetic. High-precision comparisons at sizes 96 and 128 passed
the fixed tolerances, but these finite checks are not universal rounding bounds.
The exact frozen source and its hashes are included in `data/`.

Glauber dynamics and coupling from the past remain tested alternatives in the
Python package; they were not the production samplers for this campaign. A
heat-bath face update selects the horizontal pair with probability
`ac/(ac+bd)` for clockwise top/right/bottom/left weights. Uniform selection
includes nonflippable faces. CFTP certifies only coalesced chains; capped runs
are censored, not samples.

## Variance decomposition

For two conditionally independent copies X₁ and X₂ of the same observable:

- `T = Var(X₁ − X₂) / 2 = Eω[Var(X | ω)]`.
- `D = Cov(X₁, X₂) = Varω(E[X | ω])`.
- `V = T + D` is the total variance.

The published finite-sample V is the average of the two unbiased replica
variance estimates. These identities separate contributions; they impose no
particular growth rate on either one. The scientific target here is D.

## Fixed fit and uncertainty

For each geometry, fit `Dρ(r) = aρ + bρ log(r) + c (log(r))²`. The two fractions
have separate intercepts and ordinary-log slopes, and share c. Logarithms are
natural. Fits use equal weights. The primary size window is 16–128; fixed
sensitivity windows are 32–128 and 48–128.

There are 16,384 Gamma pairs and 1,024 pairs for each deterministic control,
per size/domain, at seven sizes. All 258,048 pairs are retained. No pilot or
previous campaign is pooled into the data. Uncertainty uses 10,000 independent
whole-environment bootstrap draws within each domain/model/size cell. Each
resampled vector contains both probes and replicas. Intervals are pointwise
percentile intervals. The original additional domain contrasts are retained.

The joint control diagnostic uses 9,999 permutations of complete second-replica
probe vectors. It tests the maximum standardized absolute curvature contrast
across two controls and three fixed windows. All individual control results are
reported, including the all-one primary interval narrowly excluding zero.

## Reproduction and verification

`make analysis` verifies file hashes and dataset coverage, independently
recomputes all moments and coefficients, and uses the recorded intervals.
`make bootstrap` reruns the original resampling and control diagnostics.
See [the analysis guide](../analysis/README.md) and [data schema](../data/README.md).

Workers verified saved matching constraints and recomputed observables before
issuing completion receipts. All sampling tasks and final analysis succeeded.
Local checks reproduced 252 point estimates, 18 coefficients and nine
contrasts from the downloaded compact data. The source snapshot and batch
verification metadata are published; full tiling arrays remain on Hamilton.
