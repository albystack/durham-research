# Methods and reproducibility

## Aztec diamonds

`aztec/src/AztecDiamond.jl` implements weighted domino shuffling, height
observables, paired replicas, and SVG rendering. The package uses Julia
standard libraries only. The [implementation notes](../aztec/docs/IMPLEMENTATION.md)
record the reduction and creation rules.

An order `L` diamond is encoded by a `2L × 2L` binary tiling matrix. The
Gamma-disordered model draws its two nontrivial weight families with shape
parameters `alpha` and `beta`; the other families are gauge-fixed to one.
The illustrated quick start uses `alpha=0.2`, `beta=0.25`, and seed `20260920`.

The README diamond instead uses uniform creation probabilities, order 28,
and seed `20260920`. This keeps the introductory figure easy to read.

## Square-grid Glauber dynamics

The Python model uses an even `L × L` grid of vertices. Every undirected edge
has one fixed positive weight. Each attempted update chooses uniformly from
all `(L−1)²` bounded faces, including nonflippable faces.

For a flippable face with clockwise top/right/bottom/left weights `a,b,c,d`,
the probability of the horizontal orientation is `ac/(ac+bd)`. One sweep is
exactly `(L−1)²` attempted updates. Matching edges are the state; heights are
derived when measured. The [Python reference](../square_glauber_python/PROJECT_GUIDE.md)
specifies the exterior reference, signs, geometry, and diagnostic commands.

Coupling from the past (CFTP) returns a perfect sample only after the bounding
chains coalesce. A capped run is censored, not a sample. Keep its environment,
seed, and nested random-map history when extending it. See the
[perfect-sampling contract](../square_glauber_python/PERFECT_SAMPLING.md) and
[pre-specified analysis plan](../square_glauber_python/PAIRED_COVARIANCE_ANALYSIS_PLAN.md).

The Julia `GlauberSquareGrid` module is retained for independent validation
and Kasteleyn calculations. Its `L` describes a `2L × 2L` region, whereas the
Python `L` is the vertex side length. Do not pool their data by the size label
alone.

The README square-grid figure uses all-one weights, an all-horizontal start,
update seed `20260921`, and 500 sweeps. It is a dynamics snapshot, with no
claim of equilibrium or perfect sampling.

## What is estimated

Let `H₁` and `H₂` be conditionally independent heights (or height increments)
in the same environment `ω`. Then

\[
\tfrac12\operatorname{Var}(H_1-H_2)
=\mathbb E_\omega[\operatorname{Var}(H\mid\omega)],
\qquad
\operatorname{Cov}(H_1,H_2)
=\operatorname{Var}_\omega(\mathbb E[H\mid\omega]).
\]

The first quantity measures conditional fluctuations; the second measures
disorder-induced fluctuations. Bootstrap whole environments, keeping both
replicas and all spatial measurements together. Compare ordinary-log and
log-squared fits with recorded fit windows and uniform-weight controls.

## Reproducing a run

Record the code revision, language and dependency versions, full command,
configuration, seed, model parameters, and input checksums. Keep environment
and sample identifiers unique when extending or merging campaigns. Preserve
failed/censored tasks and mixing diagnostics alongside successful results.

`make test` runs the active Python and Julia suites. `make validate-python`
checks exact tiny-grid identities. CI runs both projects. These checks do not
certify that an arbitrary production MCMC run has mixed.

HPC scripts remain with their projects: [Julia](../aztec/hpc/) and
[Python](../square_glauber_python/hpc/). Julia Slurm wrappers are submitted
from the repository root after `bash aztec/hpc/setup_hamilton_environment.sh`;
they retain the existing Hamilton-specific Julia path. Python launch details
are in its command reference. Local setup and examples do not submit jobs.
