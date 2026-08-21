# Gamma-disordered Aztec-diamond simulations

This Julia project samples weighted Aztec-diamond tilings, records a central
face height, runs shared-environment double-dimer experiments, and compares
finite-size variance-growth models. It is self-contained and uses only Julia
standard libraries.

## Model and observables

For an order-`L` sample, draw independent weights

```text
a[i,j] ~ Gamma(alpha, 1),   alpha = 0.2,
b[i,j] ~ Gamma(beta, 1),    beta  = 0.25,
```

with the other two edge-weight families gauge-fixed to one. This is Definition
1.1 of Duits and Van Peski,
[*The Gamma-disordered Aztec diamond*](https://arxiv.org/abs/2512.03033).
The implementation uses their reduction and creation rules (1.22)-(1.23).

The binary tiling matrix has size `2L x 2L`, with one `true` entry per domino.
The staggered face-height table has size `(2L + 1) x (L + 1)`. Its top boundary
is

```text
A[1,j] = 2j - 2,
```

and moving down a face column adds `+3` across a recorded dimer edge and `-1`
otherwise. The observable is the face nearest the geometric centre:

```julia
(row = L + 1, column = fld(L, 2) + 1)
```

Production code integrates along that one column and never allocates the full
height table.

### Spatial increments and the publication experiment

The super-roughness literature predicts spatial height-difference covariances.
For two replicas in one Gamma environment, this project measures

```text
delta_a(r) = H_a(x+r) - H_a(x),  a = 1,2,
T(r) = Var(delta_1(r)-delta_2(r))/2,
D(r) = Cov(delta_1(r), delta_2(r)).
```

`T(r)` is the connected/conditional contribution and `D(r)` is the disorder
contribution. Symmetric pairs on the central row are recorded at `r/L` equal
to `1/32`, `1/16`, `1/8`, and `1/4`. A separate uniform campaign applies the
identical geometry with no shared environment.

### Single height

One observation uses a fresh environment, one conditional tiling, and its
central height `H`.

### Shared-environment double dimer

One paired observation uses a fresh environment and two conditionally
independent tilings with central heights `H1` and `H2`. Conditional
independence gives

```text
Var(H1-H2)/2 = E[Var(H | environment)],
Cov(H1,H2)   = Var(E[H | environment]).
```

The difference isolates conditional tiling noise; the paired covariance
directly measures the disorder contribution.

## Installation and tests

Julia 1.10 or newer is supported.

```bash
julia --project=aztec -e 'using Pkg; Pkg.test()'
```

The tests include hand-computed recurrence values, bit-for-bit agreement
between the probability and pre-drawn-coin implementations, full geometric
tiling validation at small orders, direct-vs-full height agreement,
reproducibility, double-dimer pairing, and the finite-sample variance identity.

## Run a small end-to-end experiment

Single height:

```bash
JULIA_NUM_THREADS=4 julia --project=aztec \
  aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_smoke.csv \
  --output-dir aztec/output/height_smoke

julia --project=aztec aztec/scripts/merge_height_batches.jl \
  --inputs aztec/output/height_smoke \
  --output aztec/output/height_smoke.csv

julia --project=aztec aztec/scripts/analyze_height_campaign.jl \
  --results-dir aztec/output/height_smoke.csv \
  --output-dir aztec/output/height_smoke_analysis \
  --bootstrap-reps 200 --min-order 4
```

Double dimer:

```bash
JULIA_NUM_THREADS=4 julia --project=aztec \
  aztec/scripts/run_double_dimer_campaign.jl \
  --config aztec/configs/double_dimer_smoke.csv \
  --output-dir aztec/output/double_smoke

julia --project=aztec aztec/scripts/merge_double_dimer_batches.jl \
  --inputs aztec/output/double_smoke \
  --output aztec/output/double_smoke.csv
```

Campaign rows specify an order, first sample ID, sample count, and atomic batch
size. Every sample seed is derived only from the campaign seed, order, and
sample ID. Complete batches are validated and skipped on rerun; malformed or
seed-inconsistent batches stop the campaign rather than being overwritten.

## Reproduce the retained analyses

```bash
# Single height.
julia --project=aztec aztec/scripts/analyze_height_campaign.jl \
  --results-dir aztec/data/height/center_height_samples.csv \
  --output-dir aztec/output/height_analysis \
  --bootstrap-reps 10000 --bootstrap-seed 20260802 --min-order 24

julia --project=aztec aztec/scripts/plot_height_campaign.jl \
  --analysis-dir aztec/output/height_analysis \
  --output-dir aztec/output/height_analysis

# Shared-environment pairs and variance decomposition.
julia --project=aztec aztec/scripts/analyze_double_dimer_campaign.jl \
  --paired-results aztec/data/double_dimer/pairs.csv \
  --single-results aztec/data/height/center_height_samples.csv \
  --output-dir aztec/output/double_analysis \
  --bootstrap-reps 10000 --bootstrap-seed 20260802 --min-order 24

julia --project=aztec aztec/scripts/plot_double_dimer_campaign.jl \
  --analysis-dir aztec/output/double_analysis \
  --output-dir aztec/output/double_analysis
```

## Run the publication-focused spatial experiment

The complete Gamma campaign, uniform control, bootstrap analysis, and figures
are resumable through one command:

```bash
JULIA_NUM_THREADS=8 sh aztec/scripts/run_spatial_publication_campaign.sh
```

The primary nested comparison is

```text
H0: a + b log(r)
H1: a + b log(r) + c (log(r))^2.
```

The analysis reports paired-environment bootstrap intervals for `c`, BIC for
the nested models, cutoff sensitivity across four relative separations, and
prediction error on the two largest held-out orders. A joint analysis preserves
the shared environments across separations and estimates one common quadratic
coefficient with fraction-specific intercepts and log slopes. Both unweighted
and inverse-bootstrap-variance-weighted fits are reported. The uniform control
checks that the workflow remains consistent with ordinary logarithmic growth.

## Fitted curves and statistical scope

The scripts compare

```text
a + b log(L),
a + b (log(L))^2,
C (log(L))^p.
```

The first two are unweighted ordinary-least-squares fits to one estimated
variance or covariance per size. Their reported BIC is a heuristic Gaussian
residual comparison; it is not a likelihood fitted to the raw heights. The
third curve is fitted after taking logarithms and has no additive constant, so
its effective exponent `p` is descriptive and is not equivalent to either
affine model.

Bootstrap resampling is stratified by order. For paired data, `(H1,H2)` stays
together in every resample. Results should always be reported with bootstrap
intervals and checked across lower-order cutoffs. None of these finite-size
fits is an asymptotic proof.

## Implementation notes

The sampler follows deletion, sliding, and creation. For the Gamma model it
retains only two rolling `L x L` floating-point weight buffers. Potential
creation decisions are packed as bits and replayed during the forward shuffle.
For a single order-`L` tiling:

- runtime is `O(L^3)`;
- rolling weight and final-tiling storage are `O(L^2)`;
- stored creation decisions are `O(L^3)` bits.

For two conditional tilings, the expensive weight reduction is performed once
and two independent bit sequences are drawn from the same reduced environment.
Comments in [`src/AztecDiamond.jl`](src/AztecDiamond.jl) connect each code stage
to the corresponding mathematical operation. A longer, file-by-file explanation
is available in [`docs/IMPLEMENTATION.md`](docs/IMPLEMENTATION.md).

## Repository layout

```text
aztec/
├── src/AztecDiamond.jl       sampler, height, validation, rendering
├── test/runtests.jl          unit and mathematical-reference tests
├── docs/IMPLEMENTATION.md    detailed code and data-flow walkthrough
├── scripts/                  campaign, merge, analysis, and plot CLIs
├── configs/                  documented sample schedules
├── data/                     retained observations and checksums
├── results/                  final tables, figures, and analysis report
├── output/                   ignored resumable/scratch output
└── mathematica/              optional legacy matrix renderer
```

See [`configs/README.md`](configs/README.md), [`data/README.md`](data/README.md),
and [`results/README.md`](results/README.md) for artifact-level documentation.

## Square-grid paired Temperley experiment

The repository now also contains an exact square-grid spanning-tree/dimer
pipeline. It preserves the original directed site-i.i.d. random environment,
adds a separate undirected-conductance comparison, samples complete weighted
spanning trees with Wilson's algorithm, constructs the complementary dual tree,
and obtains the generalized Temperley perfect matching.

The spatial observable is an exact dimer-height increment across a fixed
central horizontal cut. See [`docs/SQUARE_GRID_MODEL.md`](docs/SQUARE_GRID_MODEL.md)
for the graph, outer-face, matching, height, and randomness conventions.

Tiny local smoke run:

```bash
JULIA_NUM_THREADS=4 julia --project=aztec \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model baseline \
  --config aztec/configs/square_grid_smoke.csv \
  --output-dir aztec/output/square_grid_baseline_smoke \
  --base-seed 20260810

JULIA_NUM_THREADS=4 julia --project=aztec \
  aztec/scripts/run_square_grid_campaign.jl \
  --environment-model directed_site_iid \
  --distribution gamma --parameter 0.5 \
  --config aztec/configs/square_grid_smoke.csv \
  --output-dir aztec/output/square_grid_directed_smoke \
  --base-seed 20260811
```

The smoke schedule has only two sizes and is intentionally not fitted. The
batch schema is compatible with `analyze_spatial_campaign.jl`; analyse the
five-size pilot with the Hamilton wrapper `hpc/analyze_square_grid.slurm` after
both baseline and directed arrays finish.

Hamilton Slurm wrappers and transfer instructions are in `hpc/README.md`.
