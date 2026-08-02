# Implementation walkthrough

This document is a code-oriented companion to the comments in
`src/AztecDiamond.jl`. It explains what each stage stores, which randomness it
uses, and which checks protect the retained data.

## 1. Matrix encoding

An order-`L` tiling is stored as a `2L x 2L` Boolean matrix. A `true` entry is
the centre/parity label of one domino, not a board cell. Row and column parity
determine one of four orientation classes. This encoding is compact, makes the
shuffle local on `2 x 2` blocks, and matches the retained example matrices.

`validate_tiling` converts every true entry to the two geometric cells covered
by its domino. It then checks all four invariants:

1. exactly `L(L+1)` dominoes;
2. no cell covered twice;
3. no missing cell of the Aztec diamond;
4. no cell outside the Aztec diamond.

Full geometric validation uses sets and is intended for tests and illustrated
examples. Production samples use the same tested shuffle path and perform the
cheaper exact domino-count check.

## 2. Gamma environment

`gamma_disordered_weights` draws two independent `L x L` arrays:

```text
a[i,j] ~ Gamma(alpha, 1),
b[i,j] ~ Gamma(beta, 1).
```

The dependency-free Gamma generator is the Marsaglia-Tsang sampler. Shapes
below one use the standard shape-boosting identity before the rejection step.
The tests check its mean and variance under a nontrivial shape/scale choice.

## 3. Weight reduction

`reduce_gamma_weights!` implements equation (1.22) of Duits-Van Peski. An
order-`k` environment becomes order `k-1`. Production code owns two full-sized
buffers and alternates their roles, writing only the leading `(k-1) x (k-1)`
block at each step. Stale values outside that block are never read.

There are two front ends:

- `gamma_disordered_probabilities` retains every reduced Float64 probability
  table. It is transparent and useful for small examples and reference tests.
- `gamma_disordered_creation_choices` immediately draws Bernoulli decisions
  and stores them as packed bits. It keeps only the two rolling weight buffers
  and is used for large campaigns.

The tests compare these front ends bit for bit with the same RNG stream, as
well as checking an order-two recurrence by hand.

## 4. Creation decisions

At a reduced face with weights `a,b,1,1`, equation (1.23) chooses the two
possible domino pairs with probabilities `a/(a+b)` and `b/(a+b)`. In this
matrix encoding, `true` denotes the latter pair.

Decisions are drawn for every potential creation site before shuffling. During
the forward shuffle, some sites are not holes and their bits are ignored. This
does not change the distribution: independent unused Bernoulli variables can
be sampled and discarded without conditioning the variables that are used.

For `copies=2`, each reduced environment is visited once and two distinct RNG
draws are made per site. The copies therefore share the environment but no
creation coin.

## 5. Domino shuffling

For each level from 2 through `L`:

1. Embed the old tiling with a one-cell border.
2. Delete opposite domino pairs that would collide.
3. Slide every surviving domino according to its parity class.
4. Fill every empty face using its pre-drawn creation bit.

The live tiling is a `Matrix{Bool}` rather than a packed `BitMatrix`. Scalar
updates dominate this stage, and byte-addressed Booleans are faster. Creation
tables remain packed because their total size is cubic.

## 6. Height observable

`height_function` constructs the full staggered face table and is used as a
reference implementation. `center_height` starts from the known top boundary
and integrates down only the central column. Tests compare both values at many
small orders.

The direct computation still requires the final tiling, but avoids allocating
the additional `(2L+1) x (L+1)` integer height table.

## 7. Random streams and reproducibility

Campaign scripts convert `(base_seed, order, sample_id)` to a 64-bit seed with
SplitMix64. SplitMix64 is only a deterministic hash here; Xoshiro is the RNG
used by the simulation.

This design has three consequences:

- thread scheduling cannot change a sample;
- changing batch size cannot change a sample;
- an interrupted campaign can resume without regenerating earlier batches.

Before skipping an existing batch, the runner checks its row count, schema,
orders, sample IDs, deterministic seeds, centre indices, and (for pairs) stored
height differences. A mismatch stops the run.

## 8. Atomic batches and merging

A batch is written to `batch_XXXX.csv.tmp` and renamed only after the file is
closed. A final batch filename therefore never represents a partial write.

Merge scripts recursively select only `batch_*.csv`, validate every row, reject
duplicate `(order,sample_id)` keys and duplicate seeds, sort canonically, and
atomically replace the merged output.

## 9. Bootstrap analysis

Single observations are resampled independently within each order. Double
observations are resampled as intact `(H1,H2)` pairs. Each bootstrap replicate
recomputes the size-level variance or covariance before refitting the curves.

The affine log and squared-log fits use unweighted ordinary least squares on
the size-level estimates. Their BIC values are therefore curve-comparison
diagnostics, not a full likelihood for raw height observations. The power fit
uses

```text
log variance = log C + p log(log L)
```

and is evaluated only when every fitted bootstrap estimate is positive. The
fit report records how many bootstrap replicates met that requirement.

## 10. File map

- `src/AztecDiamond.jl`: mathematical and sampling kernels.
- `scripts/run_*campaign.jl`: resumable threaded data generation.
- `scripts/merge_*batches.jl`: strict canonical data assembly.
- `scripts/analyze_*campaign.jl`: bootstrap summaries and fits.
- `scripts/plot_*campaign.jl`: dependency-free SVG rendering.
- `test/runtests.jl`: kernel, recurrence, geometry, and identity tests.
- `test/smoke_workflows.sh`: every public CLI connected end to end.
