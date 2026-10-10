# Numerical remediation and sampler feasibility, 1 October 2026

**Updated 2 October:** The direct sampler has completed the fixed small-square
cohort using the original archived Hamilton weights. This is a new algorithm
and seed namespace; it does not retrospectively complete or alter the CFTP
campaign. Hamilton validation job **18872074 completed successfully** in one
minute: 120 Python tests, 21,970 Julia assertions, the numerical audit and the
four-pair benchmark passed. The approved corrected Aztec refinement was
submitted on 2 October as array **18897605**. All six tasks completed; the
3 October frozen analysis gives c=0.3568, 95% interval [0.0723, 0.6327].
See [the completed refinement report](AZTEC_REFINEMENT_20261003.md) for cutoff
sensitivity and the remaining square-domain feasibility gap.

## Exact arithmetic for the stored square environments

Each stored positive Float64 edge weight is an exact dyadic rational n/2^k.
Let D be a common power-of-two denominator. Multiplying every edge weight by
D makes integer weights. Every perfect matching has L^2/2 edges, so all matching
weights acquire the same factor. A further common integer factor can be
cancelled without changing the distribution.

`integer_sampling.py` uses arbitrary-size integer rejection sampling: draw k
fair bits where 2^k covers the integer range, reject integers outside that range,
and compare the accepted integer to the exact numerator. This preserves rare
outcomes even when their probabilities are below the resolution of an ordinary
53-bit floating uniform. Python's documented `getrandbits` supplies arbitrary
bit counts ([official documentation](https://docs.python.org/3/library/random.html#random.getrandbits)).
The proof is conditional on ideal fair bits; the implementation uses seeded
pseudorandom bits. It does not turn the initial Gamma generator into an exact
continuous-distribution sampler.

The new face ratio is exactly ac/(ac+bd) for the represented inputs, computed
as a rational without floating products or division. A new common heat-bath
map uses a uniformly selected face and one shared rational coin, including on
null faces. This is a correctness reference. It is not an accelerated replacement
inserted into the historical CFTP checkpoint format.

## Direct frontier sampler

Scan vertices in row-major order. A bit mask records vertices already matched
from earlier positions. At a free vertex, only its right and lower neighbours
can be new partners. Define the integer partition function of completions by

    Z(position, mask) = sum over allowed edges e of W_e * Z(next state).

At an already occupied vertex there is only the skip transition. At the end,
Z=1 for an empty mask and 0 otherwise. Cache these values, then choose each
branch with probability W_e Z(next)/Z(current) using integer random draws.
Multiplying these conditional probabilities telescopes to the matching's
integer weight divided by the total partition function. This gives the desired
Gibbs law conditional on stored edge weights, without burn-in or coalescence.

The implementation is restricted to even squares L<=12 and a two-million-state
budget. Width dependence is exponential; no scalability to L=32 is claimed.
It provides an independent small-width oracle and a possible route for modest
width experiments, not yet the large-square production solution.

## Fixed-cohort outcome

Original environment master seed: 2026092902. New sampler seeds are derived
from `sha256(integer-frontier-v1:original_chain_seed)`, with separate replica
seeds. Original environment seeds are retained, but new trajectories do not continue
old CFTP streams. A later byte-hash audit qualified environment identity, as
described below.

| Model | Side | Pairs completed | Total partition setup | Cached states/environment |
|---|---:|---:|---:|---:|
| Two-family random Gamma | 8 | 16/16 | 0.0504 s | 3,361 |
| Two-family random Gamma | 12 | 4/4 | 0.3310 s | 104,284 |
| Deterministic means and uniform | 4, 8, 12 | 6/6 total | See saved summary | Up to 104,284 |

**Original-array resolution, 2 October:** All 240 original Hamilton weight
archives were downloaded and verified against the pilot's recorded hashes.
Direct comparison found 176 locally regenerated environments identical and 64
with actual array-value differences, all in the random Gamma arm. Differences
were at most 9 ULP, with maximum relative difference 1.225e-15. This establishes
that the discrepancy is not solely archive packaging; the underlying numerical
cause has not been isolated. Earlier seeded-replay diagnostics retain this
qualification.

A separate rerun loaded the original archives directly and completed all 26
pairs: 20 random pairs, four overlapping controls, and two new side-4 controls.
All 24 overlapping input archives match their original hashes. All 52 saved
tilings passed matching validation, all file checksums passed, and all 224 saved
observable values were independently recomputed and matched. Original-input
partition setup took 0.0495 s in total for the 16 side-8 random cases and
0.3572 s for the four side-12 cases. These remain local timings.

The table above records the earlier regenerated-input run, retained under
`square_glauber_python/outputs/integer_feasibility_20261001/`. The authenticated
original-input rerun is separate, under
`square_glauber_python/outputs/integer_original_feasibility_20261001/`.
Its compact records are `integer_original_feasibility_summary.json` and
`integer_original_verification.json` in `docs/validation/numerical_20261001/`.
These small-cohort measurements have not been used to claim scaling behaviour.
The runner requires either `--original-cases` with checksum verification or an
explicit diagnostic-only `--allow-seed-regeneration` flag.

Tests include all matching masses on 2x2 and 4x4 grids versus independently
computed rational Gibbs products; exact detailed balance on every flippable
4x4 matching; scripted rare outcomes; integer rejection boundaries; common-map
height order; and deterministic replay. Timing excludes interpreter startup
and is local-machine evidence, not a Hamilton scaling estimate.

## Aztec coin repair and reduction audit

`robust_gamma_coins.jl` introduces `gamma-interval-coins-v1`. For each represented
reduced a,b pair, neighbouring Float64 values give outward bounds for b/(a+b).
A 53-bit uniform prefix determines an interval for U. If that interval is wholly
on one side of the probability bounds, the outcome is decided. Otherwise the
routine converts a,b to exact rationals and extends the SAME binary prefix in
64-bit blocks until integer comparisons decide U<p. It never restarts a uniform
conditionally on being close to the threshold.

The new routine uses the corrected array orientation and fresh versioned RNG
consumption. Scripted tests force rare outcomes near both zero and one; tests
also cover subnormals and sums overflowing Float64. The old sampling functions
and existing outputs remain intact for reproducibility.

This repairs Bernoulli rounding for the reduced weights. The reduction itself
still uses Float64. Three fixed environments were independently reduced in
256- and 512-bit arithmetic:

| Order | Maximum probability discrepancy: Float64 reduction vs 512-bit | Twice summed discrepancy | Old probabilities rounded to 0 or 1 |
|---|---:|---:|---:|
| 8 | 1.13e-16 | 3.40e-15 | 0 |
| 32 | 1.13e-15 | 6.10e-13 | 1 |
| 128 | 5.52e-15 | 1.30e-10 | 277 |

The 256-vs-512-bit maximum discrepancies are below 4.1e-76. The sum is a
coupling-error proxy across potential creation coins, not a rigorous global
interval bound or a guarantee for every environment. The audit tolerance was
max discrepancy <1e-10, twice summed discrepancy <1e-6 and reference discrepancy
<1e-50. These three environments pass. Reduction positivity/finite-value checks
remain enforced in each run.

## More precise corrected Aztec experiment: submitted 2 October

The fixed plan requests 2,048 paired environments at each order 16, 32, 48, 64,
96 and 128, using master seed 2026100102 and the new robust coins. It uses only
rho=1/8 and 1/4, which have complete distinct-probe coverage at every size.
Windows are 16–128, 32–128 and 48–128. Analysis uses 10,000 whole-environment
bootstrap resamples with seed 2026100103. No optional stopping or automatic
sample-count expansion is implemented.

A four-pair benchmark at each size passed the complete output/analysis pipeline.
Those 24 pairs are integration checks, not the refinement result. The order-128
invocation took 0.566 seconds inside the runner for four pairs on the local
machine; scaling this to Hamilton is uncertain. The planned Slurm array allows
two concurrent tasks, each one CPU, 4 GB and at most two hours. All source and
case-file hashes are checked on resume. A stale lock from a killed process must
be checked against scheduler state before manual removal.

The planned half-width 0.25 is a precision objective suggested by the earlier
pilot, not a guarantee. The refined data would improve the corrected benchmark;
it would not alone answer the matched boundary question. The direct square
sampler must first be assessed at wider grids or another validated exact
method developed. Generalized shuffling for weighted planar regions is one
literature route to investigate, not implemented or claimed here
([Janvresse, de la Rue and Velenik paper](https://www.unige.ch/math/folks/velenik/papers/JRV06.pdf)).

## Reproduction and remaining action

From the repository root:

```bash
MPLCONFIGDIR=tmp/matplotlib .venv/bin/python -m pytest square_glauber_python/tests -q
JULIA_DEPOT_PATH=.runtime/julia-depot .runtime/julia-1.10.10/bin/julia --startup-file=no --project=.runtime/aztec-checks aztec/test/runtests.jl
JULIA_DEPOT_PATH=.runtime/julia-depot .runtime/julia-1.10.10/bin/julia --startup-file=no --project=.runtime/aztec-checks aztec/scripts/audit_robust_gamma.jl
```

The integer feasibility runner refuses an existing output directory. For the
next authenticated replay supply `--original-cases ORIGINAL_CASES_DIRECTORY`,
`--pilot-status ORIGINAL_STREAM_STATUS.csv`, and `--output NEW_DIRECTORY`. It
loads the archived arrays after checking their original file hashes. The
isolated source snapshot has been uploaded and all 135 source hashes verified
on Hamilton. Validation job 18872074 passed, and
`aztec/hpc/robust_refinement.slurm` was submitted as array 18897605.
Download the six completed size directories
(or their compact tables, manifests and summaries) and run:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=square_glauber_python/src .venv/bin/python square_glauber_python/scripts/analyse_robust_refinement.py --input LOCAL_REFINEMENT_ROOT --output LOCAL_ANALYSIS_ROOT
```

The analysis refuses incomplete sizes, incorrect counts, hash mismatches,
geometry inconsistencies, wrong sampler/seed namespaces and repeated seeds.
Do not use `--validation-only` for scientific results.
