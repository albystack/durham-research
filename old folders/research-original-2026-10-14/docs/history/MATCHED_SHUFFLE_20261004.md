# Matched-weight square/Aztec comparison, prepared 4–5 October 2026

**Status, 5 October: production and analysis completed successfully. All 30,720 pairs are verified. The primary curvature difference is inconclusive and a control-curvature discrepancy needs follow-up. See [the completed results](MATCHED_RESULTS_20261005.md).**

Sunil's preferred comparison uses the original Aztec law on the ordinary
square: independent north-from-black weights Gamma(shape 0.2, scale 1),
east-from-black weights Gamma(shape 0.25, scale 1), and all other edges one.
The two random families therefore have means 0.2 and 0.25. They are not
renormalised to mean one. The earlier all-edge Gamma(shape 1) square data
answer a different question and are not substituted into this comparison.

## Removing the square sampling bottleneck

The exact frontier method validated the law at small widths but grew too
quickly for a useful scaling window. `embedded-square-shuffle-v1` instead
isolates the square in an Aztec host of order L−1. Original square weights are
unchanged; crossing edges are exactly zero and exterior edges are one. An
explicit exterior perfect matching is checked. Every allowed host matching
splits into a square matching and an exterior matching, so

Z_host = Z_square × Z_exterior.

Consequently, projecting a weighted host sample onto the square gives the
square's original Gibbs law. The host is an algorithmic device: it does not
replace the square's boundary conditions or change its weight distribution.

Weighted domino shuffling reduces active cells with weights a,b,c,d using
cell factor ac+bd and opposite-edge weights divided by that factor. Growth
moves lone old edges to their opposite edges, deletes opposite old pairs,
and fills old empty cells with probabilities proportional to ac and bd.
The implementation rejects a singular zero cell, arithmetic underflow or
nonfinite arithmetic. It does not implement the paper's general singular-cell
rules. For nonsingular reductions, the positive-weight identity extends by
continuity to the zero crossing weights. See [Janvresse, de la Rue and
Velenik, *A Note on Domino Shuffling*, section 2](https://www.unige.ch/math/folks/velenik/papers/JRV06err.pdf).

Rational arithmetic is an exact small-grid oracle (conditional on ideal
random bits and stored input weights). The scalable implementation uses
Float64 reductions followed by exact rational creation coins for their
represented weights. It is **not** unqualified exact sampling for the original
weights: reduction rounding is assessed separately. Neither implementation
uses a mixing-time or burn-in assumption. No singular case is repaired by
adding epsilon.

Local evidence:

- Exact equality with the independent Gibbs distribution for every matching
  at L=2 and L=4 for random Gamma and both controls; an additional arbitrary
  all-edge-weight L=4 comparison agrees with the integer frontier oracle.
- All 6,728 uniform L=6 matchings have exactly probability 1/6,728.
- A fixed four-environment audit compares all potential creation probabilities
  at L=8,16,32,64 with 80- and 160-digit reductions. At L=64, maximum Float64
  probability discrepancy is 1.79e−15; twice the sum across all potential
  coins is 7.94e−12. This is finite-cohort evidence, not a universal error bound.
- All 146 Python tests pass. A separate end-to-end local check completed four
  pairs in each of 36 domain/model/size cells (144 pairs); these are validation
  data and excluded from production. Raw verification is recorded separately.
- The fixed L=64 audit took about 0.51 s to construct a Float64 sampler and
  0.26 s for two draws locally, excluding cold geometry construction. Actual
  Hamilton timing and memory must be checked before submission of production.

## Fixed scientific design

Plan: `square_glauber_python/configs/matched_shuffle_20261004.json`.
Its SHA256 is `b1d4ee08dbe49d794558ef9badc61bc66f722c5b79594722de67300eb3108e70`.

| Choice | Fixed value |
|---|---|
| Domains | Ordinary even square and fresh corrected Aztec arm |
| Sizes | 8, 16, 24, 32, 48, 64 in each domain's stated convention |
| Random environments | 2,048 per size per domain |
| Deterministic controls | Family means 0.2/0.25 and all-one; 256 pairs per size/domain/model |
| Total | 30,720 paired environments, including controls |
| Replicas | Two independently sampled tilings per environment |
| Probes | Symmetric diagonal increments, rho=1/8 and 1/4; height scale 1 |
| Physical distances | sqrt(2) × size × rho; same in both domains |
| Production master seed | 2026100407, separate from every earlier cohort |
| Validation master seed | 2026100408; four pairs per cell, excluded from inference |
| Primary window | 8–64 |
| Sensitivity windows | 16–64 and 24–64 |
| Bootstrap | 10,000 whole-environment resamples, independent by domain/model/size |
| Scheduler ceiling | 36 tasks, two concurrent, one CPU, 4 GiB and two hours each |

For odd probe separation, the pre-existing fixed midpoint tie-break is
retained. Actual distances, rather than the requested fraction alone, enter
fits. Different geometries have different vertex counts; equal numeric size
means the documented domain convention and equal probe separation, not equal
area.

Estimate thermal T=Var(X1−X2)/2 and disorder D=Cov(X1,X2). Fit
D_p(d)=a_p+b_p log(d)+c(log(d))², with a separate intercept and slope for each
probe and a shared quadratic coefficient within each domain. The primary
comparison is **delta_c=c_Aztec−c_square**, with a bootstrap confidence
interval for that difference. Comparing separate significance decisions would
not test the domain contrast. The deterministic controls have true disorder
covariance zero; one joint replica-permutation audit covers their two probes,
both domains, models and all sizes (9,999 permutations, seed 2026100411).

Report all predetermined windows. A square interval containing zero does not
prove absence of a squared-log term, and a positive contrast on these sizes
does not prove a boundary-driven asymptotic theorem. No optional stopping,
environment replacement or outcome-dependent window selection is specified.
Every declared environment and both probes must be present before fitting.

## Commands and provenance

`run_matched_shuffle.py` dispatches one plan cell. Square inputs and both
matchings are saved as compressed arrays with checksums; Julia saves Aztec
families, matching tables and checksums. A fresh output directory is required;
interrupted work remains incomplete and is never treated as an equilibrium
cohort. Original files and input seeds are retained for investigation.

The source snapshot excludes local runtimes, outputs, caches and the locally
resolved Julia Manifest. The Hamilton environment resolves the project's
standard-library dependencies. `RESEARCH_JULIA_PROJECT` can point to the
existing `.runtime/aztec-checks` project for local smoke checks; its default
on Hamilton is the snapshot's `aztec` directory.

1. Restore SSH access using the current shared socket.
2. Upload the frozen source archive and verify every source checksum.
3. Submit `square_glauber_python/hpc/matched_shuffle_validation.slurm`. It runs
   both language validation suites, the numerical audit, all 36 four-pair
   smoke cells, complete-table checks and raw matching/observable verification.
4. After successful validation and resource review, submit
   `square_glauber_python/hpc/matched_shuffle.slurm`. This is the production
   campaign already requested by the user; no result is assumed in advance.
5. On completion, run `verify_matched_shuffle_raw.py` on Hamilton and
   `analyse_matched_shuffle.py` on verified complete tables. Preserve a compact
   archive with status, manifests, table hashes, raw verification and scheduler
   exits. Update the handbook, figures and unsent supervisor draft.

Prepared plans are not job submissions. Submission IDs, source archive digest
and actual outcomes belong in `docs/validation/matched_shuffle_20261004/`.


## Hamilton launch, 5 October 2026

Renewed access succeeded. All 155 files in the frozen source snapshot were
checksum-verified on Hamilton, both before validation and again at the launch
gate. Source archive SHA256:
`6a12fc01b420838b55517da1920f9a9f885d181e12953d4d2a042fc95d755efa`.

Validation job **18907657** completed with exit 0:0 in **2 min 56 sec**:
146 Python tests, 21,970 Julia assertions, both numerical audits, all 36
four-pair integration cells, 288 raw matching validations and 576 observable
recomputations. Its stderr was empty. The job initially waited because its
20-minute request exceeded the test partition's 15-minute maximum; the live
job time limit was reduced to 15 minutes, without changing scientific code or
production budgets. Future resubmissions of the frozen validation script must
include `--time=00:15:00` on the `sbatch` command line.

The L=64 Gamma cell took 5.14 seconds for four pairs, including sampler setup,
file writes and validation, with peak RSS 348,090,368 bytes. This provides
reasonable headroom under the two-hour, 4 GiB production limits; it is a small
benchmark, not a universal runtime guarantee. Peak RSS for the complete
validation job, including high-precision audits, was 914,252 KiB.

Production array **18907659** retains all planned 30,720 pairs and runs at most
two tasks concurrently. Post-processing job **18907660** has dependency
`afterok:18907659`, so it starts only after every array task succeeds. It first
revalidates raw inputs, matching certificates and observables, then runs the
prespecified bootstrap and control audit. A failed task prevents automatic
analysis; an incomplete cohort is never silently fitted.

Remote root: `/nobackup/fvkl37/matched_shuffle_20261005`.
The separate post-processing script is recorded here alongside its checksum;
it invokes the existing frozen verification and analysis code. Validation
records and scheduler accounting have been downloaded and checked. The next
step is to retrieve complete production results and inspect the direct
curvature contrast and all fixed sensitivity checks. SSH is not required for
submitted jobs or their dependency to continue.
