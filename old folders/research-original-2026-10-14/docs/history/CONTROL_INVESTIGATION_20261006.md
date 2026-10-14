# All-one control investigation and completed confirmation, 6–7 October 2026

**The saved-data audit is complete. No pairing or geometry defect has been
identified by these checks. All 38 original pairs replayed exactly and Hamilton validation job 18919384
passed. All 12,288 fresh confirmation pairs and their analysis are complete.
The primary contrast is 0.081098 [−0.281718, 0.437180], with permutation
p=0.6615: the earlier anomaly was not reproduced at this precision. No original observation,
coefficient or bootstrap interval has been changed.**

## What the anomaly can and cannot mean

For deterministic weights, let X1=f(U1) and X2=g(U2), where the random inputs
are independent. Then E[X1 X2]=E[X1]E[X2], hence their covariance is zero.
This is true even if f and g have imperfect marginal sampling distributions.
A deterministic probe offset, an incorrect common height scale, or a bias in
the marginal Gibbs sampler alone therefore does not generate a nonzero
expected cross-replica covariance. Finite-sample chance, replica dependence,
record association errors or analysis errors are relevant possibilities.
The finite control result is not, by itself, proof that the square sampler
samples the wrong Gibbs distribution.

The original all-one primary curvature contrast remains 1.160615, with
95% bootstrap interval [0.450981,1.853177]. The prespecified maximum-coordinate
permutation test remains p=0.5882. The exploratory six-control-curvature test
remains p=0.0093. These use different statistics; neither result is removed.
The independent confirmation below was fixed after these old results, before
any confirmation samples were generated. Its data were analysed separately.

## Completed checks on existing data

- Rechecked the complete measurement tables with their hashes, environment
  identities, two-replica associations and seed derivations.
- Checked all **108,000 role seeds** across the original sampling/validation
  and new confirmation/validation cohorts: no seed is reused. Distinct seeds
  do not constitute a mathematical proof of pseudorandom-stream independence.
- Reconstructed both physical probe endpoints relative to each domain centre.
  They agree exactly at all six sizes, for both probes. The one-eighth probe
  has odd separation at sizes 8 and 24 and midpoint (−0.5,+0.5) in both
  geometries. All other chosen probes have midpoint (0,0). Parity remains a
  possible finite-size feature of Gamma fits, not a geometry mismatch or an
  explanation for nonzero expected deterministic-control covariance.
- Reconstructed the curvature contrast as a linear combination of all 24
  uniform-control covariance coordinates. Every size contributes positively
  to this particular realised contrast; the largest contribution is L=64.
- Deleted each uniform-control pair in turn, moving its two probes together,
  purely as an influence diagnostic. The resulting contrast ranges from
  **1.111913 to 1.194972**. The largest individual change is **0.048702**.
  No pair was actually deleted from the experiment or the reported analysis.
- Exhausted all 24 permutations of a four-pair, two-probe example to check the
  exact permutation covariance formula used in the exploratory statistic.
  Agreement is within 1.78e−15, including cross-probe covariance terms.
- Reviewed sampler state and RNG use. Square replicas use separate seeded
  local `random.Random` objects with read-only use of the reduction tables.
  Aztec replicas use disjoint successive random draws and distinct choice
  arrays within one seeded Xoshiro stream. Under ideal independent random
  draws this construction gives independent replicas. It is not a proof
  about all behaviour of the finite pseudorandom generators.

| Size | Contribution to the original uniform domain contrast |
|---|---:|
| 8 | 0.170899 |
| 16 | 0.162342 |
| 24 | 0.295944 |
| 32 | 0.040279 |
| 48 | 0.027393 |
| 64 | 0.463758 |

The source, raw diagnostic tables and explicit endpoint coordinates are in
`docs/validation/control_investigation_20261006/`. Lag summaries there are
exploratory descriptions, not a new significance-testing programme.

## Replay before confirmation

Replay the original saved uniform pairs with indices 1,128,256 at every
size in both domains, plus the two most influential size-64 pairs (square
183 and Aztec 46): **38 pairs** in total. The square replay generates replica
2 before replica 1 and verifies exact matching-array equality against the
saved originals; reduction tables must remain unchanged. The Julia replay
uses the original paired-coin seed and checks both matching matrices and
distinct choice arrays. This tests record provenance and possible state
contamination, not statistical correctness of all random streams.

The compute-node validation also generates a separate four-pair smoke cell
at each of 12 confirmation cells, revalidates its raw files and exercises the
confirmation analysis. Local smoke validation already passed for all 48
pairs, 96 matchings and 192 observable values. Validation data are excluded
from scientific inference. All 155 original source files remain unchanged
and hash-verified; the nine added investigation files have their own manifest.

## Prospectively fixed independent confirmation

| Choice | Fixed value |
|---|---|
| Model | All edges one; no new random Gamma environments |
| Geometries | Square and Aztec |
| Sizes | 8,16,24,32,48,64 |
| Pair count | 1,024 per size/domain; 12,288 pairs total |
| Sampler | The unchanged, validated original campaign implementation |
| Production master seed | 2026100601 |
| Validation master seed | 2026100602 |
| Primary statistic | c_Aztec−c_square over sizes 8–64 |
| Primary test | Two-sided whole-vector pairing permutation against zero |
| Permutations | 9,999; seed 2026100605 |
| Bootstrap | 10,000; square seed 2026100603, Aztec seed 2026100604 |
| Sensitivity fits | 16–64 and 24–64, descriptive |
| Resource ceiling | 12 tasks, at most two concurrent; one CPU, 4 GiB, 45 min/task |
| Stop rule | Complete the fixed cohort once; no replacement or extension based on results |

The randomisation null permutes each second-replica two-probe vector within
its own size/domain, preserving marginal distributions and within-replica
spatial dependence. The primary contrast is the single confirmatory test.
The old data will not be pooled into it. A fresh p<0.05 flags a reproducible
issue needing deeper investigation. Otherwise report that the anomaly was
not reproduced at this precision; do not claim that independence or sampler
correctness has been proved. Both outcomes will be retained and reported.

The plan SHA256 is
`0ced6e1eefa0ccd541570804d183a0f2e1cebd6c603ef1bf6753d58d2dd70867`.
The tools archive SHA256 is
`3b8cd3c80b7eec956ccc6c846ec7bd78ac7fdeb4d4c4f46babc205604c3e88f9`.
Remote root: `/nobackup/fvkl37/control_confirmation_20261006`.
The original root remains `/nobackup/fvkl37/matched_shuffle_20261005`.

No new Gamma scaling campaign is part of this confirmation. After the result,
review the control evidence and finite-size window with Sunil. Historical AI
attribution and independent human mathematical review remain outstanding.


## Validation and launch record

Hamilton validation **18919384** completed in 57 seconds with exit 0:0 and
empty stderr. It replayed all 19 square and 19 Aztec pairs exactly, including
reverse-order square generation and unchanged reduction tables. All 48
separate smoke pairs passed raw verification and the confirmation-analysis
pipeline. The compact validation archive was downloaded and checked (SHA256
`05d09eaf0ec0dbfd66e01b867b60aa0c6036094f4dbf7e37f298a7865a3a1de1`).

Fresh confirmation array **18919405** completed all 12 tasks with exit 0:0.
Dependent job **18919413** completed raw verification and analysis in 24 minutes
31 seconds, also with exit 0:0. All 13 error logs are empty. The slowest
sampling cell (square side 64) took 21 minutes 36 seconds. The scheduler
records analysis completion at 2026-10-06 14:03:21; records were retrieved and
checked locally on 7 October. No pair was dropped, replaced or added.

## Confirmation results

The predeclared primary contrast is the difference between the fitted common
squared-log coefficients, Aztec minus square. All intervals below are the
original pointwise 95% whole-pair bootstrap intervals (10,000 draws).

| Window | Aztec c [95% CI] | Square c [95% CI] | Difference [95% CI] |
|---|---|---|---|
| 8–64, primary | 0.039 [−0.213, 0.290] | −0.042 [−0.297, 0.208] | 0.081 [−0.282, 0.437] |
| 16–64, descriptive | −0.026 [−0.649, 0.593] | −0.006 [−0.592, 0.579] | −0.020 [−0.867, 0.829] |
| 24–64, descriptive | 0.934 [−0.605, 2.448] | −1.082 [−2.557, 0.378] | 2.016 [−0.108, 4.104] |

The sole confirmatory, two-sided permutation test gives **p=0.6615**, using
9,999 permutations and the predeclared seed. Its exact conditional null
standard deviation is 0.182312. The earlier 1.161 [0.451,1.853] contrast was
not reproduced at this precision. This is compatible with sampling variation
but does not identify the original cause or prove independence or sampler
correctness. The broad 24–64 interval remains visible; its descriptive role
has not been promoted to a new primary test. Neither cohort has been pooled,
reweighted or extended on the basis of its result.

Hamilton reloaded **24,576 matchings**, recomputed **49,152 observable values**
and checked **55,296 files** before analysis. Local verification checks the
archive and table hashes, the fixed plan, complete pair/probe coverage and
seed identities, then independently recomputes 72 moments, six fitted
coefficients and all three contrasts. It also reconstructs the primary
permutation statistic and repeats its fixed 9,999 permutations. Bootstrap
intervals are retained from Hamilton rather than independently rerun.

Compact local evidence is under
`square_glauber_python/outputs/control_confirmation_20261006/`; full raw
tilings remain at the remote root. The archive SHA256 is
`b3f453ebcd9d2514b85fd5c1142af15aaf4fc3912770fc86850cb02b9d29cd23`.
The independent check is reproducible with
`OPENBLAS_NUM_THREADS=1 .venv/bin/python docs/validation/control_investigation_20261006/verify_confirmation.py`.
Its output is `confirmation_verification.json` alongside the script.

## Decision after confirmation

The fixed control follow-up is finished. Preserve both cohorts and report the
non-reproduction alongside the original anomaly. The matched Gamma result
remains **−0.046 [−0.366,0.282]**; the control result does not turn it into
evidence of different growth laws.

The next scientific step is review with Sunil of the completed comparison,
its physical size window and the mathematical sampler validation. A larger
matched window could be a future experiment, but should have a fresh fixed
design and a precision/runtime justification before submission. No further
campaign has been submitted. The supervisor email remains an unsent draft.
Alberto's historical AI-attribution account and independent human mathematical
review remain outstanding; they cannot be supplied by an automated audit.
