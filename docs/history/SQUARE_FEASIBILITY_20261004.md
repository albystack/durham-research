# Wider-square feasibility, 4 October 2026

**Both stages complete and verified: all 18 retained environments now have paired samples.** The next step towards
Sunil's comparison is to establish whether the same two-family Gamma square
model can be sampled at larger sizes. The 12,288-pair corrected Aztec refinement
is complete. This new square experiment is a fixed, bounded feasibility test,
not a scaling campaign. All 143 source files and six retained input archives
passed Hamilton checksum checks. Validation job 18902368 passed all 132 tests
in 14.61 seconds. Array 18902369 completed all 18 tasks: 12 complete pairs
and six cases that reached the declared state-count limit. All weights,
complete tilings and measurements were downloaded and verified.

## Mathematical target and implementation

The input law remains independent Gamma(0.2, scale 1) on north edges from black
vertices and Gamma(0.25, scale 1) on east edges, with other weights one. The
sampler's integer arithmetic is exact conditional on these stored Float64
weights and ideal random bits. The pseudorandom generator and initial Gamma
generation retain their usual numerical qualifications.

The new `integer-frontier-budgeted-v1` implementation preserves the reference
conditional partition sampler's law while reducing memory use:

1. At each black vertex b, choose a positive rational factor g_b so all its
   incident stored weights become integers, then divide by their common gcd.
   Every matching occupies exactly one edge incident to b. Thus
   `product(w'_e over M) = product(g_b over black b) * product(w_e over M)`.
   The multiplier is independent of the matching and cancels in its Gibbs
   probability. This is a computational representation of the unchanged
   environment, not mean-normalisation of the Gamma families.
2. Skip already occupied vertices without storing redundant partition states.
   Their edge weight was counted at the earlier endpoint. Pack the remaining
   position and frontier mask into one integer key.
3. Evaluate the same acyclic recursion using an explicit postorder stack.
   Conditional sampling still chooses branches proportionally to edge weight
   times the exact remaining partition function.

The old sampler and CFTP routines are unchanged. Gauge factors change integer
bit lengths and therefore random-bit consumption; trajectories from equal
seeds need not match across versions. A new sampler-seed namespace is used.
Width dependence remains exponential, and this implementation rejects L>18.

## Validation completed locally

All **132 Python tests passed**. The new tests compare every 4x4 matching mass
and all reachable conditional probabilities against the reference and rational
Gibbs products for all three models. They check subnormal positive weights,
exact partition equality at side 8, deterministic replay within the new
version, all three resource limits, original-input requirements, refusal to
overwrite output, result-hash tampering and missing-case reporting.

A separate same-environment side-12 check found identical exact partition
functions after accounting for gauge factors. Cached states fell from 104,284
to 52,341 and partition bit length from 6,430 to 4,311. Local setup times were
0.1245 s and 0.1139 s respectively. This single check is not a claim about
Hamilton memory use or wider-grid performance.

## Frozen cohort and stopping

| Square side | Random Gamma pairs | Mean-weight control pairs | Uniform control pairs | Inputs |
|---|---:|---:|---:|---|
| 14 | 4 | 1 | 1 | New, saved before sampling |
| 16 | 4 | 1 | 1 | Original Hamilton archives, checksum verified |
| 18 | 4 | 1 | 1 | New, saved before sampling |

The side-16 random cases are the first four retained original environments,
without replacing those that failed CFTP. New environment master seed:
2026100301. The retained side-16 environment seed namespace remains 2026092902.
All replica seeds are separately derived in the new sampler namespace.
The full case identities and original archive hashes are in
`square_glauber_python/configs/frontier_scaling_20261003.json`.

Each fresh worker stops partition construction at 120 seconds, 6 GiB peak
process RSS, or 3,000,000 cached states. RSS and time are checked periodically;
these are soft limits. Slurm independently limits the worker to four minutes,
one CPU and 8 GiB. The 18-task array permits two concurrent workers, for at
most 72 allocated CPU-minutes. A separate validation job has a ten-minute
limit. There is no automatic cap increase, replacement environment or retry.

All cases write status and retain their input weights. Budget exhaustion is
reported separately from software errors and scheduler interruption. Only
complete pairs receive sample files and measurements. Measurements use the
existing symmetric physical diagonal probes, rho=1/8 and 1/4, storing actual
rounded endpoints and distances. This narrow size range and tiny cohort do
not support curvature fitting; repeated rounded separations are expected.

## Running and checking

The isolated snapshot must first pass
`square_glauber_python/hpc/frontier_validation.slurm` on Hamilton. Then submit
`square_glauber_python/hpc/frontier_scaling.slurm` with `FRONTIER_OUTPUT`,
`ORIGINAL_CASES` and `RESEARCH_PYTHON` set to the recorded Hamilton locations.
The original case root is
`/nobackup/fvkl37/matched_comparison_20260929/feasibility/cases`.
Use a new output root; the worker refuses to overwrite a case directory.
Submission state and source hashes are recorded in
`docs/validation/frontier_scaling_20261003/`.

After retrieving all case files, run:

```bash
PYTHONPATH=square_glauber_python/src .venv/bin/python square_glauber_python/scripts/summarize_frontier_scaling.py --plan square_glauber_python/configs/frontier_scaling_20261003.json --input LOCAL_CASE_ROOT --output docs/validation/frontier_scaling_20261003/results.json
```

The verifier checks plan identity, original input hashes, all complete-pair
file hashes, matching constraints and recomputed observables. It retains
missing and interrupted cases in the denominator. Scheduler outcomes must
also be checked: a scheduler success can accompany a recorded resource cap.

## Decision after this bounded test

Use completion, time and memory for every retained case to assess a useful
square size range. Completion of this small cohort does not establish runtime
tails or authorise a large production run. Keep the matched production hold
until feasible sizes, distinct physical separations, controls, sample counts
and stopping rules support the direct domain comparison.

If frontier width is the limiting factor, a separate route worth validating
is generalized shuffling with forbidden edges. Janvresse, de la Rue and
Velenik describe weighted square-lattice subgraph sampling by embedding in an
Aztec diamond with a tileable complement and handling zero weights explicitly.
Using an Aztec host would be an algorithmic representation; proving that the
extracted target keeps the square Gibbs law is essential. This route has not
been implemented or substituted into the experiment.
[Primary paper](https://www.unige.ch/math/folks/velenik/papers/JRV06err.pdf).


## First-stage outcomes and explicit follow-up decision

| Side | Complete pairs | Random-case worker time | Largest process peak RSS | Cached states |
|---|---:|---:|---:|---:|
| 14 | 6/6 | 1.47–1.73 s | 0.236 GiB | 268,759 |
| 16 | 6/6 | 7.53–7.66 s | 0.798 GiB | 1,333,766 |
| 18 | 0/6 | 16.87–17.52 s before cap | 1.247 GiB | 3,000,000 cap |

The random-case times include worker setup, sampling and output, not scheduler
startup. All six side-18 outcomes were `BUDGET_EXHAUSTED: state_limit`; no case
hit the time or memory soft limits. This is not evidence that side 18 cannot
complete within those limits. The first-stage result, with six censored cases,
is preserved in `results.json` and is not rewritten as a fully completed run.

After reviewing these outcomes, a separate six-case follow-up was frozen in
`side18_followup_plan.json`. It uses the exact saved side-18 archives, the same
sampler seeds and unchanged source code. Only the state cap rises to eight
million; the 120-second partition budget, 6 GiB RSS soft limit and four-minute,
8 GiB scheduler ceiling are unchanged. Array 18902387 runs at most two workers
at once, with at most 24 allocated CPU-minutes. This is an explicit second-stage
feasibility decision, not an automatic retry rule or a replacement of hard
inputs. All six follow-up pairs completed and passed full local verification.
Their outcomes remain separate from the original plan.


## Verified side-18 outcome and present conclusion

All six follow-up pairs completed at **6,448,885 cached states** per environment.
The four Gamma cases took 39.31–39.75 seconds of worker time; maximum peak RSS
across all six cases was **4.356 GiB**. The two controls also completed. The
first-stage and follow-up scheduler records all show exit 0:0, with empty
error logs. The follow-up result archive SHA256 is
`bb7bffe65b6002a8ced9d760cde0a9e6e34d001e4430a9456e4582588e35bfb8`.

Across the two stages, 18 distinct retained environments have 36 verified
sampled matchings and 72 recomputed observable values. All input identities,
file certificates and recomputed measurements passed. No incomplete case was
replaced, and the original six capped outcomes remain in the first-stage
report. Records are in `side18_followup_results.json` and
`completion_verification_20261004.json` beside the first-stage records.

This establishes practical exact sampling for this small cohort through side
18, including the original side-16 Gamma inputs. It does not establish tails
across a large disorder ensemble. Cached-state count increases by about a
factor of five for each increment of two in side length in this range; the
frontier algorithm therefore remains unsuitable as an assumed route to large
squares. No further cap increase or larger campaign was submitted.

The next development priority is a scalable exact square sampler, with
explicit geometry and Gibbs-law checks against the frontier oracle. The
zero-weight generalized-shuffling route described above is a concrete
candidate; it needs implementation, tiny exact-distribution checks and
numerical validation before any experiment. CFTP remains a valid alternative
if a correctly coupled acceleration can resolve the demonstrated bottleneck.
The direct matched-domain curvature comparison is still outstanding.
