# Next Steps

> This is the **current action file**. Keep it short and update it whenever the research direction changes.
> The latest supplied supervisor reply takes precedence over older plans.

## P0 — Validate the square-grid model/observable before doing more large simulations

### Goal

Determine whether the negative square-grid result is genuine or is caused by the way additional randomness is transferred through the spanning-tree/Temperley representation.

### Required work

1. **Audit the current square-grid estimator.**
   - Identify exactly which random object is sampled.
   - Identify what is fixed for a "shared environment".
   - Identify which two samples are conditionally independent.
   - Trace how the height/spatial-increment observable is obtained from the spanning-tree/LERW representation.
   - Verify that the implemented quantity matches the intended \(T(r)\)/\(D(r)\) decomposition.

2. **Revisit the refreshed-weight random-walk simulation.**
   Sunil's latest instruction is to retry the random-walk/random-environment construction with a new random weight at each step rather than a fixed environment, but now using the newer observable/analysis that motivated the Aztec result.

3. **Do not guess the pairing convention.**
   The disorder covariance \(D(r)\) requires a meaningful shared source of disorder, while fully refreshed weights remove a fixed spatial environment. Before a production run, make the estimator mathematically explicit:
   - what random variables are shared between the pair?
   - what are conditionally independent?
   - what expectation/covariance is being estimated?
   - what is the null/control?

   If the existing code and correspondence do not determine this unambiguously, surface the ambiguity before spending Hamilton time.

4. **Build a reference implementation first.**
   - smallest possible sizes;
   - slow/clear code is acceptable;
   - deterministic seeds;
   - direct checks against identities/decompositions;
   - compare optimized and reference observables exactly on small instances.

5. **Pilot before scale.**
   Only after the estimator is verified:
   - use a small set of representative sizes;
   - include the no-disorder control;
   - include at least the original Gamma setting;
   - inspect effect size and diagnostics before launching a full campaign.

## P1 — Statistical validation of the P0 pilot

For any paired observable:
- preserve environment-level blocks;
- jointly resample all \(r/L\) values from the same environment;
- compare \(a+b\log L\) with \(a+b\log L+c(\log L)^2\);
- report \(c\), uncertainty, BIC difference, and prediction error;
- repeat across several lower-size cutoffs;
- report all windows, not only the best one;
- verify the sign convention for BIC;
- keep conditional \(T(r)\), disorder \(D(r)\), and total variance separate.

Decision gate:

**If a robust positive square-grid disorder-covariance signal appears**
→ reproduce with independent seeds, expand sizes/disorder laws, and compare directly against Aztec using the same analysis code.

**If it remains null**
→ stop the current approach and prepare the evidence for Sunil. Do not silently move to a different dimer dynamics.

## P2 — Active: square-grid dimer Glauber dynamics reference and pilot

The latest supplied supervisor email authorizes this independent route and
specifies the initial local heat-bath rule.  The refreshed-weight check was
reported null, so this is now the active square-grid direction.

The model is:
- start from a valid dimer configuration;
- draw and freeze independent positive mean-one weights on dimer edges for one
  environment;
- choose square faces uniformly at random;
- a flip is possible when opposite edges of the face are occupied;
- rotate the occupied pair to the other opposite pair;
- for clockwise top-first weights `(a,b,c,d)`, resample the top/bottom pair
  with probability `ac/(ac+bd)` and the right/left pair with probability
  `bd/(ac+bd)`;
- record the central-face height after an empirically validated burn-in and
  spacing schedule.

This is intended as an independent check that does not rely on the
spanning-tree correspondence.  The reference implementation uses the supplied
tileable height boundary.  Its event-driven accelerator exactly skips
self-loops while preserving literal attempted-update time; it has been checked
against the literal kernel on small cases.

Completed first milestone:
- deterministic reference sampler in `aztec/src/GlauberSquareGrid.jl`;
- exact enumeration, dimer-height invariants, and detailed-balance tests at
  tiny sizes;
- paired same-environment chains and trace diagnostics;
- exact parallel tempering in the unchanged frozen environment, with the
  physical observable retained only at `beta=1`;
- a persistent exchange clock and alternating-pair schedule, including a
  regression test that ensures the `beta=1` neighbour is attempted when trace
  spacing is shorter than the exchange interval;
- local calibrated pilots at `L=2,3,4,6`: maximum standardized extremal-start
  gaps `1.53--2.02`, minimum chain ESS `28.9`, and minimum adjacent-pair swap
  acceptance at least `0.227` for the Gamma(shape=0.5) ladder;
- direct `L=2` comparison to exact enumeration, with total-variation distance
  at most `0.00551` across eight frozen Gamma environments;
- all-one control ladder, with maximum standardized extremal-start gap at most
  `2.19` and minimum chain ESS `178`.

A fresh-seed local repetition showed that the first calibrated schedule was not
uniformly reliable: one `L=2` environment had exact total-variation distance
`0.0682`, and one `L=6` environment had standardized extremal-start gap `4.11`.
A focused 10x stress run of those same environments reduced these diagnostics
to `0.00997`/`0.44` at `L=2` and `0.11` absolute/`0.82` standardized start gap
at `L=6`.  The longer schedule is therefore the current pilot candidate, not
a validated production schedule.

Hamilton job `18546215` then ran the longer schedule with independent base seed
`2026082107`.  All eight tasks completed.  Across `L=2,3,4,6`, maximum
standardized extremal-start gaps were `1.95, 1.30, 1.57, 1.54`; minimum chain
ESS values were `1163, 131, 265, 147`; and the largest `L=2` exact
total-variation distance was `0.00533`.  All adjacent beta pairs were attempted,
with acceptance `0.2275--0.9594`; the target-adjacent pair acceptance was
`0.7485--0.9594`.

Next gate before any production size-scaling:
- run a four-environment `L=8` pilot with the same frozen-environment pairing
  and a conservative longer schedule;
- retain raw traces, seeds, per-pair exchange diagnostics, and resource use;
- require agreement from extremal starts and acceptable ESS/exchange flow at
  `L=8`, lengthening or densifying the beta ladder if diagnostics deteriorate;
- do not treat the current eight-environment-per-size Hamilton results as an
  equilibration proof or a production certificate.

Owner authorization update, 21 August 2026:
- Alberto explicitly authorized production without further supervisor review
  because Sunil was leaving on vacation and marked the former `AGENTS.md` gate
  outdated;
- the `L=8` pilot remains running as an additional diagnostic rather than a
  prerequisite that silently changes the production estimand;
- production must retain the central-height trace, two independent chains per
  frozen environment, deterministic seeds, full exchange diagnostics, atomic
  batches, and environment-blocked analysis;
- unvalidated larger sizes remain a scientific risk and must be staged with
  conservative schedules and monitored rather than represented as proven mixed.

Production submission, 21 August 2026:
- core Gamma(shape=0.5, scale=2) job `18546281`, base seed `2026082201`:
  `L=2,4,6,8,10,12`, 64--192 environments per size, 576 restart-safe tasks;
- frontier Gamma job `18546369`, the same base seed: `L=16,20`, 64 and 32
  environments, 96 isolated one-environment tasks;
- all-one control job `18546370`, base seed `2026082202`: `L=4,6,8,10,12,16`,
  32--64 environments, 240 tasks;
- production Gamma uses 21 inverse temperatures from 0 to 1, while the control
  uses the exact accelerated target kernel; every environment retains two
  independent extremal-start chains and 2,000 central-height values per chain;
- total concurrent allocation is capped at 256 one-core tasks and every task
  has a 13.5-hour limit;
- production analysis is environment-blocked in
  `aztec/scripts/analyze_glauber_square_grid_production.jl`.

Production completion and analysis, 22 August 2026:
- all 912 production tasks completed with exit code `0:0`; the output roots
  contain the expected 2,739 atomic files, no temporary files, and no non-empty
  error logs;
- the final data contain 960 Gamma environments and 352 all-one controls, with
  two 2,000-draw chains per environment;
- deterministic environment-block bootstraps and nested
  `a+b log(L)` versus `a+b log(L)+c(log(L))^2` comparisons are implemented in
  `aztec/scripts/analyze_glauber_square_grid_scaling.jl`;
- over the full Gamma range `L=2--20`, unweighted estimates of `c` are
  `-0.048 [-0.320,0.227]` for the conditional component,
  `0.108 [-0.545,0.794]` for the disorder component, and
  `0.060 [-0.721,0.822]` for the total component; all point-estimate BIC
  differences favor ordinary log, and ordinary log has lower leave-one-size-out
  RMSE;
- inverse-bootstrap-variance fits, lower-size cutoffs, and removal of
  environments with standardized extremal-start gap above 4 do not produce a
  robust positive coefficient;
- the all-one disorder covariance remains numerically near zero, providing the
  intended negative control.

Current decision:
- do not interpret the direct Glauber campaign as positive super-rough
  evidence and do not launch more production automatically;
- retain the all-environment analysis as primary and the start-gap-filtered
  result only as a diagnostic sensitivity analysis;
- before any larger-size replication, add replica-label flow or round-trip
  diagnostics and address the `beta=0--0.05` bottleneck (pooled acceptance
  `0.062` at `L=20`); the present target-adjacent exchange rate is healthy, but
  it does not by itself certify traversal of the full ladder;
- present the current result as a finite-size null with a mixing caveat when
  the project is next reviewed.

## P3 — Cross-model consistency

Once the square-grid sampler is validated:

1. Put Aztec and square-grid outputs into one common schema:
   - model;
   - geometry;
   - \(L\);
   - disorder law/parameters;
   - environment seed;
   - sample seed(s);
   - \(r/L\);
   - \(T(r)\);
   - \(D(r)\);
   - run/version metadata.

2. Run the **same fitting implementation** for both geometries.

3. Compare:
   - sign and magnitude of \(c\);
   - cutoff stability;
   - disorder-strength dependence;
   - negative controls;
   - conditional vs disorder components.

4. Separate three possible explanations for the current discrepancy:
   - genuine geometry/model difference;
   - finite-size effect;
   - mismatch in random-environment/dimer correspondence or estimator.

## P4 — Reproducibility and code release preparation

Sunil has said the code and commentary would need to be available for a possible numerics paper.

Before any release:
- deterministic/replayable seeds;
- one config file per campaign;
- exact package/environment lockfile;
- restart-safe jobs;
- raw-data manifest and checksums;
- scripts that regenerate every table/figure from raw outputs;
- concise README for reproduction;
- tests for mathematical invariants and optimized observables;
- no hidden manual data edits.

## P5 — Manuscript boundary

Do **not** use Codex/AI to draft journal manuscript prose. Sunil explicitly said AI systems should not be used for the write-up stage.

Codex may be used here for:
- repository inspection;
- implementation;
- tests;
- numerical diagnostics;
- reproducibility tooling;
- code documentation;
- internal research planning.

## Recommended first Codex task after repository context is added

Do **not** ask Codex to "analyse everything".

Give it a scoped audit:

> Read `AGENTS.md`, `context/RESEARCH_BRIEF.md`, `context/RESULTS.md`, and `context/NEXT_STEPS.md`.
> Locate the current square-grid paired/Temperley experiment and the earlier refreshed-weight LERW implementation.
> Do not edit anything yet.
> Trace the random variables and estimator used in each implementation, and determine whether the latest proposed refreshed-weight experiment has a mathematically well-defined analogue of the paired disorder covariance \(D(r)\).
> Return only:
> 1. relevant files/functions;
> 2. random-variable dependency diagram in text;
> 3. estimator currently computed;
> 4. any mismatch with `RESEARCH_BRIEF.md`;
> 5. the smallest implementation plan needed to resolve it.
> Do not inspect unrelated repositories/files.

That audit should happen before spending API budget on a large implementation.
