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
tileable height boundary and literal random-face updates; it does **not** yet
use an active-site acceleration, whose invariant distribution would need a
separate rejection-free argument.

Completed first milestone:
- deterministic reference sampler in `aztec/src/GlauberSquareGrid.jl`;
- exact enumeration, dimer-height invariants, and detailed-balance tests at
  tiny sizes;
- paired same-environment chains and trace diagnostics.

Next gate before any size-scaling or Hamilton work:
- calibrate burn-in, thinning, autocorrelation/ESS, and agreement from the
  extremal starts over a small sequence of sizes;
- retain raw traces and all seeds;
- include the all-one control and at least one mean-one Gamma environment;
- compare sampled tiny-grid distributions directly to exact enumeration.

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
