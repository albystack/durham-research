# Research roadmap

This file records the ordered scientific plan. Completed experiments and
qualifications belong in [`RESULTS.md`](RESULTS.md).

## Current position

The Aztec-diamond spatial experiment produced positive finite-size
quadratic-log curvature in the disorder-induced covariance, while the
connected component and uniform control did not. Earlier structured
square-grid Temperley experiments and direct central-height calculations did
not show a robust comparable signal.

The active question is therefore a like-for-like direct random-bond
square-grid experiment using the same shared-environment decomposition and
spatial height increments. The current Python implementation now avoids an
ordinary-Glauber burn-in judgement by using monotone coupling from the past
with the true extremal height tilings.

## P0 — Complete the frozen perfect-sampling campaign (active)

Frozen design:

- `L = 8,12,16,20,24,28,32`;
- 18,200 Gamma(shape=1, scale=1) environments;
- 5,600 all-one negative-control blocks;
- two independent CFTP-certified samples per environment;
- centre height, global height sum, and pre-specified spatial increments;
- master seed `2026090901`;
- whole-environment bootstrap seed `2026090902`.

Completion gate:

1. all primary Slurm tasks finish;
2. every capped stream is retained and replayed with the same weights, seed,
   and nested random-map history at a larger horizon;
3. the campaign audit reports `all_pairs_cftp_certified=true`;
4. manifests, hashes, task counts, and scheduler outcomes agree;
5. no environment is discarded because it is computationally hard.

The live run is external to Git under the recorded Hamilton `/nobackup`
campaign root. The local launch record is
[`CAMPAIGN_LAUNCH_READINESS.txt`](../square_glauber_python/CAMPAIGN_LAUNCH_READINESS.txt).

## P1 — Run the frozen analysis, once only

After P0 passes, run the pre-specified analysis in
[`PAIRED_COVARIANCE_ANALYSIS_PLAN.md`](../square_glauber_python/PAIRED_COVARIANCE_ANALYSIS_PLAN.md):

- disorder component `Cov(H1,H2)`;
- connected component `0.5 Var(H1-H2)`;
- the same quantities for every retained spatial separation;
- ordinary-log versus ordinary-plus-quadratic-log fits;
- 10,000 whole-environment bootstrap repetitions;
- fixed cutoff sensitivities only;
- uniform-control analysis through the identical pipeline.

Do not change the fit window, environment set, or primary observable after
looking at the result. A positive, zero, or negative quadratic coefficient is
reportable.

## P2 — Scientific audit and interpretation

Before any claim:

- check the all-one covariance is statistically consistent with zero;
- inspect residuals and size leverage;
- compare centre and spatial observables without replacing the primary one;
- report connected and disorder components separately;
- verify bootstrap resampling preserves complete environments;
- distinguish finite-size evidence from an asymptotic proof.

## P3 — Decide whether a larger campaign is warranted

Only after P2, choose one of:

- extend to larger `L` with a separately frozen design if certification cost
  and the signal-to-noise calculation are defensible;
- investigate a different exact sampler or analytic route if CFTP tails become
  prohibitive;
- stop and report a finite-size null if the direct square-grid observable does
  not support the proposed effect.

Do not select new sizes or disorder laws solely because they increase the
quadratic coefficient.

## P4 — Mathematical follow-up

The numerical project cannot itself prove \((\log L)^2\) growth. A proof would
require a separate analytic programme: identify the precise random-bond Gibbs
measure, derive disorder correlations of the height field, control the finite
boundary geometry, and establish the asymptotic covariance. Numerical results
can motivate and test such arguments but cannot substitute for them.

## Release checklist

- one public root README and working internal links;
- reproducible Python and Julia setup commands;
- exact campaign configuration and seed contracts;
- compact reviewed results under `results/`;
- no raw private correspondence or unlicensed papers;
- no credentials; historical machine usernames/paths reviewed and either
  deliberately retained as provenance or redacted before publication;
- chosen software license and citation metadata;
- full Python and Julia tests passing from a clean checkout.
