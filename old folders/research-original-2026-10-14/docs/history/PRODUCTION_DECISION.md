# Production decision, 30 September 2026

**Hold production.** This decision follows the bounded assessments in
`FOLLOWUP_20260930.md`. It does not authorise a new scheduler submission.

## Fixed scientific target

Compare disorder covariance of symmetric diagonal height increments for the
same independent two-family Gamma(0.2)/Gamma(0.25), scale-one edge law, remaining
edges one, and the documented physical orientation. Use unit height increments
and actual Euclidean separations. The old all-edge Gamma(1) data are exploratory
context, not the matched square arm. Preserve control arms: deterministic
family means and all-one weights.

The proposed primary fit uses probe-specific intercepts and log slopes, with
a common squared-log term; compare the two domain coefficients with independent
whole-environment bootstraps. Size labels are domain-specific. Fix physical
windows, probe coverage and domain-to-domain separation matching before future
outcome inspection. Do not equate Aztec order numerically with square side
without checking physical geometry and boundaries.

## Gates before committing resources

1. **Numerical validity:** audit rare-transition handling in both implementations
   over the realised weight range. Any new kernel must retain the stated Gibbs
   law, pass tiny exact-distribution checks, preserve monotonicity where needed,
   and have a new versioned random-stream/checkpoint identity. Silent clipping,
   weight normalisation, ordinary burn-in samples labelled perfect, or skipping
   slow environments are not acceptable substitutions.
2. **Square feasibility:** complete a fixed diagnostic cohort without replacing
   difficult cases, including all first-four side-8 original environments. Then
   test a separately fixed cohort at proposed production sizes and report all
   censored observations and runtime tails. With 0/80 pairs complete now, no
   reliable cost for completed production pairs is available.
3. **Corrected benchmark precision:** independently reassess corrected Aztec
   curvature with adequate precision and cutoff checks. A significant positive
   pilot should not itself trigger publication; the planned independent dataset
   must carry the inferential claim. An inconclusive or negative result changes
   the scientific question rather than justifying selective windows.
4. **Analysis coverage:** every declared probe must be available across the
   frozen size window, with at least four distinct actual separations and full
   design rank. Rounding-induced duplicated distances are explicit. Include
   fixed cutoff sensitivities and joint control diagnostics.
5. **Budget and stopping:** specify numbers of environments, paired samples,
   wall time, memory and a hard allocation before production. Extend the same
   environments/streams under a frozen cap policy; never replace hard cases.
   If incompleteness persists, report it and withhold the complete-campaign
   curvature contrast. A cap is a budget limit, not a guarantee of certification.

## Precision calculation available now

For the exploratory two-probe corrected Aztec fit over orders 16–128, 64
environments per size gave a bootstrap 95% interval [-0.7506,2.0800], approximate
half-width 1.4153. A planning-only
1/sqrt(N) extrapolation for a half-width h is

    N per size ≈ 64 × (1.4153 / h)^2.

For illustration h=0.25 suggests about 2,052 environments per size in the
same six-size design. This is not a recommendation or approved allocation:
the pilot is small, non-Gaussian finite-sample effects and larger-size behaviour
may change the variance, and the direct two-domain contrast also has square
uncertainty. The precision target must be tied to a scientifically meaningful
coefficient difference. Power for a specified alternative has not been
established. Extrapolating to h=0.25 does not prove feasibility or positive signal.

## What is deliberately not frozen

A production sample count and square size range cannot be justified from the
current censored pilot. These remain conditional on the gates above; presenting
arbitrary values as a final design would conceal the principal uncertainty.
The next implementation task is a validated rare-event/numerical strategy and
an exact-sampling feasibility investigation on the retained small cases.


## Update, 1 October 2026

The small-width feasibility gate has advanced: the integer frontier sampler
completed all 16 retained L=8 random pairs and four retained L=12 random pairs.
Its exact arithmetic targets the stored dyadic weights and avoids the CFTP
bottleneck. The implementation is deliberately limited to L<=12; large-width
feasibility remains open. Robust Aztec coins and a three-environment
high-precision reduction audit now support a separately planned corrected
benchmark refinement. This refinement is prepared, not submitted. The hold
on the large matched two-domain production comparison remains in force.
See `NUMERICAL_REMEDIATION_20261001.md` for evidence and limitations.


**1 October provenance qualification:** Seed regeneration alone did not prove
identity with all archived Hamilton arrays. Of 24 overlapping environments in
the new integer-sampler cohort, 14 archives match the recorded byte hashes and
10 do not. Differences in archive packaging versus array values remain to be
determined. Earlier local seed-replay statements should be read subject to
this limitation; the Hamilton pilot's own completion records are unchanged.
See `docs/validation/numerical_20261001/environment_identity_audit.json` and the
updated runner's mandatory source-selection and archive-checksum checks.


## Original-array resolution and job status, 2 October 2026

All 240 original Hamilton archives were downloaded and checksum-verified.
Compared with local seeded regeneration, 64 random environments differ in
actual stored values (maximum 9 ULP; relative difference 1.225e-15); the other
176 agree. The precise numerical cause remains unisolated. This supersedes the
packaging-versus-array uncertainty above, while qualifying earlier regenerated
replays.

The integer sampler was rerun on the original files: all 20 random pairs at
L=8 and L=12 and four archived control pairs completed; two new L=4 control
pairs also completed. The 24 overlapping input archives all match the pilot
hashes. Verification rechecked 52 perfect matchings and 224 observable values.
Large-width feasibility remains unresolved.

Source upload and verification of 135 files are complete. Hamilton validation
job 18872074 completed successfully in one minute: all 120 Python tests,
21,970 Julia assertions, the numerical audit and benchmark passed. The approved
corrected Aztec refinement was submitted on 2 October as array 18897605.
It has started; scientific results await completion and the frozen analysis.
The hold on large matched two-domain production remains in force.


## Completed corrected Aztec refinement, 3 October 2026

All 12,288 pairs completed in Hamilton array 18897605. The frozen two-probe
analysis gives c=0.3568 [0.0723, 0.6327] over orders 16–128. The 32–128 and
48–128 cutoff intervals also exclude zero, but estimates drift and uncertainty
widens substantially. This supports finite-window positive curvature, without
establishing a stable asymptotic coefficient or a boundary contrast. No new
campaign was submitted. Wider-square feasibility remains the limiting step;
the hold on large matched two-domain production continues. See
[the refinement report](AZTEC_REFINEMENT_20261003.md).


## Bounded wider-square test prepared, 4 October 2026

The next feasibility plan fixes 18 cases at sides 14, 16 and 18, including
original side-16 archives, with explicit time, memory and state budgets.
All 132 tests passed both locally and on Hamilton (validation job 18902368).
Array 18902369 returned 12 complete pairs and six state-capped side-18 cases.
After an explicit resource review, separate array 18902387 reused those six
inputs and seeds with an eight-million-state limit; all six completed. Time
and memory limits were unchanged. All 18 retained environments now have
verified sample pairs, with original first-stage censoring preserved. This bounded assessment
falls within the user's instruction to continue the comparison; the hold on
large matched production remains. See [the plan](SQUARE_FEASIBILITY_20261004.md).


Side-18 Gamma worker times were 39.31–39.75 seconds; maximum peak RSS was
4.356 GiB. This is a successful small-cohort feasibility assessment, not an
adequate square curvature dataset or large-width cost guarantee. Prioritise
validation of a scalable exact sampling route before large production. No
additional jobs were submitted after these two bounded stages.

## Scalable matched sampler prepared, 4–5 October 2026

Embedded square shuffling now passes exact small-grid distribution checks,
the fixed high-precision audit through L=64, all 146 local Python tests and
a 144-pair end-to-end smoke cohort covering both geometries and controls.
The fixed matched design has 30,720 pairs across six sizes, with a direct
Aztec-minus-square curvature contrast and two prescribed cutoff windows.
See [the construction, evidence and plan](MATCHED_SHUFFLE_20261004.md).

The local computational bottleneck is resolved for this scope, subject to
Float64-reduction qualifications. The remaining launch gate is successful
Hamilton validation of the frozen source and resource review. SSH expired
before upload/submission; no new campaign is running. This is an execution
gate, not a request for fresh user approval of already authorised work.

## Launch gate passed and production submitted, 5 October 2026

Renewed SSH access allowed source verification and compute-node validation.
Job 18907657 passed in 2 min 56 sec (0:0; empty stderr): 146 Python tests,
21,970 Julia assertions, both numerical audits, all 144 smoke pairs, 288
reloaded matchings and 576 recomputed observables. The largest Gamma smoke
cell used about 348 MB and took 5.14 seconds for four pairs. The fixed
production limits provide ample measured headroom for this scope.

The launch gate is satisfied. Production array **18907659** is submitted;
post-processing job **18907660** depends on the success of the entire array
and revalidates raw files before the fixed statistical analysis. The earlier
hold is superseded by these measured checks, within the user's already
explicit instruction to perform the matched comparison. No numerical
comparison conclusion is available yet. See the status and evidence in
`docs/validation/matched_shuffle_20261004/`.

## Completed comparison and interpretation, 5 October 2026

All 30,720 pairs and the dependent raw verification/analysis completed (0:0).
The primary Aztec-minus-square curvature difference is −0.046 with 95%
interval [−0.366,0.282]. Both individual primary coefficients are positive;
the cutoff contrasts remain inconclusive. The proposed asymptotic distinction
is not established. An all-one control curvature contrast is unexpectedly
positive; a labelled post hoc six-contrast permutation audit gives p=0.0093,
whereas the prespecified maximum-coordinate control test gives p=0.5882.
These are different statistics. Investigate this discrepancy and finite-size
sensitivity before another scaling campaign; do not extend samples until a
desired conclusion appears. No new campaign was submitted. See
[MATCHED_RESULTS_20261005.md](MATCHED_RESULTS_20261005.md).

## Independent all-one confirmation, 6 October 2026

The user instructed continuation of the next steps. Read-only diagnostics
found matching physical endpoints, no seed collisions, and no single pair
explaining the uniform curvature contrast. All 38 selected saved pairs
replayed exactly on Hamilton; the 48-pair fresh confirmation pipeline passed
validation job 18919384. These checks do not establish independence of all
pseudorandom outputs or explain the original anomaly.

A separately frozen, finite confirmation uses only all-one weights, 1,024
pairs at each of six sizes in each geometry, unchanged samplers and fresh
seeds. Array 18919405 and raw-verification/analysis job 18919413 were submitted.
No outcome-dependent extension or
pooling into the old experiment is specified. No new Gamma scaling campaign
has been submitted. See [the control investigation](CONTROL_INVESTIGATION_20261006.md).

## Confirmation completed and checked, 7 October 2026

All 12,288 pairs completed, with successful raw verification and analysis.
The fixed primary contrast is 0.081098 [−0.281718,0.437180], permutation
p=0.6615. Thus the earlier anomaly was not reproduced at this precision.
The original anomalous control cohort remains reported; no data were pooled
or replaced. The inconclusive matched Gamma contrast is unchanged.

The finite confirmation is closed under its original stop rule. No new jobs
were submitted. Next: review the completed evidence and mathematical
validation with Sunil before fixing a wider matched window, and complete
human-supplied historical AI attribution. A non-rejection in this confirmation
is not proof of correct sampling or a reason to assert the desired growth-law
contrast.

## Larger campaign requested and frozen, 7–8 October 2026

Alberto requested a larger campaign without time pressure and then instructed
continuation. The bounded square feasibility array 18960990 and verification
18960991 passed, with all 24 pairs checked and the numerical audit extended
through side 128. A new separate 258,048-pair design is frozen at seven sizes
16–128, with 16,384 Gamma pairs and 1,024 pairs of each deterministic control
per size/domain. Earlier samples are excluded from the new inference.

The source snapshot retains all 155 old files unchanged and adds nine batch,
validation and analysis files. Local checks passed, including 152 Python
tests and both-language integration. Source, plan and seed records are in
`docs/validation/larger_matched_20261007/`. Production requires the new
168-pair Hamilton validation and full tests to pass first. SSH expired after
pilot retrieval; validation and production remain unsubmitted. The user's
scientific authorization persists; the outstanding requirement is connection
and technical validation, not another permission request.


## Larger campaign launch, 8 October 2026

The new gate 18965117 passed with exit zero and empty stderr in 9m10s:
152 Python tests, 21,970 Julia assertions and 168 raw-verified pairs.
All 164 frozen source hashes passed. Under Alberto's existing authorization,
arrays 18965699 and 18965700 and dependent analysis 18965701 were submitted.
The second array waits for the first to succeed, keeping sampling concurrency
at two tasks. The fixed 258,048-pair cohort is now running; no scaling result
is available yet. See `LARGER_CAMPAIGN_20261007.md` and the launch evidence.
