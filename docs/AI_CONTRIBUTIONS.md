# AI contribution record

Updated 8 October 2026. This is an evidence-based working record, not a
claim that earlier project history has been fully reconstructed. No email or
paper disclosure has been sent or published.

## Confirmed for this conversation

The assistant was OpenAI Codex. The exact deployed model build and provider
retention settings are not established by the repository and are not asserted.
Alberto supplied the supervisor correspondence, requested the work, selected
the research objective, provided Hamilton access and authorised the bounded
follow-up. Those are observed human contributions; they do not establish which
earlier mathematical or coding steps were performed without AI.

| Area | Confirmed AI assistance | Evidence and verification |
|---|---|---|
| Reference document | Synthesised repository material and supplied correspondence; drafted/edited Typst handbook and reports | `docs/research/handbook.typ`, dated validation records; compiled PDF reviewed for layout by the assistant |
| Model and coordinate mapping | Inspected conventions; implemented matched weight placement, physical probes and orientation adapter | `comparison.py`, `comparison_geometry.jl`; tiny enumeration and independent height-path checks |
| Sampler audit | Wrote exact creation-path versus Gibbs-product audit; identified the finite-environment legacy mismatch | `audit_comparison_mapping.jl`, saved order-2/order-3 results; independent mathematical review by a human is not recorded |
| Python/Julia code | Generated and revised comparison runners, recovery, audits, analysis code and tests | Current working-tree changes and source manifests; 108 original Python tests and 21,926 Julia assertions passed on Hamilton; follow-up tests recorded separately |
| HPC execution | Prepared isolated source snapshot, checked hashes, submitted the authorised small jobs and downloaded outputs | Jobs 18852141, 18852143, 18852144; scheduler states and checksums |
| Statistical analysis | Implemented joint control permutations and environment-block bootstrap; selected explicit exploratory windows, examined cutoff sensitivity | Dated plans, seeds, scripts, machine-readable estimates and source hashes |
| Follow-up diagnostics | Regenerated original pilot weights; replayed bounded CFTP maps; audited rare probabilities; ran corrected Aztec pilot | `docs/validation/followup_20260930/`, saved pilot weights and tilings |
| Interpretation and writing | Drafted explanations, limitations, no-go production decision and this disclosure | `FOLLOWUP_20260930.md`, `PRODUCTION_DECISION.md`; these are AI-assisted interpretations requiring author review |

The tests and statistical checks were themselves partly AI-written. Passing
them is evidence about the tested properties, not independent human approval
or a proof of every mathematical claim. Distinct computational checks are not
being presented as independent human authorship.

## Historical work: attribution not established

Git records show project stages for random walks, UST/LERW, Aztec experiments,
square-grid dynamics, tempering, Kasteleyn replay and later CFTP campaigns.
Examples include initial project commit `f58df997`, Aztec expansion `6382063a`,
spatial Aztec experiment `7b3f355e`, Glauber setup `6fe76b54`, Kasteleyn replay
`83018864`, and September status `719ddefc`. Their presence and commit author
fields establish repository history, not who originated each idea or whether
AI assisted. No AI/co-author trailer was found in the reviewed commit messages;
absence is not evidence that AI was not used.

The supervisor's supplied email explicitly accepts AI assistance and asks for
honest disclosure. It does not authenticate individual code contributions or
certify the new audit findings. Earlier emails and historical notes can record
what Alberto reported, but do not by themselves determine AI involvement.

Pending Alberto's account: tools used before this conversation; which proofs,
derivations and code he personally checked; portions written unaided; and any
other collaborators' contributions. Unknown entries must remain unknown rather
than being assigned to a human or to AI by inference.

## Suggested disclosure limited to confirmed work

> OpenAI Codex assisted with repository review, software implementation and
> documentation, computational validation, exploratory statistical analysis,
> preparation of bounded HPC experiments, and drafting of research notes.
> Its contributions included an audit of sampler orientation and matched-probe
> conventions. Validation used automated tests, finite exact-enumeration
> comparisons and reproducible computational records. These checks do not
> replace the authors' mathematical assessment. Historical AI involvement
> before this documented phase is still being reconstructed.

This wording is suitable as a transparent working note. Before submission of a
paper, reconcile the historical record and describe the actual human review;
do not add a claim that the authors verified all results until that has happened.

## Ongoing record

For each new substantial change retain: date, tool identity where known,
request/scope, affected files, source snapshot hash, tests and outputs, human
review status, and limitations. Do not store passwords, SSH credentials or
private account data in the contribution ledger. Research provenance should
remain separable from authentication material.


## 1 October 2026 continuation

Codex generated the exact-integer Bernoulli/reference heat-bath code, the
bounded-width frontier sampler and its rational-enumeration/detailed-balance
tests. It ran the fixed 26-environment cohort, preserving original disorder
identities with a new sampler-seed namespace. It also implemented and tested
the interval-filtered rational Aztec coins, compared Float64 reductions with
256/512-bit references on three fixed environments, and prepared the revised
Aztec runner, scheduler scripts and bootstrap analysis. The larger refinement
has not been submitted because Hamilton authentication was unavailable.
A supervisor update was drafted but not sent. These contributions remain
AI-assisted; no human mathematical approval is inferred from their tests.


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

Codex checked the completed compute-node validation, submitted the authorized
refinement array and recorded the job identifiers on 2 October. Human
mathematical review and historical AI attribution remain to be confirmed.


## 3 October 2026 analysis

Codex checked completion of all six tasks in array 18897605, downloaded the
compact archive, verified hashes, source versions, coverage and seed
identities, and independently recomputed the moment estimates. It ran the
previously prepared 10,000-resample analysis without changing its windows or
sample count, produced the result figure, and updated the handbook and unsent
supervisor draft. Raw tilings were validated by the sampler at generation;
this follow-up did not download and revalidate every raw tiling. The new
positive finite-window result and its cutoff sensitivity are documented in
`AZTEC_REFINEMENT_20261003.md`. No human mathematical review or historical
AI-attribution information is inferred from these checks. No email was sent.


## 3–4 October 2026 wider-square preparation

Codex implemented the separately versioned budgeted frontier sampler,
vertex-wise exact gauge conversion, compact partition cache and explicit
stack. It designed the fixed 18-case resource-limited assessment, wrote worker,
scheduler and outcome-verification tools, and ran all 132 Python tests. Tests
compare matching masses and conditional probabilities, rather than assuming
gauge-equivalent integer representations consume the same random bits. The
new work passed Hamilton validation job 18902368 (132 tests), after source
and original-input checksum checks. Codex submitted the fixed feasibility
array as 18902369 under the user's instruction to continue. The handbook and run guide record the proof argument,
resource limits, censoring policy and limitations. Generalized zero-weight
shuffling was investigated in the primary literature but not implemented.
No human mathematical review is inferred; no email was sent.


Codex retrieved and verified all first-stage outcomes from array 18902369:
12 complete pairs, six state-capped side-18 cases. It then documented an
explicit second-stage decision based on measured resource use: preserve the
first stage, reuse all six saved side-18 inputs and seeds, retain the time and
memory limits, and raise only the state-count cap to eight million. The
separate plan and array 18902387 are recorded; they do not alter the original
18-case result or constitute a production curvature experiment.


The six-case follow-up completed. Codex downloaded all weights, samples and
records, verified checksums and retained-input identity, revalidated all 36
matchings across the two stages, and recomputed all 72 observable values.
It checked final scheduler exits and empty error logs and updated the handbook,
run guide, production decision and unsent supervisor draft. No curvature fit
was made from the feasibility cohort and no larger campaign was submitted.

## 4–5 October 2026 matched-law shuffling preparation

Codex implemented the square-in-diamond weighted shuffling construction,
its rational small-grid oracle, exact matching-distribution comparisons and
Float64/80/160-digit numerical audit. It implemented the two-domain campaign
runners, saved-input and matching certificates, independent raw observable
recomputation, fixed bootstrap contrast and joint control audit, and prepared
bounded Slurm scripts. All 146 local Python tests and the 144-pair local smoke
cohort passed. It documented the distinction between exact rational sampling
and numerically audited Float64 reductions, and updated the Typst handbook.
This is confirmed AI assistance; no independent human mathematical review is
inferred. Local smoke samples are excluded from scientific results. Hamilton
access expired before this new campaign was uploaded or submitted. No email
was sent, and no new comparison result is claimed.

## 5 October 2026 Hamilton launch

Following the user's renewed SSH authentication, Codex verified all 155
source hashes on Hamilton and submitted validation job 18907657. It corrected
a pending partition-time-limit mismatch by reducing the live validation
allocation from 20 to 15 minutes. Validation passed with 146 Python tests,
21,970 Julia assertions, both numerical audits and raw verification of all
144 smoke pairs. Codex downloaded and verified the compact validation records,
reviewed measured runtime/memory, and submitted the already-authorised fixed
production array 18907659. It queued raw verification and the frozen analysis
as job 18907660, dependent on successful completion of the entire array.
The separate scheduler script and its checksum are recorded. The handbook,
run status and unsent supervisor draft were updated. No production conclusion
or independent human review is claimed; no email was sent.

## 5 October 2026 completed matched results

Codex checked all 36 successful sampling exits and the successful dependent
analysis, downloaded a checksummed compact results archive and confirmed full
coverage of 30,720 pairs. Raw matching and observable verification had run on
Hamilton for all 61,440 matchings and 122,880 values; the follow-up downloaded
compact records rather than every raw file. It independently recomputed all
216 moment estimates and 18 coefficients, prepared the fixed-window figure,
and documented the inconclusive primary contrast without selecting a new
window or altering the original analysis.

After seeing an unexpected all-one control curvature contrast, Codex added an
explicitly post hoc whole-vector permutation diagnostic of all six control
curvature contrasts (seed 2026100514, 9,999 permutations). Its p=0.0093 and the
prespecified coordinate-test p=0.5882 are both retained, with the distinction
between their statistics and the limitations of post hoc selection explained.
This does not prove a software bug or establish the proposed scaling contrast.
The handbook and unsent email were updated. No additional campaign or email
was sent; no independent human mathematical review is inferred.

## 6 October 2026 control investigation

Codex investigated the saved control records, reconstructed physical probe
endpoints, checked 108,000 old/new role seeds for reuse, decomposed the
uniform curvature contrast by size, and performed labelled leave-one-pair-out
and lag diagnostics without altering the primary data or estimates. It
independently checked the permutation covariance formula by exhausting all
24 pairings of a four-pair two-probe example. It distinguished what a
zero-disorder control tests (pairing independence/analysis) from validation
of the marginal Gibbs distribution.

Codex designed and froze one independent all-one confirmation, 1,024 pairs
per size/domain (12,288 pairs total), with a single primary test and no
adaptive extensions. It reused the unchanged, hash-verified 155-file sampler
snapshot and added separately checksummed replay/analysis tools. Local and
Hamilton 48-pair smoke pipelines passed; Hamilton job 18919384 also replayed
38 original pairs exactly, including reverse-order generation of square
replicas. Under the user's instruction to perform the next steps, Codex
submitted confirmation array 18919405 and dependent raw-verification/analysis
job 18919413. Results are pending. It updated the handbook and unsent draft;
no email was sent and no human mathematical review is inferred.

## 7 October 2026 completed control confirmation

Codex reconnected using Alberto's authenticated Hamilton session, checked all
12 successful sampling tasks and the successful analysis job, and downloaded
a checksummed compact evidence archive. Hamilton had verified all 24,576 raw
matchings and recomputed 49,152 observable values. Codex added a separate
local verification script that checks counts, identities and hashes, then
independently recomputes 72 moments, six coefficients, three contrasts and
the fixed 9,999-permutation primary test. This reproduces p=0.6615. Original
bootstrap intervals were retained rather than independently rerun.

It documented the fresh contrast 0.081 [−0.282,0.437] as non-reproduction of
the old control anomaly at this precision, preserving both cohorts and the
inconclusive Gamma result. It updated the handbook, report, run status and
unsent supervisor draft. No new sampling job or email was submitted. These
are AI-assisted computational checks, not independent human mathematical
review; Alberto's historical attribution account remains outstanding.

## 7 October 2026 newly supplied historical correspondence

Alberto supplied a longer July–September email chain. His 2 August message
explicitly reports ChatGPT assistance with inverse-bootstrap-variance
weighting. His reply to Sunil's 29 September message reports LLM assistance
with understanding papers and notes, debugging, Python/Julia porting and
SLURM. Record these as Alberto's historical self-report; no original prompt
logs or code-to-author mapping were supplied. They resolve some historical
uncertainty but do not establish unaided authorship of the remaining work.

Sunil's 7 August email expressed a preference against AI help in the eventual
write-up. His later 29 September email accepts AI assistance with honest
disclosure, without explicitly resolving the earlier manuscript-specific
preference. No claim about a journal's actual policy is inferred. The notes,
reports and email drafts here remain AI-assisted working documents. Any
eventual manuscript workflow should be agreed by its human authors.

Codex reconciled the newly supplied correspondence against the recorded
experiments in `docs/correspondence/reconciliation_20261007.md`, distinguishing
the completed small-window matched comparison from full replication of the
old large-size campaign and noting the change from CFTP to shuffling. This
interpretation is AI-assisted; it is not an external review or a new experiment.

## 7 October 2026 larger-campaign feasibility

After Alberto asked about a larger campaign without time pressure, Codex
checked Hamilton access and previous runtime records, identified the current
square implementation's side-128 limit, and prepared a fixed 24-pair square
feasibility gate at sides 96 and 128. It reused the unchanged 155-file source
snapshot with separately checksummed wrappers, verified seed non-reuse across
the old matched, confirmation and new pilot plans, and submitted array
18960990 with dependent raw verification 18960991. No full production run was
submitted. It implemented a labelled planning approximation for curvature
precision from old centred-product covariance estimates, with explicit
larger-size extrapolation assumptions, and documented candidate designs in
`docs/history/LARGER_CAMPAIGN_20261007.md`. These assumptions and design suggestions
are AI-assisted, not supervisor-approved theoretical predictions.

## 7–8 October 2026 fixed larger campaign and batch implementation

Codex retrieved the completed feasibility gate (all 24 pairs, both numerical
audits and raw verification), downloaded the complete checksummed archive,
and rechecked raw files locally. It froze a separate 258,048-pair matched
comparison through size 128, including 28,672 deterministic control pairs,
with fixed analysis windows and a predeclared joint control-curvature test.
It implemented independent 512-pair batches, file locking, persistent plan
identity, explicit same-seed retries preserving old attempts, raw-verification
receipts and full-cohort analysis. It added restart/provenance tests and a
two-stage validation/production launcher. All original 155 source files are
unchanged; nine new files extend the snapshot.

All 152 Python tests passed locally, an eight-pair two-language integration
check passed, the control test reproduced the earlier independent diagnostic,
and 753,648 role seeds were checked for collisions across the relevant plans.
Hamilton access expired after pilot retrieval; neither new validation nor
full production has yet been submitted. Scientific sampling algorithms are
unchanged, but the orchestration and new analysis require their recorded
compute-node validation before launch. These are AI-assisted code, design,
analysis and documentation contributions, not human mathematical approval.


## 8 October 2026 larger campaign submission

After Alberto restored SSH access, Codex uploaded the frozen archive, checked
all 164 source hashes and submitted Hamilton validation 18965117. It inspected
the successful accounting, empty error log, 152 Python tests, 21,970 Julia
assertions, 168 raw-verified paired samples and validation analysis. It then
submitted the authorized fixed production arrays 18965699 and 18965700,
with analysis 18965701 dependent on their successful completion. It saved
launch/accounting evidence and updated the campaign notes and handbook.
Production outcomes remain pending; this records AI-assisted orchestration
and documentation, not a human review or a scientific conclusion.


## 10 October 2026 completed larger-campaign retrieval

Codex checked successful completion of all 504 batches and analysis job
18965701, packaged and downloaded compact measurement/provenance records,
and verified the archive checksum. It prepared an independent point-estimate
check alongside the existing full-cohort provenance checker and drafted
`docs/RESULTS.md`, retaining all fixed cutoff and control
results. Original bootstrap intervals and the joint permutation result were
retained from Hamilton, not independently regenerated. No new sampling,
selective exclusions or supervisor email were initiated. Interpretation and
report drafting in this step are AI-assisted; human scientific review remains.


## 10 October 2026 repository cleanup and publication workflow

At Alberto's request, Codex rewrote the README around the individual Aztec
and ordinary-square coefficients, organised dated reports under `docs/history/`,
and removed legacy archives/data from Git tracking while preserving them on disk
and in earlier revisions. It exported all 258,048 paired measurements with
seeds, checksums, original batch verification provenance and the frozen Hamilton
source snapshot into the public `data/` bundle. It added readable export and
analysis entry points, a short analysis guide, a covariance figure and CI checks.

The new public-data workflow independently reran the complete 10,000-resample
bootstrap and 9,999-permutation control test. All coefficients and intervals
matched within 3.83e-14; moment intervals agreed exactly and the joint p-value
was again 0.2039. The 152 Python tests passed. Code cleanup, documentation,
figure generation, data packaging and reproduction are AI-assisted. Alberto
reports that he has sent Sunil an email; Codex did not send it. Human mathematical
and manuscript review remain distinct from these computational checks.

The Julia suite also passed all 21,970 assertions and all command-line smoke
workflows during cleanup. The local Makefile runs the standard-library test
entry point directly, avoiding an unnecessary registry download.
