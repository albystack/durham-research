# Larger matched campaign: fixed design and validation, 7–8 October 2026

**The 24-pair larger-square feasibility gate passed, including both numerical
audits and raw verification. The new 258,048-pair campaign is frozen and
validated on Hamilton. Production arrays 18965699 and 18965700 are submitted,
with dependent analysis 18965701.**

Alberto asked about running a larger campaign and said there is no rush.
A larger physical size window is scientifically motivated by the cutoff
sensitivity of the completed comparison. More samples alone shrink random
error; they do not remove finite-size bias or guarantee a domain contrast.

## Work started

The renewed Hamilton session allowed pilot retrieval before expiring. Original source checksums were
verified for all 155 files. No scientific sampler code was changed.

The current square prototype has an explicit maximum host order of 127,
hence square side 128. Its cached geometry and stored reduction levels grow
approximately cubically. Larger scheduler limits alone cannot support 256:
that needs a deliberate implementation extension and fresh validation.

One bounded feasibility gate is frozen at square sides **96 and 128**:
eight Gamma pairs, two deterministic-family-mean pairs and two all-one pairs
at each size, **24 pairs total**. It also compares every potential creation
probability at each size using Float64 and 80/160 decimal-digit reduction.
The original tolerances are unchanged. Failures are retained, without replacing
inputs. Pilot observations are excluded from future scientific production,
and no curvature coefficients or significance tests are calculated for them.

Hamilton array **18960990** runs four one-CPU tasks, at most two concurrently,
with 16 GiB and one hour per task. Dependent job **18960991** verifies every
saved matching and observable after all four tasks succeed. No full production
campaign is submitted automatically.

The two numerical audits have passed:

| Square side | Potential coins per replica | Maximum probability discrepancy | Twice sum of discrepancies | Audit peak memory |
|---|---:|---:|---:|---:|
| 96 | 290,320 | 3.32e−15 | 4.02e−11 | 1.35 GiB |
| 128 | 690,880 | 4.91e−15 | 1.17e−10 | 3.09 GiB |

These are fixed-environment numerical comparisons, not universal error bounds.
All four tasks and dependent verification completed with exit 0:0 and empty
error logs. Sampling took 97 seconds at side 96 and 189 seconds at side 128,
including the three models; raw verification took 29 seconds. The eight-pair
Gamma cells took 51.33 and 115.87 seconds respectively. The side-128 Gamma
process peak was 1.86 GB (about 1.73 GiB). The complete pilot archive was
downloaded and checksum-verified; local rechecking also verified all 24 pairs,
48 matchings, 96 observables and 96 certified files. Archive SHA256:
`2a33721ec24e4ea67f42e695735784c0dfbfa11de1158103d41204346424cdab`.

## Precision planning from the completed data

`docs/validation/larger_feasibility_20261007/precision_planning.py` uses only
the old completed matched cohort. It estimates covariance of the two
sample-covariance estimates from centred-product influence vectors, preserving
dependence between probes. Independent size/domain blocks are combined using
the least-squares curvature weights. No pilot coefficient is used to choose
a design.

For unobserved sizes above 64, the size-64 influence covariance is extrapolated
using squared-log growth of the influence magnitude. A separate scenario
multiplies that larger-size covariance by four. These scenarios are neither
confidence bounds nor reliable forecasts of unknown asymptotic behaviour.
The approximate normal interval widths do not replace a production bootstrap.

| Sizes used in primary fit | Gamma pairs per size/domain | Approximate contrast interval half-width | Square Gamma CPU-hours estimate |
|---|---:|---:|---:|
| 16,24,32,48,64,96,128 | 4,096 | 0.319–0.526 | 19.2 |
| Same seven sizes | 8,192 | 0.225–0.372 | 38.4 |
| Same seven sizes | 16,384 | 0.159–0.263 | 76.7 |
| 32,48,64,96,128 | 8,192 | 0.586–0.946 | 38.1 |
| 32,48,64,96,128,192,256 | 8,192 | 0.320–0.604 | 306.8 |

CPU estimates extrapolate the previous side-64 Gamma job cubically. They
exclude controls, validation, raw verification and scheduler waits, and must
be revised using the pilot timings. They are aggregate CPU-hours, not wall
time promises. The size-256 option requires additional implementation work.

The tradeoff is visible: larger lower cutoffs reduce the span of log distances
and can make curvature much harder to estimate. Increasing the largest size
and adding independent environments should therefore be planned together.
In the old 8–64 window this approximation predicts a 95% half-width about
0.326 for 2,048 pairs, close to the observed roughly 0.324; that agreement
does not validate extrapolation to unknown sizes.

## Frozen production design

| Choice | Fixed value |
|---|---|
| Sizes, both domains | 16,24,32,48,64,96,128 |
| Gamma pairs per size/domain | 16,384 |
| Each deterministic control per size/domain | 1,024 (family means and all-one) |
| Total | 229,376 Gamma + 28,672 control = **258,048 pairs** |
| Probes | Physical diagonal fractions 1/8 and 1/4, height scale one |
| Primary estimand | Gamma c_Aztec−c_square, window 16–128 |
| Sensitivity windows | 32–128 and 48–128; report all, no cutoff selection |
| Bootstrap | 10,000 whole-environment resamples, seeds 2026100713/2026100714 by domain |
| Control diagnostic | Joint maximum absolute standardised curvature contrast over 2 controls × 3 windows; 9,999 whole-vector permutations, seed 2026100715 |
| Batch size | 512 pairs; 504 independent batches |
| Master seeds | 20261007000000 + batch index, then original role-seed hash derivation |
| Small-size array | 360 tasks through side 64, one CPU, 4 GiB, 30 min/task |
| Large-size array | 144 tasks at 96/128, one CPU, 4 GiB, 4 h/task |
| Concurrency | At most two sampling tasks total; large array follows successful small array |
| Analysis ceiling | One CPU, 8 GiB, 8 h |
| Stop rule | Fixed complete cohort; preserve failures; no replacement or outcome-dependent extension |

The two smaller historical probe fractions are not added to this frozen stage.
The original two-probe convention is retained so the old validation and raw
measurement checks apply. Neither earlier production nor pilot observations
are pooled into this new inference. Controls are diagnostics, not a rule for
discarding data or withholding an unwanted Gamma result.

The production plan is `square_glauber_python/configs/larger_matched_20261007.json`,
SHA256 `cf831391da6b9b02d19c893f2973cb44063b463ff64eafca2a85ae1cf2fa81cc`.
The separately seeded compute-node validation has 42 four-pair batches,
168 pairs, and reduced 100-bootstrap/99-permutation settings marked
`VALIDATION_ONLY`. Its plan SHA256 is
`43a9a3334fda90631abdb183f2395856b82366f7de661bf8b1e863e6786e87c0`.

## Implementation, restart and validation

`square_glauber_python/scripts/larger_campaign.py` wraps the existing square
sampler and raw verifier. The new Julia runner differs from its original only
in the comment and allowed size/count gate, extending through 128 with at
most 512 pairs per batch. All 155 original source files have unchanged hashes;
nine added files give a 164-file frozen source snapshot.

Each batch has a persistent identity and an operating-system file lock. A
verified completion is reused after its evidence is checked. An interrupted
attempt is preserved; an explicit `--retry` creates a separate attempt with
the same plan and seeds. Previously complete raw cases must match byte for
byte. No failed environment is replaced. A batch gets a completion receipt
only after every matching and observable passes raw verification. Analysis
requires all receipts and rechecks table/plan hashes, full coverage and role
seeds. Batch prefixes make environment identifiers unique; resampling moves
whole paired vectors across batches within each size/model/domain, rather
than treating a batch as an independent bootstrap unit.

Local validation passed **152 Python tests**, including restart, interrupted
attempt preservation, tamper rejection, incomplete-coverage rejection and
fixed plan totals. Both-language integration checked eight actual paired
samples with unique merged identities. The new control test independently
reproduced all six old diagnostic contrasts and p=0.0093 using the old seed;
this is a regression check, not a second scientific experiment. All **753,648
role seeds** across old matched, confirmation, feasibility, and new production/
validation plans are distinct.

The frozen source archive SHA256 is
`4c533d765419550fe3fcee005ee8ff58ed670f3ac3d08f7c2a387412a6fad475`.
The launcher `docs/validation/larger_matched_20261007/launch_hamilton.py`
first verifies the archive, source and completed numerical gate, then submits
Hamilton validation. Its separate production stage requires a successful
validation exit, empty stderr, all Python/Julia tests, all 168 raw-verified
pairs and successful validation analysis before submitting either array.
Production analysis is queued after successful completion of the arrays.
The gate passed on 8 October; the full fixed campaign is now submitted.

Revised runtime planning using the pilot suggests roughly **105 square Gamma
CPU-hours**, plus about 12 square control CPU-hours, the Aztec arm, batch
startup and verification. With two concurrent tasks, several days is a
reasonable working expectation, excluding queue delays. This remains an
estimate from a small cohort; scheduler ceilings are larger. The original
cubic extrapolation table above is retained as the earlier planning model.

A study through 128 is an extension beyond the completed square window, not
a recreation of historical Aztec orders 128–1300 or a guarantee that boundary
effects will be resolved. Exact theoretical-observable clarification and
human mathematical review remain useful in parallel with this planning.

Evidence: `docs/validation/larger_feasibility_20261007/`.
Frozen pilot plan SHA256:
`4d29b9c4ea6ade6031be0cc6a13bf6399c766e4cceac864f88f691b4d245e684`.
Tools archive SHA256:
`f4846ec4adc29c631c4b9cee24c21ac757b60ca153721244e3dc940be4353275`.
Remote root: `/nobackup/fvkl37/larger_feasibility_20261007`.
New production root: `/nobackup/fvkl37/larger_matched_20261007` (created and source-verified).
Production status: `docs/validation/larger_matched_20261007/hamilton_status.json`.

Hamilton validation job **18965117** was submitted on 8 October after all 164
source hashes passed. It completed in 9 minutes 10 seconds with exit zero and
empty stderr: 152 Python tests, 21,970 Julia assertions, and all 168 pairs
(336 matchings, 672 observable values) passed. Storage usage before production
was about 8.8 GiB of a 600 GiB quota.

Production array **18965699** covers batches 0–359 (sizes 16–64); array
**18965700** covers batches 360–503 (sizes 96 and 128), with an after-success
dependency on the first array. Both have a two-task concurrency limit, so at
most two sampling workers run at once. Analysis job **18965701** depends on
successful completion of the second array. First sampling tasks were observed
running on 8 October. These are submitted jobs, not completed scientific results.
The saved `launch.json` and `validation-accounting.txt` record the launch evidence.


## Progress check, 9 October 2026 at 16:51 BST

Hamilton access was restored via the oct9 control socket. Scheduler accounting
and completion receipts agree on **438/504 completed batches**, representing
**224,256/258,048 raw-verified pairs (86.9%)**. Every domain/model cell through
size 96 is complete. Six size-128 square Gamma batches are complete; batches
438 and 439 are running. No failed jobs or non-empty campaign error logs were
found. Analysis 18965701 remains pending on sampling completion.

The six completed size-128 batches took 104–131 minutes each, with maximum
recorded resident memory about 1.79 GiB, below the 4 GiB request. Storage is
about 33.5 GiB of the 600 GiB quota. A two-worker estimate using observed cell
medians and cubic 96-to-128 extrapolation for remaining controls/Aztec cells
gives about 24.8 hours of sampling still required. Allow roughly 25–35 hours
plus queue delays and final analysis; runtime variation makes this provisional.
Completion percentage counts pairs, not computational work. No partial scaling
fits were inspected or used to alter the fixed campaign.

Evidence: `docs/validation/larger_matched_20261007/progress_20261009.json`
and `status-20261009-accounting.txt`. No jobs or scientific settings changed.


## Progress check, 10 October 2026 at 12:06 BST

Scheduler accounting and receipts agree: **462/504 batches**, or
**236,544/258,048 pairs (91.7%)**, complete. No failed jobs or non-empty error
logs were found. Thirty of the 32 size-128 square Gamma batches are complete;
the last two are running. Four square control batches and 36 Aztec batches
remain queued. Analysis 18965701 is waiting for sampling completion.

The size-128 Gamma median batch time is now about 99 minutes. Using this and
size-96 cell medians extrapolated cubically for unobserved size-128 cells gives
about **4.8 hours of sampling remaining** at two workers. Allow roughly
4–6 hours plus queue delays and final analysis; remaining control-cell times
are extrapolations, not measurements. No settings changed and no partial
scientific fits were inspected. Evidence is in `progress_20261010.json` and
`status-20261010-accounting.txt` under the campaign validation directory.


## Completed, 10 October 2026

All 504 batches and 258,048 pairs completed; analysis 18965701 finished at
16:56:17 with exit zero. See [the completed results](../RESULTS.md)
for all fixed Gamma windows, controls, provenance and limitations.
