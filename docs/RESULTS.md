# Completed larger matched campaign: 10 October 2026

All 504 batches and 258,048 paired environments completed. Scheduler accounting
records successful exits for all 504 sampling tasks and analysis job 18965701,
which finished at 16:56:17 on 10 October after eight minutes two seconds. No
non-empty campaign error logs were found at the 19:15 check.

## Fixed Gamma comparison

The primary fit is D_p(d) = a_p + b_p log(d) + c (log(d))^2.
Intervals below are the fixed 10,000 whole-environment bootstrap intervals.
The main finding is the positive squared-log coefficient in **each geometry**.
The additional frozen domain contrast is retained as a secondary diagnostic.

| Window | Aztec c [95% CI] | Square c [95% CI] | Difference [95% CI] |
|---|---|---|---|
| 16–128 | 0.180 [0.078, 0.277] | 0.159 [0.064, 0.255] | 0.021 [-0.118, 0.159] |
| 32–128 | 0.248 [-0.019, 0.512] | 0.076 [-0.175, 0.328] | 0.171 [-0.200, 0.530] |
| 48–128 | 0.318 [-0.322, 0.945] | 0.296 [-0.305, 0.902] | 0.022 [-0.866, 0.887] |

The primary coefficient is positive in both geometries, while the direct
difference includes zero. Both individual intervals include zero at the higher
cutoffs, and every fixed domain-contrast interval includes zero. This campaign
does not establish different growth laws, equivalence, an asymptotic log-squared
law in either geometry, or absence of such a term on the square. It is evidence
of positive fitted finite-window curvature in both geometries over 16–128.

## Controls and limitations

The predeclared joint six-contrast control permutation test gives p=0.2039.
The all-one 16–128 contrast is -0.396 [−0.778, −0.013] in its pointwise bootstrap
interval; retain this result even though the joint diagnostic does not reject.
A non-rejection does not prove sampler correctness.

The original Gamma shapes 0.2 and 0.25, scale one, and remaining weights one
were matched across geometries. Two physical diagonal probe fractions were
used. Sampling used the documented shuffling implementations and numerical
qualifications, not CFTP. This is not a replication of the historical Aztec
128–1300 window. No samples were added or removed in response to the result.

## Evidence

The compact archive SHA256 is
`90e30c0f97e0e22132110dfb35a0a6e5e7007d5bbe4c8270be3a10c9977c0254`.
It contains measurement tables, manifests, completion certificates, analysis
outputs and scheduler accounting. Raw tilings remain on Hamilton, where each
worker verified them before receiving a completion certificate. The local
verification record is [`data/verification.json`](../data/verification.json).
Bootstrap intervals and the permutation result are retained from Hamilton;
local verification independently checks the point estimates, not those resamples.

Alberto reports that he has emailed Sunil about the completed experiment.
The review should cover both positive finite-window coefficients, the corrected
sampler history, cutoff sensitivity and the precise theoretical observable.
No new sampling is scheduled by this report.

Local verification passed all 504 batch records and reproduced 252 moment
estimates, 18 coefficients and nine contrasts. Maximum moment discrepancy
was 1.78e-15 and maximum coefficient discrepancy 2.44e-14. Every fitted input
also agreed with the saved moment table.

The published measurements, checksums and original source snapshot are in
[`data/`](../data/README.md). Run `make analysis` for the point estimates and
figure, or `make bootstrap` for the full fixed resampling.

During repository cleanup, the public measurement bundle independently
reproduced the complete 10,000-resample bootstrap and 9,999-permutation control
test. All 18 coefficients and intervals and nine contrast intervals matched
to within 1e-10; the joint p-value was again 0.2039. This additional check
is recorded in `data/verification.json`.
