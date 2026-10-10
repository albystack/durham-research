# Completed matched Aztec–square comparison, 5 October 2026

**All 30,720 paired environments completed, and raw verification and the
prespecified analysis passed. The primary comparison does not establish an
Aztec–square difference in squared-log curvature over sizes 8–64. Both
individual primary coefficients are positive, and a control-curvature anomaly
requires follow-up. We cannot report “Aztec is log² while squares are not.”**

**7 October follow-up:** [saved-control investigation and independent confirmation](CONTROL_INVESTIGATION_20261006.md) found matching physical endpoints and exact replay of 38 original pairs. All 12,288 fresh all-one pairs are complete: contrast 0.081 [−0.282,0.437], permutation p=0.6615. The anomaly was not reproduced at this precision. Original Gamma results and control findings below remain unchanged.

## What was completed

Hamilton array 18907659 completed all 36 cells with exit 0:0; dependent job
18907660 completed raw verification and analysis with exit 0:0. All error logs
are empty. The slowest sampling cell, random square L=64, took 44 min 17 sec;
post-processing took 22 min 9 sec and finished at scheduler time 19:07:04 on
5 October. No cases were dropped or replaced.

Each geometry has 2,048 random environment pairs at each size 8,16,24,32,48,64,
plus 256 pairs per size for each of the deterministic-mean and all-one
controls. Both use independent Gamma(shape 0.2, scale 1) north-from-black and
Gamma(shape 0.25, scale 1) east-from-black weights, with other edges one.
Two fixed symmetric diagonal probes have fractions 1/8 and 1/4, matched
physical distances and height scale one. The complete design remains in
[MATCHED_SHUFFLE_20261004.md](MATCHED_SHUFFLE_20261004.md).

The square is sampled by the validated embedded shuffling construction;
Float64 reduction rounding retains its documented numerical-audit
qualification. The fresh Aztec arm uses corrected orientation and robust
creation coins. Neither result is an unqualified numerical proof.

## Frozen primary result and all cutoff windows

The fit is D_p(d)=a_p+b_p log(d)+c(log(d))², with probe-specific intercepts
and slopes. Intervals are 95% pointwise percentile intervals from 10,000
whole-environment bootstrap draws. The contrast is c_Aztec−c_square.

| Size window | Aztec c [95% interval] | Square c [95% interval] | Difference [95% interval] |
|---|---|---|---|
| **8–64, primary** | **0.238 [0.010, 0.467]** | **0.284 [0.050, 0.517]** | **−0.046 [−0.366, 0.282]** |
| 16–64 | 0.227 [−0.387, 0.837] | −0.197 [−0.806, 0.413] | 0.424 [−0.434, 1.294] |
| 24–64 | −1.391 [−2.847, 0.045] | 0.532 [−0.909, 1.978] | −1.922 [−4.022, 0.129] |

![Completed matched comparison](../figures/matched_shuffle_20261005.png)

Both primary individual intervals exclude zero, but **the direct difference
does not**. None of the predetermined contrast intervals excludes zero.
Estimates vary materially with the lower-size cutoff. These data neither
establish different asymptotic laws nor establish that the laws are identical.
The prior corrected Aztec result over 16–128 remains a separate result; its
larger size window cannot be silently substituted into this matched contrast.

## Control discrepancy: retain both diagnostics

In either deterministic environment the true disorder covariance is zero.
The prespecified joint permutation test of the individual control covariance
coordinates gives p=0.5882. It does not flag those coordinates collectively
under its maximum-coordinate statistic.

However, the primary-window **all-one curvature contrast is 1.161**, with
bootstrap interval **[0.451, 1.853]**. Its components are Aztec c=0.642
[0.147,1.131] and square c=−0.519 [−1.031,−0.005]. This is unexpected under the
known zero-disorder null and should not be dismissed because a different
test has a large p-value.

A clearly labelled **post hoc** diagnostic permuted complete second-replica
probe vectors within each control/domain/size group. It used the same linear
curvature contrasts, studentised by their exact permutation-null standard
deviations, and took the maximum over all six control contrasts (two control
models × three windows). With 9,999 permutations and seed 2026100514, the
maximum-statistic p-value is **0.0093**; the uniform primary contrast has
pointwise permutation p=0.0021 and six-contrast adjusted p=0.0093.

This follow-up was chosen after inspecting the results. Its six-contrast
adjustment does not account for the broader choice of diagnostic, so it is
exploratory evidence of a discrepancy, not a new prespecified hypothesis test
or proof of a sampler bug. It does not change or replace the primary Gamma
analysis. The two control tests use different statistics and are not
logically contradictory. Chance, dependence/provenance errors and numerical
or analysis problems must be distinguished with focused checks.

## Verification and provenance

Hamilton raw verification covered all **61,440 matchings**, **122,880 observable
values** and **138,240 checksummed files**. The compact archive contains
manifests, full measurement tables, status records, analysis, verification
certificates and scheduler logs. Its SHA256 is
`1e6e4b3a7b2b32d51ae0dffa50c294263829c3f37bf1339fc102885957346988`.
Raw arrays and matching files remain on Hamilton; they were not all downloaded
and checked a second time locally.

Local verification checked archive/table hashes, all declared counts,
seed identities and uniqueness, geometry and complete replica/probe vectors.
Independent direct variance/covariance calculations reproduce all 216 moment
estimates (72 disorder coordinates); independent least-squares calculations
reproduce all 18 fitted coefficients. Maximum discrepancies were 1.78e−15 and
2.25e−14, respectively. Original bootstrap draws were not rerun unnecessarily.
The additional control diagnostic is a separate file, not a modification of
the frozen analysis.

- Raw verified data root: `/nobackup/fvkl37/matched_shuffle_20261005/production`.
- Local compact data: `square_glauber_python/outputs/matched_shuffle_20261005/`.
- Frozen analysis: `production/analysis/analysis.json` under that local root.
- Independent checks: `docs/validation/matched_shuffle_20261005_results/`.
- Run status: `docs/validation/matched_shuffle_20261004/hamilton_status.json`.

## Next scientific steps

1. Report the completed matched-law result to Sunil with its uncertainty,
   cutoff sensitivity and control discrepancy. The email remains a draft.
2. Report the completed control audit and independent confirmation alongside
   the original discrepancy. Its primary p=0.6615 does not prove sampler
   correctness; no outcome-dependent extension is warranted.
3. Review physical probe positions, parity and finite-size sensitivity with
   Sunil before selecting a broader asymptotic window. If more computation is
   justified, freeze a new design and keep it distinct from this completed one.
4. Complete the human mathematical review and historical AI-attribution record.

One separately fixed all-one confirmation was subsequently completed; no
further Gamma scaling campaign was launched. The requested
comparison has been executed; the hoped-for scientific conclusion has not
been established.
