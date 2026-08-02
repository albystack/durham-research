# Numerical analysis report

This report summarizes the retained Gamma-disordered Aztec-diamond campaigns
completed on 2 August 2026. It records finite-size numerical evidence, not an
asymptotic theorem.

## Data

The single-height dataset contains 35,536 independent environments through
`L = 1300`. Orders 900 and 1,000 have 1,000 observations each; orders 1,100,
1,200, and 1,300 have 512 each.

The double-dimer dataset contains 28,304 independent environments and 56,608
conditional tilings. Each environment supplies two conditionally independent
tilings with centre heights `H1` and `H2`. There are 256 pairs per order from
512 through 1,300; smaller orders have 1,000-4,000 pairs.

## Full-range comparisons

All intervals below use 10,000 deterministic within-size bootstrap resamples.
`Delta BIC` is `BIC(log) - BIC(log^2)`, so positive values favour the affine
squared-log curve.

| Observable, fitted over `L >= 24` | Effective `p` in `C(log L)^p` | Delta BIC | Interpretation |
|---|---:|---:|---|
| Single height `Var(H)` | 0.967 `[0.894, 1.032]` | 4.434 `[-1.659, 6.479]` | Affine curves unresolved |
| Double difference `Var(H1-H2)` | 0.814 `[0.700, 0.902]` | -1.928 `[-2.890, -0.277]` | Ordinary log preferred |
| Paired covariance `Cov(H1,H2)` | 1.393 `[1.156, 1.558]` | 2.025 `[0.061, 3.378]` | Upward disorder curvature |

The paired experiment uses the exact identities

```text
Var(H1-H2)/2 = E[Var(H | environment)],
Cov(H1,H2)   = Var(E[H | environment]).
```

Conditional tiling noise is consistent with ordinary logarithmic growth. The
disorder component behaves differently: its effective exponent is
`1.859 [1.391, 2.226]` for `L >= 96` and `2.103 [1.341, 2.742]` for
`L >= 192`. At those higher cutoffs, however, the affine log and squared-log
curves are indistinguishable by BIC and intervals are wide.

## Interpretation

1. The double-dimer height difference does not display squared-log growth over
   the simulated range; ordinary log is preferred in the full-range fit.
2. The paired covariance identifies the environment as the source of the
   upward curvature in total single-height variance.
3. The covariance signal is suggestive, but only 256 paired environments are
   available at each large order. More high-order pairs are needed before
   treating the cutoff trend as stable evidence.

The BIC calculations are unweighted least-squares diagnostics on noisy
size-level estimates, not raw-data likelihood comparisons. The free-power fit
has no additive constant. These qualifications are important when comparing
the reported exponent with a theoretical squared-log prediction.

## Suggested continuation

The most informative next computation is additional precision at the existing
large orders rather than a new, still larger order. The schedule
[`../configs/double_dimer_large_continuation.csv`](../configs/double_dimer_large_continuation.csv)
adds 256 independent pairs at each order from 512 through 1,300 without
reusing sample IDs.

Principal figures:

- [`height/center_height_variance.svg`](height/center_height_variance.svg)
- [`double_dimer/double_dimer_variance_fits.svg`](double_dimer/double_dimer_variance_fits.svg)
- [`double_dimer/disorder_covariance_fits.svg`](double_dimer/disorder_covariance_fits.svg)
- [`double_dimer/variance_decomposition.svg`](double_dimer/variance_decomposition.svg)
