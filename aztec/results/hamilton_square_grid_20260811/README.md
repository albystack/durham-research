# Hamilton square-grid analysis, 11 August 2026

## Executive conclusion

The completed square-grid Temperley experiment does **not** show a robust
positive `log^2` contribution to the disorder covariance for any of the five
structured disorder laws tested.

The primary comparison is the covariance-aware pooled block-GLS fit

```text
a + b log(r)
```

against

```text
a + b log(r) + c (log(r))^2.
```

Across five disorder laws and three prespecified lower fit cutoffs, 14 of the
15 bootstrap intervals for `c` contain zero, and 14 of the 15 point-estimate
BIC comparisons favour the simpler affine-log model. The sole exception is a
**negative**, cutoff-specific coefficient for directed Gamma `k=0.5` at
`L >= 512`; it disappears at the adjacent `L >= 128` and `L >= 1024` windows.
It is therefore not evidence for positive super-rough growth.

This is a result about the structured spanning-tree/Temperley disorder laws in
this repository. It does not test the canonical model with independent random
energies on every square-lattice dimer edge.

## Experiment and primary method

- Geometry: square-grid generalized Temperley dimer model with wired boundary.
- Observable: paired spatial height increments at separations `1/32`, `1/16`,
  `1/8`, and `1/4` of the full side.
- Disorder component: covariance of two conditionally independent replicas
  sharing one environment.
- Orders: `L = 32, 64, 128, 256, 512, 768, 1024, 1536, 2048, 2560, 3072,
  4096, 5120, 6144`.
- Sample count: 39,570 independent environments per model (38,600 through
  `L=2048`, plus 970 through `L=6144`), with the no-disorder baseline reused as
  the control.
- Uncertainty: 10,000 joint environment bootstraps, preserving both replicas
  and all four separations.
- Primary fit: block GLS with one bootstrap-estimated `4 x 4` cross-separation
  covariance matrix per order and independent covariance blocks across orders.
- Sensitivity windows: minimum orders 128, 512, and 1024.
- Prediction check: the two largest orders, 5120 and 6144, are held out.

## Primary block-GLS results

Positive `delta BIC = BIC(log) - BIC(log+log^2)` favours the quadratic
extension.

| Disorder law | Minimum L | c | 95% bootstrap interval | P(c>0) | delta BIC | Held-out RMSE log / log+log^2 |
|---|---:|---:|---:|---:|---:|---:|
| Directed Gamma k=0.5 | 128 | -0.002513 | [-0.007742, 0.002546] | 0.1657 | -2.977 | 0.096 / 0.092 |
| Directed Gamma k=0.5 | 512 | -0.02480 | [-0.04424, -0.006096] | 0.0050 | 2.851 | 0.095 / 0.086 |
| Directed Gamma k=0.5 | 1024 | -0.004031 | [-0.05530, 0.04446] | 0.4162 | -3.441 | 0.079 / 0.088 |
| Directed Gamma k=1 | 128 | 0.001144 | [-0.003937, 0.005915] | 0.6590 | -3.670 | 0.096 / 0.096 |
| Directed Gamma k=1 | 512 | 0.01189 | [-0.006760, 0.03014] | 0.8850 | -2.091 | 0.103 / 0.098 |
| Directed Gamma k=1 | 1024 | 0.02968 | [-0.02061, 0.07804] | 0.8738 | -2.118 | 0.104 / 0.102 |
| Directed lognormal, variance 2 | 128 | 0.001218 | [-0.003915, 0.006297] | 0.6655 | -3.649 | 0.148 / 0.147 |
| Directed lognormal, variance 2 | 512 | -0.007155 | [-0.02596, 0.01052] | 0.2054 | -3.109 | 0.155 / 0.155 |
| Directed lognormal, variance 2 | 1024 | -0.03973 | [-0.09073, 0.01054] | 0.0581 | -1.074 | 0.197 / 0.208 |
| Directed Uniform(0,2) | 128 | 0.003370 | [-0.001555, 0.008088] | 0.9090 | -1.992 | 0.165 / 0.161 |
| Directed Uniform(0,2) | 512 | -0.004514 | [-0.02134, 0.01212] | 0.2896 | -3.413 | 0.152 / 0.167 |
| Directed Uniform(0,2) | 1024 | -0.001267 | [-0.04652, 0.04222] | 0.4529 | -3.463 | 0.163 / 0.200 |
| Undirected Gamma k=0.5 | 128 | 0.003652 | [-0.001268, 0.008608] | 0.9253 | -1.778 | 0.117 / 0.115 |
| Undirected Gamma k=0.5 | 512 | 0.006979 | [-0.01159, 0.02509] | 0.7703 | -3.114 | 0.109 / 0.110 |
| Undirected Gamma k=0.5 | 1024 | 0.02425 | [-0.02417, 0.07264] | 0.8337 | -2.502 | 0.108 / 0.122 |

The full machine-readable table, including conditional and no-disorder
components, is in `summary/spatial_gls_campaign_comparison.csv`.

## Interpretation and robustness

1. **No positive coefficient is stable across fit windows.** The small positive
   point estimates for some laws all have intervals spanning zero and negative
   BIC differences.
2. **The one nonzero interval is negative and unstable.** Directed Gamma
   `k=0.5` at `L >= 512` gives `c=-0.02480`, but its bootstrap delta-BIC
   interval is `[-3.279, 17.114]`, the quadratic wins only 75.34% of bootstrap
   draws, and the coefficient returns to zero at `L >= 1024`.
3. **Held-out prediction is mixed.** The quadratic has lower held-out RMSE in
   8 of 15 comparisons and the affine-log model in 7 of 15, without a stable
   pattern across laws or cutoffs. Most differences are small.
4. **The no-disorder control is clean.** Its marginal block-GLS coefficient
   interval contains zero at every cutoff in every analysis.
5. **The result differs from the finite-size Aztec signal.** The earlier Aztec
   analysis found positive `c` at the `L >= 128` window for the original and
   stronger Gamma laws, but not for the weak or symmetric laws, and the Aztec
   intervals widened substantially at higher cutoffs. The two geometries and
   disorder constructions should be reported separately.

The statistically defensible wording is therefore: **within the simulated
range and for these structured square-grid laws, the data do not support a
robust positive log-squared disorder contribution; ordinary logarithmic growth
is the more economical description.** This is finite-size numerical evidence,
not an asymptotic proof.

## Completion and integrity checks

- All 19 packed recovery arrays, 8 release jobs, 12 original validation
  arrays, and 15 analysis jobs completed with exit code zero.
- No failed, out-of-memory, preempted, or timed-out final jobs were found.
- The high-order production tree contains exactly 5,820 batch CSVs, 5,820
  diagnostic CSVs, and 5,820 execution-provenance files.
- The standard robustness tree contains exactly 3,596 files of each type.
- Neither production tree contains a leftover `.tmp` file.
- The downloaded analysis contains 180 nonempty files: 15 reports, 15
  block-GLS tables, the supporting tables, and SVG figures.
- Recursive SHA-256 of the downloaded analysis tree and the Hamilton source:
  `6239646814d5622fa241b5a326013bb96a7f328289e6d39fed4fae725c31b117`.
- The downloaded Hamilton job maps and task manifests retain their original
  SHA-256 records in `manifests/`.

## Package layout

- `analysis/`: all 15 Hamilton analysis directories, including sensitivity
  tables, reports, and SVG figures.
- `summary/spatial_gls_campaign_comparison.csv`: consolidated 90-row block-GLS
  table (five laws x three cutoffs x two datasets x three components).
- `summary/spatial_gls_campaign_comparison.md`: compact primary 15-row table.
- `EMAIL_TO_SUNIL.txt`: the brief, copy-ready email draft.
- `EMAIL_TO_SUNIL.md`: the longer email draft with additional numerical detail.
- `manifests/`: task lists, scheduler job maps, and checksums from the recovery.
- `CHECKSUMS.sha256`: SHA-256 inventory for the 248 files in this local
  handoff package (excluding the checksum inventory itself).
