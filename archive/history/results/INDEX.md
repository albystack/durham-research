# Reviewed results

This directory is the public index of compact, reviewed results. Large raw
traces and resumable batches remain in ignored local output directories or in
the recorded Hamilton `/nobackup` campaign root.

## Aztec diamond

The retained Aztec analyses are under [`aztec/`](aztec/INDEX.md). The primary
finite-size signal is in the shared-environment covariance of spatial height
increments:

- [`aztec/spatial/`](aztec/spatial/) contains the pooled spatial fits, cutoff
  comparisons, uniform negative control, and variance-decomposition figures;
- [`aztec/double_dimer/`](aztec/double_dimer/) contains the earlier paired
  central-height analysis;
- [`aztec/height/`](aztec/height/) contains the single central-height analysis;
- [`aztec/hamilton_square_grid_20260811/`](aztec/hamilton_square_grid_20260811/)
  contains the large structured square-grid comparison.

The positive Aztec quadratic-log coefficients are numerical finite-size
evidence, not a proof. Exact values and qualifications are recorded in the
[results ledger](../docs/RESULTS.md).

## Direct square-grid dimers

Earlier Julia Glauber and finite-volume Kasteleyn results are retained under:

- [`aztec/glauber_square_grid_20260822/`](aztec/glauber_square_grid_20260822/)
- [`aztec/glauber_square_grid_kasteleyn_20260822/`](aztec/glauber_square_grid_kasteleyn_20260822/)

The active square-grid CFTP campaign has not yet produced a reviewed scaling
result. Its code, frozen schedule, and analysis plan are under
[`../square_glauber_python/`](../square_glauber_python/PROJECT_GUIDE.md). A new
result directory will be added here only after certification, aggregation,
control checks, and the pre-specified analysis complete.

## Historical LERW results

Superseded loop-erased-random-walk results remain in
[`../archive/random_walk/results/`](../archive/random_walk/results/RESULTS_INDEX.md)
for scientific provenance.
