# Archived simulation results

This directory preserves the raw batch files and validated analyses from the
two completed LERW campaigns. It contains results only; superseded code, tests,
logs, caches, temporary report tools, and an incomplete pilot have been
removed.

## Contents

- `datasets/strict-annealed/`: 379,000 walks in 379,000 independent
  environments.
- `datasets/double-dimer/`: 758,000 walks arranged as 379,000 pairs in shared
  environments. The directory name is historical; the measured observable is
  the paired LERW winding difference
  \(\Delta W_L=W_L^{(1)}-W_L^{(2)}\).

Each campaign contains:

- `batch-results/`: the atomic simulation CSV files;
- `analysis/`: combined data, validation, summaries, exponent fits, and model
  comparisons.

The compact paired tables and publication-ready figures are also available in
[`../reports/`](../reports/).

## Checksums

```text
strict-annealed/analysis/combined_raw.csv
14ff496dc5bc3849b3006791440a2d278149a4e1e0a63038fff8573fadce425a

double-dimer/analysis/combined_raw.csv
aee0bdaafe08866b3c23ac1421fca0735217c44eefb50de6f690aada9a0da31b

double-dimer/analysis/double_dimer_pairs.csv
82a0c9b38d3149ac8c7e96f1041e2495badc164d857abd4fd648891c5e4d623d
```

Both campaigns contain 3,790 completed batch files. Do not mix either
production analysis with the removed incomplete pilot.
