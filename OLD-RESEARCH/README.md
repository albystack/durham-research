# Archived research artifacts

This directory contains every non-code artifact moved out of the active
`research-julia/` package on 24 July 2026. Nothing was deleted during the
cleanup.

The active implementation remains in `../research-julia/`. Frozen experiment
configurations remain beside the code in `../research-julia/configs/`.

## Contents

### `datasets/strict-annealed/`

- `batch-results/`: 3,790 atomic result CSVs containing 379,000 walks from
  379,000 independently seeded environments.
- `analysis/`: validated combined data, task audit, 137 cell summaries,
  15 exponent fits, model comparisons, and diagnostics.

### `datasets/double-dimer/`

- `batch-results/`: 3,790 atomic result CSVs containing 758,000 raw walks.
- `analysis/`: validated raw data, 379,000 paired winding differences, task
  audit, summaries, fits, model comparisons, and diagnostics.

### `datasets/double-dimer-pilot/`

- `batch-results-incomplete/`: 55 completed pilot tasks from a 137-task design.
  This pilot is incomplete and superseded by the full production campaign.

### `reports/`

- `strict-annealed/`: aggregate tables, PNG/PDF figures, email draft, short
  scientific report, and Professor Chhita handoff ZIP.
- `double-dimer/`: aggregate tables, figures, email draft, and the complete
  379,000-row pair table.

### `documentation/`

- `julia-technical-report/`: the 15 July 2026 LaTeX source, compiled 40-page
  technical PDF, and its build directory. It documents the strict campaign and
  predates completion of the double-dimer production run.

### `tests/`

- `julia/`: the former Julia test suite. It passed all 88 assertions immediately
  before archival on 24 July 2026.

### `legacy-code/`

- `old-python/`: archived first Python/HPC implementation. It is not part of
  the active workflow.

### `logs/`

- `strict-campaign.log`: historical strict campaign attempts and resumptions.
- `strict-final-pipeline.log`: final strict analysis and plotting pipeline.
- `double-dimer-campaign.log`: completed paired production campaign.

Intermediate errors in historical logs do not invalidate later completed
outputs; use each campaign's `analysis/validation.csv` as the final audit.

### `temporary/`

- `report-generation/`: local scripts/assets used to produce report files.
- `scripts-python-cache/`: an archived Python bytecode cache found under the
  Julia scripts directory.

## Authoritative checksums

```text
Strict frozen config
8a8deba4422603bc4315d193bb914db2a8a92dc2f86870379f810e30508d76d0

Strict combined raw data
14ff496dc5bc3849b3006791440a2d278149a4e1e0a63038fff8573fadce425a

Double frozen config
ed8d79e9cd52a785200d3df188d16d0fc353ef85200555c7ff1a7b194a93a0b0

Double combined raw data
aee0bdaafe08866b3c23ac1421fca0735217c44eefb50de6f690aada9a0da31b

Double paired data
82a0c9b38d3149ac8c7e96f1041e2495badc164d857abd4fd648891c5e4d623d
```

## Preservation notes

- Do not mix the incomplete pilot with either production analysis.
- The raw result directories are large and mostly absent from normal Git
  history. Do not run broad cleanup commands against this archive.
- Archived documents may contain historical paths from before this
  reorganization. The current locations in this index are authoritative.
- `../INFORMATION.txt` contains the complete scientific and technical handoff.
