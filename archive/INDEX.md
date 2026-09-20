# Archive

Historical material is retained for provenance. The two maintained projects
are [Aztec](../aztec/) and [Glauber dynamics](../square_glauber_python/).
Archived reports are dated records, not statements of current campaign status.

| Location | Contents |
| :--- | :--- |
| [experiments/square_grid/](experiments/square_grid/) | Frozen Julia source snapshot for spanning-tree, Temperley, and temporal experiments, including their dependencies, tests, configs, and old HPC wrappers. |
| [random_walk/](random_walk/) | Earlier loop-erased random-walk project. |
| [history/results/](history/results/INDEX.md) | Previous figures, compact result tables, checksums, and campaign manifests. |
| [history/docs/](history/docs/) | Earlier research overview, results ledger, roadmap, and project guide. |
| [history/aztec/](history/aztec/) | Earlier Julia guide, campaign notes, and Mathematica rendering code. |
| [history/square_glauber_python/](history/square_glauber_python/) | Dated launch/status notes, the provisional September report, and its one-off analysis scripts. |
| [history/research_materials/](history/research_materials/) | Local correspondence, papers, notes, and bibliography; existing private-file ignore rules are retained. |
| [hpc_recovery_202608/](hpc_recovery_202608/) | One-off August scheduler recovery scripts. |
| [previous_aztec_result_26050/](previous_aztec_result_26050/) | Superseded 26,050-sample plot. |
| `local_runs/` | Previous local Julia and Python outputs; kept on disk and ignored by Git. |

The [provisional report](history/square_glauber_python/reports/provisional_20260915/REPORT.md)
and [earlier results ledger](history/docs/RESULTS.md) retain their original
limitations. No result was promoted or reinterpreted during reorganisation.

## Finding and restoring files

[RELOCATION.csv](RELOCATION.csv) maps original paths to their new locations
and records SHA-256 at relocation. The separate local-only
`local_runs/RELOCATION.csv` inventories ignored run products. A `copy` entry
denotes a frozen dependency snapshot; active files may subsequently evolve.
Archive contents are historical originals, so embedded commands and some
cross-links still use their original repository paths. Consult the mapping
when reproducing an old analysis.

The tree-model snapshot retains an `aztec/` project inside its own root.
For its historical unit tests, explicitly opt in from the repository root:

```bash
cd archive/experiments/square_grid
julia --project=aztec -e 'using Pkg; Pkg.test()'
```

Historical analysis commands may additionally need datasets from the active
`aztec/data/`, archived results, or external campaign storage. The snapshot
does not duplicate those datasets. It is excluded from active test discovery
and CI. Avoid writing new results into the archive.
