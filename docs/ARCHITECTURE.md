# Architecture

The repository uses one dependency-light Julia package and thin command-line
drivers. Scientific kernels are kept separate from campaign orchestration and
from statistical analysis.

## Package modules

### `AztecDiamond`

[`aztec/src/AztecDiamond.jl`](../aztec/src/AztecDiamond.jl) is the package
entry point. It implements weighted Aztec-diamond domino shuffling, central and
spatial height observables, validation, and lightweight SVG output.

### `SquareGrid`

[`aztec/src/SquareGrid.jl`](../aztec/src/SquareGrid.jl) implements square-grid
random environments, Wilson spanning trees, the complementary dual tree,
Temperley matchings, and exact height increments. It is included as a
submodule of `AztecDiamond`.

### `GlauberSquareGrid`

[`aztec/src/GlauberSquareGrid.jl`](../aztec/src/GlauberSquareGrid.jl)
implements direct weighted-dimer height-function dynamics. It contains the
literal random-face heat bath, an exact self-loop-skipping accelerator, exact
small-system enumeration, and parallel tempering in a frozen edge environment.

## Command-line layers

```text
configs/*.csv
     |
     v
run_*_campaign.jl  --> atomic batch CSVs + metadata + diagnostics
     |
     v
merge_*_batches.jl --> validated retained table
     |
     v
analyze_*.jl       --> environment-blocked summaries and model comparisons
     |
     v
plot_*.jl          --> deterministic SVG/PNG figures
```

Campaign runners own seed derivation and restart validation. Analysis scripts
never infer missing pairing information: the environment identifier and both
replicas must be present in the input schema.

## Output layers

- `aztec/data/` contains compact retained observations suitable for Git.
- `aztec/results/` contains reviewed derived tables, reports, and figures.
- `aztec/output/` is ignored and holds local batches or scratch analyses.
- large HPC traces live outside Git and are represented by configs, metadata,
  manifests, checksums, and compact derived tables.

## Verification layers

1. exact enumeration and detailed balance on tiny state spaces;
2. structural invariants for tilings, trees, matchings, and heights;
3. reference-versus-optimized kernel comparisons;
4. deterministic replay and seed-collision checks;
5. command-line smoke workflows;
6. production diagnostics, environment-blocked bootstrap, and controls.
