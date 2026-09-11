# Architecture

The repository has one active Python perfect-sampling package and one
established dependency-light Julia package. Scientific kernels are kept
separate from campaign orchestration and statistical analysis.

## Active Python package

[`square_glauber_python/`](../square_glauber_python/PROJECT_GUIDE.md) contains
the authoritative random-scan heat-bath reference backend, an exactly
equivalent optional Numba backend, dimer heights, exact enumeration, monotone
CFTP, frozen-environment campaign tooling, and environment-blocked analysis.
The perfect sampler publishes a state only after the true minimum and maximum
height tilings coalesce under the same replayable random maps.

## Julia package modules

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
small-system enumeration, parallel tempering in a frozen edge environment,
and finite-volume Kasteleyn moments for arbitrary valid dual height paths.

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

The determinantal square-grid runner follows the same seed and atomic-batch
contract, but records one conditional mean and variance per environment. Its
analysis combines the mean-of-conditional-variances and the variance of
conditional means, so no synthetic replica samples are introduced.

## Output layers

- `aztec/data/` contains compact retained observations suitable for Git.
- `results/aztec/` contains reviewed derived tables, reports, and figures.
- `aztec/output/` is ignored and holds local batches or scratch analyses.
- `square_glauber_python/outputs/` is ignored and holds local validation and
  diagnostic artifacts.
- large HPC traces live outside Git and are represented by configs, metadata,
  manifests, checksums, and compact derived tables.

## Verification layers

1. exact enumeration and detailed balance on tiny state spaces;
2. structural invariants for tilings, trees, matchings, and heights;
3. reference-versus-optimized kernel comparisons;
4. deterministic replay and seed-collision checks;
5. command-line smoke workflows;
6. production diagnostics, environment-blocked bootstrap, and controls.
