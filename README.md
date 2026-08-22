# Random-environment dimers and spanning trees

Reproducible Julia experiments for height fluctuations in random tilings,
spanning trees, and loop-erased random walks. The project tests whether
variance grows like \(\log L\) or develops a super-rough \((\log L)^2\)
component, with particular attention to separating conditional sampling noise
from fluctuations induced by a shared random environment.

The repository combines exact small-system checks, deterministic Monte Carlo
campaigns, environment-blocked bootstrap inference, and restart-safe Slurm
workflows. The active Julia package has no third-party runtime dependencies.

## Research highlights

- **Aztec diamonds:** paired spatial-height increments show positive
  finite-size quadratic-log curvature in the disorder covariance for the
  original and stronger Gamma laws.
- **Structured square-grid disorder:** spanning-tree/Temperley experiments up
  to \(L=6144\) do not show a stable positive quadratic-log contribution.
- **Direct weighted-dimer dynamics:** a 1,312-environment square-grid Glauber
  campaign separates conditional, disorder, and total central-height
  variance. The fitted quadratic coefficient is compatible with zero, with a
  documented mixing caveat at \(L=16,20\).
- **Negative controls:** uniform or all-one environments do not create a
  spurious disorder component.

These are finite-size numerical findings, not asymptotic proofs. Exact
estimates, uncertainty intervals, diagnostics, and limitations are recorded in
the [results ledger](docs/RESULTS.md).

## Quick start

Requirements: Julia 1.10 or newer.

```bash
cd research

# Run mathematical-reference and workflow tests.
julia --project=aztec -e 'using Pkg; Pkg.test()'
sh aztec/test/smoke_workflows.sh
```

Run a small deterministic Aztec-height campaign:

```bash
JULIA_NUM_THREADS=4 julia --project=aztec \
  aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_smoke.csv \
  --output-dir aztec/output/height_smoke
```

Generated batches and scratch analyses belong under the ignored
`aztec/output/` directory. Retained observations and reviewed results live in
`aztec/data/` and `aztec/results/`.

## Core estimators

For two conditionally independent replicas \(H_1,H_2\) in the same frozen
environment \(\omega\),

\[
\frac12\operatorname{Var}(H_1-H_2)
=\mathbb E_\omega[\operatorname{Var}(H\mid\omega)],
\]

and

\[
\operatorname{Cov}(H_1,H_2)
=\operatorname{Var}_\omega(\mathbb E[H\mid\omega]).
\]

Every bootstrap therefore resamples whole environments. Paired replicas and
all observables from one environment remain in the same resampling block.

## Repository map

```text
.
├── aztec/                 active Julia package, data, results, and CLIs
│   ├── src/               samplers, graph constructions, and observables
│   ├── scripts/           campaign, merge, analysis, and plotting commands
│   ├── configs/           deterministic smoke, pilot, and production schedules
│   ├── test/              exact, statistical, and end-to-end checks
│   ├── data/              compact retained input datasets
│   ├── results/           reviewed tables, reports, and vector figures
│   └── reference/         historical prototype retained for provenance
├── docs/                  research overview, results, roadmap, reproducibility
├── hpc/                   Slurm wrappers and cluster workflow documentation
├── archive/               self-contained earlier LERW experiments
├── .github/workflows/     continuous integration
└── CONTRIBUTING.md        scientific and engineering contribution rules
```

Start with the [documentation index](docs/README.md), then use the
[package guide](aztec/README.md) for model-specific commands.

## Reproducibility guarantees

- deterministic seed derivation from campaign identifiers;
- independent random streams where required by the estimand;
- environment-level pairing preserved end to end;
- atomic, restart-safe output batches;
- explicit campaign metadata and execution provenance;
- exact enumeration and detailed-balance checks at small sizes;
- reference-vs-optimized observable tests;
- uncertainty and model comparisons repeated across fit windows;
- controls analyzed with the same pipeline as disordered models.

See [Reproducibility](docs/REPRODUCIBILITY.md) for data flow, validation levels,
and the commands used to regenerate retained analyses.

## Documentation

- [Research overview](docs/RESEARCH_OVERVIEW.md) — models, observables, and
  statistical definitions.
- [Results ledger](docs/RESULTS.md) — chronological numerical evidence and
  caveats.
- [Roadmap](docs/ROADMAP.md) — current validation work and decision gates.
- [Implementation guide](aztec/docs/IMPLEMENTATION.md) — sampler and data-flow
  details.
- [Square-grid model contract](aztec/docs/SQUARE_GRID_MODEL.md) — graph,
  matching, height, and environment conventions.
- [HPC workflow](hpc/README.md) — generic Hamilton/Slurm setup and submission.

## Scope

The repository supports numerical research and reproducible analysis. Large
raw HPC traces are intentionally kept outside Git; compact retained datasets,
derived tables, figures, checksums, and exact campaign configurations are kept
here when practical.

Citation metadata for the software is available in [`CITATION.cff`](CITATION.cff).
