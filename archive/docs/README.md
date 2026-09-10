# Documentation

This directory is the public research record for the repository. It separates
stable scientific definitions, numerical outcomes, current work, and
reproduction instructions so that historical decisions do not obscure the
active code.

## Start here

1. [Research overview](RESEARCH_OVERVIEW.md) — the question, models,
   observables, and statistical estimands.
2. [Results ledger](RESULTS.md) — chronological experiments, numerical
   conclusions, and limitations.
3. [Roadmap](ROADMAP.md) — current validation priorities and decision gates.
4. [Reproducibility](REPRODUCIBILITY.md) — seeds, environments, output layers,
   tests, and analysis commands.
5. [Architecture](ARCHITECTURE.md) — how the Julia modules and command-line
   workflows fit together.

## Model-specific documentation

- [Aztec implementation](../aztec/docs/IMPLEMENTATION.md)
- [Square-grid model contract](../aztec/docs/SQUARE_GRID_MODEL.md)
- [Hamilton campaign record](../aztec/docs/HAMILTON_CAMPAIGN_20260808.md)
- [HPC workflow](../hpc/README.md)

## Documentation policy

- Mathematical conventions belong in model contracts, not only in code.
- Numerical claims must link to retained tables or identify the external raw
  campaign root from which they were derived.
- Results are labelled as finite-size evidence unless a theorem is cited.
- Failed pilots and mixing limitations remain documented when they affect the
  interpretation of later runs.
- Personal correspondence, credentials, and machine-specific private paths do
  not belong in the repository.
