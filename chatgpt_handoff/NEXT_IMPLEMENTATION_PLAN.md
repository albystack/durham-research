# Next Implementation Plan

## Stage A - preserve and validate current work

Goal: keep the current Aztec pipeline frozen while the square-grid work is designed.

Expected files:

- `aztec/src/AztecDiamond.jl`
- `aztec/test/runtests.jl`
- `aztec/test/smoke_workflows.sh`
- `aztec/README.md`
- `aztec/docs/IMPLEMENTATION.md`

Checks:

- `julia --project=aztec -e 'using Pkg; Pkg.test()'`
- `sh aztec/test/smoke_workflows.sh`

## Stage B - formalise the square-grid mathematical model

Goal: write down the exact square-grid dimer/spanning-tree correspondence before coding.

Expected files:

- `aztec/docs/SQUARE_GRID_MODEL.md` or an equivalent design note
- `aztec/README.md` update, if needed

Decisions to lock:

- directed site weights versus shared undirected conductances
- Temperley correspondence and boundary conditions
- height convention and parity convention
- paired replica design and observable definitions

Checks:

- a short design review against the existing Aztec conventions
- a written checklist of unresolved assumptions

## Stage C - implement paired square-grid configurations

Goal: add a minimal square-grid sampler that can produce paired shared-environment samples.

Expected files:

- `aztec/src/SquareGrid.jl` or another clearly named module
- `aztec/scripts/run_square_grid_campaign.jl`
- `aztec/scripts/merge_square_grid_batches.jl`

Checks:

- one small exact-size sample
- seed determinism across reruns
- paired-stream independence within one environment

## Stage D - implement spatial heights and tests

Goal: compute the square-grid spatial increment observable and validate it.

Expected files:

- `aztec/src/SquareGrid.jl`
- `aztec/test/square_grid_runtests.jl`

Checks:

- direct-vs-full height agreement
- finite-sample variance identity for paired replicas
- parity/location consistency on a tiny grid

## Stage E - create a small pilot

Goal: run a tiny reproducible campaign before any scale-up.

Expected files:

- `aztec/configs/square_grid_smoke.csv`
- `aztec/output/square_grid_smoke/` only as scratch output

Checks:

- batch resume safety
- duplicate-seed rejection
- validated merged CSV schema

## Stage F - analyse the pilot

Goal: reproduce the Aztec-style bootstrap and fitting workflow on the square-grid pilot.

Expected files:

- `aztec/scripts/analyze_square_grid_campaign.jl`
- `aztec/scripts/plot_square_grid_campaign.jl`
- `aztec/results/square_grid/` or a pilot analysis folder

Checks:

- within-size bootstrap
- log and squared-log fits
- holdout prediction if multiple sizes are available

## Stage G - prepare Hamilton Slurm execution

Goal: make the run portable to Durham Hamilton.

Expected files:

- `scripts/` or `hpc/` job wrappers if the repository adopts them
- a short environment note describing Julia, modules, and thread count

Checks:

- dry-run command lines
- batch size and memory estimates
- restart semantics documented clearly

## Stage H - design the production campaign

Goal: define the large square-grid schedule only after the pilot is stable.

Expected files:

- `aztec/configs/square_grid_campaign.csv`
- `aztec/configs/square_grid_large_continuation.csv`
- `aztec/results/square_grid/` for the retained analysis

Checks:

- no sample ID reuse across stages
- clear separation between pilot, continuation, and production
- a written note explaining which results are current and which are superseded
