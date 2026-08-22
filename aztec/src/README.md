# Source modules

The active package is intentionally split by mathematical construction.

- [`AztecDiamond.jl`](AztecDiamond.jl) is the package entry point and contains
  weighted domino shuffling, height observables, validation, and rendering.
- [`SquareGrid.jl`](SquareGrid.jl) contains random-environment spanning trees,
  complementary dual trees, generalized Temperley matchings, and spatial
  height increments.
- [`GlauberSquareGrid.jl`](GlauberSquareGrid.jl) contains direct square-grid
  dimer heat-bath dynamics, exact enumeration, event-driven acceleration, and
  replica exchange.

Public functions use explicit random-number generators. Reference and
optimized paths are kept separate wherever exact comparison is possible.
Command-line parsing, files, and scheduler concerns belong in `../scripts/`,
not in these modules.
