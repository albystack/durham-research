# Open Questions

## Critical

1. Directed site weights versus shared undirected conductances. The user explicitly flagged that the earlier square-grid implementation may be using the wrong weight semantics. This must be resolved before any square-grid extension is trusted.
2. Exact square-grid Temperley correspondence. The repository does not yet contain a square-grid counterpart, so the correct bijection and observable convention still need to be pinned down.
3. Square-grid dimer height convention. The active Aztec convention is explicit, but the square-grid version is still undefined.

## Important

1. Spatial location and parity conventions. The current Aztec spatial experiment uses symmetric central-row separations; the square-grid counterpart must match this carefully or the observables will not be comparable.
2. Independence of replica streams. The current code carefully separates environment draws from creation coins; the square-grid implementation must preserve the same independence structure.
3. Finite-size fitting sensitivity. The retained squared-log evidence is finite-size numerical evidence only; it may shift under additional high-order samples or a different cutoff window.
4. Current manuscript status. No manuscript file, submission draft, or formal project outline is present in the repository.
5. Missing validation for square-grid work. The current test suite covers only the Aztec kernel, not the requested square-grid extension.

## Minor

1. Stale or superseded documentation. Some output directories and historical summaries remain on disk and should be interpreted as archival, not current authority.
2. Output tree sprawl. `aztec/output/` contains many historical scratch campaigns; they are useful for provenance but not for upload.
3. Julia version wording. `README.md` says Julia 1.10 or newer, while the manifest records Julia 1.12.6 and CI tests 1.10 and 1.
4. Archive clarity. `old/random_walk/` is clearly historical, but it is easy to confuse it with the active Aztec code if the directory is not read carefully.
