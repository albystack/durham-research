# Historical reference code

[`glauber_reference.jl`](glauber_reference.jl) is the original exploratory
height-function prototype used to recover the tileable boundary and local
update conventions for the direct square-grid sampler.

It is retained as provenance only. It is not imported by the package, is not a
production runner, and does not provide deterministic seeding or restart-safe
output. The tested implementation is
[`../src/GlauberSquareGrid.jl`](../src/GlauberSquareGrid.jl).
