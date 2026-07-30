# Random-weight Aztec diamond

This directory reproduces the workflow in Sunil Chhita's supplied 2022
`simulatorfinal.jl` and `drawingtools.nb`:

1. make a `2n × 2n` table of edge weights;
2. reduce the weights to domino-creation probabilities;
3. sample an order-`n` Aztec-diamond tiling by deletion, sliding, and creation;
4. write the binary matrix expected by the old Mathematica notebook; and
5. render the same four-colour, `-45°`-rotated picture.

The requested random-weight run uses the precise example in the old Julia
file:

```julia
rand(Float64, (2 * n, 2 * n))
```

Thus all `4n²` entries are independent `Uniform(0,1)` weights. The current
reproducible run uses `n = 200`, seed `20260728`, and Julia's `Xoshiro`
generator.

## Reproduce

No packages beyond the Julia standard library are needed:

```bash
julia --startup-file=no aztec/test/runtests.jl
julia --startup-file=no aztec/scripts/run_random_weights.jl
julia --startup-file=no aztec/scripts/run_gamma_disordered.jl
```

The default output directory is
`aztec/output/random_uniform_n200_seed20260728/`. It contains:

- `weights_uniform_0_1.txt`: the saved `400 × 400` input table;
- `tiling_matrix.txt`: the saved `400 × 400` binary output matrix;
- `tiling_mathematica_style.svg`: a directly viewable rendering; and
- `run_metadata.txt`: the parameters, validation results, timings, and SHA-256
  checksums.

The original notebook is not required to view the SVG. To use Mathematica,
open and evaluate `mathematica/draw_tiling.wl`; it imports the same output
matrix and exports a PNG.

For the completed macOS run, the SVG was also rendered to
`tiling_mathematica_style.png` as a convenient email attachment.

## Matrix convention

There is one `1` (or `true` internally) per domino. If row and column have the
same parity, the domino is horizontal in the pre-rotation coordinates;
otherwise it is vertical. The parity classes reproduce the notebook colours:

| Row | Column | Colour |
|---|---|---|
| odd | odd | green |
| odd | even | blue |
| even | even | yellow |
| even | odd | red |

For order `n`, a valid matrix therefore has exactly `n(n+1)` ones and covers
all `2n(n+1)` unit cells of the Aztec diamond without gaps or overlaps. The
runner checks all of these conditions before writing output.

## Duits–Van Peski Gamma model

The Gamma follow-up implements Definition 1.1 of Duits and Van Peski,
[*The Gamma-disordered Aztec diamond*](https://arxiv.org/abs/2512.03033)
(arXiv:2512.03033). It samples two
independent `n × n` families

```math
a_{ij}\sim\Gamma(\alpha,1),\qquad
b_{ij}\sim\Gamma(\beta,1),
```

with the other two edge families gauge-fixed to 1, and computes the shuffled
weights using their equations (1.22)–(1.23). The default run uses
`n=200`, `α=0.2`, and `β=0.25`, matching the parameters illustrated in their
Figure 2. Its output is written under
`aztec/output/gamma_disordered_alpha_0p2_beta_0p25_n200_seed20260728/`.

## Verification against the supplied files

The source archive was `/Users/alberto/Desktop/SimulatorJulia.zip` (SHA-256
`1cc0cd228428e34fa10858b2c137e7d0869e15d846ba99afb5d2bb187d743aa3`).
The reimplementation was compared directly with its `simulatorfinal.jl` at
orders 1, 2, 5, 12, and 30. Given the same random stream, every creation
probability and every entry of the final tiling matrix was identical. The
automated test suite separately checks orders 1 through 12.
