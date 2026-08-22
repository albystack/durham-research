# Aztec datasets

## `examples/`

The two order-200 examples contain the exact generated inputs, binary tiling
matrix and run metadata:

- `uniform_L200/`: one `400 × 400` i.i.d. `Uniform(0,1)` weight table;
- `gamma_L200/`: the two `200 × 200` Gamma weight arrays used by the
  Duits--Van Peski model.

Each `tiling_matrix.txt` contains one `1` per domino.

## `height/center_height_samples.csv`

This is the combined production dataset for the center-height experiment.
The columns are:

| Column | Meaning |
|---|---|
| `order` | Aztec-diamond order \(L\) |
| `sample_id` | deterministic within-order observation number |
| `seed` | unsigned 64-bit Xoshiro seed |
| `center_row` | \(L+1\) in the staggered height table |
| `center_column` | \(\lfloor L/2\rfloor+1\) |
| `center_height` | sampled height at that face |

There are 35,536 data rows, covering orders 16 through 1,300, with no repeated
`(order, sample_id)` keys and no repeated seeds. Orders 900 and 1,000 have
1,000 observations each; orders 1,100, 1,200, and 1,300 have 512 each.
`campaign_metadata.txt` records the model
parameters, seed convention, source campaigns, counts, and dataset checksum.

## `double_dimer/pairs.csv`

Each of the 28,304 rows represents one fresh Gamma environment and two tilings
sampled independently conditional on it. The columns are `order`,
`sample_id`, `seed`, center-face indices, `height_1`, `height_2`, and the
validated difference `height_1-height_2`. The data cover orders 16 through
1,300. They permit two complementary estimates:

```text
Var(height_1-height_2)/2   conditional tiling variance
Cov(height_1,height_2)     disorder-induced variance
```

The metadata file records the exact schedule and checksum.

## `glauber_square_grid_20260822/`

This directory retains 960 Gamma and 352 all-one frozen-environment blocks
from the direct weighted-dimer production campaign. It is the compact input
needed to reproduce the environment-blocked scaling analysis without storing
the full MCMC traces in Git. See its [`README.md`](glauber_square_grid_20260822/README.md)
for schemas, seeds, checksums, and the exact command.
