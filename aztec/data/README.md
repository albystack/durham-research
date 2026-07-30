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

There are 26,050 data rows and no repeated `(order, sample_id)` keys.
`campaign_metadata.txt` records the model parameters, seed convention, Julia
version and planned sample count.
