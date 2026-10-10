# Completed paired-height data

**258,048 environments**, two conditionally independent tilings per environment,
two probes per pair: **516,096 rows**. The compressed files are included in Git.

| File | Contents |
|---|---|
| `measurements.csv.gz` | Every paired spatial measurement, identifiers and seeds |
| `plan.json` | Frozen allocation, 504 batches, laws, fit windows and resampling seeds |
| `intervals.csv` | Recorded moments and environment-bootstrap intervals |
| `analysis.json` | Recorded 18 fits, nine additional contrasts and joint control test |
| `provenance.json.gz` | Original batch receipts, raw-verification certificates, source hashes and scheduler accounting |
| `source-snapshot.tar.gz` | Exact 164-file source archive used on Hamilton |
| `verification.json` | Independent local checks of the completed compact results |
| `manifest.json` | SHA256 checksums of the public scientific files |

## Allocation

Sizes: **16, 24, 32, 48, 64, 96, 128**, in each geometry.

- `aztec_gamma`: 16,384 pairs per size/domain; two independent Gamma families
  with shapes 0.2 and 0.25, scale one, remaining edge weights one.
- `aztec_mean`: 1,024 pairs per size/domain; deterministic family means.
- `uniform`: 1,024 pairs per size/domain; all edge weights one.

The Gamma family means are **0.2 and 0.25**, not one. The total is 229,376
Gamma pairs plus 28,672 control pairs. Pilot and earlier campaigns are excluded.

## Measurement columns

`domain`, `model` and `size` identify the cell. `environment_id` has a batch
prefix and uniquely identifies a pair across the whole campaign. `probe`
selects the fraction 1/8 or 1/4. `r` is the coordinate separation and `distance`
is the Euclidean separation `sqrt(2) × r`, in primal-edge units. `axis` is
`diagonal`; `height_scale` is one. `value1` and `value2` are integer height
increments. `law_id` records the weight convention; `batch` records the source
batch number. Seed columns are decimal strings, preserving full integer values.

Square rows have `environment_seed`, `replica1_seed` and `replica2_seed`.
Aztec rows have `environment_seed` and `coin_seed`, which identifies the random
stream drawing the two independent choice arrays. Inapplicable seed fields
are empty. Never pool the two replicas as independent disorder environments.

## Provenance and raw arrays

Hamilton workers verified every saved matching and recomputed its observables
before issuing completion receipts. The compact archive SHA256 is
`90e30c0f97e0e22132110dfb35a0a6e5e7007d5bbe4c8270be3a10c9977c0254`.
The frozen source archive SHA256 is
`4c533d765419550fe3fcee005ee8ff58ed670f3ac3d08f7c2a387412a6fad475`.

Full weights and tiling arrays remain under
`/nobackup/fvkl37/larger_matched_20261007/production/batches/` on Hamilton;
they are not included in Git. These measurement tables reproduce the complete
statistical analysis. The raw-verification certificates document the separate
matching checks. No raw arrays were deleted by this repository cleanup.
