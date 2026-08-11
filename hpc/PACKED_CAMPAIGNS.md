# Packed Hamilton campaign recovery

The packed wrappers accelerate one-sample Slurm arrays without changing the
scientific task identity, seed, output path, or analysis schema. They are
intended for campaigns whose original one-task-per-allocation layout is limited
by a node-count QOS while leaving most CPUs and memory on each node unused.

## Safety model

1. Hold the pending elements of the original array. Do not cancel the array:
   its `afterok` dependants must remain intact.
2. Save the exact held task IDs to an immutable, timestamped manifest.
3. Submit packed jobs against that manifest with the original campaign
   parameters.
4. Check every packed job and per-task log for success.
5. Release the original array. Its resumable runner validates the files written
   by the packed jobs and exits quickly, allowing the original dependency chain
   to complete normally.

Never run a packed task concurrently with the same unheld original task. Both
runners use deterministic seeds and atomic final writes, but avoiding duplicate
writers keeps provenance and failure recovery unambiguous.

## Wrappers

- `square_grid_packed.slurm` launches original square-grid task IDs.
- `aztec_spatial_packed.slurm` launches original Aztec spatial task IDs.

Both wrappers require an absolute `TASK_FILE`, a positive `GROUP_SIZE`, an
explicit per-step `WORKER_MEMORY`, and the same campaign environment variables
used by the original submission. Setting per-step memory is essential: without
it, Slurm may assign the entire allocation memory to the first child step and
serialize workers that were intended to run concurrently. The
one-based packed array index selects consecutive groups of lines from the task
manifest. A single multi-task `srun` launches all ranks concurrently; each rank
has a separate log containing its original task ID. The packed group fails if
any rank fails.

The submission must request at least `GROUP_SIZE` Slurm tasks. Memory must be
sized for simultaneous peak use rather than for a single child.

## Measured Hamilton memory envelope

Measurements from the 8 August 2026 high-order pilots give these approximate
per-process maxima:

| Order | Maximum RSS | Conservative packing on a 250 GB node |
|---:|---:|---:|
| 6,144 | 40 GB | 5 workers, 220 GB allocation |
| 5,120 | 28 GB | 5 workers under the same ultra-size allocation |
| 4,096 | 17.41 GB production maximum | 11 workers, 242 GB allocation |
| 3,072 and 2,560 | lower | 11 workers under the same high-size allocation |

Aztec spatial tasks used less than 1 GB in the completed production run. A
24-worker, 48 GB allocation retains a wide safety margin.

These values are campaign-specific evidence, not general defaults. Re-measure
after changing the sampler, lattice representation, Julia version, or node
type.

The 9 August recovery initially submitted the high-size manifests in groups of
eight.  Before any of those groups started, `repack_hamilton_high_20260809.sh`
superseded them with groups of eleven after production job `18243074` measured
a 17.41 GB maximum RSS.  Its separate job map retains both generations of job
IDs and the derived manifest excluding the two tasks completed by the smoke
test.

## 11 August shared-partition, size-aware recovery

Hamilton's `long` partition exposes only one standard node and is intended for
jobs that need more than three days.  The packed square-grid allocations take
minutes, not days.  On 11 August, the still-pending recovery work was therefore
moved to the 119-node `shared` pool with command-line walltime overrides.  The
running lognormal group `18249837_14` was left untouched; only held pending
elements and their obsolete release jobs were replaced.

The replacement launcher is
`repack_hamilton_size_aware_20260811.sh`.  It derives and hashes disjoint task
manifests by lattice size, verifies that their union exactly matches each
superseded manifest, and rebuilds the original-array release dependencies.  Its
remote job map and checksum are:

```text
/home/fvkl37/hamilton-booster-20260809/manifests/size_aware_shared_job_map_20260811.txt
/home/fvkl37/hamilton-booster-20260809/manifests/size_aware_shared_job_map_20260811.txt.sha256
```

The measured packing envelope on Hamilton's 246 GB standard nodes is:

| Order | Measured or conservative peak RSS | Workers | Per-worker limit | Allocation | Walltime |
|---:|---:|---:|---:|---:|---:|
| 6,144 | 40 GB | 5 | 44 GB | 220 GB | 1 hour |
| 5,120 | 28 GB | 7 | 34 GB | 238 GB | 45 minutes |
| 4,096 | 17.41 GB | 12 | 20 GB | 240 GB | 30 minutes |
| 3,072 | 9.38 GB | 20 | 12 GB | 240 GB | 20 minutes |
| 2,560 | 5.90 GB baseline; under 7 GB scaled disorder estimate | 30 | 8 GB | 240 GB | 20 minutes |

Do not reuse the 2,560 estimate after changing the sampler, lattice
representation, Julia version, or environment construction.  The replacement
keeps at least about 20% per-worker headroom against the campaign measurements
or conservative size scaling.

The held original validation arrays were also updated in place from `long` to
`shared`.  Their high-size elements retain 30 minutes and 48 GB each; their
ultra-size elements retain one hour and 96 GB each.  Existing outputs normally
validate in 5--10 seconds, while these limits still allow an unexpectedly
missing task to be recomputed instead of silently weakening the dependency
chain.  The exact post-update scheduler records and checksum are stored in
`size_aware_shared_original_updates_20260811.txt` and its `.sha256` file beside
the size-aware job map.

The first four seven-worker `L=5120` allocations began concurrently on four
different shared nodes.  The first completed in 18 minutes 13 seconds with a
27,768,828 KB peak RSS for its largest rank, below the 34 GB worker limit; all
seven task logs ended with completion markers and the Slurm allocation exited
zero.  This is the retained real-data validation of the new packing.

## Audit requirements

Retain the following together:

- original array job ID and `sacct` submission line;
- exact task manifest and its SHA-256 checksum;
- packed submission command and returned job ID;
- group and per-task logs;
- final `sacct` state for original, packed, and analysis jobs;
- counts of valid batch, diagnostic, and execution-provenance files.

This makes the acceleration reproducible while leaving the published dataset
identical to a serial completion of the original arrays.
