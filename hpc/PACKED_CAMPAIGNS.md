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
