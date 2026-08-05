# Validation status

## Completed here

- Preserved the existing Aztec source and retained result directories.
- Verified shell syntax with `bash -n` / `sh -n` for every new Slurm and smoke script.
- Verified delimiter and string balance in all new Julia files.
- Mirrored the square-grid graph, Wilson tree, dual complement, perfect matching, and height-cut indexing independently in Python.
- Checked 100 independently sampled trees for every `L=1,...,9` in the mirror implementation.
- For all checked trees, the complement had exactly `4L^2-1` edges, was connected, every primal edge received exactly one Temperley owner, and every tested height increment stayed within its deterministic cut bound.

## First required execution on Hamilton

Julia itself is not installed in the artifact-building container. The first authoritative language-level validation is therefore:

```bash
bash hpc/setup_hamilton_environment.sh
```

This runs the complete Julia package tests under Hamilton's `julia/1.10.4`. The pilot remains blocked until that succeeds and `hpc/square_grid_smoke.slurm` completes without errors.
