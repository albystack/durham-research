# Square-grid implementation installation and Hamilton workflow

## Apply to the VS Code repository

The patch archive is rooted at the existing repository root. Before applying it:

```bash
cd /Users/alberto/Desktop/research
cp -a . ../research_before_square_grid
```

Extract the patch into `/Users/alberto/Desktop/research`, allowing matching files to be replaced. Then inspect:

```bash
cd /Users/alberto/Desktop/research
git status --short
git diff --stat
code .
```

The patch does not alter retained CSV results or anything under `aztec/output/`.

## Files added

```text
aztec/src/SquareGrid.jl
aztec/test/square_grid_runtests.jl
aztec/scripts/run_square_grid_campaign.jl
aztec/scripts/merge_square_grid_batches.jl
aztec/configs/square_grid_smoke.csv
aztec/configs/square_grid_pilot.csv
aztec/docs/SQUARE_GRID_MODEL.md
hpc/setup_hamilton_environment.sh
hpc/square_grid_smoke.slurm
hpc/square_grid_pilot.slurm
hpc/analyze_square_grid.slurm
hpc/README.md
```

Existing files updated:

```text
.gitignore
aztec/src/AztecDiamond.jl
aztec/test/runtests.jl
aztec/test/smoke_workflows.sh
aztec/README.md
aztec/configs/README.md
```

## Upload with WinSCP through AppsAnywhere

1. Launch WinSCP in Durham AppsAnywhere.
2. Use SFTP, host `hamilton8.dur.ac.uk`, port `22`, username `fvkl37`.
3. Complete password and MFA.
4. In the right pane create `/home/fvkl37/durham-research-square-grid`.
5. In the left pane open `/Users/alberto/Desktop/research` if the local drive is exposed. If AppsAnywhere exposes a redirected drive, locate the repository there.
6. Upload the repository source, excluding `.git`, `aztec/output`, large raw data, `.DS_Store`, and `__MACOSX`.

The simplest upload is the compact repository ZIP. Upload it to `/home/fvkl37`, then extract it in PuTTY.

## Log in with PuTTY

1. Launch PuTTY in AppsAnywhere.
2. Host: `hamilton8.dur.ac.uk`; port `22`; connection type SSH.
3. Log in as `fvkl37`, then complete password and MFA.

## Extract an uploaded ZIP

```bash
cd /home/fvkl37
mkdir -p durham-research-square-grid
unzip -q research_square_grid_ready.zip -d durham-research-square-grid
cd durham-research-square-grid/research
ls
```

The supplied compact ZIP contains exactly one top-level `research/` folder.

## Prepare Julia 1.10.4

```bash
cd /home/fvkl37/durham-research-square-grid
bash hpc/setup_hamilton_environment.sh
```

This creates `.hamilton_env`, develops the local package into it, and runs the complete tests without modifying the Julia 1.12-generated `aztec/Manifest.toml`.

## Submit the first smoke job

```bash
mkdir -p logs
sbatch hpc/square_grid_smoke.slurm
squeue -u "$USER"
```

After it finishes:

```bash
cat logs/sqgrid-smoke-*.out
cat logs/sqgrid-smoke-*.err
```

Do not submit the pilot until this job exits successfully.
