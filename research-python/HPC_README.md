# Hamilton HPC Workflow

This project now has a small Slurm-ready batch path for loop-erased random-walk
experiments in fixed site-directional environments.

The Hamilton-specific partition and account are deliberately not hard-coded.
Check them on the cluster before submitting large jobs.

## Copy Or Clone

From your Mac, copy the current project:

```bash
scp -r /local/project/path \
  <username>@hamilton8.dur.ac.uk:~/lerw-project
```

Then log in:

```bash
ssh <username>@hamilton8.dur.ac.uk
cd ~/lerw-project
```

Useful inspection commands:

```bash
find . -maxdepth 2 -type f | sort
python3 -m py_compile main.py simulation.py run_hpc_batch.py generate_hpc_config.py
```

## Python Environment

Create and populate a virtual environment:

```bash
module avail python
module load python/3.12.6
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

Hamilton's default `python3` may be too old. The code needs Python 3.10 or
newer because it uses modern type annotations.

The batch runner itself writes CSV and uses only the standard library plus the
project modules. The plotting stack in `requirements.txt` is for the existing
local analysis code and `analyze_hpc_results.py`.

## Login-Node-Safe Test

Generate the first tiny Hamilton milestone configuration:

```bash
.venv/bin/python generate_hpc_config.py \
  --output configs/hpc_test.csv \
  --sizes 64 128 \
  --distributions baseline gamma:1.0 \
  --batches 1 \
  --num-environments 2 \
  --walks-per-environment 3
```

Run one very small task locally on the login node only as a functionality
check:

```bash
.venv/bin/python run_hpc_batch.py \
  --task-id 0 \
  --config configs/hpc_test.csv \
  --output-dir results_hpc_smoke
```

Do not run large `L`, many environments, or many walks on the login node.

## Submit The First Array

Inspect available partitions/account conventions:

```bash
sinfo
```

Edit `slurm/run_lerw_array.slurm` if Hamilton requires `--partition` or
`--account`, and set the array range to match the generator output. For the
default first milestone, the range is `0-3`.

Create the log directory before submission:

```bash
mkdir -p logs
sbatch slurm/run_lerw_array.slurm
```

Monitor:

```bash
squeue -u "$USER"
sacct -j JOB_ID
```

Inspect failures:

```bash
ls logs
tail -n 80 logs/lerw_JOBID_TASKID.err
tail -n 80 logs/lerw_JOBID_TASKID.out
```

Rerun selected failed tasks by submitting a narrower array:

```bash
sbatch --array=2 slurm/run_lerw_array.slurm
```

If a task wrote a failed or partial CSV and you intentionally want to replace
it:

```bash
.venv/bin/python run_hpc_batch.py \
  --task-id 2 \
  --config configs/hpc_test.csv \
  --output-dir results_hpc \
  --rerun-failed
```

## Result Layout

Batch files are written atomically as CSV:

```text
results_hpc/
  baseline/
    L_0064/
      batch_0000.csv
  gamma_shape_1/
    L_0064/
      batch_0000.csv
```

Each row includes separate `environment_seed` and `walk_seed`, plus
`environment_id` and `walk_id`, so annealed and quenched statistics can be
computed later without mixing them.

## Analyse Completed Results

After all tasks in a config are complete, combine and analyse them:

```bash
.venv/bin/python analyze_hpc_results.py \
  --config configs/hpc_test.csv \
  --results-dir results_hpc \
  --output-dir analysis_hpc_smoke
```

This writes:

```text
analysis_hpc_smoke/
  validation.csv
  combined_raw.csv
  summary.csv
  loglog_fits.csv
  pointwise_ratios.csv
  local_effective_exponents.csv
```

Submit the same analysis as a Slurm job:

```bash
sbatch slurm/run_analysis.slurm
```

For a large-L tail fit, set `FIT_MIN_L`:

```bash
FIT_MIN_L=256 OUTPUT_DIR=analysis_hpc_tail sbatch slurm/run_analysis.slurm
```

If plotting libraries are unavailable, the analysis still writes CSV files and
records that plots were skipped.

## Generate A Larger Config

The full template covers the current distribution list and:

```text
L = 32, 64, 128, 256, 512, 1024
```

Start with modest batch counts and scale after checking runtime at `L=512` and
`L=1024`:

```bash
.venv/bin/python generate_hpc_config.py \
  --preset hpc_full \
  --output configs/hpc_full.csv \
  --batches 4 \
  --num-environments 5 \
  --walks-per-environment 5
```

That command creates `15 * 6 * 4 = 360` array tasks. Edit
`slurm/run_lerw_array.slurm` so `#SBATCH --array=0-359`, or override at submit
time:

```bash
CONFIG=configs/hpc_full.csv OUTPUT_DIR=results_hpc_full \
  sbatch --array=0-359 slurm/run_lerw_array.slurm
```

Do not run this full config on the login node.

## Copy Results Back

From your Mac:

```bash
scp -r \
  <username>@hamilton8.dur.ac.uk:~/lerw-project/results_hpc \
  ./results_hpc
```
