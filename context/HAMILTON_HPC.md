# Hamilton8 + NCC HPC Guide for the Research Repository

> **Last verified:** 16 August 2026  
> **Primary target for this project:** Durham University **Hamilton8**  
> **Secondary/alternative system:** Durham Computer Science **NCC**  
> **Purpose:** one operational source of truth for running the LERW / spanning-tree / dimer numerical campaigns reproducibly and efficiently.
>
> This file is deliberately operational. Scientific definitions and research priorities live in
> `RESEARCH_BRIEF.md`, `RESULTS.md`, and `NEXT_STEPS.md`.

---

# 1. Which cluster should this project use?

## Hamilton8 — default

Use Hamilton8 for the current Julia Monte Carlo / combinatorial simulation work.

Why:
- large CPU capacity;
- 128 physical CPU cores per standard node;
- large shared filesystem and per-user `/nobackup` space;
- node-local scratch via `$TMPDIR`;
- Slurm job arrays;
- `shared`, `multi`, `long`, `bigmem`, `test`, and `cuda` partitions;
- well suited to many independent CPU simulations.

For the current research, **prefer CPU job arrays on Hamilton8** unless an implementation is genuinely GPU-enabled.

## NCC — only when there is a reason

NCC is a separate Department of Computer Science system. It is primarily a GPU research cluster, although it also has a CPU partition.

Use NCC only when:
- the code has a real GPU implementation;
- a CS-specific workflow requires NCC;
- Hamilton is unsuitable and the NCC allocation/partition matches the job.

Do **not** copy Hamilton partition names into NCC scripts or vice versa. They are different Slurm installations with different partitions, limits, storage assumptions, and hardware.

---

# 2. Hamilton8 system facts

Hamilton8 currently runs Rocky Linux 8.

Overall service:
- 15,616 CPU cores;
- 36 TB RAM;
- about 1.9 PB disk;
- CPU and GPU nodes.

Standard CPU nodes:
- 120 total standard compute nodes;
- 128 CPU cores per node;
- 2 × AMD EPYC 7702;
- 256 GB physical RAM, about 246 GB available to users;
- 400 GB local SSD;
- 2 CPU sockets;
- each socket has 4 NUMA domains;
- each NUMA domain has 4 chiplets;
- each chiplet contains 4 cores and shares 16 MB L3 cache;
- HDR InfiniBand interconnect.

High-memory nodes:
- 2 nodes;
- 128 CPU cores each;
- about 2 TB RAM each;
- 400 GB local SSD.

GPU node:
- 128 CPU cores;
- about 2.2 TB RAM;
- 3 TB local NVMe;
- 8 × NVIDIA H200 NVL, 144 GB each;
- Slurm exposes fractional H200 resources as well as whole GPUs.

**Research consequence:** the present simulation campaign is naturally CPU-heavy and embarrassingly parallel. Do not request GPU resources merely because they exist.

---

# 3. Login and Access VPN

## Personal connection workflow used for this project

When working off campus, **connect to Durham University Access VPN first, then SSH to Hamilton**.

This is the normal working sequence for this repository:

```text
Mac / off-campus network
        ↓
Durham Access VPN connected
        ↓
Durham network access available
        ↓
ssh <CIS_USERNAME>@hamilton8.dur.ac.uk
        ↓
Hamilton8 login node
        ↓
Slurm (`srun` / `sbatch`)
        ↓
compute node
```

The Durham Access VPN secure-logon service is available through the University's Access service and is protected by the University's authentication/MFA policy.

Practical workflow:

1. Connect the Mac to **Durham Access VPN**.
2. Confirm the VPN connection is active.
3. Open Terminal.
4. SSH to Hamilton:

```bash
ssh <CIS_USERNAME>@hamilton8.dur.ac.uk
```

5. Authenticate with the normal Durham credentials/MFA when requested.
6. Once logged in, use the login node only for light work and submit computational work through Slurm.

### If SSH does not work from off campus

Check in this order:

```text
1. Is Access VPN connected?
2. Does normal internet access work?
3. Is the hostname exactly hamilton8.dur.ac.uk?
4. Is the CIS username correct?
5. Has the Hamilton account/access been enabled?
6. Is MFA completing successfully?
7. Is Hamilton in a maintenance/at-risk period?
```

Useful local diagnostic commands:

```bash
# Verify DNS resolution.
nslookup hamilton8.dur.ac.uk

# Verbose SSH diagnostics.
ssh -v <CIS_USERNAME>@hamilton8.dur.ac.uk
```

Do not put passwords, MFA codes, VPN credentials, private SSH keys, or access tokens into this repository or into Codex prompts.

## Hamilton address

```bash
ssh <CIS_USERNAME>@hamilton8.dur.ac.uk
```

From the University network, `hamilton8` may also resolve directly.

Authentication:
- normal Durham University username/password;
- MFA is normally required off campus unless an approved University network/VPN route applies;
- for this project's normal off-campus workflow, **Access VPN is connected before SSH**.

Exit:

```bash
exit
```

## Login-node policy

Hamilton has shared login nodes.

Allowed on login nodes:
- editing;
- Git;
- compiling/installing light software;
- package/environment setup;
- file movement;
- short/light tests.

Do **not** run production simulations, large analyses, or meaningful benchmarks directly on login nodes. Request compute resources through Slurm instead.

Durham lists the login nodes as being **at risk Tuesdays 08:00–10:00** for routine maintenance. Login sessions may be interrupted, but compute-node batch jobs continue normally.

---

# 4. Recommended research directory layout on Hamilton

Keep source/configuration separate from large generated data.

```text
$HOME/
└── research/
    ├── repo/                  # git checkout, code, configs, scripts
    ├── env/                   # optional small environment metadata
    └── manifests/             # run manifests/checksums worth backing up

$NOBACKUP/
└── research/
    ├── raw/                   # large raw simulation outputs
    ├── derived/               # derived tables / intermediate data
    ├── checkpoints/           # restart files
    ├── logs/                  # campaign logs if large
    └── julia_depot/           # Julia package cache/artifacts

$TMPDIR/                       # job-local scratch; exists only for the job
```

Recommended variables:

```bash
export RESEARCH_HOME="$HOME/research"
export RESEARCH_DATA="$NOBACKUP/research"
export JULIA_DEPOT_PATH="$NOBACKUP/research/julia_depot"
```

Create once:

```bash
mkdir -p "$HOME/research/repo"
mkdir -p "$HOME/research/manifests"
mkdir -p "$NOBACKUP/research"/{raw,derived,checkpoints,logs,julia_depot}
```

---

# 5. Hamilton storage

## `$HOME`

Location:

```text
/home/<username>
```

Properties:
- default quota: 10 GB;
- backed up;
- backups normally retained for 30 days.

Use for:
- source code;
- small configs;
- job scripts;
- manifests;
- irreplaceable small metadata.

Check:

```bash
echo "$HOME"
quota
du -sh "$HOME"/* 2>/dev/null | sort -h
```

## `$NOBACKUP`

Location:

```text
/nobackup/<username>
```

Properties:
- initial quota: 600 GB;
- accessible from login and compute nodes;
- **not backed up**.

Use for:
- raw Monte Carlo output;
- large derived files;
- Julia package depot;
- restart/checkpoint files that are reproducible or separately mirrored.

Check:

```bash
echo "$NOBACKUP"
quota
du -sh "$NOBACKUP"/* 2>/dev/null | sort -h
```

Important scientific rule:

> A dataset existing only in `$NOBACKUP` is not safely archived.

Copy irreplaceable data/results elsewhere.

## `$TMPDIR`

Each standard compute node has local SSD scratch. Request scratch from Slurm:

```bash
#SBATCH --gres=tmp:50G
```

Within the job:

```bash
echo "$TMPDIR"
df -h "$TMPDIR"
```

Properties:
- node-local;
- fast for temporary I/O;
- removed automatically when the job ends.

Use for:
- temporary simulation chunks;
- repeated intermediate reads/writes;
- staging data to avoid hammering shared storage.

**Anything required after the job must be copied out before exit.**

## Project directories

Shared project directories are under:

```text
/projects/<project-name>
```

Durham's documented default project quota is 20 GB and project directories are backed up. They are group-owned University project storage, not private personal storage.

Useful commands:

```bash
df -h /projects/<project-name>
getent group <project-group>
ls -l /projects/<project-name>
chgrp <project-group> <file>
chmod g+w <file>
```

---

# 6. Transfer files to/from Hamilton

Run these commands on the **local Mac**, not inside Hamilton, unless intentionally transferring from Hamilton to another host.

## `rsync` — preferred for repository/data synchronization

Local → Hamilton:

```bash
rsync -av --progress ~/Desktop/research/ \
  <CIS_USERNAME>@hamilton8.dur.ac.uk:~/research/repo/
```

Code only, excluding bulky local files:

```bash
rsync -av --progress \
  --exclude '.git/' \
  --exclude 'raw/' \
  --exclude 'results/' \
  ~/Desktop/research/ \
  <CIS_USERNAME>@hamilton8.dur.ac.uk:~/research/repo/
```

Hamilton → local results:

```bash
rsync -av --progress \
  <CIS_USERNAME>@hamilton8.dur.ac.uk:/nobackup/<CIS_USERNAME>/research/derived/ \
  ~/Desktop/research-results/
```

## `scp`

Single file local → Hamilton:

```bash
scp my_file <CIS_USERNAME>@hamilton8.dur.ac.uk:
```

Directory local → `$NOBACKUP`:

```bash
scp -r my_directory \
  <CIS_USERNAME>@hamilton8.dur.ac.uk:/nobackup/<CIS_USERNAME>/
```

Hamilton → local:

```bash
scp \
  <CIS_USERNAME>@hamilton8.dur.ac.uk:/nobackup/<CIS_USERNAME>/file.csv \
  .
```

For repeated synchronization, prefer `rsync`.

---

# 7. Git workflow on Hamilton

If the research code has a remote:

```bash
cd "$HOME/research"
git clone <REPOSITORY_URL> repo
cd repo
git status
git log --oneline -5
```

Before every campaign, record:

```bash
git rev-parse HEAD
git status --short
```

For reproducibility, production runs should ideally use a clean commit.

Recommended:

```bash
test -z "$(git status --porcelain)" || echo "WARNING: working tree is dirty"
git rev-parse HEAD > "$NOBACKUP/research/logs/git_commit.txt"
```

Do not make untracked edits directly on compute nodes during a production campaign.

---

# 8. Discovering software on Hamilton

Do not hard-code a Julia module version into permanent documentation because installed module versions can change.

Discover current availability:

```bash
module avail
module avail julia
module spider julia
```

After loading the version you intend to use:

```bash
module load <JULIA_MODULE_SHOWN_BY_MODULE_SYSTEM>
which julia
julia --version
```

Inspect loaded modules:

```bash
module list
```

If Julia is not available as a suitable module, use an approved user-local Julia installation under `$HOME`/`$NOBACKUP`; do not assume a version.

---

# 9. Julia project setup

From the repository root:

```bash
cd "$HOME/research/repo"
```

Set the package depot somewhere larger than the 10 GB home directory:

```bash
export JULIA_DEPOT_PATH="$NOBACKUP/research/julia_depot"
mkdir -p "$JULIA_DEPOT_PATH"
```

Instantiate/precompile the project:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Verify:

```bash
julia --project=. -e 'using InteractiveUtils; versioninfo()'
```

Run tests:

```bash
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

For package-style repositories:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Prefer committed:

```text
Project.toml
Manifest.toml
```

for reproducible campaigns.

---

# 10. Slurm essentials on Hamilton

Hamilton production compute must be requested through Slurm.

Submit:

```bash
sbatch job.sh
```

Queue:

```bash
squeue --me
squeue -u "$USER"
```

System/availability:

```bash
sinfo
sfree
```

Inspect a job:

```bash
scontrol show job <JOBID>
```

Cancel:

```bash
scancel <JOBID>
```

Finished-job accounting:

```bash
sacct -j <JOBID>
```

More useful accounting fields:

```bash
sacct -j <JOBID> \
  --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

For a running job:

```bash
sstat -j <JOBID>.batch
```

Follow the default Slurm output:

```bash
tail -f slurm-<JOBID>.out
```

Durham's current Hamilton documentation states that once a batch job starts, its printed output is written to:

```text
slurm-<jobid>.out
```

unless the job script overrides the output path with `#SBATCH --output=...`.

If Slurm kills a job because it exceeded requested time or memory, inspect the end of the Slurm output and the accounting record before resubmitting:

```bash
tail -100 slurm-<JOBID>.out
sacct -j <JOBID>   --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

## Interactive compute work

For work that needs to be interactive but is too intensive for a login node, request an interactive compute allocation with `srun`.

Minimal interactive shell:

```bash
srun --pty bash
```

Small test allocation:

```bash
srun --pty --mem=2G -c 2 -p test bash
```

Run a command directly through Slurm:

```bash
srun --mem=2G -c 2 -p test <command>
```

Interactive `srun` jobs are subject to the same partition/resource controls as batch jobs.

## Reading `squeue`

Useful:

```bash
squeue --me
squeue -u "$USER"
```

Common states:
- `R` = running;
- `PD` = pending.

Common pending reasons:
- `Resources` = waiting for suitable resources to become free;
- `Priority` = waiting behind higher-priority jobs;
- `PartitionNodeLimit` = the request is incompatible with the partition and will not run as submitted.

If a job is pending, inspect:

```bash
scontrol show job <JOBID>
```

Do not repeatedly cancel/resubmit a normally pending `Resources` or `Priority` job without a reason.

## Job notifications

For batch jobs, Durham documents:

```bash
#SBATCH --mail-type=END,FAIL
```

as the useful notification pattern.

A separate `--mail-user` is not normally required to send notifications to the submitting user's Durham address. Specify it only when a different supported destination is intentionally required.

## Hamilton portal

The Hamilton portal can also be used to:
- monitor queue/job status;
- inspect graphical resource-use information for running jobs;
- help tune future CPU/RAM/GPU requests.

Use Slurm accounting (`sacct`) as the durable command-line record after jobs finish.

---

# 11. Hamilton partitions and limits

Current documented Hamilton CPU/GPU partitions:

| Partition | Intended use | Node type | Current job limit |
|---|---|---|---|
| `shared` | jobs that can share nodes | standard | 3 days |
| `multi` | one or more whole nodes | standard | 3 days |
| `long` | jobs needing >3 days | standard | 7 days |
| `bigmem` | >250 GB RAM | high-memory | 3 days |
| `test` | short tests | standard | 15 minutes |
| `cuda` | GPU jobs | GPU | 3 days |

The `shared`, `multi`, and `long` partitions draw from the same standard-node pool; current Durham documentation lists 119 nodes in that schedulable pool.

Default CPU-job request if nothing is specified:
- 1 hour;
- 1 CPU core;
- 1 GB RAM;
- 1 GB `$TMPDIR`.

Useful directives:

```bash
#SBATCH -p shared
#SBATCH -t 00-04:00:00
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH --gres=tmp:20G
#SBATCH --array=0-31
#SBATCH --mail-type=END,FAIL
```

Meanings:
- `-p`: partition;
- `-t`: maximum walltime (`dd-hh:mm:ss`);
- `-c`: CPUs for one multithreaded task;
- `-n`: number of tasks, normally MPI/distributed jobs;
- `-N`: number of nodes;
- `--mem`: RAM per node;
- `--gres=tmp`: local temporary disk;
- `--array`: job array.

**Use `-c`, not `-n`, for one multithreaded Julia process.**

---

# 12. Choosing the right Hamilton partition for this research

## `test`

Use for:
- syntax checks;
- dependency checks;
- tiny simulations;
- Slurm environment checks;
- verifying paths and outputs.

Do not use it for actual statistical pilots that exceed 15 minutes.

## `shared` — default production choice

Use when one job does not require an entire node.

For this project, this is usually the best choice for:
- array shards;
- moderate Julia thread counts;
- analysis jobs;
- bootstrap jobs;
- parameter sweeps.

Example resource starting point:

```bash
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=16G
```

Then measure actual CPU/RAM and adjust.

## `multi`

Use only when the code genuinely needs:
- a whole 128-core node; or
- multiple whole nodes.

A 128-core request is not automatically faster. Benchmark scaling first.

For independent Monte Carlo workloads, **multiple smaller array jobs often outperform one huge process operationally** because:
- easier scheduling;
- lower failure cost;
- simpler restart;
- cleaner RNG partitioning;
- independent outputs;
- easier load balancing.

## `long`

Use only when one indivisible job genuinely needs more than 3 days.

Maximum documented walltime: 7 days.

Prefer checkpointed/split campaigns over one giant 7-day process whenever scientifically equivalent.

## `bigmem`

Use only when memory need exceeds about 250 GB.

Current simulation designs should not enter `bigmem` merely to avoid fixing avoidable memory use.

## `cuda`

Use only after an actual GPU implementation exists and has been benchmarked.

---

# 13. Hamilton Julia threading: important

Hamilton sets several native-library thread variables from the value requested with `#SBATCH -c`, including OpenMP/OpenBLAS/MKL/BLIS variables.

Julia threads are separate. Explicitly set:

```bash
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
```

Check in Julia:

```bash
julia -e 'using Base.Threads; @show nthreads()'
```

## Avoid nested BLAS oversubscription

If the simulation parallelizes at the Julia-thread/environment level and does not need multithreaded BLAS, override native BLAS threading:

```bash
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export OMP_NUM_THREADS=1
```

Then:

```bash
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
```

This avoids a 32-thread Julia process accidentally spawning 32 BLAS threads per Julia worker path.

If a specific analysis is BLAS-dominated, benchmark the opposite strategy separately rather than mixing both kinds of parallelism blindly.

---

# 14. CPU binding / NUMA

Hamilton standard nodes are NUMA machines.

Start simple:

```bash
srun --cpu-bind=cores julia ...
```

Do not invent complicated NUMA binding before benchmarking.

For each major kernel, benchmark:
- 1 thread;
- 2;
- 4;
- 8;
- 16;
- 32;
- 64;
- 128 only if justified.

Record:
- environments/second;
- walltime;
- CPU efficiency;
- memory;
- output I/O.

Choose the smallest core count near the performance knee rather than the maximum core count.

---

# 15. Research RNG rules on HPC

This is non-negotiable for the Monte Carlo campaign.

## Every independent environment must have a reproducible identity

A run should be recoverable from something like:

```text
campaign_id
model_id
L
disorder_id
array_task_id
environment_index
master_seed
code_commit
```

## Do not use scheduler timing as randomness

Bad:

```julia
seed = time_ns()
```

for production science.

Better concept:

```text
seed = deterministic_function(
    master_seed,
    model,
    L,
    disorder,
    array_task_id,
    local_environment_index
)
```

## Do not let thread scheduling define the random stream

If changing `JULIA_NUM_THREADS` changes which random numbers correspond to an environment, reproduction becomes difficult.

Prefer environment-specific RNG streams/seeds.

## Paired same-environment experiments

For \(H_1,H_2\) or \(\Delta H_1,\Delta H_2\):
- one deterministic environment seed;
- two separate deterministic conditional-sample seeds;
- record all three;
- never accidentally regenerate the environment between the paired samples.

For example conceptually:

```text
environment_seed = f(master, env_id)
sample1_seed     = g(master, env_id, 1)
sample2_seed     = g(master, env_id, 2)
```

---

# 16. Preferred parallelization pattern for the current research

The scientific samples are largely independent across environments.

Preferred hierarchy:

```text
Slurm job array
    ├── task 0
    │   └── Julia process with N threads
    ├── task 1
    │   └── Julia process with N threads
    └── ...
```

Each array task gets:
- non-overlapping environment IDs;
- deterministic seeds;
- its own output shard;
- its own checkpoint;
- no writes to the same mutable output file as another task.

Later:

```text
raw shards
   ↓
validation
   ↓
merge
   ↓
analysis/bootstrap
   ↓
figures/tables
```

Do not have 50 jobs append to one CSV concurrently.

---

# 17. File naming for campaigns

Recommended immutable shard naming:

```text
raw/<campaign>/
  model=<model>/
    law=<law>/
      L=<L>/
        part-00000.arrow
        part-00001.arrow
        ...
```

or:

```text
part_<SLURM_ARRAY_TASK_ID>.csv
```

Each shard should contain metadata or have a sidecar manifest containing:
- job ID;
- array task ID;
- hostname;
- start/end UTC timestamps;
- Git commit;
- Julia version;
- module list;
- master seed;
- environment-ID range;
- model parameters;
- row count;
- checksum.

---

# 18. Minimal Hamilton test job — Julia

Save as `slurm/test_julia.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=research_test
#SBATCH --partition=test
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --gres=tmp:2G
#SBATCH --output=logs/%x_%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

# Load the Julia module discovered with `module avail julia`.
module load <JULIA_MODULE>

export JULIA_DEPOT_PATH="$NOBACKUP/research/julia_depot"
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"

# Prevent nested native-library oversubscription.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

echo "job=$SLURM_JOB_ID"
echo "host=$(hostname)"
echo "cpus=$SLURM_CPUS_PER_TASK"
echo "tmp=$TMPDIR"
julia --version

julia --project=. --startup-file=no -e '
    using Base.Threads
    @show nthreads()
    println("Julia test OK")
'
```

Create logs directory before submission:

```bash
mkdir -p logs
sbatch slurm/test_julia.sh
```

---

# 19. Recommended Hamilton production array template

Save as `slurm/sim_array.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=dimer_sim
#SBATCH --partition=shared
#SBATCH --time=00-12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --gres=tmp:20G
#SBATCH --array=0-31%8
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

module load <JULIA_MODULE>

export JULIA_DEPOT_PATH="$NOBACKUP/research/julia_depot"
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

CAMPAIGN="squaregrid_pilot_v1"
OUTROOT="$NOBACKUP/research/raw/$CAMPAIGN"
CHECKROOT="$NOBACKUP/research/checkpoints/$CAMPAIGN"

mkdir -p "$OUTROOT" "$CHECKROOT" logs

TASK="${SLURM_ARRAY_TASK_ID}"
OUTFILE="$OUTROOT/part_$(printf "%05d" "$TASK").csv"
CHECKPOINT="$CHECKROOT/part_$(printf "%05d" "$TASK").jld2"

echo "campaign=$CAMPAIGN"
echo "job=$SLURM_JOB_ID"
echo "array_job=$SLURM_ARRAY_JOB_ID"
echo "task=$TASK"
echo "host=$(hostname)"
echo "commit=$(git rev-parse HEAD)"
echo "threads=$JULIA_NUM_THREADS"
date -Is

srun --cpu-bind=cores \
  julia --project=. --startup-file=no \
  scripts/run_experiment.jl \
  --campaign "$CAMPAIGN" \
  --task-id "$TASK" \
  --output "$OUTFILE" \
  --checkpoint "$CHECKPOINT"

date -Is
```

Notes:
- `0-31%8` means 32 shards with at most 8 active simultaneously;
- choose shard count and concurrency from the actual experiment, not this example;
- every shard must cover a disjoint environment-ID range;
- make the script restart-safe before scaling.

Submit:

```bash
sbatch slurm/sim_array.sh
```

---

# 20. Staging heavy I/O through `$TMPDIR`

If each shard writes frequently, use node-local scratch.

Pattern:

```bash
WORK="$TMPDIR/$SLURM_JOB_ID"
mkdir -p "$WORK"

LOCAL_OUT="$WORK/result.csv"
FINAL_OUT="$NOBACKUP/research/raw/$CAMPAIGN/part_${SLURM_ARRAY_TASK_ID}.csv"

julia --project=. scripts/run_experiment.jl --output "$LOCAL_OUT"

# Only publish a completed shard.
tmp_final="${FINAL_OUT}.tmp.$$"
cp "$LOCAL_OUT" "$tmp_final"
mv "$tmp_final" "$FINAL_OUT"
```

The final `mv` avoids presenting a partially copied file under the final filename.

For important output, compute checksum:

```bash
sha256sum "$FINAL_OUT" > "${FINAL_OUT}.sha256"
```

---

# 21. Restart-safe production pattern

Every long simulation should support:
- chunked environment ranges;
- atomic checkpoint writes;
- skip-if-complete logic;
- deterministic replay.

Pseudo-flow:

```text
if final shard exists and validates:
    exit success

load checkpoint if present

for next chunk:
    simulate deterministic environment IDs
    append/update local state
    atomically save checkpoint

validate completed shard
atomically publish final output
```

Use temporary filename + rename for checkpoint updates:

```text
checkpoint.tmp -> checkpoint.jld2
```

Do not overwrite the only valid checkpoint in-place.

---

# 22. Slurm job dependencies

Useful for separating simulation from analysis.

Submit simulation:

```bash
SIM_JOB=$(sbatch --parsable slurm/sim_array.sh)
```

Run analysis only if simulation succeeds:

```bash
sbatch --dependency=afterok:$SIM_JOB slurm/analyse.sh
```

Run a cleanup/summary after termination regardless of success if desired:

```bash
sbatch --dependency=afterany:$SIM_JOB slurm/summarise.sh
```

For research reproducibility, keep simulation and statistical analysis as separate jobs when practical.

---

# 23. Pilot → scale workflow

Never jump directly from local code to a large Hamilton campaign.

## Stage A — local

Run tiny deterministic cases.

```bash
julia --project=. test/runtests.jl
```

## Stage B — Hamilton `test`

Validate:
- module;
- paths;
- package environment;
- RNG;
- Slurm variables;
- output;
- checkpoint;
- `$TMPDIR`.

## Stage C — small `shared` pilot

Example:
- 2–4 sizes;
- one/no-disorder control;
- one disorder law;
- enough samples to estimate runtime and memory;
- a few array shards.

Measure:

```bash
sacct -j <JOBID> \
  --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

## Stage D — scaling benchmark

Measure environments/sec for several `-c` values.

Do not assume 128 threads is optimal.

## Stage E — full production

Freeze:
- code commit;
- Project/Manifest;
- config;
- seed scheme;
- shard layout;
- requested resources;
- analysis version.

Then submit.

---

# 24. Resource-sizing rule

Estimate resources from pilots.

## Time

If a pilot produces \(n\) environments in time \(t\):

```text
throughput = n / t
```

For a target shard size \(N\):

```text
estimated_time = N / throughput
```

Request walltime with a reasonable margin, not an arbitrary 3-day maximum.

## Memory

Measure `MaxRSS`.

Then request somewhat above observed peak, accounting for larger \(L\) if memory scales with system size.

Do not request 200 GB when the pilot used 8 GB.

## CPU

Use throughput per core as well as total throughput.

If:
- 16 threads gives 15× speedup;
- 32 gives 20×;
- 64 gives 22×;

then 16–32 cores may be the better scheduler/efficiency choice.

---

# 25. Monitoring and debugging

## Pending

```bash
squeue --me
scontrol show job <JOBID>
```

Common reasons:
- `Resources`: waiting for suitable resources;
- `Priority`: waiting behind higher-priority work;
- impossible partition/resource combinations must be corrected rather than waited out.

## Running

```bash
sstat -j <JOBID>.batch
tail -f logs/<file>.out
```

## Finished

```bash
sacct -j <JOBID>
```

If Slurm kills the job for time or memory, the output/accounting should be inspected before resubmitting.

Do not simply double all resources without diagnosing why it failed.

---

# 26. Useful shell commands on Hamilton

Navigation:

```bash
pwd
ls -lah
cd
tree -L 2
```

Disk:

```bash
quota
df -h
du -sh .
du -sh * | sort -h
```

Processes:

```bash
ps -u "$USER"
top
```

Files:

```bash
find . -type f | head
grep -R "pattern" .
sha256sum file
```

Archives:

```bash
tar -czf results.tar.gz results/
tar -xzf results.tar.gz
```

Git:

```bash
git status
git log --oneline -10
git rev-parse HEAD
git diff
```

Slurm:

```bash
squeue --me
sinfo
sfree
sbatch job.sh
scancel <JOBID>
sacct -j <JOBID>
scontrol show job <JOBID>
```

---

# 27. Hamilton GPU jobs

For the current project, this section is reference only.

Hamilton GPU jobs:
- partition: `cuda`;
- current GPU: H200 NVL;
- GPU requests determine associated CPU/RAM/scratch allocation;
- Durham documentation says not to separately request CPU, memory, or tmp disk for GPU jobs because these are allocated in proportion to the GPU resource.

Fractional GPU:

```bash
#SBATCH -p cuda
#SBATCH --gres=gpu:1
```

or explicit fractional type:

```bash
#SBATCH --gres=gpu:h200_nvl_1g.18gb:1
```

Whole H200:

```bash
#SBATCH --gres=gpu:h200_nvl
```

Inside:

```bash
nvidia-smi
```

Do not port the Monte Carlo code to GPU solely to use this hardware. Benchmark algorithmic suitability first.

---

# 28. NCC: separate cluster reference

NCC is the Durham Computer Science NVIDIA CUDA Centre.

Current official NCC documentation describes:
- two head nodes;
- twelve GPU compute servers;
- six CPU blades;
- 82 physical GPUs overall;
- 246 CPU threads on CPU blades;
- Ubuntu-based environment;
- Slurm for scheduling.

Connection:

```bash
ssh <CIS_USERNAME>@ncc1.clients.dur.ac.uk
```

NCC documentation is explicit:

> all computational activity must go through Slurm.

Do not run simulations directly on NCC head nodes.

## NCC storage

Current NCC docs:
- newer accounts: typically `/home3`, default 250 GB;
- older `/home2` accounts: default 100 GB;
- filesystem is not described as a durable research archive;
- NCC recommends copying important data out;
- long jobs should checkpoint regularly.

Check:

```bash
quota
df -h "$HOME"
```

## NCC current partitions

Current NCC documentation lists access-class-specific GPU partitions plus CPU:

```text
ug-gpu-small
tpg-gpu-small
res-gpu-small
res-gpu-large
gpu-bigmem
cpu
```

Access to some partitions depends on undergraduate / taught postgraduate / research/staff status.

The NCC CPU partition currently limits a job to **32 CPU cores and 60 GB RAM per node**.

For this project's CPU-heavy simulations, this is another reason Hamilton8 is normally preferable.

## NCC QOS

Current documented QOS include:

```text
debug
short
long-high-prio
long-low-prio
long-cpu
```

Maximums currently documented:
- `debug`: 2 hours;
- `short`: 2 days;
- `long-high-prio`: 7 days;
- `long-low-prio`: 7 days and may be preempted;
- `long-cpu`: 14 days, CPU partition only.

NCC limits may change. Check current docs/Slurm before submission.

## NCC preemption

For `long-low-prio`, NCC can preempt/requeue work. Jobs using it must be checkpoint/restart capable.

Official NCC guidance uses:

```bash
#SBATCH --requeue
```

and notes that a preempted job receives termination notice before it is killed/requeued.

For this research, checkpointing should exist anyway.

## NCC module discovery

```bash
source /etc/profile
module avail
module list
```

Do not assume Hamilton module names or versions exist on NCC.

---

# 29. NCC acknowledgement

Current NCC documentation requests the following acknowledgement for research made possible through NCC:

> This work has used Durham University’s NCC cluster. NCC has been purchased through Durham University’s strategic investment funds, and is installed and maintained by the Department of Computer Science.

Use only if NCC actually contributed to the research.

---

# 30. Hamilton acknowledgement

Durham publications commonly use wording such as:

> This work made use of the Hamilton HPC Service of Durham University.

Before submission, check the current ARC Hamilton acknowledgement page and use the wording requested there.

This research already used Hamilton for production simulations, so the eventual human-written paper/report should remember the Hamilton acknowledgement.

---

# 31. Campaign manifest

Every production campaign should have a machine-readable manifest committed or archived with results.

Example:

```toml
campaign_id = "squaregrid_validation_v1"
created_utc = "2026-08-16T00:00:00Z"

[code]
git_commit = "<40-char commit>"
julia_version = "<version>"
manifest_sha256 = "<sha256>"

[model]
geometry = "square_grid"
observable = "<exact estimator>"
disorder = "<law + parameters>"

[rng]
master_seed = 123456789
scheme = "<documented deterministic mapping>"

[slurm]
cluster = "hamilton8"
partition = "shared"
cpus_per_task = 16
memory = "16G"
walltime = "12:00:00"
array = "0-31"

[data]
raw_root = "/nobackup/<user>/research/raw/squaregrid_validation_v1"
checkpoint_root = "/nobackup/<user>/research/checkpoints/squaregrid_validation_v1"
```

Never infer the observable from a filename alone. Put the mathematical definition in the campaign config/manifest.

---

# 32. Automatic run metadata

At the start of each job, save:

```bash
META="$NOBACKUP/research/logs/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.meta.txt"

{
  echo "timestamp=$(date -Is)"
  echo "hostname=$(hostname)"
  echo "pwd=$PWD"
  echo "job_id=$SLURM_JOB_ID"
  echo "array_job_id=${SLURM_ARRAY_JOB_ID:-}"
  echo "array_task_id=${SLURM_ARRAY_TASK_ID:-}"
  echo "cpus=${SLURM_CPUS_PER_TASK:-}"
  echo "git_commit=$(git rev-parse HEAD 2>/dev/null || true)"
  echo "git_dirty=$(test -n "$(git status --porcelain 2>/dev/null)" && echo yes || echo no)"
  echo "julia=$(julia --version 2>&1 || true)"
  echo "julia_threads=${JULIA_NUM_THREADS:-}"
  echo "tmpdir=${TMPDIR:-}"
  module list 2>&1 || true
} > "$META"
```

This is cheap and prevents later ambiguity.

---

# 33. Scientific validation before Hamilton scaling

For this project, HPC correctness is part of scientific correctness.

Before a large run verify:

### Randomness
- [ ] independent environments really are independent;
- [ ] paired tilings/walks share only the intended environment;
- [ ] conditional samples have separate RNG streams;
- [ ] no accidental seed collision across array tasks;
- [ ] result is reproducible independent of thread count where mathematically expected.

### Observable
- [ ] exact \(T(r)\)/\(D(r)\)/winding definition documented;
- [ ] optimized implementation matches reference on tiny systems;
- [ ] spatial points/separations are consistent for every \(L\);
- [ ] sign/parity convention validated.

### Data
- [ ] one immutable shard per task;
- [ ] no concurrent appends to same file;
- [ ] row counts validated;
- [ ] checksums recorded for important outputs;
- [ ] no silent overwrite of previous campaigns.

### Statistics
- [ ] independent-environment block structure preserved;
- [ ] bootstrap resamples environment blocks;
- [ ] all \(r/L\) observations from one environment stay together;
- [ ] controls generated through the same pipeline;
- [ ] analysis runs on completed/validated shards only.

### Reproducibility
- [ ] clean Git commit;
- [ ] `Project.toml` and `Manifest.toml`;
- [ ] campaign config frozen;
- [ ] Slurm script archived;
- [ ] seed mapping documented;
- [ ] raw data copied/archived safely after campaign.

---

# 34. Performance checklist for the Julia simulations

Before optimization:
1. profile the current code;
2. identify the true hot loop;
3. determine memory scaling with \(L\);
4. benchmark one deterministic workload.

Common Julia/HPC priorities:
- avoid global variables in hot code;
- type-stable kernels;
- avoid unnecessary allocations;
- preallocate reusable arrays;
- avoid storing full environments if the **same stochastic model** can be generated/streamed exactly;
- never alter randomness/conditioning merely to reduce memory;
- batch output rather than printing per sample;
- avoid locks in the environment-level parallel loop;
- give threads independent result buffers and reduce afterwards.

Measure allocations locally/interactive before large scaling.

---

# 35. When to use threads vs arrays

## Use Julia threads when
- one simulation instance has useful shared-memory parallelism; or
- each thread independently processes environments with low synchronization;
- per-task setup is nontrivial and amortized across many samples.

## Use Slurm arrays when
- samples/parameter blocks are independent;
- failure/restart isolation matters;
- workload varies by \(L\);
- you want simple deterministic sharding.

## Recommended here

Use **both modestly**:

```text
many array shards × 8–32 Julia threads per shard
```

Benchmark exact values.

Do not default to:

```text
one array task × 128 threads
```

without scaling evidence.

---

# 36. Analysis jobs

Bootstrap/GLS analysis can be separate from simulation.

Example `slurm/analyse.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=dimer_analysis
#SBATCH --partition=shared
#SBATCH --time=00-02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module load <JULIA_MODULE>
export JULIA_DEPOT_PATH="$NOBACKUP/research/julia_depot"
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

julia --project=. --startup-file=no scripts/analyse.jl \
  --campaign squaregrid_pilot_v1
```

Do not mix raw simulation generation and final model fitting in one opaque job if separating them makes the pipeline auditable.

---

# 37. Failure handling

A failed shard is not a failed campaign.

Maintain a list of completed shard IDs.

Find missing array outputs conceptually:

```bash
for i in $(seq 0 31); do
  f="$NOBACKUP/research/raw/$CAMPAIGN/part_$(printf "%05d" "$i").csv"
  [[ -s "$f" ]] || echo "$i"
done
```

Then rerun only missing IDs:

```bash
sbatch --array=3,7,11,19 slurm/sim_array.sh
```

This is much cheaper and safer than rerunning everything.

---

# 38. Long-job design

Even when Hamilton's `long` partition allows up to 7 days:

- checkpoint;
- split by independent environment ranges;
- save progress periodically;
- make restart deterministic;
- avoid a monolithic single file;
- do not rely on a live SSH session.

A long batch job should survive the laptop being shut or disconnected.

---

# 39. What Codex should know about Hamilton

When asking Codex to generate/modify HPC scripts, point it to this file and constrain it.

Recommended instruction:

```text
Read context/HAMILTON_HPC.md before editing any Slurm/HPC code.

Target Hamilton8 unless the task explicitly says NCC.
Do not invent partition names, Julia module versions, storage paths, or GPU directives.
Use job arrays for independent Monte Carlo shards where appropriate.
Use deterministic environment IDs/seeds.
Keep jobs restart-safe.
Do not request a whole node without benchmark evidence.
Do not change the scientific estimator to make the job easier to parallelize.
```

For an HPC-focused task, Codex should read:
1. `AGENTS.md`;
2. `context/HAMILTON_HPC.md`;
3. only the relevant scientific section of `context/RESEARCH_BRIEF.md` / `context/NEXT_STEPS.md`;
4. the specific scripts/code being edited.

It should **not** read the entire email archive.

---

# 40. Recommended first HPC audit by Codex

Before changing production scripts:

```text
Read AGENTS.md and context/HAMILTON_HPC.md.

Inspect only the current Julia simulation entry points and Slurm scripts.

Do not edit yet.

Report:
1. current parallelism model;
2. how JULIA_NUM_THREADS is set;
3. whether BLAS/native thread oversubscription is possible;
4. how seeds are partitioned across jobs/tasks/threads;
5. whether outputs are race-free and restart-safe;
6. current memory/I/O pattern and whether $TMPDIR would help;
7. resource requests vs actual algorithm requirements;
8. the smallest set of changes that would improve Hamilton8 reliability and throughput.

Do not change the stochastic model or estimator.
Keep the response concise.
```

Only after reviewing that audit should Codex implement HPC changes.

---

# 41. Quick command sheet

## Connect

```bash
ssh <user>@hamilton8.dur.ac.uk
```

## Sync code

```bash
rsync -av --progress ~/Desktop/research/ <user>@hamilton8.dur.ac.uk:~/research/repo/
```

## Module discovery

```bash
module avail julia
module spider julia
module load <julia-module>
julia --version
```

## Julia

```bash
export JULIA_DEPOT_PATH="$NOBACKUP/research/julia_depot"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Interactive Slurm test

```bash
srun --pty --mem=2G -c 2 -p test bash
```

## Submit

```bash
sbatch slurm/job.sh
```

## Monitor

```bash
squeue --me
tail -f logs/<log>
```

## Inspect

```bash
scontrol show job <jobid>
sstat -j <jobid>.batch
sacct -j <jobid> --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

## Cancel

```bash
scancel <jobid>
```

## Storage

```bash
quota
df -h "$NOBACKUP"
du -sh "$NOBACKUP/research"/*
```

## Pull results back to Mac

```bash
rsync -av --progress \
  <user>@hamilton8.dur.ac.uk:/nobackup/<user>/research/derived/ \
  ~/Desktop/research-results/
```

---

# 42. Current official/reference sources

Verified on 16 August 2026.

Durham remote access:
- Access VPN secure logon: https://www.access.durham.ac.uk/

Hamilton:
- Systems: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/systems/
- Login: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/login/
- Storage: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/storage/
- Internal storage: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/storage/filestorage/
- Data transfer: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/storage/datatransfer/
- Running jobs: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/jobs/
- CPU examples: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/jobs/examplecpujobs/
- Long jobs: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/jobs/long/
- Software: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/software/
- Portal: https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/usage/portal/

NCC:
- Main/current documentation: https://nccadmin.webspace.durham.ac.uk/
- Slurm cheat sheet: https://nccadmin.webspace.durham.ac.uk/cheatsheet/
- Preemption: https://nccadmin.webspace.durham.ac.uk/preemption/

Because HPC configurations change, **live Slurm/module output wins over this document** if the two disagree:

```bash
sinfo
sfree
module avail
module spider julia
```

Update this file when Durham changes hardware, partitions, quotas, or modules.
