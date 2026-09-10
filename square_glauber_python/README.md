# Square-grid random-weight dimers: local Python pilot

This isolated Python 3.11+ package samples perfect matchings of an even
`L x L` square grid in one frozen undirected edge environment. It is a local
proof-of-correctness and mixing/pilot implementation. It contains no Slurm,
Hamilton, active-face jump-chain, checkerboard-parallel, or production-campaign
code.

## Mathematical contract

Vertices have array coordinates `(row, column)`, with row increasing downward.
A vertex is **white** when `(row + column) % 2 == 0`, and black otherwise.
Every primal edge is oriented white-to-black. Crossing that oriented edge on
the dual lattice from its geometrical left side to its right side is clockwise
around white and changes height by `+1` if unoccupied or `-3` if occupied.
The reverse crossing has the opposite sign.

The finite induced grid has boundary vertices of degree two or three, so the
degree-four `+1/-3` rule cannot consistently identify every boundary sector as
one exterior dual vertex. The additive convention therefore uses a fixed dual
cut: exterior height is zero and the cut crosses the top-left horizontal edge
into bounded face `(0,0)`. Its occupancy-dependent increment is retained; no
interior face is reset separately per configuration. All remaining faces are
integrated through interior edges, with all alternate dual paths checked.

The matching is authoritative. Heights are derived only when measured.
Horizontal occupancy has shape `(L, L-1)`, vertical occupancy `(L-1, L)`, and
bounded-face heights `(L-1, L-1)`. The unique central face is
`((L-2)//2, (L-2)//2)`.

At every update the code chooses uniformly from **all** `(L-1)^2` faces. A
nonflippable face is a null move. On a flippable face with top/right/bottom/left
weights `(a,b,c,d)`, it redraws horizontal with probability
`a*c/(a*c+b*d)` and vertical with probability `b*d/(a*c+b*d)`. A redraw of the
existing orientation is also an unchanged heat-bath move. One sweep means
exactly `(L-1)^2` attempted updates.

## Install and test

```bash
cd square_glauber_python
python3 -m venv .venv
.venv/bin/python -m pip install '.[test]'
.venv/bin/python -m pytest
```

## Commands

Exact enumeration, stochasticity, stationarity, and detailed balance for
2x2, 2x4, and 4x4 grids:

```bash
.venv/bin/python -m square_glauber.cli validate-exact
```

Write full `L=6,8` extremal height matrices and labelled dimer/height figures:

```bash
.venv/bin/python -m square_glauber.cli height-audit \
  --sizes 6,8 --output outputs/height_convention_audit
```

Compute exact deterministic-weight centre-height PMFs on 2x2, 2x4, and 4x4,
then compare a long 4x4 Glauber diagnostic with the exact law:

```bash
.venv/bin/python -m square_glauber.cli exact-height-audit \
  --burnin-sweeps 10000 --measurement-gap-sweeps 2 \
  --num-measurements 100000 --seed 20260908 \
  --output outputs/exact_height_audit
```

Mixing traces from independent all-horizontal and all-vertical starts in one
frozen environment:

```bash
.venv/bin/python -m square_glauber.cli mixing \
  --L 16 --model gamma --gamma-shape 1.0 \
  --sweeps 5000 --record-every 5 --seed 20260907 \
  --output outputs/mixing_L16_gamma
```

Each mixing case records centre height, total bounded-face height, spatial
increments, ACF tables, IACT estimates in record and sweep units, effective
sample sizes, late-window comparisons, and a lightweight split-R-hat.

Run a multi-size local diagnostic campaign without starting the paired scaling
experiment:

```bash
.venv/bin/python -m square_glauber.cli mixing-campaign \
  --sizes 8,12,16 --gamma-environments 3 --uniform-replicates 1 \
  --gamma-shape 1.0 --sweeps 20000 --record-every 10 \
  --seed 20260908 --output outputs/mixing_validation_20260908_20000_sweeps
```

Extend only the saved Gamma environments, verifying the recorded seed tuple,
bitwise equality with seed regeneration, SHA-256 checksums, and exact agreement
of the new trace prefix with the saved 20,000-sweep trace:

```bash
.venv/bin/python -m square_glauber.cli mixing-extend-frozen \
  --source outputs/mixing_validation_20260908_20000_sweeps \
  --sweeps-by-size 8:100000,12:200000,16:200000 \
  --record-every 10 --acf-max-lag 5000 \
  --output outputs/mixing_validation_20260908_extended
```

Each extended case also records edge-weight and face log-odds hardness metrics,
four consecutive post-transient window summaries, per-chain MCSE/IACT/ESS,
signed between-start `z_gap`, and the existing centre/global/spatial traces and
ACFs. Hardness metrics are descriptive and are never used to exclude an
environment. The command refuses to overwrite a nonempty output directory.

Small paired local pilot (sample counts and schedules are always explicit):

```bash
.venv/bin/python -m square_glauber.cli pilot \
  --sizes 8,12,16,24,32 --models uniform,gamma --gamma-shape 1.0 \
  --environments 50 --burnin-sweeps 1000 \
  --measurement-gap-sweeps 100 --num-measurements 1 \
  --seed 20260907 --output outputs/pilot_001
```

Exploratory, whole-environment bootstrap analysis:

```bash
.venv/bin/python -m square_glauber.cli analyze \
  --input outputs/pilot_001 --output outputs/pilot_001/analysis
```

The analysis first averages repeated measurements within each environment.
It never counts them as independent environments. It reports
`0.5*Var(H1-H2)`, `Cov(H1,H2)`, their spatial analogues, ordinary-log versus
quadratic-log descriptive fits, a pooled spatial shared-quadratic fit,
environment-blocked bootstrap intervals, diagnostics, and simple matplotlib
figures. These outputs are exploratory and do not certify mixing.

## Output tables

- `center_measurements.csv`: one paired measurement row, retaining all seeds
  and `environment_id`.
- `spatial_measurements.csv`: long-form paired increments with requested
  fraction and actual integer `r`.
- `chain_diagnostics.csv`: attempted, flippable, changed, null, and runtime
  diagnostics for each chain/environment.
- `metadata.json`: schedule, seed derivation, runtime, and scientific caveats.

Generated `outputs/` and `figures/` contents are ignored except for placeholders.

## Optional Numba acceleration

The pure-Python implementation remains the authoritative reference backend.
Install the separate optional Numba backend with a Python version supported by
your Numba release:

```bash
python -m pip install '.[test,accel]'
```

Both backends accept deterministic row-major `face_indices` and `uniforms`
arrays for exact checkpoint-by-checkpoint equivalence tests. The ordinary
accelerated RNG API still chooses uniformly from all bounded faces, allows null
moves, and applies the identical heat-bath orientation probability. Persistent
independent PCG64 substreams draw one face index and one uniform per attempted
update; the uniform is ignored for a nonflippable face. This makes trajectories
independent of update chunking and recording cadence. Heights are derived only
at requested recording times.

Run the equivalence audit and non-scientific kernel benchmark:

```bash
python -m square_glauber.cli backend-equivalence \
  --sizes 4,6,8,12,16 --attempts 20000 --checkpoint-attempts 257 \
  --output outputs/performance_validation/equivalence

python -m square_glauber.cli benchmark-backends \
  --sizes 8,12,16,24,32 --target-attempts 2000000 --repeats 3 \
  --output outputs/performance_validation/benchmark
```

The common-randomness coupling is deliberately separate from the independent
sampler and requires a successful finite ordering audit before it will run:

```bash
python -m square_glauber.cli coupling-ordering-audit \
  --sizes 4,6,8,12,16 --attempts 20000 --streams 2 \
  --output outputs/performance_validation/ordering

python -m square_glauber.cli accelerated-stress \
  --source-case outputs/mixing_validation_20260908_extended/gamma-L16-diagnostic-02 \
  --sweeps 1000000 --record-every 10 --acf-max-lag 20000 \
  --output outputs/performance_validation/gamma-L16-diagnostic-02-numba-1m

python -m square_glauber.cli common-randomness-coupling \
  --source-case outputs/mixing_validation_20260908_extended/gamma-L16-diagnostic-02 \
  --ordering-audit outputs/performance_validation/ordering/ordering_metadata.json \
  --sweeps 1000000 --record-every 10 \
  --output outputs/performance_validation/gamma-L16-diagnostic-02-coupling
```

Finite ordering and coupling checks are exploratory; they are not proofs of
attractiveness, mixing or equilibrium, and coupled chains are not independent
samples.

## Final frozen-environment calibration and order audit

The final calibration commands keep the reference/Numba kernels unchanged and
refuse nonempty output targets:

```bash
python -m square_glauber.cli exact-extremality-audit \
  --output outputs/final_mixing_and_coupling_audit_20260908/extremality

python -m square_glauber.cli exhaustive-monotonicity-audit \
  --output outputs/final_mixing_and_coupling_audit_20260908/monotonicity

python -m square_glauber.cli final-hard-calibration \
  --source-case outputs/mixing_validation_20260908_extended/gamma-L16-diagnostic-02 \
  --one-million-case outputs/performance_validation_20260908/gamma-L16-diagnostic-02-numba-1m-stream-safe \
  --record-every 10 --acf-max-lag 100000 --windows 16 \
  --output outputs/final_mixing_and_coupling_audit_20260908/hard_L16_3m
```

Exact enumeration shows that all-horizontal/all-vertical are useful distinct
diagnostic starts but are **not** the general pointwise minimum/maximum height
tilings (already false on 4x4; all-vertical is not maximal on 2x4). Therefore
their earlier common-randomness coalescence is not an extremal sandwich or a
CFTP certificate.

## Certified perfect sampling

The optional accelerated backend now also supports monotone coupling from the
past using the **true** minimum and maximum height tilings. It replays the same
ordinary all-face heat-bath random maps from increasingly remote past times and
returns a state only when the true bounds agree at time zero. No burn-in is
chosen and no uncertified state is returned.

The implementation includes a domain-independent exhaustive local
monotonicity truth table, exact tiny-grid verification of the constructed
bounds, reference-versus-Numba replay equivalence, an exact weighted 4x4 law
diagnostic, and a bounded frozen-hard-environment horizon calibration. See
[`PERFECT_SAMPLING.md`](PERFECT_SAMPLING.md) for the mathematical argument,
limitations, references, and commands. These tools do not launch a scaling
experiment or fit a logarithmic growth model.

Run a bounded multi-environment certification-cost calibration and, if any
streams hit the declared cap, extend exactly those saved histories:

```bash
python -m square_glauber.cli cftp-cost-calibration \
  --standard-sizes 8,12,16 \
  --standard-gamma-environments 30 --standard-uniform-environments 1 \
  --standard-samples-per-environment 2 \
  --probe-sizes 24,32 --probe-gamma-environments 5 \
  --probe-uniform-environments 1 --probe-samples-per-environment 1 \
  --max-horizon-sweeps 1048576 --backend numba --seed 20260909 \
  --output outputs/cftp_cost_calibration_20260909

python -m square_glauber.cli cftp-extend-censored \
  --source outputs/cftp_cost_calibration_20260909 \
  --max-horizon-sweeps 8388608 --backend numba \
  --output outputs/cftp_cost_calibration_20260909_extensions_to_8m

python -m square_glauber.cli cftp-consolidate-calibration \
  --source outputs/cftp_cost_calibration_20260909 \
  --extension outputs/cftp_cost_calibration_20260909_extensions_to_8m \
  --output outputs/cftp_cost_calibration_20260909_final
```

Every environment array, environment/CFTP seed, SHA-256, hardness diagnostic,
certification horizon and runtime is retained. Censored streams are not treated
as samples and are never replaced with easier environments. These commands
measure computational cost only; they do not calculate covariance or fit
`log L`/`(log L)^2` models.

## Hamilton paired perfect-sampling campaign

The production workflow uses two independent certified CFTP samples in every
frozen environment. It therefore has no burn-in or measurement-spacing choice.
It retains the centre height and all configured spatial increments for both
samples, along with the global height sum as a supplementary diagnostic. An
environment contributes no scientific row until both samples certify.

The frozen launch schedule is
`configs/perfect_covariance_campaign_20260909.csv`; the authoritative expanded
manifest is `configs/perfect_covariance_campaign_20260909_launch_manifest.csv`.
It contains 18,200 Gamma(shape=1) environments and 5,600 uniform-control
blocks over `L=8,12,16,20,24,28,32`, for 47,600 perfect samples in total.
The master seed is `2026090901`. The large horizon caps are censoring limits,
not burn-in lengths: most samples certify far earlier, while every capped case
is saved for identical-history extension.

On Hamilton, after inspecting the available Python executable, prepare and
test the isolated environment:

```bash
bash hpc/setup_hamilton_python.sh
```

Then launch the smoke-gated production array into a new `/nobackup` root:

```bash
OUTPUT_DIR=/nobackup/$USER/square_glauber_perfect_covariance_20260909 \
  CONCURRENCY=32 bash hpc/launch_perfect_covariance_campaign.sh
```

The final audit job writes `campaign_status.csv` even if jobs timed out or a
CFTP stream was censored. Resubmit only the affected manifest task IDs against
the same output root, optionally with a larger cap:

```bash
sbatch --array=TASK_IDS%8 \
  --export=ALL,MANIFEST=/nobackup/$USER/square_glauber_perfect_covariance_20260909/campaign_manifest.csv,OUTPUT_DIR=/nobackup/$USER/square_glauber_perfect_covariance_20260909,MAX_HORIZON_SWEEPS=134217728 \
  hpc/perfect_campaign_task.slurm
```

Only after `all_pairs_cftp_certified` is true should the existing
environment-blocked analysis run:

```bash
python -m square_glauber.cli analyze \
  --input /path/to/completed/campaign \
  --output /path/to/completed/campaign/analysis \
  --bootstrap-replicates 10000 --bootstrap-seed 2026090902
```

This analysis compares ordinary-log and log-squared descriptions and
bootstraps whole environments jointly. The campaign is designed to test the
super-roughness prediction; neither the sampler nor the fitting code assumes
that the quadratic coefficient is positive.
