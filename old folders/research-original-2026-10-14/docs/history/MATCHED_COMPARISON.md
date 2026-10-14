# Matched Aztec/square comparison — 29 September 2026

This is a new workflow. September's manifest, rightward observables, seeds,
measurements and archived reports are unchanged. No Hamilton job has been
submitted. AI contribution reconstruction remains a subsequent task.

## Most important finding: the historical Aztec benchmark needs reassessment

The new distribution-level audit finds an orientation mismatch between the
paper-specific Gamma reduction and the tiling matrix used by the legacy sampler.
The original tests checked recurrence values, reproducibility and valid tilings;
those checks did not compare the resulting weighted distribution with the
occupied-edge Gibbs law on a nonuniform tiny diamond.

Figure 1 of [Duits–Van Peski](https://arxiv.org/abs/2512.03033) assigns `a` to red
edges and `b` to yellow edges. In this repository's one-based domino matrix:

```text
weight[2i, 2j-1] = a[i,j]     red
weight[2i, 2j]   = b[i,j]     yellow
weight[2i-1, *]  = 1         blue/green
```

For fixed arbitrary positive arrays, exact enumeration of all creation-coin
paths gives the following comparison with normalized occupied-edge products:

| Order | Tilings | Legacy maximum probability error | Corrected maximum error |
|---:|---:|---:|---:|
| 2 | 8 | 0.09305996 | 5.55e-17 |
| 3 | 64 | 0.04069997 | 4.16e-17 |

The total variation discrepancies are 0.15661 and 0.25658 respectively. The
generic weighted-shuffling route also agrees with the Gibbs probabilities to
roundoff. These are deterministic exact finite experiments, not noisy sampling
comparisons. See [the machine-readable audit](../validation/aztec_mapping_audit.json).

The corrected adapter in `aztec/scripts/comparison_geometry.jl`:

1. reverses the first axis of both input family arrays into the paper's convention;
2. applies the existing Gamma reduction;
3. reverses the first axis of every creation-probability/coin table back into the
   tiling matrix convention.

For example, write `s_ij = a_ij + b_ij`. At order two, the generic
red/yellow law gives the reduced creation probability
`b22*s12*s21 / (b22*s12*s21 + a21*s11*s22)`. The legacy Gamma path
instead gives `b12*s11*s22 / (b12*s11*s22 + a11*s12*s21)`. These coincide
for spatially constant family weights, which explains why uniform or constant
controls do not expose the mismatch. The reflections reconcile the formulas.

Its probability tables agree with the independent generic weighted reduction
through orders 1–8 in the regression suite. The corrected packed-coin route is
used by the new small-order bridge. Legacy functions remain unchanged so that
historical datasets can be reproduced. **Use the corrected bridge for new
matched work. Do not use the old positive Aztec fit as a validated benchmark for
the stated quenched model until it has been regenerated/reassessed.**

This finding does not determine what the corrected large-size scaling will be.
It also does not imply that valid historical tilings can be relabelled as draws
from the intended frozen environment. A corrected large campaign is a separate,
prospectively specified task after feasibility and Hamilton access.

## Exact square weight law

Python vertices are `(row, column)`, with rows increasing downward and white
vertices having even parity. From a black (odd-parity) vertex:

- the north edge has independent `Gamma(alpha=0.2, scale=1)` weight `a`;
- the east edge has independent `Gamma(beta=0.25, scale=1)` weight `b`;
- south and west edges have weight 1.

Each undirected edge is assigned exactly once. Equivalently, horizontal edges
whose left coordinate has odd parity carry `b`; vertical edges whose top
coordinate has even parity carry `a`. This is the red/yellow law after undoing
the paper's display rotation, up to a lattice translation choosing the colour
phase. The square uses a fixed phase across sizes. Different domains have
different boundary cuts; no claim of identical finite graphs is made.

On the 2×2 graph the horizontal pair has weight `b`, the vertical pair has weight
`a`, and the heat bath selects horizontal with `b/(a+b)`. This independent
hand calculation checks the orientation convention.

Three model names are explicit:

- `aztec_gamma`: the random law above;
- `aztec_mean`: the same placement with deterministic weights 0.2 and 0.25;
- `uniform`: every edge one.

The second control detects effects of the family bias separately from disorder.
The two Gamma families are **not mean-one normalised**. Normalising them
separately changes local matching ratios; it is not silently treated as gauge
invariance. The existing all-edge mean-one Gamma sampler is untouched.

## Matched probe geometry

Aztec's central row is a row of the *staggered height table*, not a horizontal
line of the unrotated underlying square lattice. Its one-based face `(R,C)` has
physical coordinates, in units of one primal edge,

```text
x = C - 1 - floor((R-1)/2)
y = n + 1 - C - ceil((R-1)/2).
```

Changing `C` by `r` changes `(x,y)` by `(r,-r)`. Thus the primary square probe
is a symmetric **diagonal** increment with row/column displacement `(r,r)`;
its Euclidean distance is `sqrt(2)*r`. The exact tiny-diamond audit also checks
the height against an independent two-edge dual path: the fluctuation scale is
one, with a global sign depending on order parity. Covariance is unaffected by
that sign. The new bridge uses even orders.

The square also records a symmetric horizontal probe as a diagnostic. These
are never pooled with the diagonal primary probes. For side `L`, `r` is
`floor(L*rho+0.5)`, and the endpoints straddle the central bounded face as closely
as parity permits. Odd separations use a fixed negative-half-step midpoint
shift. Zero and duplicate offsets are omitted and the exact endpoint list is
saved in `probes-L*.json`. These omissions make very small sizes unsuitable for
a four-probe curvature fit.

Size labels remain domain-specific. The analysis fits recorded physical
distance, not the assumption that Aztec order equals Python vertex side. Probe
direction and height scale are matched; local slopes, liquid-region location and
boundary distance remain scientific considerations when choosing production
windows. Matching the law alone does not isolate boundary conditions as the
unique cause of a difference.

## Local work performed

### September controls and pairing

`control_audit.py` read the frozen centre/spatial tables, verified deterministic
seed derivation, global seed uniqueness and all centre/spatial joins. The 23,799
pair count is unchanged. It then jointly tested 30 all-one covariance observables.

The test permutes entire second-replica vectors within each size, keeping centre
and every probe together. It uses the largest absolute covariance divided by
its exact permutation-null standard deviation. With 9,999 permutations and
seed 2026092903, the joint statistic is **2.26643** and the plus-one randomisation
p-value is **0.4827**. The six individually nonzero intervals therefore do not
constitute joint evidence against the all-one pairing null in this audit.
This is not a proof that the single-replica Gibbs law or all uncertainty
calculations are correct.

Report: [september_control_audit.json](../validation/september_control_audit.json).

### Historical state availability

The remeasurement tool found **200** complete certified state pairs matching the
snapshot's weights, seeds and original central heights. It recorded new
symmetric probes without sampling. **23,599** snapshot pairs lack available
verified state pairs in the searched local archive. The new table is explicitly
exploratory; its available subset is not a complete-campaign scaling dataset.

Availability report: [september_state_availability.json](../validation/september_state_availability.json).
Full missing identifiers and remeasured values are in the ignored output folder
`square_glauber_python/outputs/september_symmetric_remeasurement/`. Once logged
into Hamilton, locate the original state files and use this same read-only
workflow. The unresolved missing *campaign* pair is an additional issue beyond
local state availability; no attempt was made to replace or restart it here.

### Matched square feasibility

The fixed local design requested 4 environments at each of sides 4, 8, 12 in
each of the three models: **36 environments / 72 streams**, capped at 4096
sweeps. This is only a software and cost pilot.

| Model | Side 4 complete pairs | Side 8 | Side 12 |
|---|---:|---:|---:|
| Gamma two-family | 1/4 | 0/4 | 0/4 |
| Deterministic family means | 4/4 | 4/4 | 4/4 |
| All one | 4/4 | 4/4 | 4/4 |

There are **21 censored Gamma streams**, retained with checkpoints; all 48
control streams certified. Completed-only runtime/horizon quantiles are labelled
as such. Censored horizons are lower bounds, not observed coalescence times.
No scaling fit was run on this selectively complete Gamma subset.

Compact outputs: [matched_local_summary.json](../validation/matched_local_summary.json)
and [matched_local_summary.csv](../validation/matched_local_summary.csv). Weights,
states, certificates, checkpoints and all stream statuses remain in
`square_glauber_python/outputs/matched_local_20260929/`.

### Corrected Aztec bridge

The corrected bridge generated only **20** paired environments: four each at
orders 8, 16, 24, 32 and 48. It saves both weight arrays, both tilings, seed,
coordinates and increments in the common comparison schema. This is a workflow
check, not a replacement positive scaling result. Outputs live in
`aztec/output/comparison_bridge_20260929/`.

## Reproducible commands

Run from the repository root. Use a supported Python (3.11+) and Julia (1.10+).
The local validation environment is `.venv/`; it is isolated and ignored by Git.
The pre-existing project environments have broken interpreter symlinks on this
machine and were not overwritten. Source imports are explicit below.

```bash
export PYTHONPATH="$PWD/square_glauber_python/src"
export MPLCONFIGDIR="$PWD/.runtime/matplotlib"
export NUMBA_CACHE_DIR="$PWD/.runtime/numba"

.venv/bin/python -m pytest square_glauber_python/tests

.venv/bin/python -m square_glauber.control_audit \
  --input archive/history/square_glauber_python/reports/provisional_20260915/input \
  --output docs/validation/september_control_audit.json

.venv/bin/python -m square_glauber.comparison_pilot \
  --config square_glauber_python/configs/matched_feasibility_local.json \
  --output square_glauber_python/outputs/matched_local_20260929

.venv/bin/python -m square_glauber.remeasure \
  --snapshot archive/history/square_glauber_python/reports/provisional_20260915/input/center_measurements.csv \
  --states archive/local_runs \
  --output square_glauber_python/outputs/september_symmetric_remeasurement

julia --project=aztec aztec/scripts/audit_comparison_mapping.jl
julia --project=aztec aztec/scripts/run_comparison_bridge.jl \
  --output aztec/output/a_new_comparison_bridge

typst compile --root . docs/research/handbook.typ output/pdf/research-handbook.pdf
```

The downloaded local Julia 1.10 runtime uses `.runtime/aztec-checks` (a copy of
`Project.toml` linked to the actual source) and `.runtime/julia-depot`, avoiding
changes to the pre-existing Julia 1.12 manifest. For this local machine replace
`julia --project=aztec` by:

```bash
JULIA_DEPOT_PATH=.runtime/julia-depot \
  .runtime/julia-1.10.10/bin/julia --startup-file=no --project=.runtime/aztec-checks
```

### Pilot restart contract

The output manifest fixes models, parameters, seed, sizes and environment
counts. Raising a power-of-two horizon cap is allowed; other identity changes
require a new output directory. Weights are regenerated and compared bitwise
on load, certificates and saved state hashes are checked, and a single-writer
lock protects the output root. Both replicas are attempted even when the first
is censored. The existing checkpoint implementation replays nested random maps
from true extrema at each doubled horizon. An interrupted horizon is replayed;
there is no mid-horizon continuation claim.

A timeout can leave the global progress table incomplete; the per-stream
records/checkpoints are authoritative and rerunning rebuilds the tables. A
fully completed pilot invocation writes a summary and measurement checksum.
The runner deliberately refuses sizes above 32, more than 64 environments per
group, more than 1024 total environments or caps above 2^20. It is not a
production launcher.

### Direct curvature comparison, when suitable data exist

`comparison_analysis.py` accepts the common long-form tables and their
completeness summaries. It verifies the table hashes, the same single law,
diagonal direction and common height units. The caller must declare windows
and probes; incomplete vectors, fewer than four sizes, collapsed distances and
rank-deficient designs are refused. It fits a separate intercept and log slope
per probe with shared curvature, and independently bootstraps whole environment
vectors in each domain. The output includes each coefficient and the interval
for `c_aztec - c_square`, residuals, conditioning and log/quadratic SSE.

Example syntax **for future complete, adequately sized data** (not the current
pilot, which must be refused):

```bash
.venv/bin/python -m square_glauber.comparison_analysis \
  --aztec AZTEC_MEASUREMENTS.csv --aztec-summary AZTEC_SUMMARY.json \
  --square SQUARE_MEASUREMENTS.csv --square-summary SQUARE_SUMMARY.json \
  --law-id 'aztec_gamma:alpha=0.2:beta=0.25:scale=1:a-north,b-east' \
  --aztec-window MIN_ORDER MAX_ORDER --square-window MIN_SIDE MAX_SIDE \
  --probes diagonal:rho=0.03125 diagonal:rho=0.0625 \
           diagonal:rho=0.125 diagonal:rho=0.25 \
  --bootstrap 10000 --output CONTRAST.json
```

No production window or sample budget has been declared after viewing pilot
fluctuation values. Before production, freeze windows, complete probe coverage,
fixed cutoff sensitivities, control criteria, runtime/stopping rules and planned
precision. The current comparison is unweighted OLS with environment bootstrap;
SSE is descriptive. It does not silently reuse the old BIC as a full likelihood,
or claim a power calculation has been completed.

## Hamilton execution: 30 September 2026

Authenticated access to `hamilton8.dur.ac.uk` was restored. An isolated source
snapshot is staged at `/nobackup/fvkl37/matched_comparison_20260929/code`;
all file hashes matched `docs/validation/hamilton_20260929/source_manifest.json`.
Validation job **18852141** completed with exit code `0:0` in 59 seconds on the
`test` partition. All 108 Python tests and the Julia suite passed, including
the exact distribution audit. Logs and runtime versions are retained in
`docs/validation/hamilton_20260930/`.

After validation, capped feasibility job **18852143** and historical-state
remeasurement job **18852144** were submitted to `shared`, each with one CPU,
4 GB RAM and a two-hour wall-time limit. Their outputs are respectively
`feasibility/` and `september_remeasurement/` under the isolated remote root.
Both jobs completed with exit code `0:0`: feasibility in 15 minutes 57 seconds,
remeasurement in 50 minutes 58 seconds. Compact outputs have been downloaded to
`square_glauber_python/outputs/hamilton_matched_20260930/`; retained summaries and
verification records are in `docs/validation/hamilton_20260930/`.

All 160 control pairs certified (80 uniform, 80 deterministic-mean). None of
80 random-weight pairs certified within 65,536 sweeps: only two individual
streams at side 8 certified, with 158 of 160 random streams censored. This is a
runtime/certification result, not evidence for or against squared-log scaling.
No random-model paired covariance or curvature can be estimated from this pilot.
Do not launch a larger campaign with the present sampler budget; first investigate
coalescence costs and possible exact-sampling improvements on retained cases.

All 23,799 snapshot pairs were recovered, yielding 147,392 exploratory probe
rows. Environment IDs exactly match the frozen snapshot, there are no duplicate
environment/probe rows, and the downloaded table SHA-256 matches Hamilton.
This recovers the completed September snapshot, not the unresolved 23,800th
planned pair. The historical all-edge Gamma(1) law remains distinct from the
new two-family matched model.

The September centre and spatial CSV hashes match the archived local inputs.
The original campaign records 23,799 completed pairs and one censored pair.
Its inputs and restart files remain untouched; the unresolved pair was not
restarted. Remeasurement is read-only and produces new exploratory observables.
No production scaling campaign has been submitted.

The capped design is `configs/matched_feasibility_hamilton.json`: sides 8, 12,
16, 24 and 32; 16 environments per model/size; cap 65,536 sweeps. Its master seed
is distinct from the local pilot. `hpc/matched_feasibility.slurm` is a prepared
single-task wrapper. It has no automatic submission or credential access.

Remaining execution checks:

1. Review the retained feasibility summaries and censored histories before
   deciding on further sampling.
2. Analyse recovered September symmetric probes as exploratory observables,
   retaining the original missing-pair caveat and distinct weight law.
3. Preserve the remote certificates and checkpoints; no automatic resubmission.
4. Inspect all stream outcomes and censoring by group. Extend identical capped
   histories if worthwhile; do not replace hard cases.
5. Decide whether corrected Aztec scaling can be regenerated and what square
   sizes are feasible. The orientation finding makes this reassessment essential
   before a paper-level boundary comparison.

The AI contribution discussion follows this stage. No historical attribution
has been invented and no email has been sent.


## Subsequent five-step assessment

The bounded follow-up, uncertainty analysis, corrected Aztec pilot, production
hold decision and AI record are documented in [FOLLOWUP_20260930.md](FOLLOWUP_20260930.md).
The corrected pilot is inconclusive; the September diagonal curvature is strongly
cutoff-sensitive; all eight fixed diagnostic CFTP streams remain uncoalesced
at 262,144 sweeps. Rare-probability rounding also needs numerical review.
See [PRODUCTION_DECISION.md](PRODUCTION_DECISION.md) for the conditional design
and [AI_CONTRIBUTIONS.md](../AI_CONTRIBUTIONS.md) for confirmed assistance and
unresolved historical attribution. No production run was launched.
