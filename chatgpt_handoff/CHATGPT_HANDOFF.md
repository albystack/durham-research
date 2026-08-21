# Aztec Research Handoff

## 1. Executive Summary

Workspace root: `/Users/alberto/Desktop/research`.

The active project is the Julia package in `aztec/`. It simulates Gamma-disordered Aztec-diamond domino tilings, computes center heights, shared-environment double-dimer observables, and the newer paired spatial height-increment experiment. The authoritative implementation is `aztec/src/AztecDiamond.jl` plus the CLI scripts in `aztec/scripts/`.

Current repository state:

- Git branch: `main`
- HEAD commit: `7b3f355` (`spatial-aztec diamond experiment commit`)
- Dirty tree: clean in tracked files; generated output lives under ignored `aztec/output/`

Validated status:

- `julia --project=aztec -e 'using Pkg; Pkg.test()'` passed.
- `sh aztec/test/smoke_workflows.sh` passed.

The strongest completed Aztec result is the pooled spatial Gamma covariance signal in `aztec/results/spatial/README.md` and `aztec/results/spatial/spatial_analysis_report.md`: unweighted pooled `c = 0.5843` with 95% interval `[0.1957, 0.9520]`, and weighted pooled `c = 0.4948` with interval `[0.1932, 0.7977]`. This is numerical evidence only, not proof.

## 2. Scientific Question

The project tests whether disorder changes variance growth from ordinary logarithmic scaling to squared-logarithmic scaling in three settings:

1. center-height variance in the Gamma-disordered Aztec diamond;
2. conditional-versus-disorder variance decomposition from paired double-dimer samples;
3. paired spatial height increments at symmetric central-row separations.

The same variance identities appear repeatedly:

```text
Var(H1-H2)/2 = E[Var(H | environment)]
Cov(H1,H2)   = Var(E[H | environment])
```

The first isolates conditional tiling noise; the second isolates disorder-induced fluctuation.

## 3. Full Research Chronology

1. Earlier loop-erased random walk work was archived under `old/random_walk/`.
2. The Aztec-diamond project was introduced and the active code was consolidated under `aztec/`.
3. Single-height Gamma campaigns were run and merged into `aztec/data/height/center_height_samples.csv`.
4. Shared-environment double-dimer campaigns were run and merged into `aztec/data/double_dimer/pairs.csv`.
5. Height and double-dimer analyses were retained under `aztec/results/height/` and `aztec/results/double_dimer/`.
6. The spatial Gamma-disorder experiment and uniform control were added later, with analysis and plots in `aztec/results/spatial/`.
7. The current commit mainly adds the spatial experiment on top of the earlier height and double-dimer work.

## 4. Exact Mathematical Models

### 4.1 Gamma-disordered Aztec diamond

- Domain: order-`L` Aztec diamond represented as a `2L × 2L` Boolean matrix.
- Boundary condition: top boundary of the staggered height table starts at `A[1,j] = 2j - 2`.
- Environment: two independent weight families `a[i,j] ~ Gamma(alpha,1)` and `b[i,j] ~ Gamma(beta,1)` with the other two edge families gauge-fixed to 1.
- Weights: undirected within the intended Aztec-diamond correspondence; the code stores the reduced `a/b` environment, not site-directed transition weights.
- Transition probabilities: creation probabilities are `b/(a+b)` in the matrix encoding, or `a/(a+b)` for the alternate orientation in the legacy `creation_probabilities` routine.
- Stopping rule: not a walk; the sample is produced by domino shuffling to order `L`.
- Loop-erasure: not applicable.
- Winding convention: not applicable.
- Spanning-tree measure: not implemented here.
- Dimer construction: deletion, sliding, creation on each shuffle level.
- Height convention: staggered face table of size `(2L+1) × (L+1)`; crossing a dimer edge adds `+3`, otherwise `-1`.
- Spatial increment definition: for a symmetric pair of central-row faces, `H(right) - H(left)`.
- Replica/shared-environment design: single-sample runs use one environment and one tiling; double-dimer and spatial runs use one environment and two conditionally independent tilings.
- Exact variance identities: `Var(H1-H2)/2` and `Cov(H1,H2)` as above.

### 4.2 Uniform Aztec-diamond control

- Domain and boundary: same as Gamma model.
- Environment: no disorder; creation choices are fair coins.
- Transition probabilities: `1/2` for each potential creation orientation.
- Replica design: two independent tilings, no shared environment.
- Purpose: control for ordinary logarithmic growth in the spatial experiment.

### 4.3 Archived loop-erased random walk work

- Domain: planar LERW on `[-L,L]^2`.
- Boundary: stop on the boundary of the square.
- Environment: `site_iid` and `temporal_iid` weights.
- Replicas: one-walk and paired shared-environment variants.
- This archive is historical and not part of the active Aztec code path.

## 5. Repository Architecture

- `aztec/src/AztecDiamond.jl`: authoritative simulation and geometry kernel.
- `aztec/scripts/`: resumable campaign runners, strict mergers, bootstrap analyses, and SVG plots.
- `aztec/test/`: unit tests and end-to-end smoke workflow.
- `aztec/configs/`: frozen schedules for pilot, benchmark, production, continuation, and spatial campaigns.
- `aztec/data/`: retained raw results and metadata.
- `aztec/results/`: retained aggregate tables, fit reports, and figures.
- `aztec/output/`: ignored scratch and historical generated output; not authoritative.
- `old/`: archive of earlier research, including `old/random_walk/`.

Authoritative current code is the `aztec/` package, not `output/`.

## 6. Code-Flow Explanation

The active path is:

1. sample or supply weights;
2. reduce the Gamma environment and/or pre-draw creation coins;
3. replay domino shuffling;
4. read the height observable;
5. merge batches strictly;
6. bootstrap within order;
7. fit log, squared-log, and descriptive power curves;
8. render dependency-free SVG figures.

Key functions:

- `random_uniform_weights`, `random_gamma_weights`, `gamma_disordered_weights`
- `creation_probabilities`, `gamma_disordered_probabilities`, `gamma_disordered_creation_choices`
- `sample_tiling`, `sample_tiling_from_choices`
- `sample_gamma_center_height`, `sample_gamma_center_height_pair`
- `sample_gamma_spatial_increment_pair`, `sample_uniform_spatial_increment_pair`
- `height_function`, `face_height`, `center_height`, `symmetric_height_increment`
- `validate_tiling`, `write_table`, `write_svg`

## 7. Completed Campaigns

The retained campaigns and their primary result locations are:

| campaign | model | status | retained output |
|---|---|---|---|
| `gamma_height_*` production chain | Gamma center height | complete | `aztec/data/height/center_height_samples.csv`, `aztec/results/height/` |
| `double_dimer_campaign` | Gamma paired center heights | complete | `aztec/data/double_dimer/pairs.csv`, `aztec/results/double_dimer/` |
| `spatial_publication_gamma` / `spatial_publication_uniform` | paired spatial increments | complete | `aztec/results/spatial/` |
| `gamma_height_smoke`, `double_dimer_smoke`, `spatial_smoke` | smoke validation | complete | ignored scratch output only |
| `*_benchmark.csv` schedules | timing / memory checks | attempted or historical | ignored scratch output only |
| `old/random_walk` | archived LERW research | historical | `old/random_walk/results/` |

Verified counts from files:

- `aztec/data/height/center_height_samples.csv`: 35,537 lines including header, so 35,536 observations.
- `aztec/data/double_dimer/pairs.csv`: 28,305 lines including header, so 28,304 pairs.
- `aztec/results/spatial/spatial_summary.csv`: 73 lines, spanning the four separations for Gamma and uniform.
- `aztec/results/spatial/spatial_model_comparison.csv`: 25 lines.
- `aztec/results/spatial/spatial_pooled_model_comparison.csv`: 13 lines.
- `aztec/results/spatial/spatial_weighted_model_comparison.csv`: 25 lines.
- `aztec/results/spatial/spatial_cutoff_sensitivity.csv`: 41 lines.

## 8. Current Numerical Results

### Height

- Fit range: `L >= 24`.
- `delta_bic_log_minus_log2 = 4.434344933649989`.
- Bootstrap interval for `delta_bic`: `[-1.6589750057309498, 6.478611886808162]`.
- Effective power exponent: `p = 0.9670812264264157` with 95% interval `[0.8936958188166406, 1.0319325015459806]`.
- Interpretation: affine log and squared-log curves remain unresolved.

### Double-dimer difference

- Fit range: `L >= 24`.
- `delta_bic_log_minus_log2 = -1.928473762332553`.
- Bootstrap interval for `delta_bic`: `[-2.8904248452366255, -0.27728384490180413]`.
- Effective power exponent: `p = 0.813941225618811` with interval `[0.6996408255977046, 0.9020862661131736]`.
- Interpretation: ordinary log is preferred over the simulated range.

### Double-dimer components

- Conditional tiling component: `delta_bic = -1.9284737623325512`, `p = 0.8139412256188102`.
- Disorder covariance component: `delta_bic = 2.0252417266556346`, `p = 1.3928862747999537`.
- The disorder signal is positive but still finite-size noisy.

### Spatial Gamma disorder

- Pooled unweighted disorder covariance: `c = 0.5843`, 95% interval `[0.1957, 0.9520]`.
- Pooled weighted disorder covariance: `c = 0.4948`, 95% interval `[0.1932, 0.7977]`.
- Uniform control pooled marginal coefficient: `c = 0.00882`, interval `[-0.2623, 0.2835]`.
- Primary separation table shows the positive Gamma signal is strongest at `1/32` and `1/8`; `1/4` is noisier.

## 9. Statistical Methods

- Within-size bootstrap resampling.
- Bootstrap unit for single height: one observation at fixed order.
- Bootstrap unit for double-dimer: the paired row `(H1,H2,H1-H2)` stays intact.
- Bootstrap unit for spatial analysis: one environment, with all four separations kept aligned by environment.
- Bootstrap repetitions: 10,000 in retained analyses.
- Seed conventions: `20260802` for center-height and double-dimer, `20260805` for spatial analysis, `20260803/04` for the publication campaign runners.
- Affine fits: ordinary least squares on size-level variance/covariance estimates.
- Squared-log fits: `a + b log L + c (log L)^2` or the pooled shared-`c` version for spatial data.
- Power fits: log-log regression of the variance or covariance against `log(log L)`.
- BIC sign convention: `BIC(log) - BIC(log^2)`; positive favors squared-log.
- Holdout prediction: used in spatial analysis; the largest two orders are held out by default.
- Fit-window sensitivity: explicitly reported in `aztec/results/spatial/spatial_cutoff_sensitivity.csv`.

Implementation matches the reported method: the scripts bootstrap within order, recompute size-level statistics first, and fit on those summary values, not on raw heights.

## 10. Validation Status

Passing checks:

- `Pkg.test()` passed with 96 + 1 + 34 + 2 + 1 + 144 + 3 + 5 + 3 + 3 + 2 + 19 + 4 + 1 + 17 + 6 assertions grouped across the test sets.
- `aztec/test/smoke_workflows.sh` passed end to end.

What is validated:

- recurrence values for Gamma reduction;
- exact agreement between the pre-drawn-coin and probability-based Gamma creation paths;
- geometric tiling validation at small orders;
- center height consistency with the full height table;
- double-dimer pairing and finite-sample variance identity;
- direct spatial increment calculation;
- resumable CLI behavior and seed-consistency checks.

What remains unverified:

- the large production campaigns are not rerun here;
- every file under `aztec/output/` was not individually validated;
- the square-grid extension requested by Professor Chhita is not yet implemented.

## 11. Performance Status

- Threading: campaign runners use `Base.Threads.@threads :dynamic` over samples, with deterministic per-sample seeds.
- Environment-level parallelism: the Gamma environment is shared within one sample and reused across paired tilings or paired spatial increments.
- Main memory strategy: Gamma production uses two rolling `L × L` buffers plus packed `BitMatrix` creation tables.
- Major bottlenecks: domino shuffling is cubic in `L`; repeated height-path evaluation is linear per queried face; result-writing and validation are minor.
- Growth with `L`: runtime is `O(L^3)` for a sample; stored creation bits are `O(L^3)` bits; live weights and the final tiling are `O(L^2)`.
- Slurm support: no explicit Slurm module or job script exists yet; only shell wrappers and Julia CLIs are present.
- Durham Hamilton readiness: not yet production-ready without a Slurm wrapper, environment module notes, and a square-grid implementation plan.

## 12. Known Limitations

- The project is finite-size numerical evidence, not proof.
- The publication spatial result is a pooled fit over four separations; the fractions are not independent experiments.
- `aztec/output/` contains many historical and scratch runs that are not authoritative.
- The retained data are large enough to be useful but should not be blindly copied into every upload package.
- The archive under `old/random_walk/` is historical and not part of the current Aztec pipeline.

## 13. Unresolved Mathematical Questions

1. Whether the square-grid extension should use shared undirected edge conductances or directed site-dependent outgoing weights.
2. What the exact square-grid Temperley correspondence should be for the requested paired spatial analysis.
3. What the correct square-grid dimer height convention is.
4. How to align spatial location and parity conventions between the Aztec and square-grid geometries.
5. Whether the current finite-size squared-log signal is stable under additional high-order samples.

## 14. Professor Chhita’s Requested Next Stage

The requested next stage is:

1. Organise and document the existing Aztec code and analysis.
2. Implement the same paired spatial-height analysis for the square-grid dimer/spanning-tree model.
3. Test several random-weight choices in both geometries.
4. Keep the code reproducible and manuscript-ready.

## 15. Exact Recommended Implementation Sequence

1. Freeze the current Aztec behavior with tests and a short architecture note.
2. Formalise the square-grid mathematical model in a dedicated document before coding.
3. Decide the weight semantics: directed site weights versus shared undirected conductances.
4. Implement a minimal square-grid sampler and one exact small-case validation.
5. Add paired square-grid spatial increment sampling with deterministic seeds.
6. Add the square-grid equivalent of the variance-identity checks.
7. Run a tiny pilot on a few sizes and one or two weight choices.
8. Analyse the pilot with the same bootstrap and fit machinery.
9. Only then prepare a Hamilton/Slurm campaign layout.
10. Scale to production after the pilot reproduces the existing Aztec workflow style.

## 16. First Files Another Assistant Should Inspect

Use `chatgpt_handoff/FIRST_FILES_TO_READ.md` as the ordered checklist. The most important files are:

- `aztec/src/AztecDiamond.jl`
- `aztec/scripts/run_spatial_campaign.jl`
- `aztec/scripts/analyze_spatial_campaign.jl`
- `aztec/scripts/run_double_dimer_campaign.jl`
- `aztec/scripts/analyze_double_dimer_campaign.jl`
- `aztec/scripts/run_height_campaign.jl`
- `aztec/test/runtests.jl`
- `aztec/README.md`
- `aztec/results/spatial/README.md`
- `aztec/data/README.md`

## 17. Commands for Setup, Tests, Smoke Runs, Analysis, and Plotting

```bash
julia --project=aztec -e 'using Pkg; Pkg.test()'

sh aztec/test/smoke_workflows.sh

JULIA_NUM_THREADS=4 julia --project=aztec aztec/scripts/run_height_campaign.jl \
  --config aztec/configs/gamma_height_smoke.csv \
  --output-dir aztec/output/height_smoke

julia --project=aztec aztec/scripts/merge_height_batches.jl \
  --inputs aztec/output/height_smoke \
  --output aztec/output/height_smoke.csv

julia --project=aztec aztec/scripts/analyze_height_campaign.jl \
  --results-dir aztec/data/height/center_height_samples.csv \
  --output-dir aztec/output/height_analysis \
  --bootstrap-reps 10000 --bootstrap-seed 20260802 --min-order 24

julia --project=aztec aztec/scripts/plot_height_campaign.jl \
  --analysis-dir aztec/output/height_analysis \
  --output-dir aztec/output/height_analysis

JULIA_NUM_THREADS=4 julia --project=aztec aztec/scripts/run_double_dimer_campaign.jl \
  --config aztec/configs/double_dimer_smoke.csv \
  --output-dir aztec/output/double_smoke

julia --project=aztec aztec/scripts/analyze_double_dimer_campaign.jl \
  --paired-results aztec/data/double_dimer/pairs.csv \
  --single-results aztec/data/height/center_height_samples.csv \
  --output-dir aztec/output/double_analysis \
  --bootstrap-reps 10000 --bootstrap-seed 20260802 --min-order 24

JULIA_NUM_THREADS=8 sh aztec/scripts/run_spatial_publication_campaign.sh
```

## 18. Commands That Must Not Be Run

- `git reset --hard`
- `git clean -fdx`
- destructive checkout commands that discard user work
- rerunning the large production campaigns without an explicit need
- deleting or overwriting `aztec/data/`, `aztec/results/`, `old/`, or `aztec/output/`

## 19. Honest Verification Notes

- I verified the active code path, the retained analysis tables, the current repo state, and the test suite.
- I did not rerun the production campaigns.
- I did not inspect every file under the ignored `aztec/output/` tree individually.
- I found no `INFORMATION*`, `research_project_complete_summary.txt`, or `CITATION.cff` files in the workspace.
- The square-grid extension requested for the next phase is not present yet.