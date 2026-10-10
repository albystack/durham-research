# Corrected Aztec refinement, 3 October 2026

**The corrected Aztec experiment has positive finite-window curvature.** The
primary squared-log coefficient is **0.3568**, with a 95% environment-bootstrap
interval **[0.0723, 0.6327]**. This is evidence for a positive quadratic-log
component over the tested orders. It does not establish an asymptotic law or
an Aztec-versus-square boundary effect.

## Completed experiment and checks

Hamilton array **18897605** completed all six tasks with exit code 0:0 and empty
error logs. Orders 16, 32, 48, 64, 96 and 128 each have 2,048 environments and
two replica tilings: **12,288 pairs / 24,576 tilings**, with no missing cases.
Task elapsed times were 56 s, 43 s, 65 s, 82 s, 204 s and 347 s respectively.
The sampler validated every tiling during generation. Raw weights, tilings and
case certificates remain on Hamilton; compact tables and manifests were
downloaded. This follow-up did not download and revalidate every raw tiling.

The compact archive hash and all six measurement-table hashes passed. Checks
also passed for complete probe coverage, physical distances, model and height
conventions, common source versions, unique seeds and all 24,576 seed
derivations. An independent NumPy calculation reproduced all 36 thermal,
disorder and total moment estimates. The analysis code matches the prepared
workflow; no outcome-dependent windows or sample-count changes were made.

The model uses independent Gamma shape 0.2 and 0.25 families with scale one,
the remaining weights one, and the corrected orientation. The sampler is
`gamma-interval-coins-v1`: creation coins are rational for represented reduced
weights, while reduction remains Float64. The prior three-environment
high-precision audit is numerical evidence, not a global exact-arithmetic
proof. Both physical diagonal probes use unit heights, with relative
separations 1/8 and 1/4 and actual distance sqrt(2) times the table separation.

## Frozen analysis

Fit disorder covariance as `D_p(d) = a_p + b_p log(d) + c [log(d)]^2`, with
probe-specific intercepts and slopes and a common c. The fit is ordinary least
squares on the covariance estimates. Uncertainty comes from 10,000 whole-
environment bootstrap resamples, preserving both replicas and probes together;
sizes are independently resampled. Bootstrap seed: 2026100103. Sampling master
seed: 2026100102.

| Order window | c | 95% bootstrap interval |
|---|---:|---:|
| 16–128 (primary) | 0.3568 | [0.0723, 0.6327] |
| 32–128 | 0.7362 | [0.0049, 1.4770] |
| 48–128 | 1.9256 | [0.1376, 3.7377] |

![Corrected covariance curves and preselected cutoff estimates](../figures/refinement_20261003.png)

All intervals lie above zero, but these are correlated, pointwise intervals,
not three independent confirmations or a simultaneous confidence statement.
In particular, the 32–128 interval only narrowly excludes zero. Estimates drift
upward and uncertainty expands when smaller orders are removed, so a stable
asymptotic coefficient is not established. The primary interval half-width is
0.2802, slightly larger than the planning target 0.25. The fixed sample count
was respected; no automatic extension was triggered.

## Implication for Sunil's question

The corrected positive Aztec signal now has substantially more precise
evidence than the earlier 384-pair pilot. The unresolved part is a feasible,
validated square sampler under the same Gamma law at useful comparison sizes.
The small-width integer sampler completes the retained original environments,
but is currently limited to L<=12. The matched CFTP pilot had 0/80 random pairs
complete at its cap. Neither that censoring nor the different all-edge Gamma(1)
September data establish absence of curvature for the matched square model.

The next decision is the square sampling route and a fixed feasibility budget
for wider sizes, followed by physically matched windows and a direct domain
contrast. The large matched production campaign remains on hold. The updated
supervisor draft reports this result and asks to discuss the square route;
it has not been sent. Human review of conventions, sampler reasoning and the
interpretation remains needed, alongside historical AI-attribution completion.

## Reproduction and records

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=square_glauber_python/src .venv/bin/python square_glauber_python/scripts/analyse_robust_refinement.py --input aztec/output/robust_refinement_20261003/refinement --output docs/validation/refinement_20261003
MPLCONFIGDIR=tmp/matplotlib .venv/bin/python docs/validation/refinement_20261003/plot_results.py
```

The analysis, intervals, verification and scheduler evidence are in
`docs/validation/refinement_20261003/`. The downloaded compact archive is
`tmp/hamilton/refinement-compact-20261003.tar.gz`, SHA256
`227ef3be6943a8fbaae298136f894abcd5afc096469422625ac2520d8b3c13ab`.
The full Hamilton root is
`/nobackup/fvkl37/numerical_refinement_20261001/refinement`.
