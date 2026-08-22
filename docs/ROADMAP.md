# Research roadmap

This file records the next scientifically useful work. Completed numerical
results and historical pilots belong in [`RESULTS.md`](RESULTS.md).

## Current conclusion

A finite bipartite Kasteleyn calculation has replayed all 960 frozen Gamma
environments from the direct square-grid campaign. It reproduces exact tiny
enumeration and confirms the central-height component estimates without an
MCMC mixing assumption. No component has a robust positive
\((\log L)^2\) coefficient over \(L=2\)–20.

The next objective is therefore a **like-for-like spatial-height experiment**,
not more central-height MCMC production.

## Completed — Determinantal central-height validation

- The dense Kasteleyn orientation and height-path formula match complete
  weighted `L=1,2` enumeration.
- The environment namespace matches the MCMC production seeds.
- All 960 Gamma environments have deterministic conditional means and
  variances with recorded numerical diagnostics.
- The MCMC and determinantal aggregate conclusions agree.

Retained data and results are under
`aztec/data/glauber_square_grid_kasteleyn_20260822/` and
`aztec/results/glauber_square_grid_kasteleyn_20260822/`.

## P1 — Validate spatial height increments

Use `height_difference_moments_kasteleyn` for the same relative separations as
the Aztec paired experiment:

- \(r/L = 1/32,1/16,1/8,1/4\) where the finite grid permits them;
- conditional variance from Kasteleyn covariances;
- disorder variance from environment-to-environment conditional means;
- an all-one negative control;
- several paths with the same endpoints to test height path-independence.

Completion gate: every spatial mean and covariance agrees with exhaustive
enumeration on tiny grids and is path-independent within numerical tolerance.

## P2 — Dense-method size and conditioning pilot

Benchmark geometrically spaced sizes with the same Gamma law and the all-one
control. Record factorization time, memory, residuals, high-precision fallback
rate, and the number of independent environments achievable per size.

Do not choose production sizes until this benchmark identifies where dense
factorization ceases to be practical or numerically reliable.

## P3 — Sparse or nested-dissection implementation

If the dense pilot cannot reach a useful spatial scaling window, replace the
dense factorization with a sparse selected-inverse or nested-dissection method.
Validate it against the dense implementation on shared environments before
using it for new sizes.

## P4 — Fresh-seed spatial production

Only after P1–P3 pass, predeclare sizes, environment counts, fitting windows,
bootstrap seed, and numerical-failure policy. Analyze the fresh-seed campaign
independently before combining it with any pilot.

## P5 — Disorder-strength scan

After the primary Gamma spatial experiment is frozen, compare predeclared
mean-one Gamma shapes and the all-one control. Keep each law as a separate
estimand and report the same cutoff and weighting sensitivities.

## Optional MCMC cross-validation

Replica-flow instrumentation and a repaired temperature ladder remain useful
for observables that cannot be reduced to Kasteleyn correlations. They are no
longer prerequisites for the `L<=20` central-height conclusion.

## Cross-model comparison

Compare Aztec and square-grid spatial results through one common schema:

```text
model, geometry, L, disorder law, environment seed,
observable, conditional component, disorder component, numerical diagnostics
```

Use the same environment-blocked fitting implementation and report:

- coefficient sign and uncertainty;
- BIC sign convention;
- prediction error;
- lower-size cutoff sensitivity;
- unweighted and precision-weighted fits;
- conditional, disorder, and total components separately;
- negative controls.

The comparison must distinguish genuine geometry/model differences from
finite-size effects, boundary effects, and differences in how disorder is
represented.

## Release checklist

- deterministic replayable seeds;
- one config file per retained campaign;
- restart-safe jobs and atomic outputs;
- raw-data manifests and checksums;
- scripts that regenerate retained tables and figures;
- tests for mathematical invariants and optimized observables;
- no private correspondence, credentials, or machine-specific personal paths;
- concise public documentation of negative results and limitations.
