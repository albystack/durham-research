# Random-environment dimer height fluctuations

This repository studies height fluctuations in random tilings, dimers,
spanning trees, and loop-erased random walks. The central numerical question is
whether the disorder-induced component of the height covariance contains a
positive \((\log L)^2\) term.

The project is designed to test that possibility, not to force it. Negative
controls, exact small-volume checks, frozen analysis plans, and whole-environment
resampling are retained so that an ordinary-log or null result remains a valid
scientific outcome.

## Current focus

The active experiment is the direct square-grid random-bond dimer model in
[`square_glauber_python/`](square_glauber_python/PROJECT_GUIDE.md):

- every undirected grid edge receives one fixed positive weight;
- single-face random-scan heat-bath updates preserve the weighted Gibbs law;
- the exterior-referenced dimer height uses the documented \(\pm1/\pm3\)
  convention;
- monotone coupling from the past (CFTP) returns certified perfect samples;
- two independent perfect samples share each frozen environment;
- disorder covariance is measured from paired centre heights and spatial
  height increments.

The frozen Hamilton campaign covers \(L=8,12,16,20,24,28,32\), with 18,200
Gamma(shape=1) environments and 5,600 all-one controls. Its pre-specified
analysis is in
[`PAIRED_COVARIANCE_ANALYSIS_PLAN.md`](square_glauber_python/PAIRED_COVARIANCE_ANALYSIS_PLAN.md).
No scaling conclusion is recorded until every retained environment is either
certified or explicitly reported as censored and extended with the identical
random-map history.

## Quick start

Recommended local versions are Python 3.12 or 3.13 and Julia 1.10 or newer.

```bash
# Create the Python environment with tests and optional Numba acceleration.
make setup-python PYTHON=python3.13

# Run the Python mathematical and workflow tests.
make test-python

# Run exact tiny-grid validation.
make validate-python

# Run the established Julia/Aztec tests.
make test-julia
```

Equivalent commands without `make` are documented in the two project guides:

- [square-grid Python/CFTP guide](square_glauber_python/PROJECT_GUIDE.md)
- [Aztec/Julia guide](aztec/PROJECT_GUIDE.md)

## Statistical observable

For two conditionally independent replicas \(H_1,H_2\) in the same frozen
environment \(\omega\),

\[
\frac12\operatorname{Var}(H_1-H_2)
=\mathbb E_\omega[\operatorname{Var}(H\mid\omega)]
\]

is the connected or conditional component, while

\[
\operatorname{Cov}(H_1,H_2)
=\operatorname{Var}_\omega(\mathbb E[H\mid\omega])
\]

is the disorder-induced component in which super-rough \((\log L)^2\) growth
is being tested. The independent statistical unit is the environment; both
replicas and every observable from it remain together during bootstrap
resampling.

## Repository map

```text
.
├── README.md                  this entry point and the only maintained README
├── Makefile                   short setup, test, and validation commands
├── square_glauber_python/     active Python heat-bath and perfect sampler
├── aztec/                     established Julia tiling/tree implementations
├── results/                   reviewed, compact scientific results
│   └── aztec/                 retained Aztec and prior square-grid analyses
├── docs/                      overview, ledger, roadmap, and reproducibility
├── research_materials/        local correspondence, notes, papers, bibliography
├── hpc/                       earlier Julia Hamilton wrappers
└── archive/                   superseded experiments kept for provenance
```

Generated output remains outside Git:

- `aztec/output/` contains local Julia batches and scratch analysis;
- `square_glauber_python/outputs/` contains local Python diagnostics;
- large Hamilton campaign data live under the recorded `/nobackup` campaign
  root and are distilled into compact reviewed results only after audit.

## Documentation and results

- [documentation index](docs/INDEX.md)
- [research overview](docs/RESEARCH_OVERVIEW.md)
- [results ledger](docs/RESULTS.md)
- [current roadmap](docs/ROADMAP.md)
- [reproducibility contract](docs/REPRODUCIBILITY.md)
- [reviewed results index](results/INDEX.md)
- [Hamilton workflow](hpc/HAMILTON.md)
- [research-materials policy](research_materials/INDEX.md)

## Publication notes

Raw supervisor emails and third-party papers are private/local by default and
are ignored under `research_materials/`. Before making the repository public,
review names, email addresses, machine paths, licenses, and every newly staged
file. Add a software license and final citation metadata only after choosing
the intended terms.

All reported conclusions are finite-size numerical findings unless a cited
mathematical result states otherwise.
