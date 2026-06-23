# Random Dimer Height Proxy Experiments

This project simulates loop-erased random walks in fixed random edge-weight
environments.  The winding of the loop-erased path from the origin to the
boundary of `[-L,L]^2` is used as a tree-side proxy for the dimer height
fluctuation suggested by Temperley's bijection.

The numerical question is whether the variance grows more like
`C log L` or like the conjectural super-rough law `C (log L)^2`.

## Setup

Create a local virtual environment and install the plotting stack:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

## Run

List presets and models:

```bash
.venv/bin/python main.py list
```

Run a tiny check:

```bash
.venv/bin/python main.py smoke
```

Run a pilot or the main probe:

```bash
.venv/bin/python main.py run pilot --workers 4
.venv/bin/python main.py run super_rough_probe --workers 8
```

Run the long 512/1024-only continuation:

```bash
.venv/bin/python main.py run super_rough_large --workers 8
```

Override sizes, samples, or models:

```bash
.venv/bin/python main.py run pilot --sizes 16 32 64 --samples 100 \
  --models symmetric gamma_edges:1.0 lognormal_edges:0.7
```

## Files

```text
config.py       models, distributions, and experiment presets
simulation.py   random edge environment, random walk, loop erasure, winding
analysis.py     summaries, scaling fits, bootstrap intervals, plots
main.py         command-line runner
```

Each run writes to `results/<run_name>/`:

```text
raw_trials.csv
summary.csv                         includes variance standard errors
fits.csv                            includes bootstrap slope intervals
run_report.md
variance_scaling_all_models.png     includes error bars
variance_fit_*.png                  includes error bars and fitted curves
```

## Models

All random models assign a positive quenched weight to every undirected
nearest-neighbour edge used by the walk.  The old two-random-edge setup with
fixed north/east weights has been removed.

Available random families include Gamma, Exponential, Lognormal, Pareto,
Uniform, Beta, Weibull, inverse-Gamma, Bernoulli two-point, and triangular
weights.  See `config.py` to edit or add models.

## Notes

The output fits are finite-size diagnostics.  A model looking closer to
`C (log L)^2` in `fits.csv` is evidence to investigate further, not a proof of
the conjecture.
