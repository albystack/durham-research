# LERW Random Environment Research

This is a simplified numerical research repository for loop-erased random walks
(LERW), winding variance, and fixed random environments.

The main question is whether a fixed site-dependent random environment changes
the winding variance from the symmetric baseline

```text
Var(W_L) ~ C log L
```

toward a larger finite-size behaviour such as

```text
Var(W_L) ~ C (log L)^2
```

The project explainer is kept here:

```text
research_project_explainer.tex
research_project_explainer.pdf
```

---

## The Simple File Layout

The code is intentionally small:

- `main.py`
  - Run this file for simulations and analysis.
  - Contains the simulator, loop erasure, winding calculation, CSV analysis,
    and simple SVG figure generation.

- `experiments.py`
  - Edit this file when you want to add or change an experiment.
  - Contains named setups such as `smoke`, `focused_reproduction`, and
    `balanced_large`.

- `test_main.py`
  - One test file for the project.
  - Checks the main modelling assumptions and a small trial run.

- `README.md`
  - This guide.

- `research_project_explainer.tex` / `research_project_explainer.pdf`
  - The research write-up.

Generated `results/` and `analysis_*` folders can be recreated by `main.py`.

---

## What To Run

List the available experiment setups:

```bash
python3 main.py list
```

Run a quick smoke test:

```bash
python3 main.py run smoke
```

Run a smoke test and immediately analyse it:

```bash
python3 main.py run smoke --analyse
```

Analyse an existing CSV:

```bash
python3 main.py analyse results/research_batch_focused.csv --output-dir analysis_focused
```

Run the test suite:

```bash
python3 -m unittest -v
```

---

## How To Add Or Change A Run

Open `experiments.py` and edit the `EXPERIMENTS` dictionary.

Each experiment has:

- `description`: a short label printed by `python3 main.py list`
- `sizes`: the box half-widths `L`
- `samples`: trials per model and size
- `seed`: the base random seed
- `model_configs`: strings like `symmetric` or `gamma_balanced:1`
- `output`: combined CSV path
- `checkpoint_per_config`: whether to write one CSV per model config
- `max_steps`: optional walk cutoff; normally leave as `None`

Then run:

```bash
python3 main.py run your_experiment_name --analyse
```

---

## Model Summary

Walks run on

```text
B_L = [-L, L]^2 intersect Z^2
```

starting at the origin and stopping at the first boundary hit. The raw path is
chronologically loop-erased. The winding is

```text
number of left turns - number of right turns
```

The key model families are:

- `symmetric`
  - Baseline simple random walk with equal probabilities.

- `gamma` and `exponential`
  - Original two-weight drift-diagnostic model.
  - Useful because it shows how normalisation can create effective drift.

- `gamma4`, `lognormal4`, `pareto4`, `uniform4`
  - Four independent direction weights at each site.
  - Balanced in distribution, but not necessarily drift-free at a given site.

- `gamma_balanced`, `lognormal_balanced`, `pareto_balanced`, `uniform_balanced`
  - Main locally balanced models.
  - At each site, sample a vertical weight and a horizontal weight:

```text
w_N = w_S = vertical_weight
w_E = w_W = horizontal_weight
```

These exact opposite-pair models are the main next step for the research.

---

## Current Interpretation

The earlier focused results found that the original two-weight Gamma model can
develop strong north-east drift after normalisation. That suppresses winding
variance and means those runs should be treated as drift diagnostics, not as the
main centred random-environment test.

The next useful simulations should focus on exact locally balanced models and
larger sizes. Avoid claiming that the asymptotic scaling question is settled by
the current finite-size data.
