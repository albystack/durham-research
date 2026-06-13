# LERW Random Environment Research

This folder contains the working project on **loop-erased random walks (LERW), winding variance, and fixed random environments**, based on the email discussion with Prof. Sunil Chhita.

The project began from the question:

> If the random walk used to generate a spanning-tree branch is placed in a fixed site-dependent random environment, does the winding variance grow like `C log L`, as in the symmetric case, or more like `C(log L)^2`?

The current stage of the project is **numerical exploration and paper preparation**. The code, data, figures, analysis notes, and first paper draft are all in this folder.

---

## 1. Project Summary

We simulate random walks on the square lattice box

```text
B_L = [-L, L]^2 intersect Z^2
```

starting at the origin. Each walk is stopped when it first hits the boundary. We then chronologically erase loops from the path to obtain a loop-erased random walk.

The main observable is the winding:

```text
winding = number of left turns - number of right turns
```

For the ordinary symmetric case, where each step has probability `1/4`, the expected behaviour is:

```text
Var(W_L) ~ C log L
```

The project tests random-environment variants and compares winding variance against both:

```text
C log L
C(log L)^2
```

---

## 2. Current Main Finding

The most important finding so far is a modelling issue.

Prof. Chhita suggested using edge weights of the form:

```text
w_N = 1
w_E = 1
w_S = u
w_W = v
```

where `u` and `v` are positive mean-one random variables, then normalising to get transition probabilities.

Numerically, this **two-weight model develops a strong north-east effective drift** after normalisation, especially for moderate and strong disorder.

For example, in the focused run at `L = 256`:

| model | k | mean hit x/L | mean hit y/L |
|---|---:|---:|---:|
| gamma | 20 | 0.589 | 0.591 |
| gamma | 1 | 0.894 | 0.889 |
| gamma | 0.5 | 0.916 | 0.918 |

This means the walk is usually exiting near the north-east corner. In that regime the winding variance is suppressed, which is consistent with the drifted-walk behaviour Prof. Chhita warned about.

As a diagnostic comparison, the project also includes a balanced model:

```text
w_N, w_E, w_S, w_W iid Gamma(k, 1/k)
```

This model is called `gamma4` in the code. It stays centred in the boundary-hit diagnostics. Up to `L = 256`, it shows winding variance growth broadly similar to the symmetric baseline, with no convincing evidence yet for `(log L)^2` growth.

The next scientific step is to ask Prof. Chhita whether the intended model should be the two-weight model despite the drift, or whether a more balanced/no-drift random environment is intended.

---

## 3. Folder Map

### Top-Level Code

- `lerw_random_environment.py`  
  Core simulator. Runs random walks, loop-erases paths, computes winding, and writes CSV output.

- `run_research_batch.py`  
  Batch runner. Runs many model/size/disorder configurations and supports per-configuration checkpointing.

- `analyse_winding_results.py`  
  Reads CSV results, groups by model/size/disorder, computes summary statistics, and fits variance against `log L` and `(log L)^2`.

- `make_figures.py`  
  Generates SVG figures and summary tables using only the Python standard library.

- `test_lerw_random_environment.py`  
  Unit tests for loop erasure, winding, probability normalisation, fixed environments, and small trial behaviour.

### Main Results

- `results/research_batch_main.csv`  
  Broad exploratory sweep with `19,500` trials.

- `results/research_batch_focused.csv`  
  Focused key-case run with `24,500` trials. Use this first for the paper.

- `results/research_batch_focused_*.csv`  
  Per-configuration checkpoints from the focused run.

### Analysis Outputs

- `analysis_output/`  
  Figures, tables, and analysis dossier for the broad exploratory sweep.

- `analysis_focused/`  
  Figures, tables, and analysis dossier for the focused run. This is the main analysis folder.

- `analysis_focused/results_interpretation.md`  
  Human-readable interpretation of the focused numerical results.

### Paper Materials

- `paper/paper_outline.md`  
  Proposed paper structure and defensible claims.

- `paper/draft.md`  
  First manuscript draft.

- `paper/supervisor_update_draft.md`  
  Draft message asking Prof. Chhita to clarify the intended random-environment model.

### Original Notes

The original project-planning files are stored in `old_data/`:

- `old_data/outlook.txt`  
  Email thread with Prof. Chhita.

- `old_data/ust-lerw.tex`  
  Long LaTeX primer on UST, LERW, Wilson's algorithm, winding, and random environments.

- `old_data/research_plan.md`  
  Initial research plan.

- `old_data/project_brief.md`  
  Concise brief extracted from the email thread.

- `old_data/background_outline.md`  
  Background/write-up outline.

---

## 4. Simulator Models

The simulator supports four models.

### `symmetric`

The ordinary simple symmetric random walk:

```text
p_N = p_E = p_S = p_W = 1/4
```

This is the calibration case. Its winding variance should be consistent with `C log L`.

### `gamma`

The two-random-weight model based on the email guidance:

```text
w_N = 1
w_E = 1
w_S = u
w_W = v
```

where:

```text
u, v ~ Gamma(k, 1/k)
```

The weights are normalised to get probabilities:

```text
p_direction = w_direction / (w_N + w_E + w_S + w_W)
```

Important: this model has shown strong north-east drift after normalisation.

### `exponential`

Same structure as `gamma`, but with:

```text
u, v ~ Exponential(mean = 1)
```

This is equivalent to the `gamma` model with `k = 1`, but kept as a separate option for clarity.

### `gamma4`

Balanced four-weight diagnostic model:

```text
w_N, w_E, w_S, w_W iid Gamma(k, 1/k)
```

This model is not the original email model, but it is useful because it stays approximately centred and helps separate drift effects from random-environment effects.

---

## 5. How To Run Tests

From this folder:

```bash
python3 -m unittest -v
```

Current status:

```text
8/8 tests passing
```

The tests check:

- chronological loop erasure
- winding sign convention
- symmetric transition probabilities
- positive and normalised Gamma probabilities
- fixed-in-time random environments
- small boundary-hitting trials

---

## 6. How To Reproduce The Focused Batch

The focused batch is the main dataset for writing the paper.

Command:

```bash
python3 run_research_batch.py \
  --sizes 16 32 64 128 256 \
  --samples 700 \
  --seed 20260613 \
  --models symmetric gamma gamma4 \
  --k-values 20 1 0.5 \
  --checkpoint-per-config \
  --output results/research_batch_focused.csv
```

This produces:

```text
results/research_batch_focused.csv
results/research_batch_focused_symmetric.csv
results/research_batch_focused_gamma_k20.csv
results/research_batch_focused_gamma_k1.csv
results/research_batch_focused_gamma_k0.5.csv
results/research_batch_focused_gamma4_k20.csv
results/research_batch_focused_gamma4_k1.csv
results/research_batch_focused_gamma4_k0.5.csv
```

Total focused-run size:

```text
24,500 trials
```

---

## 7. How To Regenerate Figures And Tables

To regenerate the focused analysis:

```bash
python3 make_figures.py \
  results/research_batch_focused.csv \
  --output-dir analysis_focused
```

This writes:

```text
analysis_focused/analysis_dossier.md
analysis_focused/results_interpretation.md
analysis_focused/figures/*.svg
analysis_focused/tables/group_stats.csv
analysis_focused/tables/scaling_fits.csv
```

The figures are SVG files, which can be opened in a browser or included in markdown/LaTeX workflows.

On macOS, PNG previews can be generated with:

```bash
mkdir -p analysis_focused/png
qlmanage -t -s 1400 -o analysis_focused/png analysis_focused/figures/*.svg
```

---

## 8. Main Figures

The key paper figures are in `analysis_focused/figures`.

Recommended first set:

- `variance_vs_logL_selected.svg`  
  Winding variance against `log L` for the symmetric case, two-weight model, and balanced diagnostic cases.

- `variance_vs_logL2_selected.svg`  
  Same data plotted against `(log L)^2`.

- `mean_hit_x_over_L.svg`  
  Drift diagnostic using mean boundary hit `x/L`.

- `mean_hit_y_over_L.svg`  
  Drift diagnostic using mean boundary hit `y/L`.

- `boundary_hits_L128.svg`  
  Boundary hit scatter plot showing the north-east drift in the two-weight model.

- `gamma_variance_by_k_logL.svg`  
  Winding variance by disorder level for the two-weight model.

- `gamma4_variance_by_k_logL.svg`  
  Winding variance by disorder level for the balanced diagnostic model.

- `raw_walk_length.svg`  
  Mean raw walk lengths. Useful for showing that drifted cases exit much faster.

- `lerw_path_length.svg`  
  Mean loop-erased path lengths.

- `sample_lerw_paths.svg`  
  Example loop-erased paths.

- `environment_vector_fields.svg`  
  Visual diagnostic of local mean-step directions.

---

## 9. Focused Numerical Results

### Symmetric Baseline

| L | Var(W) | mean raw steps | mean LERW steps |
|---:|---:|---:|---:|
| 16 | 1.894 | 296.9 | 40.3 |
| 32 | 2.353 | 1239.1 | 96.2 |
| 64 | 3.148 | 4861.2 | 231.8 |
| 128 | 3.572 | 20045.7 | 559.1 |
| 256 | 4.300 | 77271.2 | 1329.5 |

Fit against `log L`:

```text
Var(W) = 0.870 log L - 0.565
R^2 = 0.992
```

Fit against `(log L)^2`:

```text
Var(W) = 0.104 (log L)^2 + 1.158
R^2 = 0.987
```

Interpretation: the symmetric baseline is consistent with the expected logarithmic behaviour.

### Two-Weight Gamma Model

At `L = 256`:

| k | Var(W) | mean hit x/L | mean hit y/L | mean raw steps |
|---:|---:|---:|---:|---:|
| 20 | 3.735 | 0.589 | 0.591 | 48969.1 |
| 1 | 2.099 | 0.894 | 0.889 | 4731.9 |
| 0.5 | 1.548 | 0.916 | 0.918 | 2891.9 |

Interpretation:

- The two-weight model is strongly drifted for `k = 1` and `k = 0.5`.
- Drift makes the raw walks much shorter.
- Winding variance is suppressed.
- These cases should not be interpreted as evidence for the intended centred random-environment scaling question.

### Balanced Gamma4 Diagnostic

At `L = 256`:

| k | Var(W) | mean hit x/L | mean hit y/L | mean raw steps |
|---:|---:|---:|---:|---:|
| 20 | 4.104 | 0.017 | -0.027 | 79695.7 |
| 1 | 4.001 | 0.028 | 0.026 | 96802.5 |
| 0.5 | 4.184 | 0.017 | -0.003 | 149513.0 |

Interpretation:

- The balanced model remains centred.
- Winding variance grows with `L`.
- Over `L <= 256`, it looks closer to logarithmic growth than to clear `(log L)^2` growth.
- Larger sizes and/or more samples would be needed for stronger asymptotic claims.

---

## 10. What We Can Safely Claim

Safe claims:

- The simulation reproduces the expected qualitative symmetric baseline.
- The two-weight model develops strong effective drift after normalisation.
- Strong drift suppresses winding variance growth.
- The balanced four-weight diagnostic remains centred.
- The balanced diagnostic does not show convincing `(log L)^2` growth up to `L = 256`.

Claims to avoid:

- Do not claim the conjecture is false.
- Do not claim the random-environment problem is solved.
- Do not describe the two-weight model as drift-free just because the weights have mean one.
- Do not present `gamma4` as the official model unless Prof. Chhita agrees.

---

## 11. Paper-Writing Status

Current paper materials:

```text
paper/paper_outline.md
paper/draft.md
paper/supervisor_update_draft.md
```

The draft already contains:

- abstract
- introduction
- model description
- numerical method
- baseline results
- two-weight model results
- balanced diagnostic results
- discussion
- conclusion

Before polishing the paper, the main open question is the model clarification:

> Should the intended random-environment model be the two-weight model despite the observed drift, or should we use a more locally balanced/no-drift environment?

---

## 12. Recommended Next Steps

1. Send `paper/supervisor_update_draft.md` to Prof. Chhita.
2. Ask whether the two-weight model is still intended.
3. If yes, frame the paper around drift diagnostics and finite-variance behaviour.
4. If no, switch the main model to a balanced/no-drift environment.
5. After clarification, run a larger final batch with the agreed model.
6. Add confidence intervals or bootstrap error bars to the plots.
7. Convert `paper/draft.md` into a polished LaTeX report.

---

## 13. Reproducibility Notes

The focused run used:

```text
seed = 20260613
sizes = 16, 32, 64, 128, 256
samples = 700 per cell
models = symmetric, gamma, gamma4
k-values = 20, 1, 0.5
```

The code stores each trial with:

- `L`
- model
- `k`
- sample index
- seed
- raw walk length
- loop-erased path length
- winding
- boundary hit location
- number of sampled environment sites

The key grouped table is:

```text
analysis_focused/tables/group_stats.csv
```

The key scaling-fit table is:

```text
analysis_focused/tables/scaling_fits.csv
```

