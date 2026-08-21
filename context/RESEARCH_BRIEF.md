# Research Brief

## One-line objective

Use reproducible numerical experiments to determine **which fluctuation component in random-environment spanning-tree/dimer models exhibits super-rough \(C(\log L)^2\) growth**, and why the Aztec-diamond and square-grid experiments currently behave differently.

## Current scientific question

The project began with uniform spanning trees (USTs), loop-erased random walks (LERWs), and winding fluctuations on square grids. It evolved into a comparison between random-environment spanning-tree models and random-weight dimer models, with the main target being the conjectured/predicted \(C(\log L)^2\) contribution.

The key lesson from the experiments is that the observable matters. A single total variance or a double-dimer difference can mix or isolate the wrong component. The current focus is the decomposition of fluctuations into:

\[
\operatorname{Var}(H)
=
\mathbb E_{\omega}\!\left[\operatorname{Var}(H\mid\omega)\right]
+
\operatorname{Var}_{\omega}\!\left(\mathbb E[H\mid\omega]\right),
\]

where \(\omega\) is the random environment.

For two conditionally independent samples \(H_1,H_2\) in the **same** environment,

\[
\frac12\operatorname{Var}(H_1-H_2)
=
\mathbb E_{\omega}\!\left[\operatorname{Var}(H\mid\omega)\right],
\]

while

\[
\operatorname{Cov}(H_1,H_2)
=
\operatorname{Var}_{\omega}\!\left(\mathbb E[H\mid\omega]\right).
\]

For spatial height increments \(\Delta H_\alpha(r)\), the working estimands are

\[
T(r)=\frac12\operatorname{Var}\!\left(\Delta H_1(r)-\Delta H_2(r)\right)
\]

for the conditional/connected component, and

\[
D(r)=\operatorname{Cov}\!\left(\Delta H_1(r),\Delta H_2(r)\right)
\]

for the disorder-induced component.

The numerical evidence so far indicates that \(T(r)\) is compatible with ordinary logarithmic growth, while the Aztec-diamond \(D(r)\) can show a positive \((\log L)^2\) term.

## Models investigated

### 1. Square-grid UST / LERW

Initial implementation:
- finite \(n\times n\) square grids;
- ordinary and wired boundary conditions;
- simple random walks, LERW, Wilson's algorithm;
- branch winding = number of left turns minus number of right turns;
- seeded/reproducible runs and structural validation.

Random-environment variants:
- site-dependent transition weights fixed in time (quenched environment);
- locally balanced four-weight models;
- several disorder distributions, including Gamma, lognormal and Pareto families;
- one walk per independent environment for strict annealed sampling;
- paired walks in the same environment;
- later, weights refreshed at every step rather than fixed at a site.

### 2. Aztec-diamond random-weight dimers

Sampling is based on domino shuffling / weight reduction. The main Gamma-disordered setup in the correspondence uses independent Gamma weight families with parameters corresponding to the Duits–Van Peski example, with the remaining edge-weight families fixed to 1.

Observables evolved from:
1. central-face height;
2. two-tiling height difference;
3. paired spatial height increments at
   \[
   r\in\{L/32,L/16,L/8,L/4\},
   \]
   analysed through \(T(r)\) and \(D(r)\).

This third formulation produced the strongest super-rough signal.

### 3. Square-grid dimer / Temperley representation

The Aztec paired analysis was transferred to the square grid via the spanning-tree/Temperley correspondence. Large Hamilton-cluster runs used multiple disorder laws and reached \(L=6144\).

These runs did **not** show a robust positive \((\log L)^2\) contribution. The important caveat is that the simulated disorder was structured through the spanning-tree/Temperley construction rather than being the canonical model with independent random weights on every dimer edge.

## Current evidence

### Strongest positive evidence

Aztec-diamond paired spatial increments:
- disorder covariance \(D(r)\) shows positive finite-size quadratic-log curvature for the original Gamma law;
- a stronger-disorder Aztec law also shows a positive broad-window coefficient;
- conditional/connected component does not show the same effect;
- uniform/no-disorder controls do not manufacture the effect.

This supports the interpretation that the \((\log L)^2\) term, when present, belongs to the **environment-induced covariance/background fluctuation**, not to the conditional tiling fluctuation.

### Negative / null evidence

- fixed/random-environment LERW winding experiments: predominantly \(p\approx1\);
- strict one-walk-per-environment runs: ordinary log preferred;
- paired LERW winding differences: ordinary log preferred;
- refreshed-weight LERW pilot: ordinary simple-random-walk-like behaviour;
- Aztec central-height variance alone: not a convincing \(p=2\) diagnostic;
- two-tiling height difference alone: ordinary-log-like;
- square-grid structured Temperley paired experiments: no robust positive quadratic-log coefficient.

## Current supervisor-directed priority

The refreshed-weight random-walk check has been reported null.  Sunil's latest
supplied email therefore authorizes the next independent route: a square-grid
random-edge-weight dimer sampler based on local height-function Glauber
dynamics.  The requested first observable is the variance of the most central
face height after sufficiently separated equilibrium samples.

For one frozen mean-one edge-weight environment, select a face uniformly.  If
it has two opposite dimers, use the local heat bath with clockwise top-first
edge weights `(a,b,c,d)`: top/bottom with probability `ac/(ac+bd)`, and
right/left with probability `bd/(ac+bd)`.  Begin with exact tiny-grid checks
and a mixing pilot; do not launch a production campaign before those checks.

## Critical ambiguity to resolve before a large new run

The newer disorder-covariance observable relies on two samples sharing a well-defined random environment, whereas the "refresh weights at every step" experiment removes a fixed spatial environment. Do **not** silently invent a pairing convention.

Before implementing a production experiment, determine from the existing code/analysis exactly what quantity is meant by "the different observable that you found" in the latest instruction. If the intended shared randomness/estimand is not unambiguously recoverable, stop before a large run and surface the ambiguity.

## Statistical principles already established

- Resample at the independent-environment level, not individual within-environment observations.
- Preserve paired observations from a shared environment in the same bootstrap block.
- When pooling several separations \(r/L\), account for their within-environment dependence.
- Compare ordinary-log and quadratic-log models over several lower-size cutoffs.
- Treat fitted \((\log L)^2\) coefficients as finite-size numerical evidence, not proof.
- Use negative/no-disorder controls.
- Report sensitivity to weighting/fitting window and not only the most favourable fit.
- Keep model-selection metrics and prediction checks distinct from parameter confidence intervals.

## Computational principles

- deterministic, recorded seeds;
- independent RNG streams where scientifically required;
- restart-safe/chunked campaigns;
- save sufficient metadata to reproduce each run;
- avoid storing full environments when the statistic can be streamed, but never change the stochastic model merely to save memory;
- validate optimized observables against a slower/reference implementation on small sizes;
- pilot before scaling to Hamilton;
- keep raw simulation output separate from derived tables/figures;
- do not change an estimand during performance optimization.

## Repository role of these context files

- `context/RESEARCH_BRIEF.md`: stable scientific context and definitions.
- `context/RESULTS.md`: numerical record and caveats.
- `context/NEXT_STEPS.md`: current ordered work plan and decision gates.
- `context/SUNIL_EMAILS.md`: correspondence archive; read only when exact wording/history is needed.
- `AGENTS.md`: concise Codex operating instructions.

## AI-use boundary

Sunil explicitly stated that journal write-up rules are strict and that AI systems should not be used for the manuscript-writing stage. These files are therefore intended for **code, numerical analysis, reproducibility, and research planning**, not for drafting journal manuscript prose.
