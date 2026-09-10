# Research overview

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
- direct square-grid frozen-edge dimers: no robust positive quadratic-log
  coefficient in the conditional, disorder, or total central-height
  components over `L=2--20`. A finite-volume Kasteleyn replay of all 960 Gamma
  environments confirms this independently of MCMC mixing.

## Direct weighted-dimer direction

The refreshed-weight random-walk check was null. The next independent route is
a square-grid random-edge-weight dimer sampler based on local height-function
Glauber dynamics. Its first observable is the variance of the most central
face height after sufficiently separated equilibrium samples.

For one frozen mean-one edge-weight environment, select a face uniformly.  If
it has two opposite dimers, use the local heat bath with clockwise top-first
edge weights `(a,b,c,d)`: top/bottom with probability `ac/(ac+bd)`, and
right/left with probability `bd/(ac+bd)`.  Begin with exact tiny-grid checks
and a mixing pilot; do not launch a production campaign before those checks.

The `L<=6` calibrated pilot checks passed locally and on Hamilton. Production
was authorized on 21 August 2026 after those validation stages. It keeps the
same frozen-edge model and central-height observable. Two
conditionally independent chains per environment are retained so the primary
annealed central-height variance can also be decomposed into conditional and
disorder components without changing the estimand.

That production campaign completed on 22 August 2026 with 960 Gamma
environments and 352 all-one controls. Environment-blocked bootstrap fits do
not support a positive nested `(log L)^2` term. Some individual MCMC estimates
showed mixing concerns at `L=16,20`, so the same 960 Gamma environments were
subsequently replayed with a finite-volume Kasteleyn calculation. The exact
conditional moments preserve the aggregate null conclusion without a
Markov-chain mixing assumption.

## Critical ambiguity to resolve before a large new run

The newer disorder-covariance observable relies on two samples sharing a well-defined random environment, whereas the "refresh weights at every step" experiment removes a fixed spatial environment. Do **not** silently invent a pairing convention.

For the direct weighted-dimer route, the production estimand is the variance of the
central-face height over frozen environments and equilibrium dimer samples.
The shared-environment pairing is explicit: each environment owns two
conditionally independent chains.  Individual MCMC draws are not independent
environment blocks. For the determinantal replay, the same decomposition is
computed directly from each environment's exact conditional mean and variance.

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

## Documentation map

- [`README.md`](README.md): documentation index.
- [`RESULTS.md`](RESULTS.md): numerical record and caveats.
- [`ROADMAP.md`](ROADMAP.md): current ordered work plan and decision gates.
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md): seeds, data layers, and
  verification workflow.
- [`ARCHITECTURE.md`](ARCHITECTURE.md): package and command-line structure.
