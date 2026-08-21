# Results Ledger

> Status: numerical research record distilled from the supplied correspondence. Values below are
> evidence from the email chain, not independently re-derived from raw data. Before publication or
> formal reporting, verify every quoted number against saved outputs.

## 1. Baseline UST implementation — April 2026

Initial Python implementation:
- ordinary and wired square-grid USTs;
- random walks, loop erasure, Wilson's algorithm;
- seeded CLI and tree validation;
- automated tests: **37/37 passing**;
- ordinary \(8\times8\), seed 123: 64 vertices / 63 edges;
- wired \(8\times8\), seed 123: 37 quotient vertices / 36 edges.

This established the computational base for the later random-environment work.

## 2. Early fixed-environment random-walk model — June 2026

A two-random-weight construction with
\(w_N=w_E=1,\ w_S=u,\ w_W=v\) produced strong effective north-east drift after normalization.

At \(L=256\):
- Gamma \(k=1\): mean boundary hit about \((0.894L,0.889L)\);
- Gamma \(k=0.5\): about \((0.916L,0.918L)\).

The winding variance appeared flat, consistent with drift. A locally balanced model with all four outgoing weights independently sampled remained centred but looked closer to ordinary \(\log L\) growth than \((\log L)^2\) at the sizes then tested.

**Consequence:** move to locally balanced random environments with an independent sample at each site and test multiple distributions.

## 3. Large LERW winding campaign — July 2026

### 3.1 Strict one-walk-per-environment run

To avoid extra within-environment averaging:
- **379,000 walks**;
- **379,000 independent environments**;
- **15 distributions**;
- \(L=16\) to \(4096\);
- baseline and Gamma(shape=1) extended to \(L=8192\).

For
\[
\operatorname{Var}(W_L)=C(\log L)^p,
\]
fitted exponents were approximately **1.01–1.16** and the 95% intervals remained far below \(p=2\).

Model comparison:
- additive ordinary-log model preferred by BIC in **14/15** cases;
- Pareto \(\alpha=2\) was essentially tied rather than positive evidence for \(p=2\);
- \(L=8192\) baseline and Gamma points stayed compatible with ordinary log.

### 3.2 Paired same-environment LERWs

For each environment, two independent LERWs were generated and the winding difference \(W_1-W_2\) was analysed:
- **379,000 environments**;
- **758,000 raw walks**;
- **15 distributions**.

Effective exponents: **1.036–1.187**.
Largest reported 95% bootstrap upper bound: **1.260**.

Examples:
- baseline to \(L=8192\): \(p=1.083\);
- Gamma to \(L=8192\): \(p=1.036\).

BIC:
- \(a+b\log L\) favoured in **12/15** cases;
- Pareto(\(\alpha=3\)) and Lognormal(\(\sigma=1\)) favoured a quadratic-log additive curve at finite size;
- Lognormal(\(\sigma=0.5\)) was approximately tied;
- corresponding power exponents remained only about 1.187, 1.166 and 1.070.

At \(L=4096\), correlations between the two same-environment windings ranged roughly from **-0.08 to 0.17**, indicating modest shared-environment dependence.

**Conclusion:** no convincing winding-based \((\log L)^2\) regime.

## 4. Refreshed-weight LERW — late July 2026

Weights/transition probabilities were regenerated at each visit/time step rather than fixed spatially.

Pilot:
- about **46,000** runs;
- \(L=16\) to \(1024\);
- baseline, Gamma, lognormal, Pareto.

Reported \(p\):
- baseline: **1.064**;
- Gamma: **1.123**;
- lognormal: **1.069**;
- Pareto: **1.105**.

Ordinary logarithmic behaviour was preferred.

### Gamma high-variance extension / length observable

Gamma shape 0.5, variance 2; extra runs at \(L=2048,4096,5000,8192\).

Mean LERW length:
- \(L=5000\): **53,982**;
- \(L=8192\): **101,430**.

Fitted length exponent:
- **1.2518**;
- 95% bootstrap interval **[1.2456, 1.2573]**,

consistent with the expected \(L^{5/4}\) scaling.

Winding:
- \(p=1.134\);
- 95% interval **[0.957,1.265]**.

**Conclusion:** refreshing i.i.d. symmetric weights at each step behaves annealed-simple-random-walk-like for these observables.

## 5. Aztec diamond — central height

Random-weight Aztec simulations used domino shuffling / Gamma-specific reduction and direct computation of the central height rather than allocating the whole height table.

### Larger central-height run

Combined dataset:
- **32,000 independent samples**;
- sizes to \(L=800\).

For
\[
\operatorname{Var}(h_L)=C(\log L)^p:
\]
- \(L\ge24\): \(p=0.896\), 95% CI **[0.825,0.962]**;
- \(L\ge96\): \(p=1.018\), CI **[0.838,1.189]**;
- \(L\ge192\): \(p=1.171\), CI **[0.811,1.542]**.

The ordinary-log and quadratic-log additive fits were close, with no robust central-height \(p=2\) conclusion.

### Extended single-height run

Later email:
- **35,536 samples**;
- sizes to \(L=1300\);
- over \(L\ge24\), reported \(p=0.967\).

The email records the 95% interval as **"[90.894, 1.032]"**. This is almost certainly a transcription typo, but it is preserved here rather than silently corrected. Verify against the actual output before reuse.

## 6. Aztec diamond — two tilings, same Gamma environment

- **28,304 paired samples**;
- two conditionally independent tilings per environment.

For the height difference:
- \(p=0.814\);
- 95% CI **[0.700,0.902]**;
- only **0.87%** of bootstrap resamples favoured the squared-log fit.

**Interpretation:** the two-tiling difference isolates the conditional/connected component, so ordinary-log-like behaviour here is not evidence against a disorder-induced quadratic-log term.

## 7. Aztec diamond — paired spatial-height decomposition

For each Gamma environment, two conditionally independent tilings were sampled. Symmetric central-row increments were measured at
\[
r=L/32,\ L/16,\ L/8,\ L/4.
\]

Definitions:
\[
T(r)=\frac12\operatorname{Var}(\Delta H_1(r)-\Delta H_2(r)),
\qquad
D(r)=\operatorname{Cov}(\Delta H_1(r),\Delta H_2(r)).
\]

Dataset:
- **6,900 independent Gamma environments**;
- \(L=128,192,256,384,512,700,900,1100,1300\);
- **5,500 independent uniform-control pairs**.

Candidate fits for disorder covariance:
\[
D(r)=a_\rho+b_\rho\log r
\]
versus
\[
D(r)=a_\rho+b_\rho\log r+c(\log r)^2,
\]
with a common quadratic coefficient \(c\) across the four relative separations.

Bootstrap:
- **10,000** repetitions;
- whole environments resampled jointly across separations.

### Pooled unweighted Gamma disorder covariance

- \(c=0.584\);
- 95% interval **[0.196,0.952]**;
- bootstrap probability \(c>0\): **99.88%**;
- BIC difference (ordinary log minus quadratic-log): **10.85**, favouring quadratic log;
- held-out two-largest-order RMSE: **0.886 → 0.588**.

### Precision-weighted analysis

- \(c=0.495\);
- 95% interval **[0.193,0.798]**;
- BIC difference: **14.57**, favouring quadratic log;
- RMSE: **1.026 → 0.608**.

### Individual relative separations

Reported \(c\):
- \(L/32\): **0.566**;
- \(L/16\): **0.338**;
- \(L/8\): **1.045**;
- \(L/4\): **0.389**.

The \(L/32\) and \(L/8\) intervals excluded zero; the other two were less precise.

### Controls/components

Gamma conditional component:
- pooled \(c=-0.201\);
- weighted estimate likewise consistent with zero.

Uniform control:
- pooled marginal-variance \(c=0.009\);
- BIC favoured ordinary logarithmic growth;
- covariance fluctuated around zero, as expected without shared disorder.

**Main numerical interpretation:** the positive quadratic-log curvature appears in the **disorder-induced covariance**, not the connected/conditional component, and is absent in the uniform control.

## 8. Aztec disorder-law sensitivity — August 2026

Later campaign summary:
- original Gamma law: broad-window block-GLS \(c\approx0.502\);
- stronger-disorder law: \(c\approx0.608\);
- both reported with confidence intervals excluding zero and positive BIC support at the broad \(L\ge128\) window;
- weaker and symmetric disorder laws did not show comparable evidence;
- estimates became substantially less precise at higher lower-size cutoffs.

**Interpretation:** evidence is finite-size and parameter-sensitive, not universal across every disorder law tested.

## 9. Square-grid paired Temperley campaign — August 2026

Hamilton runs:
- **five disorder laws**;
- sizes to **\(L=6144\)**;
- paired/Temperley analysis;
- run description states **38,600 environments per model** for the submitted full campaign.

Across five models and three fitting windows:
- **14/15** confidence intervals for the quadratic-log contribution contained zero;
- **14/15** BIC comparisons favoured ordinary logarithmic growth;
- no-disorder control consistent with zero.

Important caveat:
- these are **structured spanning-tree/Temperley disorder models**;
- they are not automatically identical to a square-lattice dimer model with independent random weights on every dimer edge.

**Current status:** this null result triggered the latest supervisor request to verify that the spanning-tree/random-environment simulation is representing the intended dimer problem under additional randomness.

## 10. What is established vs not established

### Numerically well supported
- ordinary-log behaviour for the investigated LERW winding observables;
- no quadratic-log signal in the conditional/connected Aztec component;
- positive finite-size quadratic-log contribution in Aztec disorder covariance for at least the original Gamma and stronger-disorder settings;
- no comparable robust signal in the current structured square-grid Temperley implementation.

### Not established
- an asymptotic proof of \((\log L)^2\);
- universality over disorder laws;
- that the current square-grid Temperley disorder exactly represents the desired random-bond dimer model;
- the correct implementation of the proposed square-grid Glauber alternative;
- the asymptotic sign/magnitude of \(c\) beyond current finite sizes.

## 11. Verification checklist before any formal use

Recompute from raw outputs:
- sample counts by \(L\) and disorder law;
- all confidence intervals;
- BIC sign convention;
- bootstrap/block-GLS implementation;
- held-out RMSE definition;
- exact disorder-law parameterization;
- the malformed central-height confidence interval noted above.
