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

## 10. Direct square-grid weighted-dimer Glauber pilot — August 2026

The local heat-bath implementation was validated against exact tiny-grid
detailed balance and exact `L=2` enumeration.  An exact event-driven
self-loop skip preserves literal attempted-update time.  Parallel tempering
uses the same frozen edge environment at every inverse temperature and records
the physical observable only at `beta=1`.

An initial parallel-tempering pilot had a scheduling defect: the exchange
clock and odd/even adjacent-pair parity reset at every retained-sample call.
When the trace spacing was shorter than the exchange interval, this prevented
the final adjacent pair, and hence the `beta=1` replica, from exchanging during
the retained trace.  Each applied kernel still preserved the intended joint
law, so this was a mixing failure rather than a change of stationary target.
The clock and parity are now persistent and covered by a regression test.

For eight frozen Gamma(shape=0.5, scale=2) environments per size, the corrected
calibrated local ladder gave:
- sizes `L=2,3,4,6`;
- maximum absolute extremal-start gaps `0.120, 0.224, 0.236, 0.428`;
- maximum standardized gaps `1.73, 1.93, 1.53, 2.02`;
- minimum chain ESS `578.7, 65.4, 30.5, 28.9`;
- minimum adjacent-pair exchange acceptance at least `0.227` across all chains.

At `L=2`, the empirical centre-height distributions were within total-variation
distance `0.00002--0.00551` of exact enumeration.  The all-one control ladder
also showed extremal-start agreement, with maximum standardized gap at most
`2.19` and minimum chain ESS `178` across these sizes.

A fresh-seed repetition exposed two failures of that initial schedule: one
`L=2` environment had exact total-variation distance `0.0682`, and one `L=6`
environment had absolute/standardized extremal-start gaps `1.02`/`4.11`.
This ruled out treating the first successful ladder as a robust schedule.

A focused run with 10x burn-in and spacing on those same environments reduced
the `L=2` exact distance to `0.00997` and standardized start gap to `0.44`; at
`L=6`, the absolute/standardized gaps fell to `0.11`/`0.82`, with ESS about
`130` per chain.  This supports inadequate run length as the immediate cause;
the failing `L=6` chain still had minimum adjacent-pair exchange acceptance
`0.32`, so gross swap rejection was not the observed bottleneck.

These are local pilot diagnostics from small environment sets, not proof of
mixing and not authorization for production.  The longer candidate schedule
must be replicated across the full ladder on Hamilton with fresh seeds before
cautious extension to the next size.

Hamilton job `18546215` supplied that independent replication with base seed
`2026082107`.  All eight tasks completed successfully.  For `L=2,3,4,6`, the
maximum standardized extremal-start gaps were `1.95, 1.30, 1.57, 1.54`, and
the minimum chain ESS values were `1163, 131, 265, 147`.  The largest exact
`L=2` total-variation distance was `0.00533`.  Across 64 chains, every one of
the ten adjacent beta pairs was attempted; pair acceptance ranged from
`0.2275` to `0.9594`, while acceptance for the target-adjacent pair ranged from
`0.7485` to `0.9594`.  The hardest standardized environment was `L=2`,
environment 8 (`1.945` with raw gap `-0.014`); the largest raw gap was `0.29`
at `L=3`, environment 5 (`1.303` standardized).

This clears the calibrated `L<=6` pilot gate but remains too small to authorize
production.  The next controlled step is a four-environment `L=8` pilot.

On 21 August 2026 Alberto explicitly authorized production without waiting for
further supervisor review because Sunil was leaving on vacation, and identified
the prior repository prohibition as outdated.  This authorization does not
turn mixing evidence into a proof: production retains extremal-start pairs and
per-pair exchange diagnostics, and results at sizes beyond the calibrated range
must be filtered by those diagnostics before scientific interpretation.

Three production arrays were submitted on Hamilton:
- `18546281`: core Gamma sizes `L=2,4,6,8,10,12`, 576 tasks;
- `18546369`: frontier Gamma sizes `L=16,20`, 96 tasks;
- `18546370`: all-one control sizes `L=4,6,8,10,12,16`, 240 tasks.

The Gamma arrays share base seed `2026082201` across disjoint size sets; the
control uses `2026082202`.  At startup, all configured 256 concurrent
production tasks were running with no non-empty error logs.  The frontier is
kept in a separate output root so timeouts or failed mixing diagnostics cannot
invalidate the core campaign.

All three arrays subsequently completed: 576/576 core Gamma tasks, 96/96
frontier Gamma tasks, and 240/240 control tasks, all with exit code `0:0`.
Together their output roots contain the expected 2,739 files and no temporary
or non-empty error files.  The final sample comprises 960 independent frozen
Gamma environments and 352 all-one controls.  Each environment has two
independent extremal-start chains with 2,000 retained central-height values.

The environment-blocked Gamma component estimates were:

| L | environments | conditional | disorder | total |
|---:|---:|---:|---:|---:|
| 2 | 64 | 0.456 | 0.776 | 1.232 |
| 4 | 192 | 1.183 | 1.177 | 2.360 |
| 6 | 192 | 1.350 | 1.588 | 2.938 |
| 8 | 192 | 1.709 | 1.887 | 3.596 |
| 10 | 128 | 1.822 | 2.375 | 4.197 |
| 12 | 96 | 1.855 | 2.662 | 4.517 |
| 16 | 64 | 2.301 | 2.199 | 4.500 |
| 20 | 32 | 2.302 | 3.137 | 5.439 |

At every size, `conditional + disorder` agrees with the independently computed
direct annealed variance to within 0.67 percent.  In the control, disorder
covariance lies between `-0.00131` and `0.000420`, while total variance is
effectively the conditional component.

The production scaling analysis resamples whole frozen environments within
each size and compares the nested models `a+b log(L)` and
`a+b log(L)+c(log(L))^2`.  For the unweighted full range `L=2--20`:

| component | c | 95% environment bootstrap interval | P(c>0) | delta BIC | LOOCV RMSE log / quadratic |
|:--|--:|:--|--:|--:|:--|
| conditional | -0.048 | [-0.320, 0.227] | 0.357 | -1.283 | 0.114 / 0.142 |
| disorder | 0.108 | [-0.545, 0.794] | 0.574 | -1.632 | 0.329 / 0.386 |
| total | 0.060 | [-0.721, 0.822] | 0.522 | -1.845 | 0.242 / 0.279 |

Positive delta BIC favors the quadratic extension.  Here every point-estimate
BIC difference favors ordinary log, and ordinary log also predicts held-out
sizes better.  Inverse-bootstrap-variance weighting and lower cutoffs at
`L=4,6,8` do not yield a stable positive coefficient.  Excluding environments
with standardized extremal-start gap above 4 changes the full-range
conditional/disorder/total coefficients to `-0.107`, `0.007`, and `-0.101`;
it therefore does not manufacture the target signal.  A second 5,000-replicate
bootstrap seed left these conclusions unchanged.

Mixing remains the main qualification.  Maximum standardized extremal-start
gaps reach `6.28` at `L=16` and `7.76` at `L=20`, with 4/64 and 3/32
environments above 4.  Mean signed gaps are only `0.0143` and `0.0711`, with
across-environment t statistics `0.39` and `0.92`, so there is no detected
common start direction.  Minimum chain ESS remains `88.8` and `95.5`.
Adjacent-pair exchange is monotone across the ladder: at `L=20`, pooled
acceptance is `0.062` for `beta=0--0.05` but `0.637` for `beta=0.95--1.00`.
The target replica therefore exchanges locally, but full round trips were not
recorded and cannot be certified.  The production result is consequently a
finite-size null with a larger-size mixing caveat, not a proof of asymptotic
ordinary-log behaviour.

## 11. What is established vs not established

### Numerically well supported
- ordinary-log behaviour for the investigated LERW winding observables;
- no quadratic-log signal in the conditional/connected Aztec component;
- positive finite-size quadratic-log contribution in Aztec disorder covariance for at least the original Gamma and stronger-disorder settings;
- no comparable robust signal in the current structured square-grid Temperley implementation.

### Not established
- an asymptotic proof of \((\log L)^2\);
- universality over disorder laws;
- that the current square-grid Temperley disorder exactly represents the desired random-bond dimer model;
- production-scale equilibration of the proposed square-grid Glauber
  alternative beyond the calibrated `L<=6` pilot;
- the asymptotic sign/magnitude of \(c\) beyond current finite sizes.

## 12. Verification checklist before any formal use

Recompute from raw outputs:
- sample counts by \(L\) and disorder law;
- all confidence intervals;
- BIC sign convention;
- bootstrap/block-GLS implementation;
- held-out RMSE definition;
- exact disorder-law parameterization;
- the malformed central-height confidence interval noted above.
