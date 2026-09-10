# Pre-specified paired covariance analysis

This plan was frozen before the Hamilton production results were observed. Its
purpose is to test Professor Chhita's proposed super-rough behavior without
selecting environments, sizes, observables, or fit windows after seeing the
answer.

## Sampling unit and estimands

One statistical unit is one independently drawn frozen edge environment with
two conditionally independent CFTP-certified tilings. Repeated measurements do
not occur: each CFTP call returns one exact Gibbs sample. No environment is
removed because its weights are extreme or its certification horizon is long.

For paired observations `X1, X2` in the same environment, report

- disorder contribution: `Cov(X1, X2)` across environments;
- connected contribution: `0.5 Var(X1-X2)` across environments;
- marginal means, variances, and correlation.

The primary `X` values are the existing rightward central spatial height
increments at nominal separations `1/32, 1/16, 1/8, 1/4`, using the actual
deduplicated integer offset stored in each row. Centre height is the principal
secondary observable. The global height sum is supplementary and is not a
replacement scientific endpoint.

## Models

For centre-height disorder covariance, compare

`D(L) = a + b log L`

with

`D(L) = a + b log L + c (log L)^2`.

For spatial increments, the primary pooled extension gives every nominal
separation its own intercept and ordinary-log slope and gives all separations
one shared quadratic coefficient `c`. Report coefficients, SSE, AIC, BIC and
whole-environment bootstrap 95% intervals. The ordinary-log model is the null
description. The sign of `c` is unconstrained.

The full pre-specified size range is `L=8,12,16,20,24,28,32`. Sensitivity
tables may additionally report fixed lower cutoffs `L>=16` and `L>=20`, but no
other cutoff may be selected because it improves the result. Any later larger
sizes form a separately declared extension rather than silently changing this
campaign.

## Uncertainty and controls

Bootstrap whole environments, keeping both tilings and every observable from a
selected environment together. Use 10,000 bootstrap replicates with seed
`2026090902`. Never bootstrap the two configurations independently. Report
both the confidence interval and the fraction of bootstrap coefficients above
zero as descriptive quantities, alongside model-selection metrics and
residual plots.

The all-one model is a negative control: its paired covariance should be
consistent with zero because there is no environment disorder. A nonzero or
similarly curved control invalidates a straightforward disorder interpretation
and triggers an audit rather than deletion or refitting.

## Interpretation discipline

A positive quadratic point estimate is not itself evidence. A credible
numerical signal requires stability across the pooled spatial analysis,
environment bootstrap, fixed cutoff sensitivities, residuals, and the uniform
control. Conversely, failure to find a positive coefficient is a reportable
result and must not be repaired by changing the sampler or selectively adding
or removing environments.
