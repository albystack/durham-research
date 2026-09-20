# Provisional covariance analysis — 15 September 2026

The full existing analysis suite completed on **23,799 certified pairs out of 23,800 planned**, with **10,000 whole-environment bootstrap replicates**, seed `2026090902`. Results are provisional because both sample streams for `gamma-L32-env001015` remain unresolved. No further sampling was submitted.

## Findings

The gamma(shape=1) center disorder covariance rises from 0.714 at L=8 to 1.414 at L=32 (1.463 at L=28). Uniform center controls fluctuate around zero. The simple center fit is approximately `-0.473 + 0.559 log(L)`; the bootstrap 95% interval for its slope is [0.404, 0.718].

There is no robust evidence here for positive quadratic logarithmic growth:

| Gamma model fit | Quadratic coefficient | Bootstrap 95% interval |
|---|---:|---:|
| Center, all sizes | 0.0446 | [-0.3444, 0.4254] |
| Spatial, shared coefficient across separations | -0.2244 | [-0.6370, 0.1855] |

For center fits, AIC is -41.675 for the linear-log model and -40.012 for the quadratic model: a modest preference for the simpler model. Spatial AIC slightly favors the quadratic model, while BIC is essentially tied; this is not convincing evidence for a positive quadratic term. These fit comparisons are exploratory and use the existing unweighted least-squares procedure.

The center quadratic estimate changes from 0.0446 with all sizes to -0.4354 with L>=16 and -1.5176 with L>=20. Omitting L=32 entirely changes it to +0.2017 and reverses the center AIC preference. These changes demonstrate size sensitivity; four to seven sizes cannot securely identify the asymptotic growth law.

## Unresolved environment and sensitivity

The omission is 1/23,800 overall, but **1/1,200 in the gamma L=32 group**. It is related to sampling difficulty, so ordinary bootstrap intervals do not account for possible selection bias. Every included pair is certified; the full campaign is not certified complete.

As a conservative check, true extremal center heights at L=32 range from -33 to 31. Adding a missing pair anywhere within that rectangle could put the completed L=32 sample covariance between 0.5593 and 2.2684 (currently 1.4141). The corresponding all-size center quadratic coefficient could range approximately from -0.8191 to 0.9077. These are deliberately broad support bounds, not probable outcomes or imputed data. The omission cannot be declared harmless solely from its small fraction.

## Why the last samples take so long

The sampler uses coupling from the past with true lower and upper extremal states. Certification requires those chains to meet. Both remaining streams failed that check through a horizon of 1,073,741,824 sweeps. The implementation checkpoints completed failed horizons, then replays the next doubled horizon from the extrema using the assigned randomness. It does not resume halfway through a horizon. A doubled attempt is therefore roughly twice the previous 14.6–16 hour cost, and success is not guaranteed. This is computational work to certify two samples, not time spent analyzing the existing data.

Recommendation: retain the unresolved case and checkpoints, report these provisional results, and discuss the scientific value of another extension with Sunil before allocating more sampling time.

## Outputs and validation

`results/` contains center and spatial estimates, connected components, control plots, model comparisons, minimum-size sensitivity, diagnostics, and all bootstrap intervals. `provisional_center_intervals.png` includes error bars and the incompleteness label. Additional files contain fits excluding L=32 and conservative missing-pair sensitivity bounds.

All downloaded input tables match the Hamilton snapshot SHA256 checksums. The strict completeness check remains the default; a new explicit `allow_incomplete=True` option records the unchanged campaign summary and selection warning. Four targeted tests passed. The cached bootstrap reproduces the original center and spatial output on five full-data replicates with identical random draws to tolerance 1e-12; it preserves joint resampling of both chains and all spatial observables. All reported intervals contain 10,000 replicates.

The original Hamilton analysis job `18748230` was cancelled after moving to the validated cached local implementation. The local full run completed in 104.7 seconds. Frozen campaign data and sampling checkpoints were not modified. The analysis can be reproduced with `PYTHONPATH=src:hpc .venv/bin/python hpc/run_provisional_local_20260915.py` from `square_glauber_python`.
