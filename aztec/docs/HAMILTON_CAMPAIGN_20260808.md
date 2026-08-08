# Hamilton campaign: 8 August 2026

Remote repository:
`/home/fvkl37/durham-research-square-grid/research`

Large outputs are written under `/nobackup/fvkl37` rather than the home quota.
All production arrays are gated by successful pilot jobs, and all analyses use
joint environment bootstrapping with the covariance-aware block-GLS fit.

## Square-grid paired experiment

The standard schedule contains 38,600 independent environments per law at
orders 32 through 2,048. The active arrays are:

- `18242101`: directed site-iid Gamma, shape 1;
- `18242102`: directed site-iid lognormal, variance 2;
- `18242103`: directed site-iid Uniform(0,2);
- `18242104`: undirected conductance Gamma, shape 0.5.

Jobs `18242185` through `18242196` are the dependent 10,000-bootstrap analyses
at minimum orders 32, 128, and 256.

The high-order extension contains 970 environments per model at orders 2,560,
3,072, 4,096, 5,120, and 6,144, for the baseline and five disorder models.
Pilots are `18243047` through `18243052`; production arrays are `18243072`
through `18243083`. Jobs `18243126` through `18243140` are the dependent
10,000-bootstrap analyses at minimum orders 128, 512, and 1,024.

## Aztec paired experiment

Each new Gamma law uses the original 6,900-environment schedule at orders 128
through 1,300, with two conditionally independent tilings per environment:

- stronger disorder `(alpha,beta)=(0.1,0.125)`: pilot `18243419`, production
  `18243432`;
- weaker disorder `(0.8,1.0)`: pilot `18243420`, production `18243433`;
- symmetric disorder `(0.25,0.25)`: pilot `18243421`, production `18243434`.

All nine safety pilots passed, including the three order-1,300 tasks. Jobs
`18243441` through `18243449` run 20,000-bootstrap analyses at minimum orders
128, 384, and 512. Matching analyses for the original `(0.2,0.25)` law are
jobs `18243453` through `18243455` and have completed. The focused end-to-end
runner validation was job `18243620` and passed.
Job `18243729` will consolidate all twelve block-GLS tables into one CSV and
one compact comparison report after the nine new-law analyses finish.

## Status commands

```sh
ssh -S /tmp/codex-hamilton.sock fvkl37@hamilton8.dur.ac.uk
squeue -u fvkl37
sacct -j JOB_ID --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

The square-grid output roots are
`/nobackup/fvkl37/square_grid_robustness_20260808` and
`/nobackup/fvkl37/square_grid_high_l_20260808`. The Aztec output root is
`/nobackup/fvkl37/aztec_gamma_parameters_20260808`.
