# Contributing

Contributions should preserve the mathematical model and make numerical claims
reproducible.

## Scientific requirements

- Do not change the stochastic model, pairing structure, observable, or
  estimand as a performance optimization.
- Keep deterministic seeds and independent random streams explicit.
- Treat the frozen environment—not individual MCMC draws—as the resampling
  block in quenched experiments.
- Validate optimized observables against a clear reference implementation on
  small cases.
- Pilot new schedules before production and retain failed-pilot diagnostics
  when they affect interpretation.
- Report controls, fit-window sensitivity, uncertainty, and prediction checks.

## Engineering requirements

- Add tests for mathematical identities and regressions.
- Keep library code in `aztec/src/` and command-line orchestration in
  `aztec/scripts/`.
- Write generated output below `aztec/output/`; do not commit scratch files or
  raw scheduler logs.
- Document new configs in `aztec/configs/README.md`.
- Use atomic outputs and validate existing batches before resuming.
- Avoid credentials, private correspondence, and developer-specific absolute
  paths in committed files.

## Before opening a change

```bash
julia --project=aztec -e 'using Pkg; Pkg.test()'
sh aztec/test/smoke_workflows.sh
git diff --check
```

Update the relevant model contract or results ledger whenever a change affects
scientific interpretation, data schemas, or reproduction commands.
