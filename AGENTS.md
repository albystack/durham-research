# AGENTS.md

This repository contains numerical research on random-environment LERW/spanning-tree and dimer models.

Before research changes, read:
1. `context/RESEARCH_BRIEF.md`
2. `context/NEXT_STEPS.md`
3. relevant sections of `context/RESULTS.md`

Read `context/SUNIL_EMAILS.md` only when exact supervisor wording/history is required.

For Slurm/Hamilton/NCC/HPC tasks, also read `context/HAMILTON_HPC.md`.

Rules:
- Scientific correctness and reproducibility outrank speed.
- Do not change the stochastic model, pairing structure, observable, or estimand as an optimization.
- Preserve deterministic seeds, environment-level pairing, restart-safe jobs, and enough metadata to reproduce runs.
- Validate optimized observables against a clear reference implementation on small cases.
- Pilot before launching large/HPC campaigns.
- Inspect only files relevant to the current task; avoid repo-wide wandering.
- Before editing, give a plan of at most 6 bullets. Afterward report only files changed, tests run, and scientific assumptions/risks.
- Square-grid Glauber production was explicitly authorized by Alberto on
  21 August 2026 without further supervisor review. Preserve the validated
  frozen-environment model, paired independent chains, central-height
  observable, equilibration diagnostics, and staged/restart-safe execution.
- Do not use AI to draft journal manuscript prose; the supervisor explicitly ruled this out for the write-up stage.
