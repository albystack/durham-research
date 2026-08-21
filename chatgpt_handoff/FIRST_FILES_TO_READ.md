# First Files To Read

1. `aztec/src/AztecDiamond.jl` - authoritative sampler, height geometry, paired observables, validation, and SVG rendering.
2. `aztec/README.md` - high-level project summary, model definitions, commands, and current numerical conclusions.
3. `aztec/docs/IMPLEMENTATION.md` - code-flow walkthrough that links each kernel stage to the mathematical recurrence.
4. `aztec/test/runtests.jl` - verifies the recurrence, geometry, reproducibility, and variance identities.
5. `aztec/scripts/run_height_campaign.jl` - resumable center-height campaign runner and seed logic.
6. `aztec/scripts/analyze_height_campaign.jl` - bootstrap and fit methodology for the retained single-height result.
7. `aztec/scripts/run_double_dimer_campaign.jl` - shared-environment paired campaign runner and batch validation.
8. `aztec/scripts/analyze_double_dimer_campaign.jl` - paired variance decomposition and component fits.
9. `aztec/scripts/run_spatial_campaign.jl` - Gamma/uniform spatial experiment runner with fraction handling.
10. `aztec/scripts/analyze_spatial_campaign.jl` - pooled, weighted, holdout, and cutoff-sensitivity analysis.
11. `aztec/results/analysis_report.md` - concise summary of the retained numerical conclusions.
12. `aztec/results/spatial/README.md` - compact statement of the strongest current spatial result.
13. `aztec/results/spatial/spatial_analysis_report.md` - full spatial narrative and interpretation.
14. `aztec/data/README.md` - describes the retained raw datasets and their checksums.
15. `aztec/data/height/campaign_metadata.txt` - exact single-height schedule provenance and sample counts.
16. `aztec/data/double_dimer/campaign_metadata.txt` - exact paired-campaign provenance and sample counts.
17. `aztec/configs/README.md` - explains which schedules are production, continuation, smoke, and benchmark.
18. `aztec/configs/double_dimer_campaign.csv` - main paired center-height production schedule.
19. `aztec/configs/spatial_publication_gamma.csv` - publication Gamma spatial schedule.
20. `aztec/configs/spatial_publication_uniform.csv` - publication uniform control schedule.
21. `aztec/test/smoke_workflows.sh` - end-to-end CLI smoke coverage for all public commands.
22. `old/README.md` - explains the archived historical experiments and their scope.
23. `old/random_walk/README.md` - archived LERW project summary and remaining results.
24. `.github/workflows/ci.yml` - CI entry point and validation scope.
