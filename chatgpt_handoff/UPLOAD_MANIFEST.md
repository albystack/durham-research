# Upload Manifest

Estimated upload size: about 6-7 MB if the raw retained data CSVs are excluded, and about 14 MB if they are included. This manifest assumes the compact option: upload code, docs, configs, tests, current results, metadata, and small illustrative examples, but keep the raw retained data CSVs and all generated `aztec/output/` material on disk only.

## Essential

These files are needed to understand, modify, and validate the current Aztec code path.

- `README.md` - top-level project summary and current scientific framing.
- `aztec/Project.toml` - package metadata and dependency surface.
- `aztec/Manifest.toml` - current Julia dependency lockfile.
- `aztec/README.md` - authoritative project guide and current model description.
- `aztec/docs/IMPLEMENTATION.md` - detailed code-flow explanation.
- `aztec/src/AztecDiamond.jl` - core sampler, geometry, observables, and rendering.
- `aztec/test/runtests.jl` - unit and mathematical-reference tests.
- `aztec/test/smoke_workflows.sh` - end-to-end CLI smoke test.
- `aztec/scripts/*.jl` - campaign runners, mergers, analyses, and plotters.
- `aztec/configs/*.csv` - frozen schedules for production, continuation, smoke, and benchmark runs.
- `aztec/configs/README.md` - schedule provenance and usage notes.
- `aztec/data/README.md` - retained dataset descriptions.
- `aztec/data/height/campaign_metadata.txt` - exact single-height provenance.
- `aztec/data/double_dimer/campaign_metadata.txt` - exact paired provenance.
- `aztec/results/README.md` - retained result inventory.
- `aztec/results/analysis_report.md` - current numerical summary.
- `aztec/results/height/*` - height summaries, fits, curve table, and figure.
- `aztec/results/double_dimer/*` - double-dimer summaries, component fits, curve table, and figures.
- `aztec/results/spatial/*` - spatial summaries, fit tables, report, and figures.
- `.gitignore` - generated output and local cache exclusions.
- `.github/workflows/ci.yml` - CI validation scope.
- `chatgpt_handoff/CHATGPT_HANDOFF.md` - main audit and transfer summary.
- `chatgpt_handoff/FILE_INVENTORY.tsv` - machine-readable workspace inventory.
- `chatgpt_handoff/FIRST_FILES_TO_READ.md` - ordered reading sequence for the next assistant.
- `chatgpt_handoff/REPOSITORY_TREE.txt` - collapsed repository tree with counts and sizes.
- `chatgpt_handoff/RESEARCH_RESULTS_TABLE.csv` - campaign-by-campaign result table.
- `chatgpt_handoff/OPEN_QUESTIONS.md` - unresolved scientific and implementation issues.
- `chatgpt_handoff/NEXT_IMPLEMENTATION_PLAN.md` - staged follow-on plan.
- `chatgpt_handoff/UPLOAD_COMMANDS.sh` - safe packaging helper.
- `chatgpt_handoff/UPLOAD_MANIFEST.md` - this upload manifest.

## Useful

These files clarify provenance, history, or visual interpretation.

- `aztec/data/examples/*` - tiny illustrative inputs and exact example outputs.
- `old/README.md` - archive overview.
- `old/random_walk/README.md` - historical LERW project summary.
- `old/random_walk/src/*.jl` - archived LERW implementation.
- `old/random_walk/scripts/*.jl` - archived LERW command line.
- `old/random_walk/results/README.md` - archive result documentation.
- `aztec/results/tilings/*.png` - example tiling figures for the generated Aztec examples.
- `aztec/mathematica/draw_tiling.wl` - legacy renderer reference.

## Exclude

These paths are generated, historical scratch, large raw datasets, or not needed to modify the code safely.

- `aztec/data/height/center_height_samples.csv` - large retained raw observation table.
- `aztec/data/double_dimer/pairs.csv` - large retained paired raw observation table.
- `aztec/output/` - ignored scratch, benchmark, continuation, and historical batches.
- `old/previous_aztec_result_26050/center_height_variance.png` - superseded legacy figure.
- `.github/.DS_Store` - macOS metadata.
