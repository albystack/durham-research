#!/bin/sh
set -eu

# Build a compact upload archive from the curated handoff manifest.
# The script copies files into a staging directory and then zips that staging
# tree. It never moves or deletes source material.

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
staging_root=$(mktemp -d "${TMPDIR:-/tmp}/research-chatgpt-upload.XXXXXX")
upload_root="$staging_root/research_chatgpt_upload"
archive_path="$staging_root/research_chatgpt_upload.zip"

mkdir -p "$upload_root"

copy_item() {
  item="$1"
  if [ -e "$repo_root/$item" ]; then
    mkdir -p "$upload_root/$(dirname "$item")"
    rsync -a --relative "$repo_root/$item" "$upload_root/"
    printf '%s\n' "$item"
  fi
}

echo "Included files and directories:"

for item in \
  README.md \
  .gitignore \
  .github/workflows/ci.yml \
  aztec/Project.toml \
  aztec/README.md \
  aztec/docs/IMPLEMENTATION.md \
  aztec/src/AztecDiamond.jl \
  aztec/test/runtests.jl \
  aztec/test/smoke_workflows.sh \
  aztec/scripts/analyze_double_dimer_campaign.jl \
  aztec/scripts/analyze_height_campaign.jl \
  aztec/scripts/analyze_spatial_campaign.jl \
  aztec/scripts/merge_double_dimer_batches.jl \
  aztec/scripts/merge_height_batches.jl \
  aztec/scripts/plot_double_dimer_campaign.jl \
  aztec/scripts/plot_height_campaign.jl \
  aztec/scripts/plot_spatial_campaign.jl \
  aztec/scripts/run_double_dimer_campaign.jl \
  aztec/scripts/run_gamma_disordered.jl \
  aztec/scripts/run_height_campaign.jl \
  aztec/scripts/run_random_weights.jl \
  aztec/scripts/run_spatial_campaign.jl \
  aztec/scripts/run_spatial_publication_campaign.sh \
  aztec/configs \
  aztec/data/README.md \
  aztec/data/examples \
  aztec/data/height/campaign_metadata.txt \
  aztec/data/double_dimer/campaign_metadata.txt \
  aztec/results/README.md \
  aztec/results/analysis_report.md \
  aztec/results/height \
  aztec/results/double_dimer \
  aztec/results/spatial \
  old/README.md \
  old/random_walk \
  chatgpt_handoff/CHATGPT_HANDOFF.md \
  chatgpt_handoff/FILE_INVENTORY.tsv \
  chatgpt_handoff/FIRST_FILES_TO_READ.md \
  chatgpt_handoff/NEXT_IMPLEMENTATION_PLAN.md \
  chatgpt_handoff/OPEN_QUESTIONS.md \
  chatgpt_handoff/REPOSITORY_TREE.txt \
  chatgpt_handoff/RESEARCH_RESULTS_TABLE.csv \
  chatgpt_handoff/UPLOAD_MANIFEST.md \
  chatgpt_handoff/UPLOAD_COMMANDS.sh
do
  copy_item "$item"
done

printf '\nCreating zip archive...\n'
(cd "$staging_root" && zip -qr "$archive_path" research_chatgpt_upload)
du -sh "$archive_path"
printf 'Archive path: %s\n' "$archive_path"
printf 'Staging directory: %s\n' "$upload_root"
