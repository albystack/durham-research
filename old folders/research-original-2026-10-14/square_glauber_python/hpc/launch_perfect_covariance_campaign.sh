#!/bin/bash
set -euo pipefail

repository_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$repository_root"

: "${OUTPUT_DIR:?Set a new production OUTPUT_DIR under /nobackup}"
MANIFEST="${MANIFEST:-configs/perfect_covariance_campaign_20260909_launch_manifest.csv}"
# This remains user-configurable. Hamilton sustained 256 workers in the
# September 2026 campaign; 128 is a conservative shared-partition default.
CONCURRENCY="${CONCURRENCY:-128}"

[[ "$OUTPUT_DIR" == /nobackup/* ]] || {
  echo "OUTPUT_DIR must be under /nobackup on Hamilton" >&2
  exit 2
}
test -x .hamilton_python/bin/python
test -f "$MANIFEST"
test -f "${MANIFEST%.csv}.json"
test ! -e "$OUTPUT_DIR" || {
  echo "refusing to overwrite existing production root $OUTPUT_DIR" >&2
  exit 2
}

.hamilton_python/bin/python -c '
import hashlib, json, sys
manifest, metadata_path = sys.argv[1:]
expected = json.load(open(metadata_path))["manifest_sha256"]
actual = hashlib.sha256(open(manifest, "rb").read()).hexdigest()
assert actual == expected, (actual, expected)
' "$MANIFEST" "${MANIFEST%.csv}.json"

mkdir -p logs "$OUTPUT_DIR"
cp "$MANIFEST" "$OUTPUT_DIR/campaign_manifest.csv"
if [[ -f "${MANIFEST%.csv}.json" ]]; then
  cp "${MANIFEST%.csv}.json" "$OUTPUT_DIR/campaign_manifest.json"
fi
frozen_manifest="$OUTPUT_DIR/campaign_manifest.csv"
task_count=$(.hamilton_python/bin/python -c \
  'import pandas as pd,sys; print(len(pd.read_csv(sys.argv[1])))' "$frozen_manifest")
last_task=$((task_count - 1))

smoke_output="${OUTPUT_DIR}_smoke"
smoke_job=$(sbatch --parsable \
  --export=ALL,SMOKE_OUTPUT="$smoke_output" \
  hpc/perfect_campaign_smoke.slurm)
production_job=$(sbatch --parsable \
  --dependency="afterok:${smoke_job}" \
  --array="0-${last_task}%${CONCURRENCY}" \
  --export=ALL,MANIFEST="$frozen_manifest",OUTPUT_DIR="$OUTPUT_DIR" \
  hpc/perfect_campaign_task.slurm)
audit_job=$(sbatch --parsable \
  --dependency="afterany:${production_job}" \
  --export=ALL,MANIFEST="$frozen_manifest",OUTPUT_DIR="$OUTPUT_DIR" \
  hpc/perfect_campaign_aggregate.slurm)

printf '%s\n' \
  "smoke_job=${smoke_job}" \
  "production_job=${production_job}" \
  "audit_job=${audit_job}" \
  "task_count=${task_count}" \
  "manifest=${frozen_manifest}" \
  "output=${OUTPUT_DIR}" | tee "$OUTPUT_DIR/launch_jobs.txt"
