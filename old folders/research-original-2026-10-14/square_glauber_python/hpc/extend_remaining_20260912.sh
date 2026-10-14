#!/bin/bash
# One-off authorized extension of the existing frozen campaign; run on Hamilton.
set -euo pipefail
cd /home/fvkl37/durham-research-square-grid/research/square_glauber_python
export OUTPUT_DIR=/nobackup/fvkl37/square_glauber_perfect_covariance_20260910
export MANIFEST="$OUTPUT_DIR/campaign_manifest.csv"
export ANALYSIS_DIR="$OUTPUT_DIR/analysis_frozen_20260912"
export MAX_HORIZON_SWEEPS=268435456
extension_archive="$OUTPUT_DIR/extension_268435456_20260912"
task_ids=98,149,238,244,300,321,516,670,682,736,744,831,855,876,889
test ! -e "$extension_archive"
test ! -e "$ANALYSIS_DIR"
.hamilton_python/bin/python -c '
import hashlib,json,sys
from pathlib import Path
r=Path(sys.argv[1]); s=json.loads((r/"campaign_summary.json").read_text())
assert s["complete_environments"]==23785 and s["censored_environments"]==15
assert s["incomplete_slurm_array_spec"]==sys.argv[2]
assert hashlib.sha256((r/"campaign_manifest.csv").read_bytes()).hexdigest()=="2b0da28ac20a9c2155073961c5ee44b81a70a3d7545a611dfe386a0fd64bffed"
' "$OUTPUT_DIR" "$task_ids"
mkdir "$extension_archive"
mkdir "$extension_archive/tasks"
cp "$OUTPUT_DIR/campaign_summary.json" "$OUTPUT_DIR/campaign_status.csv" \
  "$OUTPUT_DIR/center_measurements.csv" "$OUTPUT_DIR/spatial_measurements.csv" \
  "$OUTPUT_DIR/chain_diagnostics.csv" "$OUTPUT_DIR/censored_environment_ids.txt" \
  "$OUTPUT_DIR/incomplete_task_ids.txt" "$OUTPUT_DIR/incomplete_slurm_array_spec.txt" \
  "$extension_archive/"
for task in 98 149 238 244 300 321 516 670 682 736 744 831 855 876 889; do
  printf -v task_file 'task-%06d.json' "$task"
  cp "$OUTPUT_DIR/tasks/$task_file" "$extension_archive/tasks/"
done
extension_job=$(sbatch --parsable --time=1-12:00:00 --mem=1G \
  --array="${task_ids}%15" \
  --export=ALL,MANIFEST,OUTPUT_DIR,MAX_HORIZON_SWEEPS \
  hpc/perfect_campaign_task.slurm)
echo "extension_job=$extension_job"
finish_job=$(sbatch --parsable --dependency="afterany:$extension_job" \
  --export=ALL,MANIFEST,OUTPUT_DIR,ANALYSIS_DIR \
  hpc/perfect_campaign_finish.slurm)
echo "finish_job=$finish_job"
echo "archive=$extension_archive"
