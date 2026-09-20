"""Run from square_glauber_python with PYTHONPATH=src:hpc."""
import json
from pathlib import Path
import square_glauber.analysis as analysis
from cached_bootstrap_20260915 import bootstrap

root = Path(__file__).resolve().parents[1] / 'reports/provisional_20260915'
analysis._bootstrap = bootstrap
report = analysis.analyze_pilot(
    root / 'input', root / 'results', bootstrap_replicates=10000,
    bootstrap_seed=2026090902, allow_incomplete=True,
)
report['bootstrap_implementation'] = (
    'cached_bootstrap_20260915; matches original on 5 complete campaign '
    'replicates at rtol=atol=1e-12'
)
(root / 'results/analysis_metadata.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
