"""Explicitly provisional analysis; never alter the frozen campaign or its summary."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import time

ROOT = Path('/nobackup/fvkl37/square_glauber_perfect_covariance_20260910')
STAGE = ROOT / 'provisional_analysis_20260915'

def main():
    STAGE.mkdir(exist_ok=False)
    snapshot = STAGE / 'input'
    snapshot.mkdir()
    hashes = {}
    for name in ('campaign_summary.json', 'campaign_status.csv', 'center_measurements.csv',
                 'spatial_measurements.csv', 'chain_diagnostics.csv', 'campaign_manifest.csv',
                 'censored_environment_ids.txt'):
        shutil.copy2(ROOT / name, snapshot / name)
        hashes[name] = hashlib.sha256((snapshot / name).read_bytes()).hexdigest()
    (STAGE / 'input_sha256.json').write_text(json.dumps(hashes, indent=2))
    summary = json.loads((snapshot / 'campaign_summary.json').read_text())
    assert summary['complete_environments'] == 23799
    assert summary['censored_environments'] == 1
    note = ('PROVISIONAL: 23,799 certified pairs of 23,800 planned. '
            'gamma-L32-env001015 is unresolved and omitted. Selection may be informative; '
            'bootstrap intervals do not account for the omission. No additional sampling submitted.\n')
    (STAGE / 'README.txt').write_text(note)
    spec = importlib.util.spec_from_file_location('provisional_analysis', Path(__file__).with_name('analysis_provisional_20260915.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    started = time.monotonic()
    print(note, flush=True)
    report = module.analyze_pilot(snapshot, STAGE / 'results', bootstrap_replicates=10000,
                                  bootstrap_seed=2026090902, allow_incomplete=True)
    print(json.dumps(report, indent=2), flush=True)
    print(f'Completed in {time.monotonic()-started:.1f}s', flush=True)

if __name__ == '__main__':
    main()
