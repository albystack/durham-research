import json

import numpy as np
import pandas as pd
import pytest

from square_glauber.accelerated import NUMBA_AVAILABLE
from square_glauber.cli import build_parser
from square_glauber.perfect_campaign import create_campaign_manifest, run_campaign_task


@pytest.mark.parametrize('backend', ['reference', 'numba'])
def test_campaign_chunking_preserves_samples_and_counters(tmp_path, backend):
    if backend == 'numba' and not NUMBA_AVAILABLE:
        pytest.skip('optional Numba backend not installed')
    schedule = tmp_path / 'schedule.csv'
    pd.DataFrame([dict(L=4, model='gamma', gamma_shape=1.0, environments=1,
                       environments_per_task=1, max_horizon_sweeps=4096)]).to_csv(schedule, index=False)
    manifest = tmp_path / 'manifest.csv'
    create_campaign_manifest(schedule, manifest, campaign_id='chunk-equivalence',
                             master_seed=99017, backend=backend)
    results = []
    for chunk in [17, 4_000_000]:
        root = tmp_path / str(chunk)
        report = run_campaign_task(manifest, 0, root, chunk_attempts=chunk)
        assert report['chunk_attempts'] == chunk
        assert report['complete_pair_count'] == 1
        case = root / 'cases/gamma/L004/environment-000000'
        results.append(json.loads((case / 'pair.json').read_text()))
    for sample in ['sample1', 'sample2']:
        a, b = results[0][sample].copy(), results[1][sample].copy()
        a.pop('elapsed_seconds'); b.pop('elapsed_seconds')
        assert a == b
    assert results[0]['environment_sha256'] == results[1]['environment_sha256']


def test_chunk_argument_and_invalid_size(tmp_path):
    args = build_parser().parse_args(['perfect-campaign-task', '--manifest', 'm',
                                     '--task-id', '0', '--output', 'o', '--chunk-attempts', '4000000'])
    assert args.chunk_attempts == 4_000_000
    with pytest.raises(ValueError, match='chunk_attempts'):
        run_campaign_task(tmp_path / 'absent', 0, tmp_path, chunk_attempts=0)
