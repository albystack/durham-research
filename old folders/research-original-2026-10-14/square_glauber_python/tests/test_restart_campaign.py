import json
from pathlib import Path

import pandas as pd
import pytest

from square_glauber.experiment import derive_seeds
from square_glauber.perfect_campaign import _prepare_environment, create_campaign_manifest
from square_glauber.restart_campaign import (
    consolidate_case, exclusive_lock, frozen_case, run_sample_worker,
)


@pytest.fixture
def frozen(tmp_path):
    schedule = tmp_path / "schedule.csv"
    pd.DataFrame([dict(L=4, model="gamma", gamma_shape=1., environments=1,
        environments_per_task=1, max_horizon_sweeps=16384)]).to_csv(schedule, index=False)
    manifest = tmp_path / "manifest.csv"
    create_campaign_manifest(schedule, manifest, campaign_id="restart-test",
                             master_seed=917, backend="reference")
    output = tmp_path / "output"
    case = output / "cases/gamma/L004/environment-000000"
    _prepare_environment(case, L=4, model="gamma", gamma_shape=1., environment_index=0,
                         seeds=derive_seeds(917, 4, "gamma", 0), campaign_id="restart-test")
    return manifest, output, case


def test_workers_and_single_coordinator_preserve_frozen_inputs(frozen):
    manifest, output, case = frozen
    before = {name: (case / name).read_bytes() for name in
              ("environment.json", "environment_weights.npz")}
    first = run_sample_worker(manifest, 0, output, 0, 1, max_horizon_sweeps=16384,
                              chunk_attempts=29)
    assert first["certified_perfect_sample"]
    assert not (case / "pair.json").exists()
    assert not consolidate_case(manifest, 0, output, 0)["complete"]
    second = run_sample_worker(manifest, 0, output, 0, 2, max_horizon_sweeps=16384,
                               chunk_attempts=41)
    assert first["cftp_seed"] != second["cftp_seed"]
    assert not (case / "pair.json").exists()
    assert consolidate_case(manifest, 0, output, 0)["complete"]
    pair = (case / "pair.json").read_bytes()
    assert run_sample_worker(manifest, 0, output, 0, 1,
                             max_horizon_sweeps=16384) == first
    assert consolidate_case(manifest, 0, output, 0)["complete"]
    assert (case / "pair.json").read_bytes() == pair
    assert all((case / name).read_bytes() == data for name, data in before.items())
    assert not (output / "tasks").exists()


def test_missing_environment_is_not_created(frozen):
    manifest, output, case = frozen
    (case / "environment_weights.npz").unlink()
    with pytest.raises(FileNotFoundError):
        run_sample_worker(manifest, 0, output, 0, 1, max_horizon_sweeps=16)
    assert not (case / "environment_weights.npz").exists()


def test_seed_and_hash_mismatch_abort(frozen):
    manifest, output, case = frozen
    path = case / "environment.json"
    original = json.loads(path.read_text())
    for key in ("chain1_seed", "environment_sha256", "weights_npz_sha256"):
        damaged = dict(original)
        damaged[key] = 1 if key == "chain1_seed" else "wrong"
        path.write_text(json.dumps(damaged))
        with pytest.raises(RuntimeError):
            frozen_case(manifest, 0, output, 0)


def test_locks_reject_duplicate_writer_and_busy_coordinator(frozen):
    manifest, output, case = frozen
    with exclusive_lock(case / ".sample1.restart.lock"):
        with pytest.raises(RuntimeError, match="another restart worker"):
            run_sample_worker(manifest, 0, output, 0, 1, max_horizon_sweeps=16)
        with pytest.raises(RuntimeError, match="another restart worker"):
            consolidate_case(manifest, 0, output, 0)
        # Independent second sample has a different lock.
        with exclusive_lock(case / ".sample2.restart.lock"):
            pass


def test_manifest_hash_and_case_membership(frozen):
    manifest, output, case = frozen
    with pytest.raises(ValueError):
        frozen_case(manifest, 0, output, 1)
    with manifest.open("a") as stream:
        stream.write("\n")
    with pytest.raises(RuntimeError, match="manifest hash"):
        frozen_case(manifest, 0, output, 0)
