"""Batch provenance and restart checks using genuine tiny square samples."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import larger_campaign as campaign


def plan_file(tmp_path):
    plan = {
        "purpose": "local restart validation",
        "sizes": [8],
        "models": {"uniform": 4},
        "validation_only": True,
        "batches": [
            {
                "index": i,
                "domain": domain,
                "model": "uniform",
                "size": 8,
                "pairs": 2,
                "master_seed": 20261007300000 + i,
            }
            for i, domain in enumerate(("square", "square", "aztec", "aztec"))
        ],
    }
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    return path


def test_square_restart_and_corruption(tmp_path):
    plan = plan_file(tmp_path)
    root = tmp_path / "data"
    campaign.run_batch(plan, 0, root)
    receipt_path = root / "batch-0000/complete.json"
    saved = receipt_path.read_bytes()
    campaign.run_batch(plan, 0, root)
    assert receipt_path.read_bytes() == saved
    assert len(list((root / "batch-0000").glob("attempt-*"))) == 1
    table = root / "batch-0000/attempt-0001/cells/cell-00/measurements.csv"
    table.write_text(table.read_text() + "\n")
    with pytest.raises(ValueError, match="evidence changed"):
        campaign.run_batch(plan, 0, root)


def test_retry_preserves_attempt_and_samples(tmp_path):
    plan = plan_file(tmp_path)
    root = tmp_path / "data"
    campaign.run_batch(plan, 0, root)
    (
        root / "batch-0000/complete.json"
    ).unlink()  # Crash after raw verification, before receipt.
    with pytest.raises(ValueError, match="explicit --retry"):
        campaign.run_batch(plan, 0, root)
    campaign.run_batch(plan, 0, root, retry=True)
    assert len(list((root / "batch-0000").glob("attempt-*"))) == 2
    assert (
        json.loads((root / "batch-0000/complete.json").read_text())["attempt"]
        == "attempt-0002"
    )


def test_plan_coverage_and_seed_identity(tmp_path):
    path = plan_file(tmp_path)
    plan = json.loads(path.read_text())
    plan["batches"][1]["master_seed"] = plan["batches"][0]["master_seed"]
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="duplicate"):
        campaign.read_plan(path)


def test_missing_batches_block_collection(tmp_path):
    plan = plan_file(tmp_path)
    root = tmp_path / "data"
    campaign.run_batch(plan, 0, root)
    with pytest.raises(FileNotFoundError):
        campaign.collect(plan, root)


def test_changed_plan_cannot_reuse_outputs(tmp_path):
    path = plan_file(tmp_path)
    root = tmp_path / "data"
    campaign.run_batch(path, 0, root)
    plan = json.loads(path.read_text())
    plan["purpose"] = "changed plan"
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="identity"):
        campaign.run_batch(path, 0, root)


def test_frozen_plan_totals():
    spec = importlib.util.spec_from_file_location(
        "freeze", SCRIPTS / "freeze_larger_campaign.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    plan = module.make_plan()
    assert len(plan["batches"]) == 504
    assert sum(b["pairs"] for b in plan["batches"]) == 258048
    assert all(b["size"] <= 64 for b in plan["batches"][:360])
    assert all(b["size"] > 64 for b in plan["batches"][360:])
