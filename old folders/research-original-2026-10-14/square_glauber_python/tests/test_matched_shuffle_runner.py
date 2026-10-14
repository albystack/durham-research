"""Integration/provenance checks for the square runner and complete-data gate."""

import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from analyse_matched_shuffle import load_complete
from run_matched_shuffle import run_cell


def test_saved_pair_and_tampered_table_gate(tmp_path):
    source = SCRIPTS.parent / "configs/matched_shuffle_20261004.json"
    plan = json.loads(source.read_text())
    plan["cells"] = [plan["cells"][0]]
    p = tmp_path / "plan.json"
    p.write_text(json.dumps(plan))
    output = tmp_path / "run"
    result = run_cell(p, 0, output, True)
    assert result["state"] == "COMPLETE"
    data, _ = load_complete(output, p, True)
    assert len(data) == 8
    with pytest.raises(FileExistsError):
        run_cell(p, 0, output, True)
    table = output / "cell-00/measurements.csv"
    table.write_text(table.read_text().replace("square", "tampered", 1))
    with pytest.raises(ValueError, match="checksum"):
        load_complete(output, p, True)


def test_production_cannot_accept_validation_or_incomplete_cell(tmp_path):
    source = SCRIPTS.parent / "configs/matched_shuffle_20261004.json"
    plan = json.loads(source.read_text())
    plan["cells"] = [plan["cells"][0]]
    p = tmp_path / "plan.json"
    p.write_text(json.dumps(plan))
    output = tmp_path / "run"
    run_cell(p, 0, output, True)
    with pytest.raises(ValueError, match="incomplete"):
        load_complete(output, p, False)
    status = output / "cell-00/status.json"
    value = json.loads(status.read_text())
    value["state"] = "ERROR"
    status.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="incomplete"):
        load_complete(output, p, True)
