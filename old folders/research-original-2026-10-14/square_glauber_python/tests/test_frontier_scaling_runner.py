"""The worker must preserve inputs and report censoring without sample files."""

import importlib.util
import json
from pathlib import Path

import pytest

from square_glauber.budgeted_frontier import VERSION

SPEC = importlib.util.spec_from_file_location(
    "frontier_runner",
    Path(__file__).parents[1] / "scripts/run_frontier_scaling.py",
)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)
SUMMARY_SPEC = importlib.util.spec_from_file_location(
    "frontier_summary",
    Path(__file__).parents[1] / "scripts/summarize_frontier_scaling.py",
)
summary = importlib.util.module_from_spec(SUMMARY_SPEC)
SUMMARY_SPEC.loader.exec_module(summary)


def plan_file(tmp_path, max_states):
    p = tmp_path / "plan.json"
    p.write_text(
        json.dumps(
            dict(
                sampler_version=VERSION,
                budgets=dict(
                    max_states=max_states, max_seconds=10, max_rss_bytes=2 * 1024**3
                ),
                cases=[
                    dict(
                        case_index=0,
                        L=4,
                        model="aztec_gamma",
                        environment_index=0,
                        environment_seed="513",
                        replica_seeds=["15", "16"],
                        input_source="generated_saved",
                    )
                ],
            )
        )
    )
    return p


def test_worker_certifies_complete_pairs_and_refuses_overwrite(tmp_path):
    p = plan_file(tmp_path, 1000)
    status = runner.run_case(p, 0, tmp_path / "out")
    assert status["state"] == "COMPLETE" and status["complete_pair"]
    case = tmp_path / "out/case-000"
    for name, digest in status["files_sha256"].items():
        assert runner.digest(case / name) == digest
    assert len(json.loads((case / "measurements.json").read_text())) == 1
    report = summary.summarize(p, tmp_path / "out")
    assert report["complete_pairs"] == 1 and report["all_pairs_complete"]
    with pytest.raises(FileExistsError):
        runner.run_case(p, 0, tmp_path / "out")
    (case / "sample1.npz").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        summary.summarize(p, tmp_path / "out")


def test_worker_keeps_censored_environment_without_partial_pair(tmp_path):
    p = plan_file(tmp_path, 1)
    status = runner.run_case(p, 0, tmp_path / "out")
    assert status["state"] == "BUDGET_EXHAUSTED"
    assert not status["complete_pair"] and status["reason"] == "state_limit"
    case = tmp_path / "out/case-000"
    assert (case / "weights.npz").is_file()
    assert not list(case.glob("sample*.npz"))
    assert not (case / "measurements.json").exists()
    report = summary.summarize(p, tmp_path / "out")
    assert report["complete_pairs"] == 0 and report["all_outcomes_recorded"]
    assert report["outcomes"] == {"BUDGET_EXHAUSTED": 1}


def test_missing_original_inputs_are_errors_not_regeneration(tmp_path):
    p = plan_file(tmp_path, 1000)
    plan = json.loads(p.read_text())
    plan["cases"][0]["input_source"] = "original_archive"
    p.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="original archive root"):
        runner.run_case(p, 0, tmp_path / "out")
    status = json.loads((tmp_path / "out/case-000/status.json").read_text())
    assert status["state"] == "ERROR" and not status["complete_pair"]


def test_summary_keeps_missing_cases_in_denominator(tmp_path):
    p = plan_file(tmp_path, 1000)
    report = summary.summarize(p, tmp_path / "out")
    assert report["expected_cases"] == 1 and not report["all_outcomes_recorded"]
    assert report["outcomes"] == {"MISSING": 1}
