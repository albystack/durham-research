"""Write immutable production/validation plans for the larger comparison."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def make_plan(validation=False):
    models = (
        {"aztec_gamma": 4, "aztec_mean": 4, "uniform": 4}
        if validation
        else {"aztec_gamma": 16384, "aztec_mean": 1024, "uniform": 1024}
    )
    plan = {
        "date": "2026-10-07",
        "version": "larger-batches-v1",
        "purpose": "Fresh independent larger matched-law comparison; no predetermined outcome",
        "validation_only": validation,
        "sizes": [16, 24, 32, 48, 64, 96, 128],
        "models": models,
        "probes": [0.125, 0.25],
        "height_scale": 1,
        "weights": {
            "a": "Gamma(shape=0.2,scale=1), north from black",
            "b": "Gamma(shape=0.25,scale=1), east from black",
            "remaining": 1,
            "normalization": "none",
        },
        "analysis": {
            "primary": "Gamma c_Aztec-c_square, fixed window 16..128",
            "primary_window": [16, 128],
            "windows": [[16, 128], [32, 128], [48, 128]],
            "bootstrap": 100 if validation else 10000,
            "bootstrap_seed": 2026100723 if validation else 2026100713,
            "permutations": 99 if validation else 9999,
            "control_seed": 2026100725 if validation else 2026100715,
            "resampling": "Whole environment vectors; pool batches within size/model/domain; independent groups",
            "control_test": "Predeclared maximum absolute standardized curvature contrast over both deterministic models and all three windows",
            "interpretation": "Finite-window contrast only; no absence/equivalence claim from non-rejection; control test is a diagnostic and never an exclusion rule",
        },
        "completion": "All batches and pairs required; preserve interrupted attempts; exact seed replay only; no replacement or sign-based extension; earlier cohorts and pilots excluded",
        "scheduler": {
            "concurrency_total": 2,
            "cpu_per_task": 1,
            "memory_GiB": 4,
            "small_batch_minutes": 30,
            "large_batch_minutes": 240,
            "small_tasks": 360,
            "large_tasks": 144,
            "maximum_sampling_cpu_hours": 756,
            "analysis_hours": 8,
            "analysis_memory_GiB": 8,
        },
        "batches": [],
    }
    for size in plan["sizes"]:
        for domain in ("square", "aztec"):
            for model, count in models.items():
                for start in range(0, count, 512):
                    index = len(plan["batches"])
                    plan["batches"].append(
                        {
                            "index": index,
                            "size": size,
                            "domain": domain,
                            "model": model,
                            "pairs": min(512, count - start),
                            "master_seed": (
                                20261007100000 if validation else 20261007000000
                            )
                            + index,
                        }
                    )
    return plan


if __name__ == "__main__":
    for validation in (False, True):
        suffix = "_validation" if validation else ""
        path = (
            ROOT / f"square_glauber_python/configs/larger_matched_20261007{suffix}.json"
        )
        text = json.dumps(make_plan(validation), indent=2) + "\n"
        if path.exists() and path.read_text() != text:
            raise ValueError("refusing to alter an existing frozen plan")
        path.write_text(text)
        print(path.name)
