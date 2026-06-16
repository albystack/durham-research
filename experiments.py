"""Experiment presets for the LERW random-environment project.

Edit this file when you want to add or change a run. The main script reads
these dictionaries and does the actual simulation.

Each experiment has:

- description: short human label printed by `python3 main.py list`
- sizes: box half-widths L for [-L, L]^2
- samples: number of independent trials for every model/size pair
- seed: base seed; main.py expands this into one seed per trial
- model_configs: models to run, written as "model" or "model:parameter"
- output: combined CSV path
- checkpoint_per_config: whether to write one temporary CSV per model config
- max_steps: optional safety cap for a walk; None uses the default cap
"""

from __future__ import annotations


EXPERIMENTS = {
    # Very small run for checking that the whole pipeline works.
    "smoke": {
        "description": "Quick balanced-model smoke test.",
        "sizes": [16, 32],
        "samples": 5,
        "seed": 20260615,
        "model_configs": [
            "symmetric",
            "gamma_balanced:1",
            "lognormal_balanced:0.5",
            "pareto_balanced:2",
            "uniform_balanced:0.5",
        ],
        "output": "results/smoke_balanced_distributions.csv",
        "checkpoint_per_config": False,
        "max_steps": None,
    },

    # Recreates the focused diagnostic dataset from the current write-up.
    "focused_reproduction": {
        "description": "Reproduce the focused drift-diagnostic dataset.",
        "sizes": [16, 32, 64, 128, 256],
        "samples": 700,
        "seed": 20260613,
        "model_configs": [
            "symmetric",
            "gamma:20",
            "gamma:1",
            "gamma:0.5",
            "gamma4:20",
            "gamma4:1",
            "gamma4:0.5",
        ],
        "output": "results/research_batch_focused.csv",
        "checkpoint_per_config": True,
        "max_steps": None,
    },

    # Main supervisor-directed run: exact local balance at every site.
    "balanced_large": {
        "description": "Large locally balanced run for the next research pass.",
        "sizes": [16, 32, 64, 128, 256, 512],
        "samples": 1000,
        "seed": 20260615,
        "model_configs": [
            "symmetric",
            "gamma_balanced:20",
            "gamma_balanced:1",
            "gamma_balanced:0.5",
            "lognormal_balanced:0.5",
            "lognormal_balanced:1",
            "pareto_balanced:3",
            "pareto_balanced:2",
            "uniform_balanced:0.5",
        ],
        "output": "results/research_batch_balanced_large.csv",
        "checkpoint_per_config": True,
        "max_steps": None,
    },
}
