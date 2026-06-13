#!/usr/bin/env python3
"""Run a moderate batch of LERW winding experiments."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from lerw_random_environment import TrialResult, run_experiment, write_csv


DEFAULT_SIZES = [16, 32, 64, 128, 256]
DEFAULT_K_VALUES = [20.0, 10.0, 5.0, 2.0, 1.0, 0.5]


def stable_config_seed(base_seed: int, model: str, k: float | None) -> int:
    label = f"{base_seed}:{model}:{k}"
    digest = hashlib.sha256(label.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")


def config_label(model: str, k: float | None) -> str:
    if k is None:
        return model
    return f"{model}_k{k:g}"


def run_config(
    model: str,
    k: float | None,
    sizes: list[int],
    samples: int,
    base_seed: int,
    max_steps: int | None,
) -> list[TrialResult]:
    print(f"\n== {config_label(model, k)} ==", flush=True)
    print(f"sizes={sizes}, samples_per_size={samples}", flush=True)
    results = run_experiment(
        sizes=sizes,
        samples=samples,
        model=model,
        k=k,
        seed=stable_config_seed(base_seed, model, k),
        max_steps=max_steps,
    )
    print(f"completed {len(results)} trials", flush=True)
    return results


def checkpoint_path(output: Path, model: str, k: float | None) -> Path:
    return output.with_name(f"{output.stem}_{config_label(model, k)}{output.suffix}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the main LERW research batch.")
    parser.add_argument("--sizes", type=int, nargs="+", default=DEFAULT_SIZES)
    parser.add_argument("--samples", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260612)
    parser.add_argument("--max-steps", type=int, default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/research_batch_main.csv"),
    )
    parser.add_argument(
        "--models",
        nargs="+",
        choices=["symmetric", "gamma", "gamma4"],
        default=["symmetric", "gamma", "gamma4"],
        help="Models to include. gamma/gamma4 run all k values.",
    )
    parser.add_argument("--k-values", type=float, nargs="+", default=DEFAULT_K_VALUES)
    parser.add_argument(
        "--checkpoint-per-config",
        action="store_true",
        help="Write one CSV per model/k as soon as that configuration finishes.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    all_results: list[TrialResult] = []

    if "symmetric" in args.models:
        results = run_config(
            model="symmetric",
            k=None,
            sizes=args.sizes,
            samples=args.samples,
            base_seed=args.seed,
            max_steps=args.max_steps,
        )
        all_results.extend(results)
        if args.checkpoint_per_config:
            path = checkpoint_path(args.output, "symmetric", None)
            write_csv(path, results)
            print(f"checkpoint wrote {path}", flush=True)

    for model in ("gamma", "gamma4"):
        if model not in args.models:
            continue
        for k in args.k_values:
            results = run_config(
                model=model,
                k=k,
                sizes=args.sizes,
                samples=args.samples,
                base_seed=args.seed,
                max_steps=args.max_steps,
            )
            all_results.extend(results)
            if args.checkpoint_per_config:
                path = checkpoint_path(args.output, model, k)
                write_csv(path, results)
                print(f"checkpoint wrote {path}", flush=True)

    write_csv(args.output, all_results)
    print(f"\nwrote {len(all_results)} total rows to {args.output}", flush=True)


if __name__ == "__main__":
    main()
