#!/usr/bin/env python3
"""Command-line entry point for the dimer-height/LERW experiments."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import time
from pathlib import Path

from analysis import make_fit_table, make_summary_table, write_csv, write_plots, write_run_report
from config import EXPERIMENTS, available_models_text, model_label, parse_model
from simulation import TrialTask, make_trial_tasks, run_trial


def list_presets() -> None:
    """Print configured experiment presets and model catalogue."""

    print("Experiment presets:\n")
    for name, preset in EXPERIMENTS.items():
        total = len(preset.sizes) * preset.samples_per_size * len(preset.models)
        workers = "auto" if preset.workers == 0 else str(preset.workers)
        print(name)
        print(f"  {preset.description}")
        print(f"  sizes: {list(preset.sizes)}")
        print(f"  samples per size: {preset.samples_per_size}")
        print(f"  models: {len(preset.models)}")
        print(f"  workers: {workers}")
        print(f"  total trials: {total}")
        print()
    print("Models:\n")
    print(available_models_text())


def _run_tasks(tasks: list[TrialTask], workers: int) -> list[dict[str, object]]:
    """Run tasks serially or through multiprocessing."""

    if workers <= 1:
        return [run_trial(task) for task in tasks]
    chunksize = max(1, len(tasks) // (workers * 8))
    with mp.Pool(processes=workers) as pool:
        return pool.map(run_trial, tasks, chunksize=chunksize)


def run_preset(
    preset_name: str,
    output_name: str | None,
    sizes_override: list[int] | None,
    samples_override: int | None,
    models_override: list[str] | None,
    workers_override: int | None,
    no_plots: bool,
) -> Path:
    """Run one preset and write raw data, summaries, fits, and plots."""

    if preset_name not in EXPERIMENTS:
        available = ", ".join(EXPERIMENTS)
        raise ValueError(f"unknown preset {preset_name!r}; choose one of: {available}")

    preset = EXPERIMENTS[preset_name]
    sizes = sizes_override if sizes_override else list(preset.sizes)
    samples = samples_override if samples_override is not None else preset.samples_per_size
    model_strings = models_override if models_override else list(preset.models)
    workers = workers_override if workers_override is not None else preset.workers
    if workers == 0:
        workers = max(1, mp.cpu_count() - 1)

    # Validate model strings before launching a long run.  The parsed values are
    # only used for printing here; task creation parses them again cleanly.
    parsed_models = [parse_model(model_string) for model_string in model_strings]

    run_name = output_name or preset_name
    output_dir = Path("results") / run_name
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Running preset: {preset_name}")
    print(f"Output folder: {output_dir}")
    print(f"Sizes: {sizes}")
    print(f"Samples per size: {samples}")
    print(f"Workers: {workers}")
    print("Models:")
    for model, parameter in parsed_models:
        print(f"  - {model_label(model, parameter)}")
    print()

    start = time.perf_counter()
    raw_rows: list[dict[str, object]] = []
    progress_rows: list[dict[str, object]] = []
    total_batches = len(model_strings) * len(sizes)
    batch_number = 0

    for model_string in model_strings:
        model, parameter = parse_model(model_string)
        for L in sizes:
            batch_number += 1
            batch_start = time.perf_counter()
            print(
                f"[{batch_number}/{total_batches}] "
                f"{model_label(model, parameter)} at L={L} ({samples} samples)",
                flush=True,
            )
            tasks = make_trial_tasks(
                sizes=[L],
                samples_per_size=samples,
                model_strings=[model_string],
                base_seed=preset.base_seed,
                max_steps=preset.max_steps,
            )
            batch_rows = _run_tasks(tasks, workers=workers)
            raw_rows.extend(batch_rows)
            batch_seconds = time.perf_counter() - batch_start
            progress_rows.append(
                {
                    "batch": batch_number,
                    "model": model,
                    "parameter": "" if parameter is None else parameter,
                    "L": L,
                    "samples": len(batch_rows),
                    "seconds": batch_seconds,
                    "cumulative_trials": len(raw_rows),
                }
            )
            write_csv(output_dir / "raw_trials.csv", raw_rows)
            write_csv(output_dir / "progress.csv", progress_rows)
            print(
                f"    wrote checkpoint: {len(raw_rows)} cumulative trials, "
                f"batch {batch_seconds:.2f}s",
                flush=True,
            )

    seconds = time.perf_counter() - start

    summary_rows = make_summary_table(raw_rows)
    fit_rows = make_fit_table(
        summary_rows=summary_rows,
        raw_rows=raw_rows,
        bootstrap_reps=preset.bootstrap_reps,
        bootstrap_seed=preset.base_seed,
    )

    write_csv(output_dir / "raw_trials.csv", raw_rows)
    write_csv(output_dir / "summary.csv", summary_rows)
    write_csv(output_dir / "fits.csv", fit_rows)
    write_csv(output_dir / "progress.csv", progress_rows)

    if no_plots:
        plot_status = "plots skipped by --no-plots"
    else:
        plot_status = write_plots(output_dir, summary_rows, fit_rows)

    write_run_report(
        output_dir / "run_report.md",
        preset_name=preset_name,
        raw_rows=raw_rows,
        summary_rows=summary_rows,
        fit_rows=fit_rows,
        seconds=seconds,
        plot_status=plot_status,
    )

    print(f"Wrote {len(raw_rows)} trials in {seconds:.2f}s")
    print(f"Summary: {output_dir / 'summary.csv'}")
    print(f"Fits: {output_dir / 'fits.csv'}")
    print(f"Report: {output_dir / 'run_report.md'}")
    print(plot_status)
    return output_dir


def build_parser() -> argparse.ArgumentParser:
    """Build the small CLI."""

    parser = argparse.ArgumentParser(
        prog="python3 main.py",
        description="Run random-edge LERW winding experiments for dimer-height scaling.",
    )
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("list", help="List presets and model names.")

    smoke = subparsers.add_parser("smoke", help="Run the tiny smoke preset.")
    smoke.add_argument("--workers", type=int, default=None)
    smoke.add_argument("--no-plots", action="store_true")

    run = subparsers.add_parser("run", help="Run a named preset.")
    run.add_argument("preset", type=str)
    run.add_argument("--output-name", type=str, default=None)
    run.add_argument("--sizes", type=int, nargs="*", default=None)
    run.add_argument("--samples", type=int, default=None)
    run.add_argument("--models", type=str, nargs="*", default=None)
    run.add_argument("--workers", type=int, default=None)
    run.add_argument("--no-plots", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> None:
    """CLI dispatcher."""

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command in (None, "list"):
        list_presets()
        return
    if args.command == "smoke":
        run_preset(
            preset_name="smoke",
            output_name=None,
            sizes_override=None,
            samples_override=None,
            models_override=None,
            workers_override=args.workers,
            no_plots=args.no_plots,
        )
        return
    if args.command == "run":
        run_preset(
            preset_name=args.preset,
            output_name=args.output_name,
            sizes_override=args.sizes,
            samples_override=args.samples,
            models_override=args.models,
            workers_override=args.workers,
            no_plots=args.no_plots,
        )
        return

    parser.print_help()


if __name__ == "__main__":
    main()
