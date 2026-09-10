"""Command-line interface for exact validation, diagnostics, and local pilots."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Sequence

import numpy as np

from .analysis import analyze_pilot
from .audit import write_exact_center_height_audit, write_height_convention_audit
from .cftp_calibration import (
    consolidate_cftp_cost_calibration,
    extend_censored_cftp_calibration,
    run_cftp_cost_calibration,
)
from .diagnostics import compare_extremal_starts
from .environment import sample_environment
from .exact import standard_exact_validations
from .experiment import _seed_uint64, run_pilot
from .final_calibration import (
    build_three_run_comparison,
    run_final_hard_calibration,
    write_final_results_summary,
)
from .geometry import DEFAULT_RELATIVE_SEPARATIONS, spatial_separations, validate_L
from .mixing_campaign import run_mixing_campaign, write_mixing_case
from .mixing_extension import build_old_new_comparison, run_frozen_gamma_extension
from .performance import (
    audit_common_randomness_ordering,
    benchmark_backends,
    run_accelerated_frozen_stress,
    run_backend_equivalence_audit,
    run_common_randomness_coupling,
)
from .order_audit import (
    write_cftp_prerequisites,
    write_exact_extremality_audit,
    write_exhaustive_monotonicity_audit,
    write_local_pattern_monotonicity_audit,
)
from .perfect_audit import (
    run_exact_perfect_sampling_diagnostic,
    run_frozen_hard_cftp_calibration,
    run_perfect_pair_smoke,
)
from .perfect_campaign import (
    aggregate_campaign,
    create_campaign_manifest,
    run_campaign_task,
)


def _comma_ints(value: str) -> tuple[int, ...]:
    try:
        result = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not result:
        raise argparse.ArgumentTypeError("expected at least one integer")
    return result


def _comma_strings(value: str) -> tuple[str, ...]:
    result = tuple(item.strip() for item in value.split(",") if item.strip())
    if not result:
        raise argparse.ArgumentTypeError("expected at least one value")
    return result


def _comma_floats(value: str) -> tuple[float, ...]:
    try:
        result = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated numbers") from error
    if not result:
        raise argparse.ArgumentTypeError("expected at least one number")
    return result


def _size_sweep_mapping(value: str) -> dict[int, int]:
    result: dict[int, int] = {}
    try:
        for item in value.split(","):
            size_text, sweeps_text = item.split(":", maxsplit=1)
            result[int(size_text.strip())] = int(sweeps_text.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "expected comma-separated L:sweeps pairs, for example 8:100000,12:200000"
        ) from error
    if not result:
        raise argparse.ArgumentTypeError("expected at least one L:sweeps pair")
    return result


def command_validate_exact(args: argparse.Namespace) -> int:
    reports = standard_exact_validations(args.seed)
    tolerance = args.tolerance
    passed = True
    for report in reports:
        row = asdict(report)
        row["passed"] = (
            report.row_sum_residual <= tolerance
            and report.stationarity_residual <= tolerance
            and report.detailed_balance_residual <= tolerance
            and report.illegal_transition_count == 0
        )
        passed &= bool(row["passed"])
        print(json.dumps(row, sort_keys=True))
    expected_counts = {(2, 2): 2, (2, 4): 5, (4, 4): 36}
    passed &= all(expected_counts[(item.rows, item.columns)] == item.matching_count for item in reports)
    print(f"exact validation: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


def command_mixing(args: argparse.Namespace) -> int:
    L = validate_L(args.L)
    output = Path(args.output or f"outputs/mixing_L{L}_{args.model}")
    output.mkdir(parents=True, exist_ok=True)
    root = np.random.SeedSequence(args.seed)
    environment_sequence, chain1_sequence, chain2_sequence = root.spawn(3)
    environment_seed = _seed_uint64(environment_sequence)
    chain1_seed = _seed_uint64(chain1_sequence)
    chain2_seed = _seed_uint64(chain2_sequence)
    environment = sample_environment(
        L,
        np.random.Generator(np.random.PCG64(environment_seed)),
        args.model,
        gamma_shape=args.gamma_shape,
    )
    separations = spatial_separations(L, args.spatial_rhos)
    comparison = compare_extremal_starts(
        environment,
        np.random.Generator(np.random.PCG64(chain1_seed)),
        np.random.Generator(np.random.PCG64(chain2_seed)),
        sweeps=args.sweeps,
        record_every=args.record_every,
        separations=separations,
    )
    metadata = {
        "environment_id": f"{args.model}-L{L}-mixing",
        "L": L,
        "model": args.model,
        "gamma_shape": args.gamma_shape if args.model == "gamma" else None,
        "master_seed": args.seed,
        "environment_seed": environment_seed,
        "chain1_seed": chain1_seed,
        "chain2_seed": chain2_seed,
        "actual_separations": [asdict(value) for value in separations],
    }
    write_mixing_case(
        comparison,
        output,
        metadata=metadata,
        max_lag_records=args.acf_max_lag,
    )
    summary = comparison.summary() | metadata
    print(json.dumps(summary, indent=2))
    print(f"wrote mixing diagnostics to {output}")
    return 0


def command_height_audit(args: argparse.Namespace) -> int:
    report = write_height_convention_audit(args.output, args.sizes)
    print(json.dumps(report, indent=2))
    print(f"wrote height-convention audit to {args.output}")
    return 0


def command_exact_height_audit(args: argparse.Namespace) -> int:
    report = write_exact_center_height_audit(
        args.output,
        seed=args.seed,
        burnin_sweeps=args.burnin_sweeps,
        measurement_gap_sweeps=args.measurement_gap_sweeps,
        num_measurements=args.num_measurements,
    )
    print(json.dumps(report, indent=2))
    print(f"wrote exact centre-height audit to {args.output}")
    return 0


def command_mixing_campaign(args: argparse.Namespace) -> int:
    report = run_mixing_campaign(
        sizes=args.sizes,
        gamma_environments=args.gamma_environments,
        uniform_replicates=args.uniform_replicates,
        gamma_shape=args.gamma_shape,
        sweeps=args.sweeps,
        record_every=args.record_every,
        master_seed=args.seed,
        output=args.output,
        relative_separations=args.spatial_rhos,
        max_lag_records=args.acf_max_lag,
    )
    print(json.dumps(report, indent=2))
    print(f"wrote mixing campaign to {args.output}")
    return 0


def command_mixing_extend_frozen(args: argparse.Namespace) -> int:
    report = run_frozen_gamma_extension(
        source=args.source,
        output=args.output,
        sweeps_by_size=args.sweeps_by_size,
        record_every=args.record_every,
        relative_separations=args.spatial_rhos,
        max_lag_records=args.acf_max_lag,
    )
    print(json.dumps(report, indent=2))
    print(f"wrote frozen-environment mixing extension to {args.output}")
    return 0


def command_mixing_finalize(args: argparse.Namespace) -> int:
    table = build_old_new_comparison(args.source, args.output)
    print(f"wrote {len(table)} assessed comparison rows to {args.output}")
    return 0


def command_backend_equivalence(args: argparse.Namespace) -> int:
    report = run_backend_equivalence_audit(
        args.output,
        sizes=args.sizes,
        attempts=args.attempts,
        checkpoint_attempts=args.checkpoint_attempts,
        seed=args.seed,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_benchmark_backends(args: argparse.Namespace) -> int:
    report = benchmark_backends(
        args.output,
        sizes=args.sizes,
        target_attempts=args.target_attempts,
        repeats=args.repeats,
        seed=args.seed,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_ordering_audit(args: argparse.Namespace) -> int:
    report = audit_common_randomness_ordering(
        args.output,
        sizes=args.sizes,
        attempts=args.attempts,
        streams_per_model_size=args.streams,
        seed=args.seed,
    )
    print(json.dumps(report, indent=2))
    return 0 if report["ordering_preserved_in_all_checked_cases"] else 1


def command_accelerated_stress(args: argparse.Namespace) -> int:
    report = run_accelerated_frozen_stress(
        args.source_case,
        args.output,
        sweeps=args.sweeps,
        record_every=args.record_every,
        relative_separations=args.spatial_rhos,
        max_lag_records=args.acf_max_lag,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_common_randomness_coupling(args: argparse.Namespace) -> int:
    report = run_common_randomness_coupling(
        args.source_case,
        args.output,
        ordering_audit=args.ordering_audit,
        sweeps=args.sweeps,
        record_every=args.record_every,
        seed=args.seed,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_exact_extremality_audit(args: argparse.Namespace) -> int:
    report = write_exact_extremality_audit(args.output)
    print(json.dumps(report, indent=2))
    return 0


def command_exhaustive_monotonicity_audit(args: argparse.Namespace) -> int:
    report = write_exhaustive_monotonicity_audit(args.output)
    print(json.dumps(report, indent=2))
    return 0 if report["all_exhaustive_small_state_checks_preserved_order"] else 1


def command_cftp_prerequisites(args: argparse.Namespace) -> int:
    extremality = json.loads(Path(args.extremality_summary).read_text())
    monotonicity = json.loads(Path(args.monotonicity_summary).read_text())
    write_cftp_prerequisites(args.output, extremality, monotonicity)
    print(f"wrote CFTP prerequisite audit to {args.output}")
    return 0


def command_local_monotonicity_audit(args: argparse.Namespace) -> int:
    report = write_local_pattern_monotonicity_audit(args.output)
    print(json.dumps(report, indent=2))
    return 0 if report["all_cases_preserve_order"] else 1


def command_exact_perfect_sampling(args: argparse.Namespace) -> int:
    report = run_exact_perfect_sampling_diagnostic(
        args.output,
        samples=args.samples,
        seed=args.seed,
        max_horizon_sweeps=args.max_horizon_sweeps,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_cftp_hard_calibration(args: argparse.Namespace) -> int:
    report = run_frozen_hard_cftp_calibration(
        args.source_case,
        args.output,
        samples=args.samples,
        seed=args.seed,
        max_horizon_sweeps=args.max_horizon_sweeps,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_perfect_pair_smoke(args: argparse.Namespace) -> int:
    report = run_perfect_pair_smoke(
        args.output,
        sizes=args.sizes,
        seed=args.seed,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_create_perfect_campaign_manifest(args: argparse.Namespace) -> int:
    frame = create_campaign_manifest(
        args.schedule,
        args.output,
        campaign_id=args.campaign_id,
        master_seed=args.seed,
        backend=args.backend,
    )
    print(
        json.dumps(
            {
                "manifest": args.output,
                "task_count": len(frame),
                "environment_count": int(frame["environment_count"].sum()),
                "perfect_sample_count": int(2 * frame["environment_count"].sum()),
            },
            indent=2,
        )
    )
    return 0


def command_perfect_campaign_task(args: argparse.Namespace) -> int:
    report = run_campaign_task(
        args.manifest,
        args.task_id,
        args.output,
        max_horizon_override=args.max_horizon_sweeps,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_aggregate_perfect_campaign(args: argparse.Namespace) -> int:
    report = aggregate_campaign(
        args.manifest,
        args.output,
        require_complete=not args.allow_incomplete,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_cftp_cost_calibration(args: argparse.Namespace) -> int:
    report = run_cftp_cost_calibration(
        args.output,
        standard_sizes=args.standard_sizes,
        standard_gamma_environments=args.standard_gamma_environments,
        standard_uniform_environments=args.standard_uniform_environments,
        standard_samples_per_environment=args.standard_samples_per_environment,
        probe_sizes=args.probe_sizes,
        probe_gamma_environments=args.probe_gamma_environments,
        probe_uniform_environments=args.probe_uniform_environments,
        probe_samples_per_environment=args.probe_samples_per_environment,
        gamma_shape=args.gamma_shape,
        master_seed=args.seed,
        max_horizon_sweeps=args.max_horizon_sweeps,
        backend=args.backend,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_cftp_extend_censored(args: argparse.Namespace) -> int:
    report = extend_censored_cftp_calibration(
        args.source,
        args.output,
        max_horizon_sweeps=args.max_horizon_sweeps,
        backend=args.backend,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_cftp_consolidate(args: argparse.Namespace) -> int:
    report = consolidate_cftp_cost_calibration(
        args.source,
        args.extension,
        args.output,
    )
    print(json.dumps(report, indent=2))
    return 0


def command_final_hard_calibration(args: argparse.Namespace) -> int:
    report = run_final_hard_calibration(
        args.source_case,
        args.one_million_case,
        args.output,
        sweeps=3_000_000,
        record_every=args.record_every,
        relative_separations=args.spatial_rhos,
        max_lag_records=args.acf_max_lag,
        window_count=args.windows,
    )
    table = build_three_run_comparison(
        args.source_case, args.one_million_case, args.output
    )
    print(json.dumps(report, indent=2))
    print(table.to_string(index=False))
    return 0


def command_final_audit_summary(args: argparse.Namespace) -> int:
    path = write_final_results_summary(
        args.root,
        classification=args.classification,
        recommendation=args.recommendation,
        late_window_assessment=args.late_window_assessment,
        visual_assessment=args.visual_assessment,
        test_results=args.test_results,
    )
    print(f"wrote final audit summary to {path}")
    return 0


def command_pilot(args: argparse.Namespace) -> int:
    metadata = run_pilot(
        sizes=args.sizes,
        models=args.models,
        gamma_shape=args.gamma_shape,
        environments=args.environments,
        burnin_sweeps=args.burnin_sweeps,
        measurement_gap_sweeps=args.measurement_gap_sweeps,
        num_measurements=args.num_measurements,
        master_seed=args.seed,
        output=args.output,
        relative_separations=args.spatial_rhos,
    )
    print(json.dumps(metadata, indent=2))
    print(f"wrote pilot to {args.output}")
    return 0


def command_analyze(args: argparse.Namespace) -> int:
    report = analyze_pilot(
        args.input,
        args.output,
        bootstrap_replicates=args.bootstrap_replicates,
        bootstrap_seed=args.bootstrap_seed,
    )
    print(json.dumps(report, indent=2))
    print(f"wrote exploratory analysis to {args.output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m square_glauber.cli",
        description="Local proof-of-correctness and pilot tools for square-grid dimers.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    exact = subparsers.add_parser("validate-exact", help="run deterministic tiny-grid proofs")
    exact.add_argument("--seed", type=int, default=20260907)
    exact.add_argument("--tolerance", type=float, default=1e-12)
    exact.set_defaults(function=command_validate_exact)

    height_audit = subparsers.add_parser(
        "height-audit", help="write extremal height matrices and labelled figures"
    )
    height_audit.add_argument("--sizes", type=_comma_ints, default=(6, 8))
    height_audit.add_argument("--output", required=True)
    height_audit.set_defaults(function=command_height_audit)

    exact_height = subparsers.add_parser(
        "exact-height-audit", help="write exact centre PMFs and a long 4x4 comparison"
    )
    exact_height.add_argument("--seed", type=int, default=20260908)
    exact_height.add_argument("--burnin-sweeps", type=int, default=10_000)
    exact_height.add_argument("--measurement-gap-sweeps", type=int, default=2)
    exact_height.add_argument("--num-measurements", type=int, default=100_000)
    exact_height.add_argument("--output", required=True)
    exact_height.set_defaults(function=command_exact_height_audit)

    mixing = subparsers.add_parser("mixing", help="compare two extremal starting states")
    mixing.add_argument("--L", type=int, required=True)
    mixing.add_argument("--model", choices=("uniform", "gamma"), default="gamma")
    mixing.add_argument("--gamma-shape", type=float, default=1.0)
    mixing.add_argument("--sweeps", type=int, required=True)
    mixing.add_argument("--record-every", type=int, default=5)
    mixing.add_argument("--seed", type=int, default=20260907)
    mixing.add_argument("--spatial-rhos", type=_comma_floats, default=DEFAULT_RELATIVE_SEPARATIONS)
    mixing.add_argument("--acf-max-lag", type=int, default=200)
    mixing.add_argument("--output")
    mixing.set_defaults(function=command_mixing)

    campaign = subparsers.add_parser(
        "mixing-campaign", help="run a local uniform/Gamma extremal-start campaign"
    )
    campaign.add_argument("--sizes", type=_comma_ints, required=True)
    campaign.add_argument("--gamma-environments", type=int, default=3)
    campaign.add_argument("--uniform-replicates", type=int, default=1)
    campaign.add_argument("--gamma-shape", type=float, default=1.0)
    campaign.add_argument("--sweeps", type=int, required=True)
    campaign.add_argument("--record-every", type=int, default=10)
    campaign.add_argument("--seed", type=int, default=20260908)
    campaign.add_argument(
        "--spatial-rhos", type=_comma_floats, default=DEFAULT_RELATIVE_SEPARATIONS
    )
    campaign.add_argument("--acf-max-lag", type=int, default=200)
    campaign.add_argument("--output", required=True)
    campaign.set_defaults(function=command_mixing_campaign)

    extension = subparsers.add_parser(
        "mixing-extend-frozen",
        help="extend chains in exact saved Gamma environments without redrawing",
    )
    extension.add_argument("--source", required=True)
    extension.add_argument("--output", required=True)
    extension.add_argument("--sweeps-by-size", type=_size_sweep_mapping, required=True)
    extension.add_argument("--record-every", type=int, default=10)
    extension.add_argument(
        "--spatial-rhos", type=_comma_floats, default=DEFAULT_RELATIVE_SEPARATIONS
    )
    extension.add_argument("--acf-max-lag", type=int, default=5_000)
    extension.set_defaults(function=command_mixing_extend_frozen)

    finalize = subparsers.add_parser(
        "mixing-finalize-assessment",
        help="merge manual visual assessments into an extended comparison table",
    )
    finalize.add_argument("--source", required=True)
    finalize.add_argument("--output", required=True)
    finalize.set_defaults(function=command_mixing_finalize)

    equivalence = subparsers.add_parser(
        "backend-equivalence",
        help="require checkpoint equality under identical supplied randomness",
    )
    equivalence.add_argument("--sizes", type=_comma_ints, default=(4, 6, 8, 12, 16))
    equivalence.add_argument("--attempts", type=int, default=20_000)
    equivalence.add_argument("--checkpoint-attempts", type=int, default=257)
    equivalence.add_argument("--seed", type=int, default=20260908)
    equivalence.add_argument("--output", required=True)
    equivalence.set_defaults(function=command_backend_equivalence)

    benchmark = subparsers.add_parser(
        "benchmark-backends",
        help="benchmark reference and Numba supplied-stream update kernels",
    )
    benchmark.add_argument("--sizes", type=_comma_ints, default=(8, 12, 16, 24, 32))
    benchmark.add_argument("--target-attempts", type=int, default=2_000_000)
    benchmark.add_argument("--repeats", type=int, default=3)
    benchmark.add_argument("--seed", type=int, default=20260908)
    benchmark.add_argument("--output", required=True)
    benchmark.set_defaults(function=command_benchmark_backends)

    ordering = subparsers.add_parser(
        "coupling-ordering-audit",
        help="test extremal face-height ordering after every common update",
    )
    ordering.add_argument("--sizes", type=_comma_ints, default=(4, 6, 8, 12, 16))
    ordering.add_argument("--attempts", type=int, default=20_000)
    ordering.add_argument("--streams", type=int, default=2)
    ordering.add_argument("--seed", type=int, default=20260908)
    ordering.add_argument("--output", required=True)
    ordering.set_defaults(function=command_ordering_audit)

    stress = subparsers.add_parser(
        "accelerated-stress",
        help="run an accelerated independent-chain stress diagnostic",
    )
    stress.add_argument("--source-case", required=True)
    stress.add_argument("--sweeps", type=int, default=1_000_000)
    stress.add_argument("--record-every", type=int, default=10)
    stress.add_argument(
        "--spatial-rhos", type=_comma_floats, default=DEFAULT_RELATIVE_SEPARATIONS
    )
    stress.add_argument("--acf-max-lag", type=int, default=20_000)
    stress.add_argument("--output", required=True)
    stress.set_defaults(function=command_accelerated_stress)

    coupling = subparsers.add_parser(
        "common-randomness-coupling",
        help="run a non-independent coupled diagnostic after an ordering audit",
    )
    coupling.add_argument("--source-case", required=True)
    coupling.add_argument("--ordering-audit", required=True)
    coupling.add_argument("--sweeps", type=int, default=1_000_000)
    coupling.add_argument("--record-every", type=int, default=10)
    coupling.add_argument("--seed", type=int, default=20260908)
    coupling.add_argument("--output", required=True)
    coupling.set_defaults(function=command_common_randomness_coupling)

    extremality = subparsers.add_parser(
        "exact-extremality-audit",
        help="enumerate tiny grids and identify true pointwise height extrema",
    )
    extremality.add_argument("--output", required=True)
    extremality.set_defaults(function=command_exact_extremality_audit)

    monotonicity = subparsers.add_parser(
        "exhaustive-monotonicity-audit",
        help="test every ordered tiny-grid pair, face, and heat-bath regime",
    )
    monotonicity.add_argument("--output", required=True)
    monotonicity.set_defaults(function=command_exhaustive_monotonicity_audit)

    prerequisites = subparsers.add_parser(
        "cftp-prerequisite-audit",
        help="write the conservative mathematical/CFTP prerequisite table",
    )
    prerequisites.add_argument("--extremality-summary", required=True)
    prerequisites.add_argument("--monotonicity-summary", required=True)
    prerequisites.add_argument("--output", required=True)
    prerequisites.set_defaults(function=command_cftp_prerequisites)

    local_monotonicity = subparsers.add_parser(
        "local-monotonicity-audit",
        help="exhaust the complete local height-order update truth table",
    )
    local_monotonicity.add_argument("--output", required=True)
    local_monotonicity.set_defaults(function=command_local_monotonicity_audit)

    exact_perfect = subparsers.add_parser(
        "exact-perfect-sampling",
        help="compare certified CFTP samples with the exact weighted 4x4 law",
    )
    exact_perfect.add_argument("--samples", type=int, default=20_000)
    exact_perfect.add_argument("--seed", type=int, default=20260909)
    exact_perfect.add_argument("--max-horizon-sweeps", type=int, default=16_384)
    exact_perfect.add_argument("--output", required=True)
    exact_perfect.set_defaults(function=command_exact_perfect_sampling)

    cftp_hard = subparsers.add_parser(
        "cftp-hard-calibration",
        help="measure certified true-bound CFTP horizons in the frozen hard L=16 case",
    )
    cftp_hard.add_argument("--source-case", required=True)
    cftp_hard.add_argument("--samples", type=int, default=50)
    cftp_hard.add_argument("--seed", type=int, default=20260909)
    cftp_hard.add_argument("--max-horizon-sweeps", type=int, default=1_048_576)
    cftp_hard.add_argument("--output", required=True)
    cftp_hard.set_defaults(function=command_cftp_hard_calibration)

    perfect_pair = subparsers.add_parser(
        "perfect-pair-smoke",
        help="generate two independent certified samples in each small fixed environment",
    )
    perfect_pair.add_argument("--sizes", type=_comma_ints, default=(4, 8, 12, 16))
    perfect_pair.add_argument("--seed", type=int, default=20260909)
    perfect_pair.add_argument("--output", required=True)
    perfect_pair.set_defaults(function=command_perfect_pair_smoke)

    perfect_manifest = subparsers.add_parser(
        "create-perfect-campaign-manifest",
        help="expand a production schedule into deterministic disjoint task rows",
    )
    perfect_manifest.add_argument("--schedule", required=True)
    perfect_manifest.add_argument("--campaign-id", required=True)
    perfect_manifest.add_argument("--seed", type=int, required=True)
    perfect_manifest.add_argument(
        "--backend", choices=("reference", "numba"), default="numba"
    )
    perfect_manifest.add_argument("--output", required=True)
    perfect_manifest.set_defaults(function=command_create_perfect_campaign_manifest)

    perfect_task = subparsers.add_parser(
        "perfect-campaign-task",
        help="run or resume one paired perfect-sampling manifest task",
    )
    perfect_task.add_argument("--manifest", required=True)
    perfect_task.add_argument("--task-id", type=int, required=True)
    perfect_task.add_argument("--max-horizon-sweeps", type=int)
    perfect_task.add_argument("--output", required=True)
    perfect_task.set_defaults(function=command_perfect_campaign_task)

    perfect_aggregate = subparsers.add_parser(
        "aggregate-perfect-campaign",
        help="audit expected perfect pairs and write analysis-ready environment tables",
    )
    perfect_aggregate.add_argument("--manifest", required=True)
    perfect_aggregate.add_argument("--output", required=True)
    perfect_aggregate.add_argument("--allow-incomplete", action="store_true")
    perfect_aggregate.set_defaults(function=command_aggregate_perfect_campaign)

    cftp_cost = subparsers.add_parser(
        "cftp-cost-calibration",
        help="run bounded multi-environment and larger-size CFTP cost probes",
    )
    cftp_cost.add_argument("--standard-sizes", type=_comma_ints, default=(8, 12, 16))
    cftp_cost.add_argument("--standard-gamma-environments", type=int, default=30)
    cftp_cost.add_argument("--standard-uniform-environments", type=int, default=1)
    cftp_cost.add_argument("--standard-samples-per-environment", type=int, default=2)
    cftp_cost.add_argument("--probe-sizes", type=_comma_ints, default=(24, 32))
    cftp_cost.add_argument("--probe-gamma-environments", type=int, default=5)
    cftp_cost.add_argument("--probe-uniform-environments", type=int, default=1)
    cftp_cost.add_argument("--probe-samples-per-environment", type=int, default=1)
    cftp_cost.add_argument("--gamma-shape", type=float, default=1.0)
    cftp_cost.add_argument("--seed", type=int, default=20260909)
    cftp_cost.add_argument("--max-horizon-sweeps", type=int, default=1_048_576)
    cftp_cost.add_argument("--backend", choices=("reference", "numba"), default="numba")
    cftp_cost.add_argument("--output", required=True)
    cftp_cost.set_defaults(function=command_cftp_cost_calibration)

    cftp_extend = subparsers.add_parser(
        "cftp-extend-censored",
        help="extend only capped CFTP streams with identical saved arrays and histories",
    )
    cftp_extend.add_argument("--source", required=True)
    cftp_extend.add_argument("--max-horizon-sweeps", type=int, required=True)
    cftp_extend.add_argument("--backend", choices=("reference", "numba"), default="numba")
    cftp_extend.add_argument("--output", required=True)
    cftp_extend.set_defaults(function=command_cftp_extend_censored)

    cftp_consolidate = subparsers.add_parser(
        "cftp-consolidate-calibration",
        help="combine capped and extended CFTP streams into final cost tables",
    )
    cftp_consolidate.add_argument("--source", required=True)
    cftp_consolidate.add_argument("--extension", required=True)
    cftp_consolidate.add_argument("--output", required=True)
    cftp_consolidate.set_defaults(function=command_cftp_consolidate)

    final_hard = subparsers.add_parser(
        "final-hard-calibration",
        help="run the exact frozen hard L=16 environment to 3m sweeps per chain",
    )
    final_hard.add_argument("--source-case", required=True)
    final_hard.add_argument("--one-million-case", required=True)
    final_hard.add_argument("--record-every", type=int, default=10)
    final_hard.add_argument("--acf-max-lag", type=int, default=100_000)
    final_hard.add_argument("--windows", type=int, default=16)
    final_hard.add_argument(
        "--spatial-rhos", type=_comma_floats, default=DEFAULT_RELATIVE_SEPARATIONS
    )
    final_hard.add_argument("--output", required=True)
    final_hard.set_defaults(function=command_final_hard_calibration)

    final_summary = subparsers.add_parser(
        "final-audit-summary",
        help="write the post-visual-review RESULTS_SUMMARY.txt",
    )
    final_summary.add_argument("--root", required=True)
    final_summary.add_argument(
        "--classification",
        choices=(
            "CREDIBLE FOR PILOT CALIBRATION",
            "BORDERLINE — EXTEND FURTHER",
            "NOT CREDIBLE — ORDINARY GLAUBER STILL TOO SLOW",
        ),
        required=True,
    )
    final_summary.add_argument(
        "--recommendation",
        choices=(
            "A. proceed to multi-environment mixing calibration with accelerated ordinary Glauber",
            "B. prototype monotone CFTP before paired sampling",
            "C. continue investigating mixing / sampling before any scaling study",
        ),
        required=True,
    )
    final_summary.add_argument("--late-window-assessment", required=True)
    final_summary.add_argument("--visual-assessment", required=True)
    final_summary.add_argument("--test-results", required=True)
    final_summary.set_defaults(function=command_final_audit_summary)

    pilot = subparsers.add_parser("pilot", help="run paired chains in shared environments")
    pilot.add_argument("--sizes", type=_comma_ints, required=True)
    pilot.add_argument("--models", type=_comma_strings, default=("uniform", "gamma"))
    pilot.add_argument("--gamma-shape", type=float, default=1.0)
    pilot.add_argument("--environments", type=int, required=True)
    pilot.add_argument("--burnin-sweeps", type=int, required=True)
    pilot.add_argument("--measurement-gap-sweeps", type=int, required=True)
    pilot.add_argument("--num-measurements", type=int, default=1)
    pilot.add_argument("--seed", type=int, default=20260907)
    pilot.add_argument("--spatial-rhos", type=_comma_floats, default=DEFAULT_RELATIVE_SEPARATIONS)
    pilot.add_argument("--output", required=True)
    pilot.set_defaults(function=command_pilot)

    analyze = subparsers.add_parser("analyze", help="analyze pilot output by environment")
    analyze.add_argument("--input", required=True)
    analyze.add_argument("--output", required=True)
    analyze.add_argument("--bootstrap-replicates", type=int, default=500)
    analyze.add_argument("--bootstrap-seed", type=int, default=20260908)
    analyze.set_defaults(function=command_analyze)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.function(args))
    except (TypeError, ValueError) as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
