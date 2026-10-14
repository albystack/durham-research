"""Fixed numerical/geometry audit before any embedded-square campaign.

Compare every potential creation probability, not only realised coin outcomes,
using Float64 versus 80/160 decimal-digit reductions of identical stored
inputs. This is evidence on a fixed cohort, not a universal interval proof.
"""

import argparse
import json
from decimal import Decimal, localcontext
from pathlib import Path
from time import perf_counter

import numpy as np
from square_glauber.budgeted_frontier import peak_rss_bytes
from square_glauber.comparison import ComparisonModel, comparison_environment
from square_glauber.embedded_shuffle import VERSION, EmbeddedSquareSampler, diamond


def probabilities(sampler):
    result = []
    with localcontext() as ctx:
        ctx.prec = 160
        for n in sorted(sampler.levels):
            weights = sampler.levels[n]
            for ids in diamond(n).face_edges:
                a, b, c, d = [Decimal(w) for w in (weights[i] for i in ids)]
                result.append(a * c / (a * c + b * d))
    return result


def audit(output, sizes=(8, 16, 32, 64)):
    report = {
        "version": VERSION,
        "cases": [],
        "tolerances": {
            "max_probability_error": 1e-10,
            "two_replica_sum_abs_error": 1e-6,
            "reference_disagreement": 1e-65,
        },
        "caveat": "Fixed environments only; not a global rounding-error proof",
    }
    for L in sizes:
        env = comparison_environment(
            L, np.random.default_rng(2026100401 + L), ComparisonModel()
        )
        start = perf_counter()
        floating = EmbeddedSquareSampler(env, cast=float)
        setup = perf_counter() - start
        start = perf_counter()
        floating.sample(2026100402)
        floating.sample(2026100403)
        sampling = perf_counter() - start
        pf = probabilities(floating)
        with localcontext() as ctx:
            ctx.prec = 80
            medium = EmbeddedSquareSampler(env, cast=Decimal)
        pm = probabilities(medium)
        with localcontext() as ctx:
            ctx.prec = 160
            high = EmbeddedSquareSampler(env, cast=Decimal)
        ph = probabilities(high)
        assert len(pf) == len(pm) == len(ph)
        error = [abs(a - b) for a, b in zip(pf, ph)]
        reference = max(abs(a - b) for a, b in zip(pm, ph))
        case = {
            "L": L,
            "seed": 2026100401 + L,
            "coins_per_replica": len(pf),
            "max_probability_error": float(max(error)),
            "two_replica_sum_abs_error": float(2 * sum(error)),
            "reference_disagreement": float(reference),
            "float_setup_seconds": setup,
            "float_pair_seconds": sampling,
            "process_peak_rss_bytes": peak_rss_bytes(),
        }
        case["passed"] = all(case[k] < v for k, v in report["tolerances"].items())
        report["cases"].append(case)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(case), flush=True)
        if not case["passed"]:
            raise ArithmeticError("numerical audit tolerance failed")
        del floating, medium, high, pf, pm, ph
    report["all_passed"] = True
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.output)
