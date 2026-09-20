"""Bounded execution benchmark, not a scientific sampling campaign.

Run on an allocated compute core. Uses the unchanged CFTP replay implementation,
checks identical endpoint occupancies and counters across chunk partitions, and
profiles the default replay. Prints JSON for archival in Slurm stdout.
"""
import argparse
import cProfile
import hashlib
import io
import json
import pstats
import socket
import time
from dataclasses import asdict

import numpy as np

from square_glauber.environment import DimerEnvironment
from square_glauber.glauber import UpdateCounts
from square_glauber.perfect import MonotoneCFTPSampler, _random_block


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--environment', required=True)
    parser.add_argument('--L', type=int, default=32)
    parser.add_argument('--sweeps', type=int, default=8192)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    with np.load(args.environment) as weights:
        env = DimerEnvironment(args.L, weights['horizontal_weights'], weights['vertical_weights'])
    baseline = None
    rows = []
    # Warm compilation; this does not enter benchmark timings.
    warm = MonotoneCFTPSampler(env)
    warm._apply_random_blocks(warm.minimum.state.copy(), warm.maximum.state.copy(),
                             [_random_block(123, args.L, 0, 1)], UpdateCounts(), UpdateCounts(), verify_order=True)
    for repeat in range(args.repeats):
        for chunk in ([65536, 262144, 1000000, 4000000] if repeat % 2 == 0 else [4000000, 1000000, 262144, 65536]):
            sampler = MonotoneCFTPSampler(env, chunk_attempts=chunk)
            lower, upper = sampler.minimum.state.copy(), sampler.maximum.state.copy()
            lc, uc = UpdateCounts(), UpdateCounts()
            started = time.perf_counter()
            sampler._apply_random_blocks(lower, upper, [_random_block(20260912, args.L, 0, args.sweeps)], lc, uc, verify_order=True)
            elapsed = time.perf_counter() - started
            arrays = [lower.horizontal_occupied, lower.vertical_occupied, upper.horizontal_occupied, upper.vertical_occupied]
            signature = (hashlib.sha256(b''.join(a.tobytes() for a in arrays)).hexdigest(), asdict(lc), asdict(uc))
            if baseline is None:
                baseline = signature
            assert signature == baseline, 'Chunk partition changed endpoint or update counts'
            lower.validate(); upper.validate()
            rows.append(dict(repeat=repeat, chunk_attempts=chunk, seconds=elapsed,
                             coupled_instructions_per_second=args.sweeps*(args.L-1)**2/elapsed))
    profile = cProfile.Profile()
    sampler = MonotoneCFTPSampler(env)
    profile.enable()
    sampler._apply_random_blocks(sampler.minimum.state.copy(), sampler.maximum.state.copy(),
                                [_random_block(20260912, args.L, 0, args.sweeps)], UpdateCounts(), UpdateCounts(), verify_order=True)
    profile.disable()
    stream = io.StringIO()
    pstats.Stats(profile, stream=stream).sort_stats('cumulative').print_stats(18)
    print(json.dumps(dict(host=socket.gethostname(), L=args.L, sweeps=args.sweeps,
                         repeats=args.repeats, endpoint_and_counts_equivalent=True, rows=rows,
                         profile=stream.getvalue()), indent=2), flush=True)


if __name__ == '__main__':
    main()
