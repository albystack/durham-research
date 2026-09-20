"""Optional CFTP restart at completed doubling-horizon boundaries.

No endpoint is composed through a longer past. Each new horizon replays its
entire unchanged random history from the true bounds. Interrupted horizons are
replayed; this deliberately does NOT promise mid-horizon checkpoint recovery.
Callers must hold a single-writer lock. Legacy censor records are trusted local
provenance, not independently verifiable certificates of noncoalescence.
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from .glauber import UpdateCounts
from .height import face_heights
from .perfect import CFTPIteration, CFTPNotCoalesced, PerfectSample, _random_block, _states_equal
from .perfect_campaign import _atomic_json
from .mixing_extension import environment_hashes


def _digest(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def sample_with_checkpoint(sampler, seed, *, checkpoint_path, max_horizon_sweeps,
                           legacy_censor_path=None):
    """Canonical 1,2,4,... CFTP, preserving completed failed horizon work."""
    cap = int(max_horizon_sweeps)
    if cap != max_horizon_sweeps or cap <= 0 or cap & (cap - 1):
        raise ValueError("restart requires a positive power-of-two horizon cap")
    path = Path(checkpoint_path)
    identity = dict(cftp_seed=int(seed), L=sampler.environment.L,
                    environment_sha256=environment_hashes(sampler.environment)["environment_sha256"])
    rows = []
    lower_counts, upper_counts = UpdateCounts(), UpdateCounts()
    elapsed = 0.0
    info = dict(checkpoint_granularity="completed_failed_horizon",
                update_count_scope="all_replayed_horizons",
                omitted_logical_attempts_per_bound=0)
    if path.exists():
        envelope = json.loads(path.read_text())
        data = envelope["payload"]
        if envelope["sha256"] != _digest(data) or data["schema"] != 1:
            raise ValueError("invalid checkpoint checksum/schema")
        if data["identity"] != identity:
            raise ValueError("checkpoint environment/seed mismatch")
        rows = [CFTPIteration(**row) for row in data["iterations"]]
        lower_counts = UpdateCounts(**data["lower_counts"])
        upper_counts = UpdateCounts(**data["upper_counts"])
        elapsed, info = data["elapsed_seconds"], data["restart_info"]
    elif legacy_censor_path is not None:
        data = json.loads(Path(legacy_censor_path).read_text())
        if (data.get("cftp_seed") != int(seed) or
            data.get("environment_sha256") != identity["environment_sha256"] or
            data.get("certified_perfect_sample") is not False or data.get("censored") is not True):
            raise ValueError("legacy censor identity/status mismatch")
        rows = [CFTPIteration(**row) for row in data["iterations"]]
        if not rows or rows[-1].horizon_sweeps != data["horizon_sweeps"]:
            raise ValueError("legacy horizon mismatch")
        info.update(update_count_scope="since_legacy_import",
                    omitted_logical_attempts_per_bound=(2 * rows[-1].horizon_sweeps - 1) * (sampler.environment.L - 1)**2,
                    omitted_flippable_and_change_counts="unknown",
                    legacy_source_sha256=hashlib.sha256(Path(legacy_censor_path).read_bytes()).hexdigest())
    for index, row in enumerate(rows):
        if (row.horizon_sweeps != 2**index or row.block_count != index+1 or
            row.lower_upper_coalesced_at_time_zero or row.maximum_pointwise_height_gap_at_time_zero <= 0):
            raise ValueError("checkpoint must contain consecutive failed canonical horizons")
    for counts in (lower_counts, upper_counts):
        if not 0 <= counts.actual_changes <= counts.flippable_proposals <= counts.attempted_updates:
            raise ValueError("invalid checkpoint counters")
    def censored():
        error = CFTPNotCoalesced(rows[-1].horizon_sweeps, tuple(rows))
        error.restart_info = info
        raise error
    if rows and rows[-1].horizon_sweeps >= cap:
        censored()
    horizon = 2 * rows[-1].horizon_sweeps if rows else 1
    started = perf_counter()
    while True:
        count = horizon.bit_length()
        blocks = [_random_block(int(seed), sampler.environment.L, i, 1 if i == 0 else 2**(i-1))
                  for i in reversed(range(count))]
        lower, upper = sampler.minimum.state.copy(), sampler.maximum.state.copy()
        sampler._apply_random_blocks(lower, upper, blocks, lower_counts, upper_counts, verify_order=True)
        lh, uh = face_heights(lower), face_heights(upper)
        if not np.all(lh <= uh):
            raise AssertionError("CFTP bounds lost height order")
        equal = _states_equal(lower, upper)
        rows.append(CFTPIteration(horizon, count, equal, int(np.max(uh-lh))))
        if equal:
            lower.validate()
            result = PerfectSample(lower, int(seed), sampler.backend, True, horizon,
                tuple(rows), tuple(blocks), lower_counts, upper_counts,
                elapsed + perf_counter()-started, sampler.minimum.directed_flips,
                sampler.maximum.directed_flips)
            result.restart_info = info
            return result
        payload = dict(schema=1, identity=identity, iterations=[asdict(row) for row in rows],
            lower_counts=asdict(lower_counts), upper_counts=asdict(upper_counts),
            elapsed_seconds=elapsed+perf_counter()-started, restart_info=info)
        _atomic_json(path, dict(payload=payload, sha256=_digest(payload)))
        if horizon >= cap:
            censored()
        horizon *= 2
