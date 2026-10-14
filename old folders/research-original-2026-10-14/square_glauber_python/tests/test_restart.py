from dataclasses import asdict
import json
import numpy as np
import pytest

from square_glauber.environment import sample_environment
from square_glauber.accelerated import NUMBA_AVAILABLE
from square_glauber.mixing_extension import environment_hashes
from square_glauber.perfect import MonotoneCFTPSampler, CFTPNotCoalesced
from square_glauber.restart import sample_with_checkpoint


@pytest.mark.parametrize("backend", ["reference", "numba"])
@pytest.mark.parametrize("model", ["uniform", "gamma"])
def test_restart_matches_original_with_different_chunks(tmp_path, backend, model):
    if backend == "numba" and not NUMBA_AVAILABLE:
        pytest.skip("optional numba absent")
    env = sample_environment(4, np.random.default_rng(53), model=model)
    sampler = MonotoneCFTPSampler(env, backend=backend, chunk_attempts=13)
    expected = sampler.sample(91, max_horizon_sweeps=16384)
    path = tmp_path / "checkpoint.json"
    with pytest.raises(CFTPNotCoalesced):
        sample_with_checkpoint(sampler, 91, checkpoint_path=path, max_horizon_sweeps=1)
    sampler.chunk_attempts = 47
    actual = sample_with_checkpoint(sampler, 91, checkpoint_path=path, max_horizon_sweeps=16384)
    np.testing.assert_array_equal(actual.state.horizontal_occupied, expected.state.horizontal_occupied)
    np.testing.assert_array_equal(actual.state.vertical_occupied, expected.state.vertical_occupied)
    assert actual.iterations == expected.iterations
    assert actual.blocks == expected.blocks
    assert actual.lower_update_counts == expected.lower_update_counts
    assert actual.upper_update_counts == expected.upper_update_counts
    with pytest.raises(ValueError, match="seed mismatch"):
        sample_with_checkpoint(sampler, 92, checkpoint_path=path, max_horizon_sweeps=16384)
    path.write_text(path.read_text().replace('"schema": 1', '"schema": 2'))
    with pytest.raises(ValueError):
        sample_with_checkpoint(sampler, 91, checkpoint_path=path, max_horizon_sweeps=16384)


def test_legacy_import_and_no_endpoint_composition(tmp_path):
    env = sample_environment(4, np.random.default_rng(53), model="gamma")
    sampler = MonotoneCFTPSampler(env, backend="reference", chunk_attempts=31)
    with pytest.raises(CFTPNotCoalesced) as caught:
        sampler.sample(91, max_horizon_sweeps=1)
    legacy = tmp_path / "legacy.json"
    legacy.write_text(json.dumps(dict(cftp_seed=91, **environment_hashes(env),
        certified_perfect_sample=False, censored=True, horizon_sweeps=1,
        iterations=[asdict(row) for row in caught.value.iterations])))
    actual = sample_with_checkpoint(sampler, 91, checkpoint_path=tmp_path/"new.json",
                                    max_horizon_sweeps=16384, legacy_censor_path=legacy)
    expected = sampler.sample(91, max_horizon_sweeps=16384)
    np.testing.assert_array_equal(actual.state.horizontal_occupied, expected.state.horizontal_occupied)
    assert actual.iterations == expected.iterations
    assert actual.lower_update_counts.attempted_updates + 9 == expected.lower_update_counts.attempted_updates
    assert actual.restart_info["update_count_scope"] == "since_legacy_import"
    with pytest.raises(ValueError, match="power-of-two"):
        sample_with_checkpoint(sampler, 91, checkpoint_path=tmp_path/"x", max_horizon_sweeps=7)
