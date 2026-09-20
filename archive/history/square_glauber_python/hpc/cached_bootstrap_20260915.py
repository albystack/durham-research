"""Equivalent environment bootstrap with observable lookups cached once."""
import numpy as np
import pandas as pd

def bootstrap(center_environment, spatial_environment, replicates, seed):
    from square_glauber.analysis import GROUP_COLUMNS, _pair_statistics
    rng = np.random.Generator(np.random.PCG64(seed))
    groups = []
    for keys, center in center_environment.groupby(GROUP_COLUMNS, dropna=False, sort=True):
        base = dict(zip(GROUP_COLUMNS, keys, strict=True))
        spatial_groups = []
        if not spatial_environment.empty:
            selected = spatial_environment[
                (spatial_environment.model == keys[0]) & (spatial_environment.L == keys[2]) &
                (spatial_environment.gamma_shape.isna() if pd.isna(keys[1]) else spatial_environment.gamma_shape == keys[1])
            ]
            for (rho, r), frame in selected.groupby(['requested_rho', 'r'], sort=True):
                assert frame.environment_id.is_unique
                assert set(frame.environment_id) == set(center.environment_id)
                aligned = frame.set_index('environment_id').loc[center.environment_id]
                spatial_groups.append((dict(requested_rho=rho, r=r), aligned[['DeltaH1_environment_mean', 'DeltaH2_environment_mean']].to_numpy()))
        groups.append((base, center[['H1_environment_mean', 'H2_environment_mean']].to_numpy(), spatial_groups))
    center_rows, spatial_rows = [], []
    for replicate in range(replicates):
        for base, values, spatial_groups in groups:
            chosen = rng.integers(0, len(values), size=len(values))
            selected = values[chosen]
            center_rows.append(base | _pair_statistics(selected[:, 0], selected[:, 1]) | dict(replicate=replicate))
            for spatial_key, spatial_values in spatial_groups:
                selected = spatial_values[chosen]
                spatial_rows.append(base | spatial_key | _pair_statistics(selected[:, 0], selected[:, 1]) | dict(replicate=replicate))
        if (replicate + 1) % 1000 == 0:
            print(f'Bootstrap {replicate + 1}/{replicates}', flush=True)
    return pd.DataFrame(center_rows), pd.DataFrame(spatial_rows)
