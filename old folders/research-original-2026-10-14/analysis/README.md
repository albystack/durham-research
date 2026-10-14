# Analysis in two commands

From the repository root, run `make setup-python`, then:

```bash
make analysis     # fast: verify data, recompute points, use recorded intervals
make bootstrap    # full: independently rerun uncertainty and control test
```

`run.py` is the entry point. It checks `data/manifest.json`, loads the complete
measurements, rejects missing or duplicated pairs/probes, and computes results
for each domain separately. No observations are filtered by their values.

## What the measurements mean

For an environment ω, the two tilings are independent conditional on ω.
`value1` and `value2` are the same spatial height increment in the two replicas.
Across independent environments we estimate:

- `D = Cov(value1, value2)`: disorder-induced fluctuations.
- `T = Var(value1 − value2) / 2`: average conditional fluctuations.
- `V = T + D`: average of the two replica variances.

`point_table()` computes these sample moments. The shared, tested routines are
in `square_glauber_python/src/square_glauber/followup_analysis.py`.

## Fitting the squared-log term

The fit is `Dρ(r) = aρ + bρ log(r) + c (log(r))²`. Each probe fraction has its
own intercept and ordinary-log slope; both share one coefficient `c`.
Distances are actual Euclidean separations, and logarithms are natural.
Fits use equal weights, as fixed in the campaign plan.

The primary window is sizes **16–128**. The separately fixed sensitivity
windows are **32–128** and **48–128**. We report all three for each geometry.
The primary scientific finding is the coefficient for **each geometry**.
The original additional domain contrast and all control diagnostics remain
in the output for completeness.

## Bootstrap and controls

The full command resamples whole paired environments within each domain,
model and size. A resampled environment carries both replicas and both probes;
these observations must not be resampled separately. Different cells have
independent resampling streams. The seeds and 10,000 repetitions are fixed
in `data/plan.json`; intervals are pointwise 95% percentile intervals.

The control diagnostic permutes complete second-replica probe vectors within
each deterministic control cell. It evaluates the maximum standardized absolute
curvature contrast over two controls and three windows, with 9,999 permutations.
This reuses the previously validated implementation in `larger_campaign.py`.
A large joint p-value does not prove sampler correctness.

## Outputs

| File in `output/analysis/` or `output/recomputed/` | Contents |
|---|---|
| `coefficients.csv` | 18 coefficients: both domains, three models, three windows |
| `moments.csv` | 84 rows of T, D and V, with recorded or recomputed intervals |
| `analysis.json` | Fits, cutoff results, original contrasts and control test |
| `gamma_covariance.png` | Both Gamma covariance curves in the same physical units |

The quick output explicitly labels intervals as recorded Hamilton values.
`--bootstrap 10000` recomputes them from the published paired observations.
Other repetition counts are exploratory and do not reproduce the frozen plan.

To regenerate the public data bundle from the original verified compact archive:

```bash
python analysis/export_data.py path/to/results-compact-20261010.tar.gz
```

The export checks the archive and original table checksums and gives pair IDs
a batch prefix. It never changes any height measurement. The original archive
is retained locally and on Hamilton.
