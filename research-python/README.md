# Archived Python prototype

This directory preserves the original Python implementation used during early
model development and for the `L=32` through `L=1024` simulations included in
the aggregate report.

The maintained implementation is the Julia package in
[`../research-julia/`](../research-julia/). New experiments should use that
package because it provides bounded memory, multithreading, restart-safe batch
execution, and the current statistical pipeline.

## Historical scope

The Python code contains:

- an undirected edge-weight prototype;
- the four-direction site-IID environment used by the historical batch runner;
- a directional drift model retained only for result reproducibility;
- command-line simulation, batch, and analysis utilities;
- unit tests for configuration, seeding, environments, and output validation.

These models are separate and should not be mixed without checking the
`environment_model` field in each result row.

## Reproduce the Python tests

Python 3.10 or newer is required.

```bash
cd research-python
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m unittest discover -s tests -v
```

## Historical smoke run

```bash
.venv/bin/python generate_hpc_config.py \
  --output configs/hpc_test.csv \
  --sizes 64 128 \
  --distributions baseline gamma:1.0 \
  --batches 1 \
  --num-environments 2 \
  --walks-per-environment 3

.venv/bin/python run_hpc_batch.py \
  --task-id 0 \
  --config configs/hpc_test.csv \
  --output-dir results_hpc_smoke
```

Generated outputs are intentionally excluded from version control. Compact
cross-language aggregate tables are published under [`../reports/`](../reports/).
