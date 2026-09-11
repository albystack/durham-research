# Local schedules

The validated v1 runner takes its complete schedule from explicit CLI flags;
no scientifically unjustified burn-in or size scaling is embedded in code.
Keep reviewed command lines or future small JSON schedule files here. Generated
measurements belong in the ignored `../outputs/` directory.

Suggested first smoke command:

```bash
.venv/bin/python -m square_glauber.cli pilot \
  --sizes 4,6,8 --models uniform,gamma --gamma-shape 1.0 \
  --environments 4 --burnin-sweeps 20 --measurement-gap-sweeps 5 \
  --num-measurements 2 --seed 20260907 --output outputs/smoke
```

This schedule is an integration smoke test, not evidence of mixing.
