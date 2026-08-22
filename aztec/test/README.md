# Test suite

Run all unit and mathematical-reference checks with:

```bash
julia --project=aztec -e 'using Pkg; Pkg.test()'
```

Run the command-line smoke workflows with:

```bash
sh aztec/test/smoke_workflows.sh
```

The suite covers exact tiny tiling counts, weighted recurrence reference
values, detailed balance, height invariants, tree/dual-tree/matching structure,
seed reproducibility, paired-environment identities, optimized-vs-reference
agreement, replica-exchange scheduling, environment-blocked estimators, and
end-to-end CSV/SVG workflows.

Tests deliberately use small systems. Production mixing is assessed from
retained diagnostics and cannot be certified by unit tests alone.
