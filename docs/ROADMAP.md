# Research roadmap

This file records the next scientifically useful work. Completed numerical
results and historical pilots belong in [`RESULTS.md`](RESULTS.md).

## Current conclusion

The direct frozen-edge square-grid Glauber campaign does not provide robust
evidence for a positive \((\log L)^2\) coefficient over \(L=2\)–20. The
all-one control behaves correctly and the variance decomposition closes, but
the largest sizes have unresolved mixing diagnostics.

The next objective is therefore **frontier validation**, not a larger general
production campaign.

## P0 — Freeze and audit the completed campaign

- Preserve the raw core, frontier, and control roots without manual edits.
- Record campaign configs, base seeds, code revision, file counts, and SHA-256
  manifests.
- Retain all-environment estimates as the primary analysis.
- Label start-gap filtering as sensitivity analysis only.
- Recompute a conservative core-only \(L=2\)–12 fit that does not depend on the
  unresolved frontier mixing.

Completion gate: a fresh checkout plus the external raw roots can regenerate
the environment-block and scaling tables without undocumented steps.

## P1 — Measure replica flow directly

Extend the parallel-tempering diagnostics without changing the target model:

- attach a persistent label to each replica;
- record temperature-index occupancy and transitions;
- count completed cold-to-hot-to-cold round trips;
- retain adjacent-pair attempts and acceptances;
- report label-flow diagnostics separately for both independent chains.

Validate the instrumentation on tiny exactly enumerable systems. The added
bookkeeping must not affect configurations, acceptance decisions, seeds, or
the \(\beta=1\) observable.

## P2 — Repair the hot-end exchange bottleneck

At \(L=20\), pooled acceptance for the \(\beta=0\leftrightarrow0.05\) pair is
approximately 0.062, while target-adjacent exchange is healthy. Use a pilot to
construct a denser ladder near \(\beta=0\), keeping:

- the same frozen i.i.d. mean-one edge environment;
- the same local heat-bath rule;
- the same \(\beta=1\) target distribution;
- the same two-chain extremal-start comparison;
- deterministic, recorded seeds.

Do not optimize the ladder by changing the disorder law, central-height
observable, pairing, or estimand.

## P3 — Fresh-seed frontier pilot

Run a small, restart-safe validation campaign at \(L=16,20\) only. Predeclare
the schedule and diagnostic gate before inspecting the estimates.

Required evidence:

- replica flow across the full ladder rather than only local target exchange;
- no common signed extremal-start bias;
- adequate chain ESS for the retained central-height trace;
- stable conditional, disorder, and total components across the two starts;
- agreement with exact/reference behaviour wherever a smaller analogue is
  available;
- no errors, incomplete batches, or missing provenance.

If the pilot fails, lengthen or redesign the ladder and repeat the pilot. Do
not compensate by pooling more poorly mixed environments.

## P4 — Independent frontier replication

Only after the P3 gate passes, repeat \(L=16,20\) with a fresh campaign seed.
Analyze it independently before combining anything with the completed run.

- If the replication agrees, update the full-range scaling comparison and
  report the improved mixing evidence.
- If it differs materially, exclude the original frontier from the primary
  analysis and retain \(L\le12\) as the defensible range.

## P5 — Cross-model comparison

Once the direct square-grid frontier is validated, compare Aztec and
square-grid results through one common schema:

```text
model, geometry, L, disorder law, environment seed, replica seeds,
observable, conditional component, disorder component, diagnostics
```

Use the same environment-blocked fitting implementation and report:

- coefficient sign and uncertainty;
- BIC sign convention;
- prediction error;
- lower-size cutoff sensitivity;
- unweighted and precision-weighted fits;
- conditional, disorder, and total components separately;
- negative controls.

The comparison must distinguish genuine geometry/model differences from
finite-size effects and from differences in how disorder is represented.

## Release checklist

- deterministic replayable seeds;
- one config file per retained campaign;
- restart-safe jobs and atomic outputs;
- raw-data manifests and checksums;
- scripts that regenerate retained tables and figures;
- tests for mathematical invariants and optimized observables;
- no private correspondence, credentials, or machine-specific personal paths;
- concise public documentation of negative results and limitations.
