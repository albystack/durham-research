STRICT-ANNEALED PRODUCTION ARCHIVE
==================================

batch-results/
    3,790 atomic CSV files.
    379,000 successful walks.
    One unique independently seeded environment per walk.

analysis/
    validation.csv       - completeness audit for all 3,790 tasks
    combined_raw.csv     - all 379,000 validated walk rows
    summary.csv          - 137 distribution-by-size cells
    loglog_fits.csv      - 15 effective-exponent fits
    scaling_model_comparison.csv
    pointwise_ratios.csv
    local_effective_exponents.csv

Frozen configuration:
    ../../../research-julia/configs/strict_annealed_reproduction.csv

Primary prose report:
    ../../reports/strict-annealed/README.md
