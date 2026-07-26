DOUBLE-DIMER PAIR PRODUCTION ARCHIVE
====================================

batch-results/
    3,790 atomic CSV files.
    758,000 successful raw walks.
    Two independent walk seeds in each of 379,000 unique environments.

analysis/
    validation.csv           - completeness audit for all 3,790 tasks
    combined_raw.csv         - all 758,000 validated raw walk rows
    double_dimer_pairs.csv   - all 379,000 winding-difference pairs
    summary.csv              - 137 distribution-by-size cells
    loglog_fits.csv          - 15 effective-exponent fits
    scaling_model_comparison.csv
    pointwise_ratios.csv
    local_effective_exponents.csv

Frozen configuration:
    ../../../research-julia/configs/double_dimer_reproduction.csv

Primary prose report:
    ../../reports/double-dimer/README.md
